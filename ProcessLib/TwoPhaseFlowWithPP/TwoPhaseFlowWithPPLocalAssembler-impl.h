/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TwoPhaseFlowWithPPLocalAssembler.cpp
 *
 * Created on October 19, 2016, 2:28 PM
 */

#ifndef OGS_TWOPHASEFLOWWITHPPLOCALASSEMBLER_IMPL_H
#define OGS_TWOPHASEFLOWWITHPPLOCALASSEMBLER_IMPL_H

#include <iostream>
#include "TwoPhaseFlowWithPPLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseFlowWithPPLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, std::vector<double> const& local_x,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();

    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
    int const n_nodes = ShapeFunction::NPOINTS;
    Eigen::MatrixXd mass_mat_coeff =
        Eigen::MatrixXd::Zero(NUM_NODAL_DOF, NUM_NODAL_DOF);
    Eigen::MatrixXd K_mat_coeff =
        Eigen::MatrixXd::Zero(NUM_NODAL_DOF, NUM_NODAL_DOF);
    Eigen::VectorXd H_vec_coeff = Eigen::VectorXd::Zero(NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<NodalVectorType>(
        local_b_data, local_matrix_size);

    typedef Eigen::Matrix<double, n_nodes, n_nodes> MatrixNN;
    MatrixNN _Mgp;
    MatrixNN _Mgpc;
    MatrixNN _Mlp;
    MatrixNN _Mlpc;
    MatrixNN _Kgp;
    MatrixNN _Kgpc;
    MatrixNN _Klp;
    MatrixNN _Klpc;

    typedef Eigen::Matrix<double, n_nodes, 1> VectorNN;
    VectorNN _Bg;
    VectorNN _Bl;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    _process_data._material->setMaterialID(pos);

    const Eigen::MatrixXd& perm = _process_data._material->getPermeability(
        t, pos, _element.getDimension());

    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to perm.rows() == _element->getDimension()
    assert(perm.rows() == GlobalDim || perm.rows() == 1);
    // x indicates the molar fraction
    MathLib::PiecewiseLinearInterpolation const& interpolated_x_co2_in_wet =
        _process_data._interpolated_co2_solubility;
    MathLib::PiecewiseLinearInterpolation const& interpolated_x_h2o_in_nonwet =
        _process_data._interpolated_h2o_solubiity;
    MathLib::PiecewiseLinearInterpolation const& interpolated_rho_co2_nonwet =
        _process_data._interpolated_co2_density;
    MathLib::PiecewiseLinearInterpolation const& interpolated_mu_co2 =
        _process_data._interpolated_co2_viscosity;
    double porosity_variable = 0.;
    double storage_variable = 0.;
    // Note: currently only isothermal case is considered, so the temperature is
    // assumed to be const
    // the variation of temperatura will be taken into account in future
    _temperature = 313.15;
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);

        double pc_int_pt = 0.;
        double pg_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pg_int_pt, pc_int_pt);

        // TODO : compute _temperature from the heat transport pcs
        double const pl = pg_int_pt - pc_int_pt;
        _pressure_wetting[ip] = pl;
        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();

        double const rho_w =
            _process_data._material->getLiquidDensity(pl, _temperature);

        double x_co2_in_wet = 0.0;
        double x_h2o_in_nonwet = 0.0;
        double x_co2_in_wet_dpn_plus = 0.0;
        double x_h2o_in_nonwet_dpn_plus = 0.0;
        double x_co2_in_wet_dpn_minus = 0.0;
        double x_h2o_in_nonwet_dpn_minus = 0.0;
		double x_co2_in_wet_henry = 0.0;
        /*! pure co2 density
        * from Span and Wagner EoS
        */
        double const rho_co2 = interpolated_rho_co2_nonwet.getValue(pg_int_pt);
        double const rho_co2_plus =
            interpolated_rho_co2_nonwet.getValue(pg_int_pt + eps*std::abs(pg_int_pt));
        double const rho_co2_minus =
            interpolated_rho_co2_nonwet.getValue(pg_int_pt - eps*std::abs(pg_int_pt));
        // comparison with reading curve
        double const rho_nonwet =
            interpolated_rho_co2_nonwet.getValue(pg_int_pt);
        /*! Mutual solubility
        * x_co2_in_wet molar fraction of CO2 dissolved in wetting phase
        * x_h2o_in_nonwet molar fraction of h2o dissolved in non-wetting phase
        * from Pruess EoS
        */
        _process_data._material->calculateMoleFractions(
            pg_int_pt, _temperature, 0.0, x_co2_in_wet, x_h2o_in_nonwet,
            rho_co2);
        // for derivative
        _process_data._material->calculateMoleFractions(
            pg_int_pt + eps*std::abs(pg_int_pt), _temperature, 0.0, x_co2_in_wet_dpn_plus,
            x_h2o_in_nonwet_dpn_plus, rho_co2_plus);
        _process_data._material->calculateMoleFractions(
            pg_int_pt - eps*std::abs(pg_int_pt), _temperature, 0.0, x_co2_in_wet_dpn_minus,
            x_h2o_in_nonwet_dpn_minus, rho_co2_minus);
		
		//---- For comparison from reading curve
		double const x_co2_in_wet_duan = _process_data._material->moleFracCO2InBrine_duan(_temperature, pg_int_pt, 0.0, rho_co2);
		double const x_co2_in_wet_test =
            interpolated_x_co2_in_wet.getValue(pg_int_pt);
        double const x_h2o_in_nonwet_test =
            interpolated_x_h2o_in_nonwet.getValue(pg_int_pt);
		//---- END of test----------------------
        double const rho_mol_w = rho_w / molar_mass_h2o / (1 - x_co2_in_wet);
        double const rho_mass_wet =
            rho_mol_w * (x_co2_in_wet * molar_mass_co2 +
                         (1 - x_co2_in_wet) * molar_mass_h2o);
        double const rho_mol_nonwet =
            rho_co2 / molar_mass_co2 / (1 - x_h2o_in_nonwet);
        double const rho_mass_nonwet =
            rho_mol_nonwet * ((1 - x_h2o_in_nonwet) * molar_mass_co2 +
                              x_h2o_in_nonwet * molar_mass_h2o);
        /* calculate derivatives based on finite diff
        * dx_co2_in_wet_dpn partial pressure
        * drho_n_dpn pure co2 density derivative w.r.t pressure nonwet
        */
        double const dx_co2_in_wet_dpn =
            (x_co2_in_wet_dpn_plus - x_co2_in_wet_dpn_minus) / 2 / eps/std::abs(pg_int_pt);
        double const dx_h2o_in_nonwet_dpn =
            (x_h2o_in_nonwet_dpn_plus - x_h2o_in_nonwet_dpn_minus) / 2 / eps/std::abs(pg_int_pt);

        //double const drho_n_dpn = (rho_co2_plus - rho_co2_minus) / 2 / eps;
		double const drho_n_dpn = interpolated_rho_co2_nonwet.getDerivative(pg_int_pt);

        double Sw = _process_data._material->getSaturation(pc_int_pt);  // pc
        
        _saturation[ip] = Sw;
        double dSwdPc = _process_data._material->getDerivSaturation(pc_int_pt);
		//if (pc_int_pt < 0)
			//dSwdPc = 0.0;
        double const poro = _process_data._material->getPorosity(
            t, pos, pg_int_pt, _temperature, porosity_variable);

        // component a -- (water component)
        // component b -- (non-wet component)
        // Assemble mass matrix, M
        // wet -- wetting phase
        // nonwetting
        mass_mat_coeff(nonwet_pressure_coeff_index,
                       nonwet_pressure_coeff_index) =
            poro * rho_mol_w * Sw * dx_co2_in_wet_dpn / (1 - x_co2_in_wet) +
            poro * (1 - Sw) * drho_n_dpn / molar_mass_co2;
        
		mass_mat_coeff(nonwet_pressure_coeff_index, cap_pressure_coeff_index) =
			poro * dSwdPc *
			(rho_mol_w*x_co2_in_wet - rho_co2 / molar_mass_co2);

        // wetting
		mass_mat_coeff(cap_pressure_coeff_index, nonwet_pressure_coeff_index) = 0.0;
            /*poro * (1 - Sw) *
            (rho_mol_nonwet * dx_h2o_in_nonwet_dpn / (1 - x_h2o_in_nonwet) +
             drho_n_dpn * x_h2o_in_nonwet / (1 - x_h2o_in_nonwet) /
                 molar_mass_co2);*/

        mass_mat_coeff(cap_pressure_coeff_index, cap_pressure_coeff_index) =
            poro * dSwdPc *
            (rho_w / molar_mass_h2o - rho_mol_nonwet * x_h2o_in_nonwet);
		//std::cout << mass_mat_coeff << std::endl;
        _Mgp.noalias() += mass_mat_coeff(nonwet_pressure_coeff_index,
                                         nonwet_pressure_coeff_index) *
                          sm.N.transpose() * sm.N * integration_factor;
        _Mgpc.noalias() += mass_mat_coeff(nonwet_pressure_coeff_index,
                                          cap_pressure_coeff_index) *
                           sm.N.transpose() * sm.N * integration_factor;
        _Mlp.noalias() += mass_mat_coeff(cap_pressure_coeff_index,
                                         nonwet_pressure_coeff_index) *
                          sm.N.transpose() * sm.N * integration_factor;
        _Mlpc.noalias() +=
            mass_mat_coeff(cap_pressure_coeff_index, cap_pressure_coeff_index) *
            sm.N.transpose() * sm.N * integration_factor;

        double const k_rel_L =
            _process_data._material->getrelativePermeability_liquid(Sw);

        double const k_rel_G =
            _process_data._material->getrelativePermeability_gas(Sw);
		
        // gas viscosity
		double const mu_gas = 
            _process_data._material->gasViscosity(pg_int_pt, _temperature, rho_co2);
        double const mu_liquid =
            _process_data._material->getLiquidViscosity(pl, _temperature);

        double const lambda_G = k_rel_G / mu_gas;
        double const lambda_L = k_rel_L / mu_liquid;

        double const diffusion_coeff_componentb =
            _process_data._diffusion_coeff_componentb(t, pos)[0];
        double const diffusion_coeff_componenta =
            _process_data._diffusion_coeff_componenta(t, pos)[0];
        /*
        *construct the K matrix
        */
        K_mat_coeff(nonwet_pressure_coeff_index, nonwet_pressure_coeff_index) =
            (rho_co2 / molar_mass_co2) * perm(0, 0) * lambda_G +
            rho_mol_w * x_co2_in_wet * perm(0, 0) * lambda_L +
			poro * Sw * diffusion_coeff_componentb * rho_mol_w *
			dx_co2_in_wet_dpn -
            poro * (1 - Sw) * diffusion_coeff_componenta *
                (rho_mol_nonwet * dx_h2o_in_nonwet_dpn / (1 - x_h2o_in_nonwet) +
                 drho_n_dpn * x_h2o_in_nonwet / (1 - x_h2o_in_nonwet) /
                     molar_mass_co2);
        K_mat_coeff(nonwet_pressure_coeff_index, cap_pressure_coeff_index) =
            -rho_mol_w * x_co2_in_wet * perm(0, 0) * lambda_L;

        // water
        K_mat_coeff(cap_pressure_coeff_index, nonwet_pressure_coeff_index) =
            rho_w * perm(0, 0) * lambda_L / molar_mass_h2o +
            rho_mol_nonwet * x_h2o_in_nonwet * perm(0, 0) * lambda_G -
            poro * Sw * diffusion_coeff_componentb * rho_mol_w *
                dx_co2_in_wet_dpn  +
			poro * (1 - Sw) * diffusion_coeff_componenta *
			(rho_mol_nonwet * dx_h2o_in_nonwet_dpn / (1 - x_h2o_in_nonwet) +
				drho_n_dpn * x_h2o_in_nonwet / (1 - x_h2o_in_nonwet) /
				molar_mass_co2);// / (1 - x_co2_in_wet) std::pow(poro * Sw,10/3)*(1/poro/poro)

        K_mat_coeff(cap_pressure_coeff_index, cap_pressure_coeff_index) =
            -rho_w * perm(0, 0) * lambda_L / molar_mass_h2o;
        //std::cout << K_mat_coeff << std::endl;
        // assembly the mass matrix

        _Kgp.noalias() += K_mat_coeff(nonwet_pressure_coeff_index,
                                      nonwet_pressure_coeff_index) *
                          sm.dNdx.transpose() * sm.dNdx * integration_factor;
        _Kgpc.noalias() +=
            K_mat_coeff(nonwet_pressure_coeff_index, cap_pressure_coeff_index) *
            sm.dNdx.transpose() * sm.dNdx * integration_factor;
        _Klp.noalias() +=
            K_mat_coeff(cap_pressure_coeff_index, nonwet_pressure_coeff_index) *
            sm.dNdx.transpose() * sm.dNdx * integration_factor;
        _Klpc.noalias() +=
            K_mat_coeff(cap_pressure_coeff_index, cap_pressure_coeff_index) *
            sm.dNdx.transpose() * sm.dNdx * integration_factor;

        // std::cout << local_K << std::endl;
        H_vec_coeff(nonwet_pressure_coeff_index) =
            (rho_co2 / molar_mass_co2) * rho_mass_nonwet * perm(0, 0) *
                lambda_G +
            rho_mol_w * x_co2_in_wet * perm(0, 0) * lambda_L * rho_mass_wet;

        H_vec_coeff(cap_pressure_coeff_index) =
            rho_w * perm(0, 0) * lambda_L * rho_mass_wet / molar_mass_h2o +
            rho_mol_nonwet * x_h2o_in_nonwet * rho_mass_nonwet * perm(0, 0) *
                lambda_G;

        // std::cout << H_vec_coeff << std::endl;
        if (_process_data._has_gravity)
        {
            auto const body_force = _process_data._specific_body_force(t, pos);
            assert(body_force.size() == GlobalDim);
            auto const b =
                MathLib::toVector<GlobalDimVectorType>(body_force, GlobalDim);
            _Bg.noalias() += sm.dNdx.transpose() *
                             H_vec_coeff(nonwet_pressure_coeff_index) * b *
                             integration_factor;
            _Bl.noalias() += sm.dNdx.transpose() *
                             H_vec_coeff(cap_pressure_coeff_index) * b *
                             integration_factor;

        }  // end of has gravity
           // std::cout << local_b << std::endl;
    }      // end of GP
    if (_process_data._has_mass_lumping)
    {
        for (unsigned row = 0; row < _Mgpc.cols(); row++)
        {
            for (unsigned column = 0; column < _Mgpc.cols(); column++)
            {
                if (row != column)
                {
                    _Mgpc(row, row) += _Mgpc(row, column);
                    _Mgpc(row, column) = 0.0;
                    _Mgp(row, row) += _Mgp(row, column);
                    _Mgp(row, column) = 0.0;
                    _Mlpc(row, row) += _Mlpc(row, column);
                    _Mlpc(row, column) = 0.0;
                    _Mlp(row, row) += _Mlp(row, column);
                    _Mlp(row, column) = 0.0;
                }
            }
        }
    }
    // assembler fully coupled mass matrix
    local_M
        .block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index)
        .noalias() += _Mgp;
    local_M
        .block<nonwet_pressure_size, cap_pressure_size>(
            nonwet_pressure_matrix_index, cap_pressure_matrix_index)
        .noalias() += _Mgpc;
    local_M
        .block<cap_pressure_size, nonwet_pressure_size>(
            cap_pressure_matrix_index, nonwet_pressure_matrix_index)
        .noalias() += _Mlp;
    local_M
        .block<cap_pressure_size, cap_pressure_size>(cap_pressure_matrix_index,
                                                     cap_pressure_matrix_index)
        .noalias() += _Mlpc;
    local_K
        .block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index)
        .noalias() += _Kgp;
    local_K
        .block<nonwet_pressure_size, cap_pressure_size>(
            nonwet_pressure_matrix_index, cap_pressure_matrix_index)
        .noalias() += _Kgpc;
    local_K
        .block<cap_pressure_size, nonwet_pressure_size>(
            cap_pressure_matrix_index, nonwet_pressure_matrix_index)
        .noalias() += _Klp;
    local_K
        .block<cap_pressure_size, cap_pressure_size>(cap_pressure_matrix_index,
                                                     cap_pressure_matrix_index)
        .noalias() += _Klpc;

    local_b.block<nonwet_pressure_size, 1>(nonwet_pressure_matrix_index, 0) +=
        _Bg;
    local_b.block<cap_pressure_size, 1>(cap_pressure_matrix_index, 0) += _Bl;
}

}  // end of namespace
}  // end of namespace

#endif
