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
// using namespace Eigen;
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
    //_material_properties.setMaterialID(pos);

    const Eigen::MatrixXd& perm = _process_data._material->getPermeability(
        t, pos, _element.getDimension());

    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to perm.rows() == _element->getDimension()
    assert(perm.rows() == GlobalDim || perm.rows() == 1);

    double porosity_variable = 0.;
    double storage_variable = 0.;
    // Note: currently only isothermal case is considered, so the temperature is
    // assumed to be const
    // the variation of temperatura will be taken into account in future
    _temperature = 293.15;
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
        double const rho_gas =
            _process_data._material->getGasDensity(pg_int_pt, _temperature);
        double const rho_w =
            _process_data._material->getLiquidDensity(pl, _temperature);
        double const mu_gas =
            _process_data._material->getGasViscosity(pg_int_pt, _temperature);
        double const mu_liquid =
            _process_data._material->getLiquidViscosity(pl, _temperature);

        double Sw = _process_data._material->getSaturation(pc_int_pt);  // pc

        if (pc_int_pt < 0)
            Sw = 1.0;
        double dSwdPc = _process_data._material->getDerivSaturation(pc_int_pt);

        double const k_rel_L =
            _process_data._material->getrelativePermeability_liquid(Sw);

        double const k_rel_G =
            _process_data._material->getrelativePermeability_gas(Sw);  // Sw
        double const drhogas_dpg = _process_data._material->getDerivGasDensity(
            pg_int_pt, _temperature);
        double const poro = _process_data._material->getPorosity(
            t, pos, pg_int_pt, _temperature, porosity_variable);

        // Assemble mass matrix, M
        _saturation[ip] = Sw;
        // nonwetting
        mass_mat_coeff(nonwet_pressure_coeff_index,
                       nonwet_pressure_coeff_index) =
            poro * (1 - Sw) * drhogas_dpg;  // dPG
        mass_mat_coeff(nonwet_pressure_coeff_index, cap_pressure_coeff_index) =
            -poro * rho_gas * dSwdPc;  // dPC
        // wetting
        mass_mat_coeff(cap_pressure_coeff_index, nonwet_pressure_coeff_index) =
            0.0;  // dPG
        mass_mat_coeff(cap_pressure_coeff_index, cap_pressure_coeff_index) =
            poro * dSwdPc * rho_w;  // dPC

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

        double const lambda_G = k_rel_G / mu_gas;
        double const lambda_L = k_rel_L / mu_liquid;
        /*
        *construct the K matrix
        */
        K_mat_coeff(nonwet_pressure_coeff_index, nonwet_pressure_coeff_index) =
            rho_gas * perm(0, 0) * lambda_G;
        K_mat_coeff(nonwet_pressure_coeff_index, cap_pressure_coeff_index) =
            0.0;

        // water
        K_mat_coeff(cap_pressure_coeff_index, nonwet_pressure_coeff_index) =
            rho_w * perm(0, 0) * lambda_L;
        K_mat_coeff(cap_pressure_coeff_index, cap_pressure_coeff_index) =
            -rho_w * perm(0, 0) * lambda_L;
        // std::cout << K_mat_coeff << std::endl;
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
            rho_gas * rho_gas * perm(0, 0) * lambda_G;

        H_vec_coeff(cap_pressure_coeff_index) =
            rho_w * rho_w * perm(0, 0) * lambda_L;
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
