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

    NodalMatrixType _Mgp;
    NodalMatrixType _Mgpc;
    NodalMatrixType _Mlp;
    NodalMatrixType _Mlpc;
    NodalMatrixType _Kgp;
    NodalMatrixType _Kgpc;
    NodalMatrixType _Klp;
    NodalMatrixType _Klpc;

    NodalVectorType _Bg;
    NodalVectorType _Bl;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    _process_data._material->setMaterialID(pos);

    const Eigen::MatrixXd& perm = _process_data._material->getPermeability(
        t, pos, _element.getDimension());

    MathLib::PiecewiseLinearInterpolation const& interpolated_Pc =
        *_process_data._curves.at("curve_PC_S");
    MathLib::PiecewiseLinearInterpolation const& interpolated_Kr_wet =
        *_process_data._curves.at("curve_S_Krel_wet");
    MathLib::PiecewiseLinearInterpolation const& interpolated_Kr_nonwet =
        *_process_data._curves.at("curve_S_Krel_nonwet");

    assert(perm.rows() == GlobalDim || perm.rows() == 1);

    // Note: currently only isothermal case is considered, so the temperature is
    // assumed to be const
    // the variation of temperatura will be taken into account in future
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);

        double pc_int_pt = 0.;
        double pg_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pg_int_pt, pc_int_pt);

        _pressure_wetting[ip] = pg_int_pt - pc_int_pt;

        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();
        double const rho_gas =
            _process_data._material->getGasDensity(pg_int_pt, _temperature);
        double const rho_w =
            _process_data._material->getLiquidDensity(_pressure_wetting[ip], _temperature);

        double const Sw =
            (pc_int_pt < 0) ? 1.0 : interpolated_Pc.getValue(pc_int_pt);

        _saturation[ip] = Sw;
        double dSwdPc = interpolated_Pc.getDerivative(pc_int_pt);
        if (pc_int_pt > interpolated_Pc.getSupportMax())
            dSwdPc =
                interpolated_Pc.getDerivative(interpolated_Pc.getSupportMax());
        else if (pc_int_pt < interpolated_Pc.getSupportMin())
            dSwdPc =
                interpolated_Pc.getDerivative(interpolated_Pc.getSupportMin());

        double const porosity = _process_data._material->getPorosity(
            t, pos, pg_int_pt, _temperature, 0);

        // Assemble M matrix
        // nonwetting
        double const drhogas_dpg = _process_data._material->getDerivGasDensity(
            pg_int_pt, _temperature);
        mass_mat_coeff(nonwet_pressure_coeff_index,
                       nonwet_pressure_coeff_index) =
            porosity * (1 - Sw) * drhogas_dpg;
        mass_mat_coeff(nonwet_pressure_coeff_index, cap_pressure_coeff_index) =
            -porosity * rho_gas * dSwdPc;
        // wetting
        mass_mat_coeff(cap_pressure_coeff_index, nonwet_pressure_coeff_index) =
            0.0;
        mass_mat_coeff(cap_pressure_coeff_index, cap_pressure_coeff_index) =
            porosity * dSwdPc * rho_w;

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


        double const k_rel_G = interpolated_Kr_nonwet.getValue(Sw);
        double const mu_gas =
            _process_data._material->getGasViscosity(pg_int_pt, _temperature);
        double const lambda_G = k_rel_G / mu_gas;

        double const k_rel_L = interpolated_Kr_wet.getValue(Sw);
        double const mu_liquid =
            _process_data._material->getLiquidViscosity(_pressure_wetting[ip], _temperature);
        double const lambda_L = k_rel_L / mu_liquid;
        // Assemble M matrix
        // nonwet
        K_mat_coeff(nonwet_pressure_coeff_index, nonwet_pressure_coeff_index) =
            rho_gas * perm(0, 0) * lambda_G;
        K_mat_coeff(nonwet_pressure_coeff_index, cap_pressure_coeff_index) =
            0.0;

        // wet
        K_mat_coeff(cap_pressure_coeff_index, nonwet_pressure_coeff_index) =
            rho_w * perm(0, 0) * lambda_L;
        K_mat_coeff(cap_pressure_coeff_index, cap_pressure_coeff_index) =
            -rho_w * perm(0, 0) * lambda_L;

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

        H_vec_coeff(nonwet_pressure_coeff_index) =
            rho_gas * rho_gas * perm(0, 0) * lambda_G;

        H_vec_coeff(cap_pressure_coeff_index) =
            rho_w * rho_w * perm(0, 0) * lambda_L;

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
