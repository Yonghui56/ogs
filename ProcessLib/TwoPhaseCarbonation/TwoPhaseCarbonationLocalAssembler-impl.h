/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
* common nomenclature
* --------------primary variable----------------------
* pn_int_pt    pressure for nonwetting phase at each integration point
* pc_int_pt    capillary pressure at each integration point
* --------------secondary variable--------------------
* temperature              capillary pressure
* Sw wetting               phase saturation
* dSw_dpc                  derivative of wetting phase saturation with respect
* to capillary pressure
* rho_mol_nonwet           molar density of nonwetting phase
* drhononwet_dpn           derivative of nonwetting phase density with respect
*to nonwetting phase pressure
* rho_wet                  density of wetting phase
* k_rel_nonwet             relative permeability of nonwetting phase
* mu_nonwet                viscosity of nonwetting phase
* lambda_nonwet            mobility of nonwetting phase
* k_rel_wet                relative permeability of wetting phase
* mu_wet                   viscosity of wetting phase
* lambda_wet               mobility of wetting phase
* _pressure_wet            output vector for wetting phase pressure with respect
* to each integration point
* _saturation              output vector for wetting phase saturation with
* respect to each integration point
*/
#pragma once

#include "TwoPhaseCarbonationLocalAssembler.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseCarbonationProcessData.h"

using MaterialLib::PhysicalConstant::IdealGasConstant;
using MaterialLib::PhysicalConstant::MolarMass::Water;
namespace ProcessLib
{
namespace TwoPhaseCarbonation
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseCarbonationLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, std::vector<double> const& local_x,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();

    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    auto Mgp =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mgpc = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);

    auto Mlpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    NodalMatrixType laplace_operator;
    laplace_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kgp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Kgpc = local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);
    auto Klp = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Klpc = local_K.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    auto Bg = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);

    auto Bl =
        local_b.template segment<cap_pressure_size>(cap_pressure_matrix_index);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id =
        _process_data._material->getMaterialID(pos.getElementID().get());

    const Eigen::MatrixXd& perm = _process_data._material->getPermeability(
        material_id, t, pos, _element.getDimension());
    assert(perm.rows() == _element.getDimension() || perm.rows() == 1);
    GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
        _element.getDimension(), _element.getDimension());
    if (perm.rows() == _element.getDimension())
        permeability = perm;
    else if (perm.rows() == 1)
        permeability.diagonal().setConstant(perm(0, 0));

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];

        double pc_int_pt = 0.;
        double pn_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pn_int_pt, pc_int_pt);

        _pressure_wet[ip] = pn_int_pt - pc_int_pt;

        auto const& wp = _integration_method.getWeightedPoint(ip);

        double const dt = _process_data._dt;//get current time step size
        const double temperature = _process_data._temperature(t, pos)[0];
        double const rho_nonwet =
            _process_data._material->getGasDensity(pn_int_pt, temperature);

        double const rho_mol_nonwet = pn_int_pt / IdealGasConstant / temperature;
        double const Henry_constant_co2 = 0.164e+4 * 101325;
        double const x_mol_co2_wet = pn_int_pt / Henry_constant_co2;//dissolved co2
        double const x_mol_co2_nonwet = 1.0;
        double const rho_water = _process_data._material->getLiquidDensity(
            _pressure_wet[ip], temperature);
        double const rho_mol_water = rho_water / Water;
        double const rho_mol_wet = rho_water / Water/(1- x_mol_co2_wet);

        double const rho_mol_total_co2 = rho_mol_nonwet*x_mol_co2_nonwet + rho_mol_wet*x_mol_co2_nonwet;
        double & rho_mol_sio2_wet = _ip_data[ip].rho_mol_sio2;
        double const Sw =
            (pc_int_pt < 0)
                ? 1.0
                : _process_data._material->getSaturation(
                      material_id, t, pos, pn_int_pt, temperature, pc_int_pt);

        _saturation[ip] = Sw;

        double dSw_dpc = _process_data._material->getSaturationDerivative(
            material_id, t, pos, pn_int_pt, temperature, Sw);

        double const porosity = _process_data._material->getPorosity(
            material_id, t, pos, pn_int_pt, temperature, 0);

        // Assemble M matrix
        // nonwetting
        double const drhononwet_dpn =
            _process_data._material->getGasDensityDerivative(pn_int_pt,
                                                             temperature);
        double const d_rho_mol_nonwet_d_pn = 1 / IdealGasConstant / temperature;
        double const d_x_mol_co2_wet_d_pn = 1 / Henry_constant_co2;

        double const d_rho_mol_wet_d_pn = rho_water * d_x_mol_co2_wet_d_pn / Water / (1 - x_mol_co2_wet) / (1 - x_mol_co2_wet);
        Mgp.noalias() +=
            porosity * ((1 - Sw) * x_mol_co2_nonwet *d_rho_mol_nonwet_d_pn +
                Sw*rho_mol_wet*d_x_mol_co2_wet_d_pn + Sw*x_mol_co2_wet*d_rho_mol_wet_d_pn) * _ip_data[ip].massOperator;
        Mgpc.noalias() +=
            -porosity * (rho_mol_wet*x_mol_co2_wet - rho_mol_nonwet*x_mol_co2_nonwet) * dSw_dpc * _ip_data[ip].massOperator;

        Mlpc.noalias() +=
            porosity * dSw_dpc * rho_mol_water * _ip_data[ip].massOperator;

        double const diffusion_coeff_component_co2 =
            _process_data._diffusion_coeff_component_b(t, pos)[0];
        // nonwet
        double const k_rel_nonwet =
            _process_data._material->getNonwetRelativePermeability(
                t, pos, pn_int_pt, temperature, Sw);
        double const mu_nonwet =
            _process_data._material->getGasViscosity(pn_int_pt, temperature);
        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;

        // wet
        double const k_rel_wet =
            _process_data._material->getWetRelativePermeability(
                t, pos, _pressure_wet[ip], temperature, Sw);
        double const mu_wet = _process_data._material->getLiquidViscosity(
            _pressure_wet[ip], temperature);
        double const lambda_wet = k_rel_wet / mu_wet;

        laplace_operator.noalias() = sm.dNdx.transpose() * permeability *
                                     sm.dNdx * _ip_data[ip].integration_weight;

        Kgp.noalias() += rho_mol_nonwet *x_mol_co2_nonwet* lambda_nonwet * laplace_operator
            + rho_mol_wet*x_mol_co2_wet*lambda_wet *laplace_operator
            +Sw*porosity*rho_mol_wet*diffusion_coeff_component_co2*d_x_mol_co2_wet_d_pn*_ip_data[ip].diffusionOperator;;

        Kgpc.noalias() += -rho_mol_nonwet*x_mol_co2_nonwet*lambda_wet*laplace_operator;

        Klp.noalias() += rho_mol_water*lambda_wet*laplace_operator 
            - Sw*porosity*rho_mol_wet*diffusion_coeff_component_co2*d_x_mol_co2_wet_d_pn*_ip_data[ip].diffusionOperator;
        Klpc.noalias() += -rho_mol_water * lambda_wet * laplace_operator;

        if (_process_data._has_gravity)
        {
            auto const& b = _process_data._specific_body_force;
            Bg.noalias() += (rho_mol_nonwet *x_mol_co2_nonwet * lambda_nonwet *  rho_nonwet 
                + rho_mol_wet*x_mol_co2_wet*lambda_wet*rho_water) *
                            sm.dNdx.transpose() * permeability * b *
                            _ip_data[ip].integration_weight;
            Bl.noalias() += rho_mol_water * rho_water * lambda_wet *
                            sm.dNdx.transpose() * permeability * b *
                            _ip_data[ip].integration_weight;
        }  // end of has gravity
        double const fluid_volume_rate(0.0);
        double const quartz_dissolute_rate(0.0);
        if (Sw > 0.3)
        {
            Bg.noalias() += 0.0;
            Bl.noalias() += sm.dNdx.transpose() * fluid_volume_rate *_ip_data[ip].integration_weight;
            rho_mol_sio2_wet = rho_mol_sio2_wet + quartz_dissolute_rate*dt;//cumulative dissolved 
        }
    }
    if (_process_data._has_mass_lumping)
    {
        for (unsigned row = 0; row < Mgpc.cols(); row++)
        {
            for (unsigned column = 0; column < Mgpc.cols(); column++)
            {
                if (row != column)
                {
                    Mgpc(row, row) += Mgpc(row, column);
                    Mgpc(row, column) = 0.0;
                    Mgp(row, row) += Mgp(row, column);
                    Mgp(row, column) = 0.0;
                    Mlpc(row, row) += Mlpc(row, column);
                    Mlpc(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // end of namespace
}  // end of namespace
