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
* rho_nonwet               density of nonwetting phase
* drhononwet_dpn           derivative of nonwetting phase density with respect
*to nonwetting phase pressure
* rho_wet                  density of wetting phase
* k_rel_nonwet             relative permeability of nonwetting phase
* mu_nonwet                viscosity of nonwetting phase
* lambda_nonwet            mobility of nonwetting phase
* k_rel_wet                relative permeability of wetting phase
* mu_wet                   viscosity of wetting phase
* lambda_wet               mobility of wetting phase
*/
#pragma once

#include "TwoPhaseFlowWithPPLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseFlowWithPPProcessData.h"

using MaterialLib::PhysicalConstant::MolarMass::Water;
using MaterialLib::PhysicalConstant::MolarMass::Air;
using MaterialLib::PhysicalConstant::IdealGasConstant;
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

    auto Mlp = local_M.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mlpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    NodalMatrixType laplace_operator;
    laplace_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kgp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    
    auto Kgpc = local_K.template block<nonwet_pressure_size, cap_pressure_size>(
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
        _process_data.material->getMaterialID(pos.getElementID().get());

    const Eigen::MatrixXd& perm = _process_data.material->getPermeability(
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

        const double temperature = _process_data.temperature(t, pos)[0];
        double const rho_nonwet =
            _process_data.material->getGasDensity(pn_int_pt, temperature);
        double const rho_wet = _process_data.material->getLiquidDensity(
            _pressure_wet[ip], temperature);

        double const Sw = (pc_int_pt<0) ? 1 : _process_data.material->getSaturation(
            material_id, t, pos, pn_int_pt, temperature, pc_int_pt);

        _saturation[ip] = Sw;

        double dSw_dpc = (pc_int_pt < 0) ? 0 : _process_data.material->getSaturationDerivative(
            material_id, t, pos, pn_int_pt, temperature, Sw);

        double const porosity = _process_data.material->getPorosity(
            material_id, t, pos, pn_int_pt, temperature, 0);

        double p_sat = get_P_sat(temperature);

        const double rho_mol_water = _rho_l_std / Water;
        double kelvin_term = exp(pc_int_pt / rho_mol_water / IdealGasConstant / temperature);
        double d_kelvin_term_d_pc = kelvin_term / rho_mol_water / IdealGasConstant / temperature;
        if (pc_int_pt < 0)
        {
            kelvin_term= exp(0 / rho_mol_water / IdealGasConstant / temperature);
            d_kelvin_term_d_pc = 0;
        }
        double const x_vapor_nonwet = get_x_nonwet_vapor(pn_int_pt, p_sat, kelvin_term);
        double const x_water_wet = pn_int_pt*x_vapor_nonwet*kelvin_term / p_sat;
        double const x_air_nonwet = 1- x_vapor_nonwet;
        double const x_air_wet = pn_int_pt*x_air_nonwet / _hen_L_air;

        double const d_x_vapor_nonwet_d_pg = get_derivative_x_nonwet_h2o_d_pg(pn_int_pt, p_sat, kelvin_term);
        double const d_x_vapor_nonwet_d_kelvin = get_derivative_x_nonwet_h2o_d_kelvin(pn_int_pt, p_sat, kelvin_term);
        double const d_x_vapor_nonwet_d_pc = d_x_vapor_nonwet_d_kelvin * d_kelvin_term_d_pc;

        double const d_x_air_nonwet_d_pg = -d_x_vapor_nonwet_d_pg;
        double const d_x_air_nonwet_d_pc = -d_x_vapor_nonwet_d_pc;
        double const d_x_air_wet_d_pg = x_air_nonwet / _hen_L_air + pn_int_pt * d_x_air_nonwet_d_pg / _hen_L_air;
        double const d_x_air_wet_d_pc= pn_int_pt * d_x_air_nonwet_d_pc / _hen_L_air;

        const double rho_mol_nonwet = pn_int_pt / IdealGasConstant / temperature;
        const double rho_mol_wet = rho_mol_water / x_water_wet;
        double const d_rho_mol_nonwet_d_pg = 1 / IdealGasConstant / temperature;
        double const d_rho_mol_wet_d_pg= -rho_mol_water * kelvin_term *(x_vapor_nonwet / p_sat +
            pn_int_pt * d_x_vapor_nonwet_d_pg / p_sat) /
            x_water_wet / x_water_wet;
        double const d_rho_mol_wet_d_pc = -rho_mol_water *(pn_int_pt * d_x_vapor_nonwet_d_pc * kelvin_term / p_sat
            + pn_int_pt * x_vapor_nonwet * d_kelvin_term_d_pc / p_sat) /
            x_water_wet / x_water_wet;
        double const rho_mass_nonwet = rho_mol_nonwet*(x_vapor_nonwet*Water + x_air_nonwet*Air);
        double const rho_mass_wet = rho_mol_wet *
            (x_water_wet* Water + x_air_wet*Air);

        // Assemble M matrix
        // nonwetting

        Mgp.noalias() +=
            porosity * (d_rho_mol_nonwet_d_pg*(1 - Sw)*x_air_nonwet + rho_mol_nonwet*(1 - Sw)*d_x_air_nonwet_d_pg
                + rho_mol_wet * Sw* d_x_air_nonwet_d_pg*pn_int_pt / _hen_L_air
                + rho_mol_wet*Sw*x_air_wet / _hen_L_air)
            * _ip_data[ip].massOperator;//dpn
        Mgpc.noalias() +=
            porosity * ((rho_mol_wet*x_air_nonwet-rho_mol_nonwet*x_air_nonwet) * dSw_dpc +
                Sw*d_rho_mol_wet_d_pc*x_air_wet
                +Sw*rho_mol_wet*d_x_air_wet_d_pc)* _ip_data[ip].massOperator;//dpc

        Mlp.noalias() += porosity *(d_rho_mol_nonwet_d_pg*(1 - Sw)
            + d_rho_mol_wet_d_pg *Sw)*_ip_data[ip].massOperator;//dpn
        Mlpc.noalias() +=
            porosity * (dSw_dpc * (rho_mol_wet-rho_mol_nonwet) +d_rho_mol_wet_d_pc*Sw)
            * _ip_data[ip].massOperator;//dpc
        /*std::cout << Mgp << std::endl;
        std::cout << Mgpc << std::endl;
        std::cout << Mlp << std::endl;
        std::cout << Mlpc << std::endl;*/
        // nonwet
        double const k_rel_nonwet =
            _process_data.material->getNonwetRelativePermeability(
                t, pos, pn_int_pt, temperature, Sw);
        double const mu_nonwet =
            _process_data.material->getGasViscosity(pn_int_pt, temperature);
        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;
        double const diffusion_coeff_component_air =
            _process_data.diffusion_coeff_component_b(t, pos)[0];
        // wet
        double const k_rel_wet =
            _process_data.material->getWetRelativePermeability(
                t, pos, _pressure_wet[ip], temperature, Sw);
        double const mu_wet = _process_data.material->getLiquidViscosity(
            _pressure_wet[ip], temperature);
        double const lambda_wet = k_rel_wet / mu_wet;
        double const diffusion_coeff_component_h2o =
            _process_data.diffusion_coeff_component_a(t, pos)[0];

        laplace_operator.noalias() = sm.dNdx.transpose() * permeability *
                                     sm.dNdx * _ip_data[ip].integration_weight;

        Kgp.noalias() += (lambda_nonwet*rho_mol_nonwet*x_air_nonwet + lambda_wet*rho_mol_wet*x_air_wet)
            *laplace_operator
            + (porosity*(1 - Sw)*rho_mol_nonwet*diffusion_coeff_component_air*d_x_air_nonwet_d_pg
                + porosity*Sw*rho_mol_wet*diffusion_coeff_component_h2o*(pn_int_pt * d_x_air_nonwet_d_pg / _hen_L_air + x_air_nonwet / _hen_L_air))
            *_ip_data[ip].diffusion_operator;

        Kgpc.noalias() +=-lambda_wet*rho_mol_wet*x_air_wet*laplace_operator
            + (porosity*(1 - Sw)*rho_mol_nonwet*diffusion_coeff_component_air*d_x_air_nonwet_d_pc
                + porosity*Sw*rho_mol_wet*diffusion_coeff_component_h2o*(pn_int_pt * d_x_air_nonwet_d_pc / _hen_L_air))
            *_ip_data[ip].diffusion_operator;
        Klp.noalias() += (lambda_nonwet*rho_mol_nonwet + lambda_wet*rho_mol_wet)*laplace_operator;
        Klpc.noalias() += -lambda_wet * rho_mol_wet* laplace_operator;
        /*std::cout << Kgp << std::endl;
        std::cout << Kgpc << std::endl;
        std::cout << Klp << std::endl;
        std::cout << Klpc << std::endl;*/

        if (_process_data.has_gravity)
        {
            auto const& b = _process_data.specific_body_force;

            NodalVectorType gravity_operator = sm.dNdx.transpose() *
                                               permeability * b *
                                               _ip_data[ip].integration_weight;
            Bg.noalias() +=
                rho_mol_nonwet*x_air_nonwet * rho_mass_nonwet * lambda_nonwet * gravity_operator
                + rho_mol_wet*x_air_wet*rho_mass_wet*lambda_wet*gravity_operator;
            Bl.noalias() += rho_mol_wet * rho_mass_wet * lambda_wet * gravity_operator
                + rho_mol_nonwet*rho_mass_nonwet*lambda_nonwet*gravity_operator;
        }  // end of has gravity
    }
    if (_process_data.has_mass_lumping)
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
