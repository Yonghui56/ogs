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

#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseNCompCarbonationLocalAssembler.h"
#include "TwoPhaseNCompCarbonationProcessData.h"

using MaterialLib::PhysicalConstant::IdealGasConstant;
using MaterialLib::PhysicalConstant::MolarMass::Water;
using MaterialLib::PhysicalConstant::MolarMass::Air;

namespace ProcessLib
{
namespace TwoPhaseNCompCarbonation
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseNCompCarbonationLocalAssembler<
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

    auto Mcp =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mcpc = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);
    auto Mcx =
        local_M.template block<nonwet_pressure_size, molar_fraction_size>(
            nonwet_pressure_matrix_index, molar_fraction_matrix_index);

    auto Map = local_M.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mapc = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);
    auto Max = local_M.template block<cap_pressure_size, molar_fraction_size>(
        cap_pressure_matrix_index, molar_fraction_matrix_index);

    auto Mlpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
        molar_fraction_matrix_index, cap_pressure_matrix_index);
    auto Mlx = local_M.template block<cap_pressure_size, molar_fraction_size>(
        molar_fraction_matrix_index, molar_fraction_matrix_index);

    NodalMatrixType laplace_operator;
    laplace_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kcp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Kcpc =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, cap_pressure_matrix_index);
    auto Kcx =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, molar_fraction_matrix_index);

    auto Kap =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            cap_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Kapc =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            cap_pressure_matrix_index, cap_pressure_matrix_index);
    auto Kax =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            cap_pressure_matrix_index, molar_fraction_matrix_index);

    auto Klp = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
        molar_fraction_matrix_index, nonwet_pressure_matrix_index);

    auto Klpc = local_K.template block<cap_pressure_size, cap_pressure_size>(
        molar_fraction_matrix_index, cap_pressure_matrix_index);
    auto Klx = local_K.template block<cap_pressure_size, cap_pressure_size>(
        molar_fraction_matrix_index, molar_fraction_matrix_index);

    auto Bc = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);
    auto Ba = local_b.template segment<nonwet_pressure_size>(
        cap_pressure_matrix_index);
    auto Bl = local_b.template segment<cap_pressure_size>(
        molar_fraction_matrix_index);

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
        double x_co2_gas_int_pt = 0.;  // co2 molar fraction in the gas phase
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pn_int_pt, pc_int_pt,
                                         x_co2_gas_int_pt);

        _pressure_wet[ip] = pn_int_pt - pc_int_pt;

        auto const& wp = _integration_method.getWeightedPoint(ip);

        double const dt = _process_data._dt;  // get current time step size
        const double temperature = _process_data._temperature(t, pos)[0];

        double const Henry_constant_co2 = 0.164e+4 * 101325;
        double const Henry_constant_air = 0.164e+4 * 101325;
        double const p_co2_nonwet = pn_int_pt * x_co2_gas_int_pt;
        double const p_air_nonwet = pn_int_pt * (1 - x_co2_gas_int_pt);

        double const x_mol_air_nonwet = p_air_nonwet / pn_int_pt;
        double const x_mol_air_wet = p_air_nonwet / Henry_constant_air;
        double const d_x_mol_air_nonwet_d_pn = 0.0;
        double const d_x_mol_air_nonwet_d_x = -1;
        double const d_x_mol_air_wet_d_pn =
            x_mol_air_nonwet / Henry_constant_air;
        double const d_x_mol_air_wet_d_x = -pn_int_pt / Henry_constant_air;

        double const x_mol_co2_wet =
            p_co2_nonwet / Henry_constant_co2;  // dissolved co2
        double const d_x_mol_co2_wet_d_pn =
            x_co2_gas_int_pt / Henry_constant_co2;
        double const d_x_mol_co2_wet_d_x = pn_int_pt / Henry_constant_co2;

        double const rho_water = _process_data._material->getLiquidDensity(
            _pressure_wet[ip], temperature);
        double const rho_mol_water = rho_water / Water;
        double const rho_mol_wet =
            rho_water / Water / (1 - x_mol_co2_wet - x_mol_air_wet);

        double const rho_mol_nonwet =
            pn_int_pt / IdealGasConstant / temperature;
        double const rho_mass_nonwet =
            rho_mol_nonwet *
            (x_mol_air_nonwet * Air + x_co2_gas_int_pt * 0.044);

        double const Sw =
            (pc_int_pt < 0)
                ? 1.0
                : _process_data._material->getSaturation(
                      material_id, t, pos, pn_int_pt, temperature, pc_int_pt);

        _saturation[ip] = Sw;

        double dSw_dpc = _process_data._material->getSaturationDerivative(
            material_id, t, pos, pn_int_pt, temperature, Sw);

        double& rho_mol_sio2_wet = _ip_data[ip].rho_mol_sio2;
        /*double const porosity = _process_data._material->getPorosity(
            material_id, t, pos, pn_int_pt, temperature, 0);*/
        double& porosity = _ip_data[ip].porosity;
        // get cumulative consumed co2
        double& rho_mol_co2_cumul_total =
            _ip_data[ip].rho_mol_co2_cumul_total;  // get cumulative co2

        // calculate the current ammount of co2
        double rho_mol_total_co2 =
            _ip_data[ip].porosity_prev *
            (rho_mol_nonwet * x_co2_gas_int_pt * (1 - Sw) +
             rho_mol_wet * x_mol_co2_wet * Sw);
        double test =
            _ip_data[ip].porosity_prev *
            (rho_mol_nonwet * 1 * 1e-6 + rho_mol_wet * x_mol_co2_wet * 1);

        double const pH = bi_interpolation(
            _ip_data[ip].rho_mol_sio2_prev,
            _ip_data[ip].rho_mol_co2_cumul_total_prev, _pH_at_supp_pnt);
        _pH_value[ip] = pH;
        // Assemble M matrix
        // nonwetting
        _co2_concentration[ip] = rho_mol_nonwet*x_co2_gas_int_pt;
        double const d_rho_mol_nonwet_d_pn = 1 / IdealGasConstant / temperature;
        double const d_rho_mol_nonwet_d_x = 0.0;
        double const x_mol_water_wet = 1 - x_mol_air_wet - x_mol_co2_wet;
        double const d_rho_mol_wet_d_pn =
            rho_water * (d_x_mol_co2_wet_d_pn + d_x_mol_air_wet_d_pn) / Water /
            x_mol_water_wet / x_mol_water_wet;
        double const d_rho_mol_wet_d_x =
            rho_water * (d_x_mol_co2_wet_d_x + d_x_mol_air_wet_d_x) / Water /
            (1 - x_mol_co2_wet - x_mol_air_wet) /
            (1 - x_mol_co2_wet - x_mol_air_wet);

        Mcp.noalias() += porosity *
                         ((1 - Sw) * x_co2_gas_int_pt * d_rho_mol_nonwet_d_pn +
                          Sw * rho_mol_wet * d_x_mol_co2_wet_d_pn +
                          Sw * x_mol_co2_wet * d_rho_mol_wet_d_pn) *
                         _ip_data[ip].massOperator;
        Mcpc.noalias() += porosity * (rho_mol_wet * x_mol_co2_wet -
                                      rho_mol_nonwet * x_co2_gas_int_pt) *
                          dSw_dpc * _ip_data[ip].massOperator;
        Mcx.noalias() += porosity *
                         ((1 - Sw) * x_co2_gas_int_pt * d_rho_mol_nonwet_d_x +
                          (1 - Sw) * rho_mol_nonwet +
                          Sw * rho_mol_wet * d_x_mol_co2_wet_d_x +
                          Sw * x_mol_co2_wet * d_rho_mol_wet_d_x) *
                         _ip_data[ip].massOperator;

        Map.noalias() += porosity *
                         ((1 - Sw) * x_mol_air_nonwet * d_rho_mol_nonwet_d_pn +
                          (1 - Sw) * d_x_mol_air_nonwet_d_pn * rho_mol_nonwet +
                          Sw * x_mol_air_wet * d_rho_mol_wet_d_pn +
                          Sw * rho_mol_wet * d_x_mol_air_wet_d_pn) *
                         _ip_data[ip].massOperator;
        Mapc.noalias() += porosity * (rho_mol_wet * x_mol_air_wet -
                                      rho_mol_nonwet * x_mol_air_nonwet) *
                          dSw_dpc * _ip_data[ip].massOperator;
        Max.noalias() += porosity *
                         ((1 - Sw) * x_mol_air_nonwet * d_rho_mol_nonwet_d_x +
                          (1 - Sw) * d_x_mol_air_nonwet_d_x * rho_mol_nonwet +
                          Sw * d_x_mol_air_wet_d_x * rho_mol_wet +
                          Sw * x_mol_air_wet * d_rho_mol_wet_d_x) *
                         _ip_data[ip].massOperator;

        Mlpc.noalias() +=
            porosity * dSw_dpc * rho_mol_water * _ip_data[ip].massOperator;

        double const diffusion_coeff_wet =
            _process_data._diffusion_coeff_component_b(t, pos)[0];
        double const diffusion_coeff_nonwet =
            _process_data._diffusion_coeff_component_a(t, pos)[0];

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
        auto num_nodes = ShapeFunction::NPOINTS;
        auto xn_nodal_values = 
            Eigen::Map<const NodalVectorType>(&local_x[2 * num_nodes], num_nodes);
        GlobalDimVectorType velocity_diff_nonwet =
            -diffusion_coeff_nonwet*porosity*(1 - Sw)*sm.dNdx*xn_nodal_values;

        Kcp.noalias() +=
            rho_mol_nonwet * x_co2_gas_int_pt * lambda_nonwet *
                laplace_operator +
            rho_mol_wet * x_mol_co2_wet * lambda_wet * laplace_operator +
            Sw * porosity * rho_mol_wet * diffusion_coeff_wet *
                d_x_mol_co2_wet_d_pn * _ip_data[ip].diffusionOperator;

        Kcpc.noalias() +=
            -rho_mol_wet * x_mol_co2_wet * lambda_wet * laplace_operator;
        Kcx.noalias() +=
            Sw * porosity * rho_mol_wet * diffusion_coeff_wet *
                d_x_mol_co2_wet_d_x * _ip_data[ip].diffusionOperator +
            (1 - Sw) * porosity * rho_mol_nonwet * diffusion_coeff_nonwet *
                _ip_data[ip].diffusionOperator;

        Kap.noalias() +=
            rho_mol_nonwet * x_mol_air_nonwet * lambda_nonwet *
                laplace_operator +
            rho_mol_wet * x_mol_air_wet * lambda_wet * laplace_operator +
            Sw * porosity * rho_mol_wet * diffusion_coeff_wet *
                d_x_mol_air_wet_d_pn * _ip_data[ip].diffusionOperator;

        Kapc.noalias() +=
            -rho_mol_wet * x_mol_air_wet * lambda_wet * laplace_operator;
        Kax.noalias() +=
            Sw * porosity * rho_mol_wet * diffusion_coeff_wet *
                d_x_mol_air_wet_d_x * _ip_data[ip].diffusionOperator -
            (1 - Sw) * porosity * rho_mol_nonwet * diffusion_coeff_nonwet *
                _ip_data[ip].diffusionOperator;

        Klp.noalias() +=
            rho_mol_water * lambda_wet * laplace_operator -
            Sw * porosity * rho_mol_wet * diffusion_coeff_wet *
                d_x_mol_co2_wet_d_pn * _ip_data[ip].diffusionOperator -
            Sw * porosity * rho_mol_wet * diffusion_coeff_wet *
                d_x_mol_air_wet_d_pn * _ip_data[ip].diffusionOperator;
        Klpc.noalias() += -rho_mol_water * lambda_wet * laplace_operator;
        Klx.noalias() +=
            -Sw * porosity * rho_mol_wet * diffusion_coeff_wet *
                d_x_mol_co2_wet_d_x * _ip_data[ip].diffusionOperator -
            Sw * porosity * rho_mol_wet * diffusion_coeff_wet *
                d_x_mol_air_wet_d_x * _ip_data[ip].diffusionOperator;
        if (_process_data._has_gravity)
        {
            auto const& b = _process_data._specific_body_force;
            Bc.noalias() +=
                (rho_mol_nonwet * x_co2_gas_int_pt * lambda_nonwet *
                     rho_mass_nonwet +
                 rho_mol_wet * x_mol_co2_wet * lambda_wet * rho_water) *
                sm.dNdx.transpose() * permeability * b *
                _ip_data[ip].integration_weight;
            Ba.noalias() +=
                (rho_mol_nonwet * x_mol_air_nonwet * lambda_nonwet *
                     rho_mass_nonwet +
                 rho_mol_wet * x_mol_air_wet * lambda_wet * rho_water) *
                sm.dNdx.transpose() * permeability * b *
                _ip_data[ip].integration_weight;
            Bl.noalias() += rho_mol_water * rho_water * lambda_wet *
                            sm.dNdx.transpose() * permeability * b *
                            _ip_data[ip].integration_weight;
        }  // end of has gravity

        double const flag_carbon = bi_interpolation(
            _ip_data[ip].rho_mol_sio2_prev,
            _ip_data[ip].rho_mol_co2_cumul_total_prev, _flag_carbon_suppt_pnt);
        double& fluid_volume = _ip_data[ip].fluid_volume;
        fluid_volume = bi_interpolation(
            _ip_data[ip].rho_mol_sio2_prev,
            _ip_data[ip].rho_mol_co2_cumul_total_prev, _fluid_volume_suppt_pnt);
        double quartz_dissolute_rate = 31557600*bi_interpolation(
            _ip_data[ip].rho_mol_sio2_prev,
            _ip_data[ip].rho_mol_co2_cumul_total_prev, _quartz_rate_suppt_pnt);
        // quartz_dissolute_rate is always nonpositive.
        if (quartz_dissolute_rate > 0)
            quartz_dissolute_rate = 0;

        if (Sw > 0.3 && dt > 0)  // threshhold value
        {
            double const fluid_volume_rate =
                (fluid_volume - _ip_data[ip].fluid_volume_prev) / dt;
            if (_ip_data[ip].rho_mol_co2_cumul_total_prev >=
                3800)  // means carbonation stops
                rho_mol_total_co2 = 0.0;
            // update the current cumulated co2 consumption
            rho_mol_co2_cumul_total =
                _ip_data[ip].rho_mol_co2_cumul_total_prev + rho_mol_total_co2;
            // co2 consumption
            Bc.noalias() -= sm.N.transpose() * (rho_mol_total_co2 / dt) *
                            _ip_data[ip].integration_weight;
            // air compensate
            Ba.noalias() += sm.N.transpose() * (rho_mol_total_co2 / dt) *
                _ip_data[ip].integration_weight;
            // water source/sink term
            Bl.noalias() += sm.N.transpose() * fluid_volume_rate *
                            _ip_data[ip].integration_weight;
            // update the amount of dissolved sio2
            rho_mol_sio2_wet =
                _ip_data[ip].rho_mol_sio2_prev -
                quartz_dissolute_rate * dt;  // cumulative dissolved sio2
            // update the porosity
            porosity =
                bi_interpolation(_ip_data[ip].rho_mol_sio2_prev,
                                 _ip_data[ip].rho_mol_co2_cumul_total_prev,
                                 _porosity_at_supp_pnts);
        }
        _porosity_value[ip] = porosity;
        _accum_co2_concentr[ip] = rho_mol_co2_cumul_total;
        _accum_diss_sio2_concentr[ip] = rho_mol_sio2_wet;
        for (unsigned d = 0; d < GlobalDim; ++d)
        {
           _darcy_velocities[d][ip] = velocity_diff_nonwet[d];
        }
    }
    if (_process_data._has_mass_lumping)
    {
        for (unsigned row = 0; row < Mcpc.cols(); row++)
        {
            for (unsigned column = 0; column < Mcpc.cols(); column++)
            {
                if (row != column)
                {
                    Mcpc(row, row) += Mcpc(row, column);
                    Mcpc(row, column) = 0.0;
                    Mcp(row, row) += Mcp(row, column);
                    Mcp(row, column) = 0.0;
                    Mcx(row, row) += Mcx(row, column);
                    Mcx(row, column) = 0.0;
                    Mapc(row, row) += Mapc(row, column);
                    Mapc(row, column) = 0.0;
                    Map(row, row) += Map(row, column);
                    Map(row, column) = 0.0;
                    Max(row, row) += Max(row, column);
                    Max(row, column) = 0.0;
                    Mlpc(row, row) += Mlpc(row, column);
                    Mlpc(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // end of namespace
}  // end of namespace
