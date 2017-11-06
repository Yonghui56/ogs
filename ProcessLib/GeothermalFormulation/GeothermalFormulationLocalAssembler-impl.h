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

#include "GeothermalFormulationLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "GeothermalFormulationProcessData.h"

namespace ProcessLib
{
namespace GeothermalFormulation
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void GeothermalFormulationLocalAssembler<
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

    auto const num_nodes = ShapeFunction::NPOINTS;

    auto Mgp =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            0, 0);
    auto Mgh = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        0, num_nodes);

    auto Mtp = local_M.template block<cap_pressure_size, cap_pressure_size>(
        num_nodes, 0);
    auto Mth = local_M.template block<cap_pressure_size, cap_pressure_size>(
        num_nodes, num_nodes);

    NodalMatrixType laplace_operator =
        NodalMatrixType::Zero(num_nodes, num_nodes);

    auto Kgp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            0, 0);
    auto Kgh =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            0, num_nodes);


    auto Ktp = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
        num_nodes, 0);

    auto Kth = local_K.template block<cap_pressure_size, cap_pressure_size>(
        num_nodes, num_nodes);

    auto Bg = local_b.template segment<nonwet_pressure_size>(0);

    auto Bt =
        local_b.template segment<cap_pressure_size>(num_nodes);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto p_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[0], num_nodes);
    auto const& b = _process_data.specific_body_force;
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
        double p_int_pt = 0.;//store the primary variable for pressure
        double h_int_pt = 0.; //store the primary variable for overall enthalpy
        NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, p_int_pt,
                                         h_int_pt);

        //calculate the secondary variable
        double& Sw = _ip_data[ip].sw;
        /// Here only consider one component in gas phase
        double const X_h2_nonwet = 1.0;
        double& hs = _ip_data[ip].h_steam;
        double& hw= _ip_data[ip].h_water;
        double& dSwdp = _ip_data[ip].dsw_dp;
        double& dSwdh = _ip_data[ip].dsw_dh;
        double& dhs_dp = _ip_data[ip].d_h_steam_d_p;
        double& dhs_dh = _ip_data[ip].d_h_steam_d_h;
        double& dhw_dp = _ip_data[ip].d_h_water_d_p;
        double& dhw_dh = _ip_data[ip].d_h_water_d_h;
        if (!_ip_data[ip].mat_property.computeConstitutiveRelation(
            t,
            pos,
            material_id,
            p_int_pt,
            h_int_pt,
            Sw,
            hs,
            hw,
            dSwdp,
            dSwdh,
            dhs_dp,
            dhs_dh,
            dhw_dp,
            dhw_dh))
            OGS_FATAL("Computation of local constitutive relation failed.");
        _pressure_wet[ip] = p_int_pt;
        auto const density_solid = _process_data.density_solid(t, pos)[0];

        //calculate the saturation pressure
        double p_sat = 0.0;
        // Determine which phases are present based on the value of z
        bool is_liquid = false;
        bool is_gas = false;
        bool is_twophase = false;
        double dTdh = 0;
        double dTdp = 0;
        double rho_water = 0;
        double rho_steam = 0;
        double d_rho_w_d_p = 0;
        double d_rho_w_d_h = 0;
        double d_rho_s_d_p = 0;
        double d_rho_s_d_h = 0;
        double temperature = 273.15;
        double const rho_solid = 2500;
        double const specific_heat_capacity_solid = 0.0;
        if (Sw > 1 - 1e-6) {//water region
            is_liquid = true;
        }
        else if (Sw < 1e-6) {//steam region
            is_gas = true;
        }
        else{ // two - phase region
            is_twophase = true;
        }

        if (is_twophase) {
            temperature = _process_data.material->getTemperatureWaterRegion(
                p_int_pt, hw);
            p_sat = get_P_sat(temperature);
            dTdp = -2 * 9.31415e-15*p_int_pt;

            dTdh = (2.56222e-4 - 2 * 2.2568e-11*hw)*dhw_dh;
            rho_water =
                _process_data.material->getDensityWater(p_sat, h_int_pt);
            rho_steam =
                _process_data.material->getDensitySteam(p_sat, h_int_pt);
            d_rho_w_d_p = _process_data.material->getDerivWaterDensity_dPressure(p_sat, h_int_pt);
            d_rho_w_d_h = _process_data.material->getDerivWaterDensity_dEnthalpy(p_sat, h_int_pt);
            d_rho_s_d_p = _process_data.material->getDerivSteamDensity_dPressure(p_sat, h_int_pt);
            d_rho_s_d_h== _process_data.material->getDerivSteamDensity_dEnthalpy(p_sat, h_int_pt);
        }
        else if (is_liquid)
        {
            temperature = _process_data.material->getTemperatureWaterRegion(
                p_int_pt, h_int_pt);
            p_sat = get_P_sat(temperature);
            dTdh = 2.56222e-4 - 2 * 2.2568e-11*h_int_pt;
            dTdp = -2 * 9.31415e-15*p_int_pt;
            rho_water= _process_data.material->getDensityWater(p_int_pt, h_int_pt);
            d_rho_w_d_p = _process_data.material->getDerivWaterDensity_dPressure(p_sat, h_int_pt);
            d_rho_w_d_h = _process_data.material->getDerivWaterDensity_dEnthalpy(p_sat, h_int_pt);
        }
        else if (is_gas)
        {
            temperature= _process_data.material->getTemperatureSteamRegion(
                p_int_pt, h_int_pt);
            p_sat = get_P_sat(temperature);
            dTdp = 4.79921e-5 - 2 * 6.33606e-13*p_int_pt
                + 2 * 3.3372e+24 * std::pow(h_int_pt, -2) * std::pow(p_int_pt, -3)
                - 3 * 3.57154e+16* std::pow(p_int_pt, -4)
                - 1.1725e-24 * std::pow(h_int_pt, 3);
            dTdh = 2 * 7.39386e-11 * h_int_pt
                + 2 * 3.3372e+24 * std::pow(h_int_pt, -3) * std::pow(p_int_pt, -2)
                - 3 * 1.1725e-24*std::pow(h_int_pt, 2) + 4 * 2.26861e+27 * std::pow(h_int_pt, -5);
            rho_steam =
                _process_data.material->getDensitySteam(p_int_pt, h_int_pt);
            d_rho_s_d_p=_process_data.material->getDerivSteamDensity_dPressure(p_int_pt, h_int_pt);
            d_rho_s_d_h = _process_data.material->getDerivSteamDensity_dEnthalpy(p_int_pt, h_int_pt);
        }

        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

        _saturation[ip] = Sw;

        double const porosity = _process_data.material->getPorosity(
            material_id, t, pos, p_int_pt, temperature, 0);

        // Assemble M matrix
        // nonwetting
        
        Mgp.noalias() +=
            porosity * (Sw*d_rho_w_d_p+(1-Sw)*d_rho_s_d_p+ dSwdp*(rho_water-rho_steam))* _ip_data[ip].massOperator;
        Mgh.noalias() +=
            porosity * (Sw*d_rho_w_d_h + (1 - Sw)* d_rho_s_d_h + dSwdh*(rho_water - rho_steam)) * _ip_data[ip].massOperator;

        Mtp.noalias()+= (density_solid * specific_heat_capacity_solid * (1 - porosity)*dTdp+
            porosity*((1-Sw)*hs*d_rho_s_d_p-rho_steam*hs*dSwdp+rho_steam*(1-Sw)*dhs_dp+
                d_rho_w_d_p*Sw*hw+rho_water*dSwdp*hw+rho_water*Sw*dhw_dp)) * _ip_data[ip].massOperator;
        Mth.noalias() +=
            (density_solid * specific_heat_capacity_solid * (1 - porosity)*dTdh +
                porosity*((1 - Sw)*hs*d_rho_s_d_h - rho_steam*hs*dSwdh + rho_steam*(1 - Sw)*dhs_dh +
                    d_rho_w_d_h*Sw*hw + rho_water*dSwdh*hw + rho_water*Sw*dhw_dh) )* _ip_data[ip].massOperator;

        auto const thermal_conductivity_solid =
            _process_data.thermal_conductivity_solid(t, pos)[0];
        auto const thermal_conductivity_fluid =
            _process_data.thermal_conductivity_fluid(t, pos)[0];

        double const thermal_conductivity =
            thermal_conductivity_solid * (1 - porosity) +
            thermal_conductivity_fluid * porosity;

        // nonwet
        double const k_rel_nonwet =
            _process_data.material->getNonwetRelativePermeability(
                t, pos, p_int_pt, temperature, Sw);
        double const mu_nonwet =
            _process_data.material->getGasViscosity(p_int_pt, temperature);
        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;

        // wet
        double const k_rel_wet =
            _process_data.material->getWetRelativePermeability(
                t, pos, p_int_pt, temperature, Sw);
        double const mu_wet = _process_data.material->getLiquidViscosity(
            _pressure_wet[ip], temperature);
        double const lambda_wet = k_rel_wet / mu_wet;

        GlobalDimVectorType const velocity_steam =
            _process_data.has_gravity
            ? GlobalDimVectorType(-permeability*lambda_nonwet *
            (_ip_data[ip].dNdx * p_nodal_values - rho_steam * b))
            : GlobalDimVectorType(-permeability*lambda_nonwet * _ip_data[ip].dNdx * p_nodal_values);//K_over_mu * dNdx * p_nodal_values

        GlobalDimVectorType const velocity_water =
            _process_data.has_gravity
            ? GlobalDimVectorType(-permeability*lambda_wet *
            (_ip_data[ip].dNdx * p_nodal_values - rho_water * b))
            : GlobalDimVectorType(-permeability*lambda_wet * _ip_data[ip].dNdx * p_nodal_values);//K_over_mu * dNdx * p_nodal_values

        laplace_operator.noalias() = _ip_data[ip].dNdx.transpose() *
                                     permeability * _ip_data[ip].dNdx *
                                     _ip_data[ip].integration_weight;

        Kgp.noalias() += (rho_steam*lambda_nonwet+rho_water*lambda_wet) * laplace_operator;
        Kgh.noalias() += 0 * laplace_operator;

        Kth.noalias() += (_ip_data[ip].dNdx.transpose() * (thermal_conductivity *dTdh)* _ip_data[ip].dNdx +
            _ip_data[ip].N.transpose() * velocity_steam.transpose() * _ip_data[ip].dNdx *
            rho_steam * dhs_dh + _ip_data[ip].N.transpose() * velocity_water.transpose() * _ip_data[ip].dNdx *
            rho_water*dhw_dh) * _ip_data[ip].integration_weight;
        Ktp.noalias() += (_ip_data[ip].dNdx.transpose() * (thermal_conductivity*dTdp*I  + 
            permeability*lambda_nonwet*rho_steam*hs + permeability*lambda_wet*rho_water*hw)* _ip_data[ip].dNdx +
            _ip_data[ip].N.transpose() * velocity_steam.transpose() * _ip_data[ip].dNdx *
            rho_steam * dhs_dp + _ip_data[ip].N.transpose() * velocity_water.transpose() * _ip_data[ip].dNdx *
            rho_water*dhw_dp) * _ip_data[ip].integration_weight;

        if (_process_data.has_gravity)
        {

            NodalVectorType gravity_operator = _ip_data[ip].dNdx.transpose() *
                                               permeability * b *
                                               _ip_data[ip].integration_weight;
            Bg.noalias() +=
                (rho_steam  * rho_steam * lambda_nonwet+rho_water*rho_water*lambda_wet) * gravity_operator;
            Bt.noalias() +=  (rho_steam*lambda_nonwet*hs*rho_steam+rho_water*lambda_wet*hw*rho_water) * gravity_operator;
        }  // end of has gravity
    }
    if (_process_data.has_mass_lumping)
    {
        for (unsigned row = 0; row < Mgp.cols(); row++)
        {
            for (unsigned column = 0; column < Mgp.cols(); column++)
            {
                if (row != column)
                {
                    Mgp(row, row) += Mgp(row, column);
                    Mgp(row, column) = 0.0;
                    Mgh(row, row) += Mgh(row, column);
                    Mgh(row, column) = 0.0;
                    Mtp(row, row) += Mtp(row, column);
                    Mtp(row, column) = 0.0;
                    Mth(row, row) += Mth(row, column);
                    Mth(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // end of namespace
}  // end of namespace
