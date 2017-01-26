/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
* Common convenitions for naming:
* X_air_nonwet mass fraction of air in nonwetting phase
* x_air_nonwet molar fraction of air in nonwetting phase
*
*/
#pragma once

#include "ThermalTwoPhaseFlowWithPPLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "ThermalTwoPhaseFlowWithPPProcessData.h"

using MaterialLib::PhysicalConstant::CelsiusZeroInKelvin;
using MaterialLib::PhysicalConstant::MolarMass::Water;
using MaterialLib::PhysicalConstant::MolarMass::Air;
using MaterialLib::PhysicalConstant::IdealGasConstant;
namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void ThermalTwoPhaseFlowWithPPLocalAssembler<
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
    auto Mgt = local_M.template block<nonwet_pressure_size, temperature_size>(
        nonwet_pressure_matrix_index, temperature_matrix_index);

    auto Mlp = local_M.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mlpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);
    auto Mlt = local_M.template block<cap_pressure_size, temperature_size>(
        cap_pressure_matrix_index, temperature_matrix_index);

    auto Mep = local_M.template block<temperature_size, nonwet_pressure_size>(
        temperature_matrix_index, nonwet_pressure_matrix_index);
    auto Mepc = local_M.template block<temperature_size, cap_pressure_size>(
        temperature_matrix_index, cap_pressure_matrix_index);
    auto Met = local_M.template block<temperature_size, temperature_size>(
        temperature_matrix_index, temperature_matrix_index);

    NodalMatrixType laplace_operator;
    laplace_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kgp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Kgpc = local_K.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);
    auto Kgt = local_K.template block<nonwet_pressure_size, temperature_size>(
        nonwet_pressure_matrix_index, temperature_matrix_index);

    auto Klp = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Klpc = local_K.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);
    auto Klt = local_K.template block<cap_pressure_size, temperature_size>(
        cap_pressure_matrix_index, temperature_matrix_index);

    auto Kep = local_K.template block<temperature_size, nonwet_pressure_size>(
        temperature_matrix_index, nonwet_pressure_matrix_index);
    auto Kepc = local_K.template block<temperature_size, cap_pressure_size>(
        temperature_matrix_index, cap_pressure_matrix_index);
    auto Ket = local_K.template block<temperature_size, temperature_size>(
        temperature_matrix_index, temperature_matrix_index);

    auto Bg = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);
    auto Bl =
        local_b.template segment<cap_pressure_size>(cap_pressure_matrix_index);
    auto Be =
        local_b.template segment<temperature_size>(temperature_matrix_index);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id =
        _process_data._material->getMaterialID(pos.getElementID().get());

    auto const num_nodes = ShapeFunction::NPOINTS;
    auto pg_nodal_values =
        Eigen::Map<const Eigen::VectorXd>(&local_x[0], num_nodes);
    auto pc_nodal_values =
        Eigen::Map<const Eigen::VectorXd>(&local_x[num_nodes], num_nodes);

    const Eigen::MatrixXd& perm = _process_data._material->getPermeability(
        t, pos, _element.getDimension());
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

        double pg_int_pt = 0.;
        double pc_int_pt = 0.;
        double T_int_pt = 0.0;  // Temperature
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pg_int_pt, pc_int_pt,
                                         T_int_pt);

        auto const& wp = _integration_method.getWeightedPoint(ip);

        double const rho_water =
            _process_data._material->getLiquidDensity(pg_int_pt, T_int_pt);

        double const Sw = (pc_int_pt < 0)
                              ? 1.0
                              : _process_data._material->getSaturation(
                                    t, pos, pg_int_pt, T_int_pt, pc_int_pt);

        _saturation[ip] = Sw;
        double dSwdPc =
            (pc_int_pt > _process_data._material->getCapillaryPressure(
                             t, pos, pg_int_pt, T_int_pt, 0.0))
                ? 0.0
                : _process_data._material->getSaturationDerivative(
                      t, pos, pg_int_pt, T_int_pt, Sw);

        /// calculate the water vapor pressure
        double const p_vapor_nonwet =
            _process_data._material->calculateVaporPressureNonwet(pc_int_pt,
                                                                  T_int_pt);
        /// partial pressure of air
        double const p_air_nonwet = pg_int_pt - p_vapor_nonwet;
        /// molar fraction of air in nonwet phase
        double const x_air_nonwet = p_air_nonwet / pg_int_pt;
        /// molar fraction of water vapor in nonwet phase
        double const x_vapor_nonwet = p_vapor_nonwet / pg_int_pt;

        /// mass fraction of air in the nonwet phase
        double const X_air_nonwet =
            x_air_nonwet / (x_air_nonwet + x_vapor_nonwet * Water / Air);
        /// mass density of air in nonwet phase
        double const rho_air_nonwet =
            p_air_nonwet * Air / IdealGasConstant / T_int_pt;
        double const rho_mol_nonwet = pg_int_pt / IdealGasConstant / T_int_pt;
        double const rho_mol_water = rho_water / Water;
        /// calculate derivatives
        double const dRhoMolNonwetdPg = 1 / IdealGasConstant / T_int_pt;
        double const dRhoMolAirNonwetdPg = dRhoMolNonwetdPg;
        double const dPvaporNonwetdT =
            _process_data._material->calculateDerivativedPgwdT(pc_int_pt,
                                                               T_int_pt);
        double const dPvaporNonwetdPc =
            _process_data._material->calculateDerivativedPgwdPC(pc_int_pt,
                                                                T_int_pt);
        /// Derivative
        double const dRhoMolNonwetdT =
            -pg_int_pt / IdealGasConstant / T_int_pt / T_int_pt;
        double const dxAirNonwetdPg = p_vapor_nonwet / pg_int_pt / pg_int_pt;
        double const dxAirNonwetdPc = -dPvaporNonwetdPc / pg_int_pt;
        double const dxAirNonwetdT = -dPvaporNonwetdT / pg_int_pt;

        /// mass density of the nonwet phase
        double const rho_nonwet_air =
            p_air_nonwet * Air / IdealGasConstant / T_int_pt;
        double const rho_nonwet_vapor =
            p_vapor_nonwet * Water / IdealGasConstant / T_int_pt;
        double const rho_nonwet = rho_nonwet_air + rho_nonwet_vapor;
        double const rho_wet = rho_water;
        double const rho_solid = _process_data._density_solid(t, pos)[0];
        /// Derivative of nonwet phase density in terms of T
        double const dRhoNonwetdT =
            _process_data._material->calculatedRhoNonwetdT(
                p_air_nonwet, p_vapor_nonwet, pc_int_pt, T_int_pt);
        /// Derivative of nonwet phase density in terms of Pc
        double const dRhoNonwetdPc =
            Water * dPvaporNonwetdPc / IdealGasConstant / T_int_pt;
        /// Derivative of mass density of air in nonwet phase
        double const dRhoNonwetAirdPg = Air / IdealGasConstant / T_int_pt;
        double const dRhoNonwetAirdPc =
            -Air * dPvaporNonwetdPc / IdealGasConstant / T_int_pt;
        double const dRhoNonwetAirdT =
            -p_air_nonwet * Air / IdealGasConstant / T_int_pt / T_int_pt -
            Air * dPvaporNonwetdT / IdealGasConstant / T_int_pt;
        /// Derivative of mass density of water vapor in nonwet phase
        double const dRhoNonwetVapordPc =
            Water * dPvaporNonwetdPc / IdealGasConstant / T_int_pt;
        double const dRhoNonwetVapordT =
            -p_vapor_nonwet * Water / IdealGasConstant / T_int_pt / T_int_pt +
            Water * dPvaporNonwetdT / IdealGasConstant / T_int_pt;

        _pressure_wetting[ip] = pg_int_pt - pc_int_pt;
        /// heat capacity of nonwet phase
        double const heat_capacity_dry_air =
            _process_data._specific_heat_capacity_dry_air(t, pos)[0];
        const double heat_capacity_water_vapor =
            _process_data._specific_heat_capacity_water_vapor(t, pos)[0];

        double const heat_capacity_water =
            _process_data._specific_heat_capacity_water(t, pos)[0];
        double const heat_capacity_solid =
            _process_data._specific_heat_capacity_solid(t, pos)[0];
        /// specific internal energy and specific enthalpy
        double const latent_heat_evaporation =
            _process_data._latent_heat_evaporation(t, pos)[0];
        double const enthalpy_nonwet_air =
            heat_capacity_dry_air * (T_int_pt - CelsiusZeroInKelvin) +
            IdealGasConstant * (T_int_pt - CelsiusZeroInKelvin) / Air;
        double const enthalpy_wet =
            heat_capacity_water * (T_int_pt - CelsiusZeroInKelvin);
        double const enthalpy_nonwet_vapor =
            heat_capacity_water_vapor * (T_int_pt - CelsiusZeroInKelvin) +
            latent_heat_evaporation;
        double const enthalpy_nonwet =
            enthalpy_nonwet_air * X_air_nonwet +
            enthalpy_nonwet_vapor * (1 - X_air_nonwet);
        double const internal_energy_nonwet =
            enthalpy_nonwet - pg_int_pt / rho_nonwet;
        double const internal_energy_wet = enthalpy_wet;
        /// Derivative
        double const dEnthalpyNonwetAirdT =
            heat_capacity_dry_air + IdealGasConstant / Air;
        double const dEnthalpyNonwetdT =
            heat_capacity_water * (1 - X_air_nonwet) +
            dEnthalpyNonwetAirdT * X_air_nonwet;
        // Assemble M matrix
        // nonwetting
        double const porosity = _process_data._material->getPorosity(
            t, pos, pg_int_pt, T_int_pt, 0);

        Mgp.noalias() += porosity *
                         ((1 - Sw) * (rho_mol_nonwet * dxAirNonwetdPg +
                                      x_air_nonwet * dRhoMolNonwetdPg)) *
                         _ip_data[ip].massOperator;
        Mgpc.noalias() += porosity *
                          ((1 - Sw) * rho_mol_nonwet * dxAirNonwetdPc -
                           rho_mol_nonwet * x_air_nonwet * dSwdPc) *
                          _ip_data[ip].massOperator;
        Mgt.noalias() += porosity *
                         ((1 - Sw) * (rho_mol_nonwet * dxAirNonwetdT +
                                      x_air_nonwet * dRhoMolNonwetdT)) *
                         _ip_data[ip].massOperator;

        Mlpc.noalias() +=
            porosity *
            ((1 - Sw) * dPvaporNonwetdPc / IdealGasConstant / T_int_pt +
             rho_mol_nonwet * x_vapor_nonwet * (-dSwdPc) +
             dSwdPc * rho_mol_water) *
            _ip_data[ip].massOperator;
        Mlt.noalias() +=
            porosity *
            ((1 - Sw) *
             (dPvaporNonwetdT / IdealGasConstant / T_int_pt -
              p_vapor_nonwet / IdealGasConstant / T_int_pt / T_int_pt)) *
            _ip_data[ip].massOperator;

        Mep.noalias() +=
            porosity *
            ((x_air_nonwet * Air + x_vapor_nonwet * Water) * dRhoMolNonwetdPg *
                 enthalpy_nonwet -
             rho_mol_nonwet * (Water - Air) * dxAirNonwetdPg * enthalpy_nonwet -
             1) *
            (1 - Sw) * _ip_data[ip].massOperator;
        Mepc.noalias() += porosity * (rho_wet * internal_energy_wet -
                                      rho_nonwet * internal_energy_nonwet) *
                              dSwdPc * _ip_data[ip].massOperator +
                          porosity * ((Water - Air) * enthalpy_nonwet /
                                      IdealGasConstant / T_int_pt) *
                              (1 - Sw) * dPvaporNonwetdPc *
                              _ip_data[ip].massOperator;
        Met.noalias() +=
            ((1 - porosity) * rho_solid * heat_capacity_solid +
             porosity * ((1 - Sw) * (dRhoNonwetdT * enthalpy_nonwet +
                                     rho_nonwet * dEnthalpyNonwetdT) +
                         Sw * rho_wet * heat_capacity_water)) *
            _ip_data[ip].massOperator;

        // nonwet
        double const k_rel_nonwet =
            _process_data._material->getNonwetRelativePermeability(
                t, pos, _pressure_wetting[ip], T_int_pt, Sw);
        double const mu_nonwet = _process_data._material->getGasViscosity(
            _pressure_wetting[ip], T_int_pt);
        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;
        double const diffusion_coeff_component_air =
            _process_data._diffusion_coeff_component_b(t, pos)[0];

        // wet
        double const k_rel_wet =
            _process_data._material->getWetRelativePermeability(
                t, pos, pg_int_pt, T_int_pt, Sw);
        double const mu_wet =
            _process_data._material->getLiquidViscosity(pg_int_pt, T_int_pt);
        double const lambda_wet = k_rel_wet / mu_wet;

        GlobalDimVectorType const velocity_nonwet =
            -lambda_nonwet * permeability * (sm.dNdx * pg_nodal_values);
        GlobalDimVectorType const velocity_wet =
            -lambda_wet * permeability *
            (sm.dNdx * (pg_nodal_values - pc_nodal_values));

        laplace_operator.noalias() = sm.dNdx.transpose() * permeability *
                                     sm.dNdx * _ip_data[ip].integration_weight;

        /*Kep.noalias() +=
            integration_factor * sm.N.transpose() *
                ((x_air_nonwet * Air +
                  x_vapor_nonwet * Water) *
                     dRhoMolNonwetdPg * enthalpy_nonwet -
                 rho_mol_nonwet * (Water - Air) *
                     dxAirNonwetdPg * enthalpy_nonwet) *
                velocity_nonwet.transpose() * sm.dNdx +
            integration_factor * sm.N.transpose() * rho_wet *
                heat_capacity_water * velocity_wet.transpose() * sm.dNdx;
        Kepc.noalias() += integration_factor * sm.N.transpose() *
                          ((Water - Air) * enthalpy_nonwet /
                          IdealGasConstant/ T_int_pt) *
                          dPvaporNonwetdPc * velocity_nonwet.transpose() *
                          sm.dNdx;*/
        Ket.noalias() += _ip_data[ip].integration_weight * sm.N.transpose() *
                             (dRhoNonwetdT * enthalpy_nonwet +
                              rho_nonwet * dEnthalpyNonwetdT) *
                             velocity_nonwet.transpose() * sm.dNdx +
                         _ip_data[ip].integration_weight * sm.N.transpose() *
                             heat_capacity_water * rho_water *
                             velocity_wet.transpose() * sm.dNdx;
        /// heat conductivity
        double const heat_conductivity_dry_solid =
            _process_data._heat_conductivity_dry_solid(t, pos)[0];
        double const heat_conductivity_wet_solid =
            _process_data._heat_conductivity_wet_solid(t, pos)[0];
        double heat_conductivity_unsaturated =
            _process_data._material->calculateUnsatHeatConductivity(
                t, pos, Sw, heat_conductivity_dry_solid,
                heat_conductivity_wet_solid);
        // Laplace
        Kgp.noalias() +=
            (rho_mol_nonwet * x_air_nonwet * lambda_nonwet) * laplace_operator +
            ((1 - Sw) * porosity * diffusion_coeff_component_air *
             rho_mol_nonwet * dxAirNonwetdPg) *
                _ip_data[ip].diffusionOperator;
        Kgpc.noalias() += ((1 - Sw) * porosity * diffusion_coeff_component_air *
                           rho_mol_nonwet * dxAirNonwetdPc) *
                          _ip_data[ip].diffusionOperator;
        Kgt.noalias() += ((1 - Sw) * porosity * diffusion_coeff_component_air *
                          rho_mol_nonwet * dxAirNonwetdT) *
                         _ip_data[ip].diffusionOperator;

        Klp.noalias() += (rho_mol_nonwet * x_vapor_nonwet * lambda_nonwet) *
                             laplace_operator +
                         rho_mol_water * lambda_wet * laplace_operator -
                         ((1 - Sw) * porosity * diffusion_coeff_component_air *
                          rho_mol_nonwet * dxAirNonwetdPg) *
                             _ip_data[ip].diffusionOperator;
        Klpc.noalias() += (-rho_mol_water * lambda_wet * laplace_operator) -
                          ((1 - Sw) * porosity * diffusion_coeff_component_air *
                           rho_mol_nonwet * dxAirNonwetdPc) *
                              _ip_data[ip].diffusionOperator;
        Klt.noalias() += -((1 - Sw) * porosity * diffusion_coeff_component_air *
                           rho_mol_nonwet * dxAirNonwetdT) *
                         _ip_data[ip].diffusionOperator;

        Kep.noalias() += (lambda_nonwet * rho_nonwet * enthalpy_nonwet +
                          lambda_wet * rho_wet * enthalpy_wet) *
                             laplace_operator +
                         (1 - Sw) * porosity * diffusion_coeff_component_air *
                             rho_mol_nonwet * (Air * enthalpy_nonwet_air -
                                               Water * enthalpy_nonwet_vapor) *
                             dxAirNonwetdPg * _ip_data[ip].diffusionOperator;
        Kepc.noalias() +=
            -lambda_wet * enthalpy_wet * rho_wet * laplace_operator +
            (1 - Sw) * porosity * diffusion_coeff_component_air *
                rho_mol_nonwet *
                (Air * enthalpy_nonwet_air - Water * enthalpy_nonwet_vapor) *
                dxAirNonwetdPc * _ip_data[ip].diffusionOperator;
        Ket.noalias() += sm.dNdx.transpose() * heat_conductivity_unsaturated *
                             sm.dNdx * _ip_data[ip].integration_weight +
                         (1 - Sw) * porosity * diffusion_coeff_component_air *
                             rho_mol_nonwet * (Air * enthalpy_nonwet_air -
                                               Water * enthalpy_nonwet_vapor) *
                             dxAirNonwetdT * _ip_data[ip].diffusionOperator;

        if (_process_data._has_gravity)
        {
            auto const& b = _process_data._specific_body_force;
            Bg.noalias() +=
                (rho_mol_nonwet * x_air_nonwet * lambda_nonwet * rho_nonwet) *
                sm.dNdx.transpose() * permeability * b * _ip_data[ip].integration_weight;
            Bl.noalias() +=
                (rho_mol_water * lambda_wet * rho_wet +
                 rho_mol_nonwet * x_vapor_nonwet * lambda_nonwet * rho_nonwet) *
                sm.dNdx.transpose() * permeability * b * _ip_data[ip].integration_weight;
            Be.noalias() +=
                (lambda_nonwet * rho_nonwet * rho_nonwet * enthalpy_nonwet +
                 lambda_wet * rho_wet * rho_wet * enthalpy_wet) *
                sm.dNdx.transpose() * permeability * b * _ip_data[ip].integration_weight;
        }  // end of has gravity
    }
    if (_process_data._has_mass_lumping)
    {
        for (unsigned row = 0; row < Mgp.cols(); row++)
        {
            for (unsigned column = 0; column < Mgp.cols(); column++)
            {
                if (row != column)
                {
                    Mgp(row, row) += Mgp(row, column);
                    Mgp(row, column) = 0.0;
                    Mgpc(row, row) += Mgpc(row, column);
                    Mgpc(row, column) = 0.0;
                    Mgt(row, row) += Mgt(row, column);
                    Mgt(row, column) = 0.0;
                    Mlp(row, row) += Mlp(row, column);
                    Mlp(row, column) = 0.0;
                    Mlpc(row, row) += Mlpc(row, column);
                    Mlpc(row, column) = 0.0;
                    Mlt(row, row) += Mlt(row, column);
                    Mlt(row, column) = 0.0;
                    Mep(row, row) += Mep(row, column);
                    Mep(row, column) = 0.0;
                    Mepc(row, row) += Mepc(row, column);
                    Mepc(row, column) = 0.0;
                    Met(row, row) += Met(row, column);
                    Met(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // end of namespace
}  // end of namespace
