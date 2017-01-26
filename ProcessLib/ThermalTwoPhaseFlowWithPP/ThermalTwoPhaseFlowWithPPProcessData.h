/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ThermalTwoPhaseFlowWithPPMaterialProperties.h"
namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace ThermalTwoPhaseFlowWithPP
{
struct ThermalTwoPhaseFlowWithPPProcessData
{
    ThermalTwoPhaseFlowWithPPProcessData(
        Eigen::VectorXd const specific_body_force_,
        bool const has_gravity_,
        bool const has_mass_lumping_,
        Parameter<double> const& diffusion_coeff_component_b_,
        Parameter<double> const& diffusion_coeff_component_a_,
        Parameter<double> const& specific_heat_capacity_solid_,
        Parameter<double> const& specific_heat_capacity_water_,
        Parameter<double> const& specific_heat_capacity_dry_air_,
        Parameter<double> const& specific_heat_capacity_water_vapor_,
        Parameter<double> const& heat_conductivity_dry_solid_,
        Parameter<double> const& heat_conductivity_wet_solid_,
        Parameter<double> const& density_solid_,
        Parameter<double> const& latent_heat_evaporation_,
        std::unique_ptr<ThermalTwoPhaseFlowWithPPMaterialProperties>&&
            material_)
        : _specific_body_force(specific_body_force_),
          _has_gravity(has_gravity_),
          _has_mass_lumping(has_mass_lumping_),
          _diffusion_coeff_component_b(diffusion_coeff_component_b_),
          _diffusion_coeff_component_a(diffusion_coeff_component_a_),
          _specific_heat_capacity_solid(specific_heat_capacity_solid_),
          _specific_heat_capacity_water(specific_heat_capacity_water_),
          _specific_heat_capacity_dry_air(specific_heat_capacity_dry_air_),
          _specific_heat_capacity_water_vapor(
              specific_heat_capacity_water_vapor_),
          _heat_conductivity_dry_solid(heat_conductivity_dry_solid_),
          _heat_conductivity_wet_solid(heat_conductivity_wet_solid_),
          _density_solid(density_solid_),
          _latent_heat_evaporation(latent_heat_evaporation_),
          _material(std::move(material_))

    {
    }

    ThermalTwoPhaseFlowWithPPProcessData(
        ThermalTwoPhaseFlowWithPPProcessData&& other)
        : _specific_body_force(other._specific_body_force),
          _has_gravity(other._has_gravity),
          _has_mass_lumping(other._has_mass_lumping),
          _diffusion_coeff_component_b(other._diffusion_coeff_component_b),
          _diffusion_coeff_component_a(other._diffusion_coeff_component_a),
          _specific_heat_capacity_solid(other._specific_heat_capacity_solid),
          _specific_heat_capacity_water(other._specific_heat_capacity_water),
          _specific_heat_capacity_dry_air(
              other._specific_heat_capacity_dry_air),
          _specific_heat_capacity_water_vapor(
              other._specific_heat_capacity_water_vapor),
          _heat_conductivity_dry_solid(other._heat_conductivity_dry_solid),
          _heat_conductivity_wet_solid(other._heat_conductivity_wet_solid),
          _density_solid(other._density_solid),
          _latent_heat_evaporation(other._latent_heat_evaporation),
          _material(std::move(other._material))
    {
    }

    //! Copies are forbidden.
    ThermalTwoPhaseFlowWithPPProcessData(
        ThermalTwoPhaseFlowWithPPProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermalTwoPhaseFlowWithPPProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermalTwoPhaseFlowWithPPProcessData&&) = delete;
    Eigen::VectorXd const _specific_body_force;

    bool const _has_gravity;
    bool const _has_mass_lumping;
    Parameter<double> const& _diffusion_coeff_component_b;
    Parameter<double> const& _diffusion_coeff_component_a;
    Parameter<double> const& _specific_heat_capacity_solid;
    Parameter<double> const& _specific_heat_capacity_water;
    Parameter<double> const& _specific_heat_capacity_dry_air;
    Parameter<double> const& _specific_heat_capacity_water_vapor;
    Parameter<double> const& _heat_conductivity_dry_solid;
    Parameter<double> const& _heat_conductivity_wet_solid;
    Parameter<double> const& _density_solid;
    Parameter<double> const& _latent_heat_evaporation;
    std::unique_ptr<ThermalTwoPhaseFlowWithPPMaterialProperties> _material;
};

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
