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
        : specific_body_force(specific_body_force_),
          has_gravity(has_gravity_),
          has_mass_lumping(has_mass_lumping_),
          diffusion_coeff_component_b(diffusion_coeff_component_b_),
          diffusion_coeff_component_a(diffusion_coeff_component_a_),
          specific_heat_capacity_solid(specific_heat_capacity_solid_),
          specific_heat_capacity_water(specific_heat_capacity_water_),
          specific_heat_capacity_dry_air(specific_heat_capacity_dry_air_),
          specific_heat_capacity_water_vapor(
              specific_heat_capacity_water_vapor_),
          heat_conductivity_dry_solid(heat_conductivity_dry_solid_),
          heat_conductivity_wet_solid(heat_conductivity_wet_solid_),
          density_solid(density_solid_),
          latent_heat_evaporation(latent_heat_evaporation_),
          material(std::move(material_))

    {
    }

    ThermalTwoPhaseFlowWithPPProcessData(
        ThermalTwoPhaseFlowWithPPProcessData&& other)
        : specific_body_force(other.specific_body_force),
          has_gravity(other.has_gravity),
          has_mass_lumping(other.has_mass_lumping),
          diffusion_coeff_component_b(other.diffusion_coeff_component_b),
          diffusion_coeff_component_a(other.diffusion_coeff_component_a),
          specific_heat_capacity_solid(other.specific_heat_capacity_solid),
          specific_heat_capacity_water(other.specific_heat_capacity_water),
          specific_heat_capacity_dry_air(
              other.specific_heat_capacity_dry_air),
          specific_heat_capacity_water_vapor(
              other.specific_heat_capacity_water_vapor),
          heat_conductivity_dry_solid(other.heat_conductivity_dry_solid),
          heat_conductivity_wet_solid(other.heat_conductivity_wet_solid),
          density_solid(other.density_solid),
          latent_heat_evaporation(other.latent_heat_evaporation),
          material(std::move(other.material))
    {
    }

    //! Copies are forbidden.
    ThermalTwoPhaseFlowWithPPProcessData(
        ThermalTwoPhaseFlowWithPPProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermalTwoPhaseFlowWithPPProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermalTwoPhaseFlowWithPPProcessData&&) = delete;
    Eigen::VectorXd const specific_body_force;

    bool const has_gravity;
    bool const has_mass_lumping;
    Parameter<double> const& diffusion_coeff_component_b;
    Parameter<double> const& diffusion_coeff_component_a;
    Parameter<double> const& specific_heat_capacity_solid;
    Parameter<double> const& specific_heat_capacity_water;
    Parameter<double> const& specific_heat_capacity_dry_air;
    Parameter<double> const& specific_heat_capacity_water_vapor;
    Parameter<double> const& heat_conductivity_dry_solid;
    Parameter<double> const& heat_conductivity_wet_solid;
    Parameter<double> const& density_solid;
    Parameter<double> const& latent_heat_evaporation;
    std::unique_ptr<ThermalTwoPhaseFlowWithPPMaterialProperties> material;
};

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
