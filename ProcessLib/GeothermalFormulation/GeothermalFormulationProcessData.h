/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once
#include "GeothermalFormulationMaterialProperties.h"

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace GeothermalFormulation
{
struct GeothermalFormulationProcessData
{
    GeothermalFormulationProcessData(
        Eigen::VectorXd const specific_body_force_,
        bool const has_gravity_,
        bool const has_mass_lumping_,
        Parameter<double> const& temperature_,
        std::unique_ptr<GeothermalFormulationMaterialProperties>&& material_,
        ProcessLib::Parameter<double> const& density_solid_,
        ProcessLib::Parameter<double> const& specific_heat_capacity_solid_,
        ProcessLib::Parameter<double> const& thermal_conductivity_solid_,
        ProcessLib::Parameter<double> const& thermal_conductivity_fluid_)
        : specific_body_force(specific_body_force_),
          has_gravity(has_gravity_),
          has_mass_lumping(has_mass_lumping_),
          temperature(temperature_),
          material(std::move(material_)),
          density_solid(density_solid_),
          specific_heat_capacity_solid(specific_heat_capacity_solid_),
          thermal_conductivity_solid(thermal_conductivity_solid_),
          thermal_conductivity_fluid(thermal_conductivity_fluid_)

    {
    }

    GeothermalFormulationProcessData(GeothermalFormulationProcessData&& other)
        : specific_body_force(other.specific_body_force),
          has_gravity(other.has_gravity),
          has_mass_lumping(other.has_mass_lumping),
          temperature(other.temperature),
          material(std::move(other.material)),
          density_solid(other.density_solid),
          specific_heat_capacity_solid(other.specific_heat_capacity_solid),
          thermal_conductivity_solid(other.thermal_conductivity_solid),
          thermal_conductivity_fluid(other.thermal_conductivity_fluid)
    {
    }

    //! Copies are forbidden.
    GeothermalFormulationProcessData(GeothermalFormulationProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(GeothermalFormulationProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(GeothermalFormulationProcessData&&) = delete;

    //! Specific body forces applied to solid and fluid.
    //! It is usually used to apply gravitational forces.
    //! A vector of displacement dimension's length.
    Eigen::VectorXd const specific_body_force;

    bool const has_gravity;

    //! Enables lumping of the mass matrix.
    bool const has_mass_lumping;
    Parameter<double> const& temperature;
    std::unique_ptr<GeothermalFormulationMaterialProperties> material;
    Parameter<double> const& density_solid;
    Parameter<double> const& specific_heat_capacity_solid;
    Parameter<double> const& thermal_conductivity_solid;
    Parameter<double> const& thermal_conductivity_fluid;
};

}  // namespace GeothermalFormulation
}  // namespace ProcessLib
