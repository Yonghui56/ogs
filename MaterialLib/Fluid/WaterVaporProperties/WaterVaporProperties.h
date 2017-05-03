/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   WaterVaporProperties.h
 */

#pragma once

#include <array>
#include <string>

namespace MaterialLib
{
namespace Fluid
{
/**
* This class provides a series of functions to calculate the water vapor
* properties, including the saturation vapor pressure calculation(Kelvin
* equation regularized vapor pressure), and corresponding derivatives
* calculations.
*/
class WaterVaporProperties
{
public:
    WaterVaporProperties() = default;
    /// water vapor saturation pressure
    double calculateSaturatedVaporPressure(const double T) const;

    /// Derivative of SaturatedVaporPressure in terms of T
    double calculateDerivativedPsatdT(const double T) const;
    /// partial water vapor pressure in nonwetting phase
    /// Kelvin equation
    double calculateVaporPressureNonwet(
        const double pc /*capillary pressure*/, const double T /*temperature*/,
        const double rho_mass_h2o /*mass density of water*/) const;
    /// Derivative of partial vapor pressure in terms of T
    double calculateDerivativedPgwdT(const double pc, const double T,
                                     const double rho_mass_h2o) const;
    /// Derivative of partial vapor pressure in terms of PC
    double calculateDerivativedPgwdPC(const double pc, const double T,
                                      const double rho_mass_h2o) const;
    /// Derivative of vapor density in terms of T
    double calculatedRhoNonwetdT(const double p_air_nonwet,
                                 const double p_vapor_nonwet, const double pc,
                                 const double T,
                                 const double rho_mass_h2o) const;
    /// Specific enthalpy of water vapor
    double getWaterVaporEnthalpySimple(const double temperature,
        const double heat_capacity_water_vapor,
        const double /*pressure*/,
        const double /*latent_heat_evaporation*/) const;
};

}  // end namespace
}  // end namespace
