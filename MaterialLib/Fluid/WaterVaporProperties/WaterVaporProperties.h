/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   VaporProperty.h
 *
 * Created on August 12, 2016, 3:34 PM
 */

#pragma once

#include <array>
#include <string>
#include "MaterialLib/PhysicalConstant.h"

using MaterialLib::PhysicalConstant::MolarMass::Water;
using MaterialLib::PhysicalConstant::MolarMass::Air;
using MaterialLib::PhysicalConstant::IdealGasConstant;
namespace MaterialLib
{
namespace Fluid
{
/*
* This class provides a series of functions to calculate the water vapor properties.
* including the saturation vapor pressure calculation(Kelvin equation regulized vapor pressure)
* and corresponding derivatives calculations
*/
class WaterVaporProperties
{
public:
    WaterVaporProperties() = default;
    // water vapor saturation pressure
    const double calculateSaturatedVaporPressure(
            const double T) const
    {
        return _p_0 * exp(((1 / _temperature_0) - (1 / T)) * Water * _h_wg / IdealGasConstant);
    }
    // Derivative of SaturatedVaporPressure in terms of T
    const double calculateDerivativedPsatdT(
        const double T) const
    {
        return _p_0 * (Water * _h_wg / IdealGasConstant) * (1. / T / T) *
            exp(((1. / _temperature_0) - (1. / T)) * Water * _h_wg / IdealGasConstant);
    }
    // partial water vapor pressure in nonwetting phase
    // Kelvin equation
    const double calculateVaporPressureNonwet(
            const double pc/*capillary pressure*/, const double T/*temperature*/, 
        const double rho_mass_h2o/*mass density of water*/) const
    {
        const double p_sat = calculateSaturatedVaporPressure(T);
        const double c_w = Water / IdealGasConstant / T;
        return p_sat * exp(-pc * c_w / rho_mass_h2o);
    }
    // Derivative of partial vapor pressure in terms of T
    const double calculateDerivativedPgwdT(
        const double pc, const double T, const double rho_mass_h2o) const
    {
        const double c_w = Water / IdealGasConstant / T;
        const double p_sat = calculateSaturatedVaporPressure(T);
        const double dPsatdT = calculateDerivativedPsatdT(T);
        return dPsatdT * exp(-pc * c_w / rho_mass_h2o) +
            p_sat * exp(-pc * c_w / rho_mass_h2o) *
            (pc * Water / rho_mass_h2o / IdealGasConstant / T / T);
    }
    // Derivative of partial vapor pressure in terms of PC
    const double calculateDerivativedPgwdPC(
        const double pc, const double T, const double rho_mass_h2o) const
    {
        const double c_w = Water / IdealGasConstant / T;
        const double p_sat = calculateSaturatedVaporPressure(T);
        return p_sat * exp(-pc * c_w / rho_mass_h2o) * (-c_w / rho_mass_h2o);
    }
    // Derivative of vapor density in terms of T
    const double calculatedRhoNonwetdT(
        const double p_air_nonwet, const double p_vapor_nonwet, const double pc,
        const double T, const double rho_mass_h2o) const
    {
        const double dPgwdT = calculateDerivativedPgwdT(pc, T, rho_mass_h2o);
        return -((p_air_nonwet * Air + p_vapor_nonwet * Water) / IdealGasConstant /
            T / T) +
            (Water - Air) * dPgwdT / IdealGasConstant / T;
    }

private:
    const double _temperature_0 = 373.15;// reference temperature in [K]
    const double _p_0 = 101325.0; // reference pressure 
    const double _h_wg = 2258000.0;// latent heat of water evaporation
};

}  // end namespace
}  // end namespace
