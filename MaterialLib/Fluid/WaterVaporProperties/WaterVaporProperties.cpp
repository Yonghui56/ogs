/**
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WaterVaporProperties.cpp
 *
 */

#include "WaterVaporProperties.h"
#include <array>
#include <cmath>
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialLib
{
using PhysicalConstant::MolarMass::Water;
using PhysicalConstant::MolarMass::Air;
using PhysicalConstant::IdealGasConstant;
using PhysicalConstant::CelsiusZeroInKelvin;
namespace Fluid
{
static const double temperature_0 = 373.15;  /// reference temperature in [K]
static const double p_0 = 101325.0;          /// reference pressure
static const double h_wg = 2258000.0;  /// latent heat of water evaporation

double WaterVaporProperties::calculateSaturatedVaporPressure(
    const double T) const
{
    return p_0 * std::exp(((1 / temperature_0) - (1 / T)) * Water * h_wg /
                           IdealGasConstant);
}

double WaterVaporProperties::calculateDerivativedPsatdT(const double T) const
{
    return p_0 * (Water * h_wg / IdealGasConstant) * (1. / T / T) *
           std::exp(((1. / temperature_0) - (1. / T)) * Water * h_wg /
                    IdealGasConstant);
}

double WaterVaporProperties::calculateVaporPressureNonwet(
    const double pc /*capillary pressure*/, const double T /*temperature*/,
    const double rho_mass_h2o /*mass density of water*/) const
{
    const double p_sat = calculateSaturatedVaporPressure(T);
    const double c_w = Water / IdealGasConstant / T;
    return p_sat * std::exp(-pc * c_w / rho_mass_h2o);
}
double WaterVaporProperties::calculateDerivativedPgwdT(
    const double pc, const double T, const double rho_mass_h2o) const
{
    const double c_w = Water / IdealGasConstant / T;
    const double p_sat = calculateSaturatedVaporPressure(T);
    const double dPsatdT = calculateDerivativedPsatdT(T);
    return dPsatdT * std::exp(-pc * c_w / rho_mass_h2o) +
           p_sat * std::exp(-pc * c_w / rho_mass_h2o) *
               (pc * Water / rho_mass_h2o / IdealGasConstant / T / T);
}
double WaterVaporProperties::calculateDerivativedPgwdPC(
    const double pc, const double T, const double rho_mass_h2o) const
{
    const double c_w = Water / IdealGasConstant / T;
    const double p_sat = calculateSaturatedVaporPressure(T);
    return p_sat * std::exp(-pc * c_w / rho_mass_h2o) * (-c_w / rho_mass_h2o);
}
double WaterVaporProperties::calculatedRhoNonwetdT(
    const double p_air_nonwet, const double p_vapor_nonwet, const double pc,
    const double T, const double rho_mass_h2o) const
{
    const double dPgwdT = calculateDerivativedPgwdT(pc, T, rho_mass_h2o);
    return -((p_air_nonwet * Air + p_vapor_nonwet * Water) / IdealGasConstant /
             T / T) +
           (Water - Air) * dPgwdT / IdealGasConstant / T;
}
double WaterVaporProperties::getWaterVaporEnthalpySimple(const double temperature,
    const double heat_capacity_water_vapor,
    const double /*pressure*/) const
{
    return heat_capacity_water_vapor * (temperature - CelsiusZeroInKelvin) +
        h_wg;
}

}  // end namespace
}  // end namespace
