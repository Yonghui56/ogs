/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   VanGenuchtenCapillaryPressureSaturation.cpp
 *
 *  Created on October 28, 2016, 6:05 PM
 */

#include "VanGenuchtenCapillaryPressureSaturation.h"

#include <cmath>

#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace PorousMedium
{
double VanGenuchtenCapillaryPressureSaturation::getCapillaryPressure(
    const double saturation) const
{
    const double S = MathLib::limitValueInInterval(
        saturation, _Sr + _minor_offset, _Smax - _minor_offset);
    const double Se = (S - _Sr) / (_Smax - _Sr);
    const double pc =
        _pb * std::pow(std::pow(Se, (-1.0 / _mm)) - 1.0, 1.0 - _mm);
    return MathLib::limitValueInInterval(pc, _minor_offset, _Pc_max);
}
double VanGenuchtenCapillaryPressureSaturation::getRegularizedCapillaryPressure(
	const double saturation) const
{
	double Sg = 1 - saturation;
	if (Sg <= 1- _Sr && Sg >= 0)
	{
		return getPc_bar_vG_Sg(Sg);
	}
	else if (Sg < 0.0)
	{
		return getPc_bar_vG_Sg(0.0) + get_dPCdS_vG_bar(0.0)*(Sg - 0.0);
	}
	else
	{
		return getPc_bar_vG_Sg(1 - _Sr) + get_dPCdS_vG_bar(1 - _Sr)*(Sg - 1 + _Sr);
	}
}

double VanGenuchtenCapillaryPressureSaturation::getSaturation(
    const double capillary_pressure) const
{
    const double pc =
                (capillary_pressure < 0.0) ? _minor_offset : capillary_pressure;
    double Se = std::pow(pc / _pb, 1.0 / (1.0 - _mm)) + 1.0;
    Se = std::pow(Se, -_mm);
    const double S = Se * (_Smax - _Sr) + _Sr;
    return MathLib::limitValueInInterval(S, _Sr + _minor_offset,
                                         _Smax - _minor_offset);
}

double VanGenuchtenCapillaryPressureSaturation::getdPcdS(
    const double saturation) const
{
    const double S = MathLib::limitValueInInterval(
        saturation, _Sr + _minor_offset, _Smax - _minor_offset);
    const double val1 = std::pow(((S - _Sr) / (_Smax - _Sr)), -1.0 / _mm);
    const double val2 = std::pow(val1 - 1.0, -_mm);
    return _pb * (_mm - 1.0) * val1 * val2 / (_mm * (S - _Sr));
}

double VanGenuchtenCapillaryPressureSaturation::getRegularizedPcdS(const double saturation) const
{
	double const Sg = 1 - saturation;
	if (Sg >= 0.0 && Sg <= 1 - _Sr)
	{
		return -get_dPCdS_vG_bar(Sg);
	}
	else if (Sg < 0.0)
	{
		return -get_dPCdS_vG_bar(0.0);
	}
	else
	{
		return -get_dPCdS_vG_bar(1 - _Sr);
	}
}


}  // end namespace
}  // end namespace
