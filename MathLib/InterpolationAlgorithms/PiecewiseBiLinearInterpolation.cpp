/**
 * \file
 * \author Thomas Fischer
 * \date   2010-09-07
 * \brief  Implementation of the PiecewiseLinearInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <cmath>
#include <limits>
#include <utility>

#include "BaseLib/Error.h"
#include "BaseLib/quicksort.h"

#include "PiecewiseBiLinearInterpolation.h"

namespace MathLib
{
	PiecewiseBiLinearInterpolation::PiecewiseBiLinearInterpolation(
    std::vector<double>&& supporting_points_x,
	std::vector<double>&& supporting_points_y,
    std::vector<double>&& values_at_supp_pnts,
    bool supp_pnts_sorted)
    : _supp_pnts_x(std::move(supporting_points_x)),
	  _supp_pnts_y(std::move(supporting_points_y)),
      _values_at_supp_pnts(std::move(values_at_supp_pnts))
{
    

    const auto it_x = std::adjacent_find(_supp_pnts_x.begin(), _supp_pnts_x.end());
    if (it_x != _supp_pnts_x.end())
    {
        const std::size_t ix = std::distance(_supp_pnts_x.begin(), it_x);
        OGS_FATAL(
            "Variable %d and variable %d are the same. "
            "Piecewise bi-linear interpolation is not possible\n",
            ix, ix + 1);
    }
	const auto it_y = std::adjacent_find(_supp_pnts_y.begin(), _supp_pnts_y.end());
	if (it_y != _supp_pnts_y.end())
	{
		const std::size_t iy = std::distance(_supp_pnts_y.begin(), it_y);
		OGS_FATAL(
			"Variable %d and variable %d are the same. "
			"Piecewise bi-linear interpolation is not possible\n",
			iy, iy + 1);
	}
}

double PiecewiseBiLinearInterpolation::getBiValue(
    double pnt_x_to_interpolate, double pnt_y_to_interpolate) const
{
    auto const& it_X(std::lower_bound(_supp_pnts_x.begin(), _supp_pnts_x.end(),
                                      pnt_x_to_interpolate));
    std::size_t const x_interval_idx =
        std::distance(_supp_pnts_x.begin(), it_X) - 1;
    auto const& it_Y(std::lower_bound(_supp_pnts_y.begin(), _supp_pnts_y.end(),
                                      pnt_y_to_interpolate));
    std::size_t const y_interval_idx =
        std::distance(_supp_pnts_y.begin(), it_Y) - 1;

	const double x_supp_pnt_size = _supp_pnts_x.size();

    const double f_r1 =
		_values_at_supp_pnts[x_interval_idx  + x_supp_pnt_size * y_interval_idx] *
            (_supp_pnts_x[x_interval_idx + 1] - pnt_x_to_interpolate) /
            (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]) +
		_values_at_supp_pnts[x_interval_idx + 1 + x_supp_pnt_size * y_interval_idx] *
            (pnt_x_to_interpolate- _supp_pnts_x[x_interval_idx]) /
            (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]);
    const double f_r2 =
		_values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * (y_interval_idx + 1)] *
            (_supp_pnts_x[x_interval_idx + 1] - pnt_x_to_interpolate) /
            (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]) +
		_values_at_supp_pnts[x_interval_idx + 1 + x_supp_pnt_size * (y_interval_idx + 1)] *
            (pnt_x_to_interpolate- _supp_pnts_x[x_interval_idx]) /
            (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]);
    const double f_p =
        f_r1 * (_supp_pnts_y[y_interval_idx + 1] - pnt_y_to_interpolate) /
            (_supp_pnts_y[y_interval_idx + 1] - _supp_pnts_y[y_interval_idx]) +
        f_r2 * (pnt_y_to_interpolate- _supp_pnts_y[y_interval_idx]) /
            (_supp_pnts_y[y_interval_idx + 1] - _supp_pnts_y[y_interval_idx]);
    const double f_val_11 =
		_values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * y_interval_idx] *
        (_supp_pnts_x[x_interval_idx + 1] - pnt_x_to_interpolate) *
        (_supp_pnts_y[y_interval_idx + 1] - pnt_y_to_interpolate) /
        (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]) /
        (_supp_pnts_y[y_interval_idx + 1] - _supp_pnts_y[y_interval_idx]);
    const double f_val_21 =
		_values_at_supp_pnts[x_interval_idx + 1 + x_supp_pnt_size * y_interval_idx] *
        (pnt_x_to_interpolate- _supp_pnts_x[x_interval_idx]) *
        (_supp_pnts_y[y_interval_idx + 1] - pnt_y_to_interpolate) /
        (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]) /
        (_supp_pnts_y[y_interval_idx + 1] - _supp_pnts_y[y_interval_idx]);
    const double f_val_12 =
		_values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * (y_interval_idx + 1)] *
        (_supp_pnts_x[x_interval_idx + 1] - pnt_x_to_interpolate) *
        (pnt_y_to_interpolate- _supp_pnts_y[y_interval_idx]) /
        (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]) /
        (_supp_pnts_y[y_interval_idx + 1] - _supp_pnts_y[y_interval_idx]);
    const double f_val_22 =
		_values_at_supp_pnts[x_interval_idx + 1 + x_supp_pnt_size * (y_interval_idx + 1)] *
        (pnt_x_to_interpolate- _supp_pnts_x[x_interval_idx]) *
        (pnt_y_to_interpolate- _supp_pnts_y[y_interval_idx]) /
        (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]) /
        (_supp_pnts_y[y_interval_idx + 1] - _supp_pnts_y[y_interval_idx]);

    return f_val_11 + f_val_21 + f_val_12 + f_val_22;
}

double PiecewiseBiLinearInterpolation::getBiDerivativeDx(double pnt_x_to_interpolate, double pnt_y_to_interpolate) const
{
	auto const& it_X(std::lower_bound(_supp_pnts_x.begin(), _supp_pnts_x.end(),
		pnt_x_to_interpolate));
	std::size_t const x_interval_idx =
		std::distance(_supp_pnts_x.begin(), it_X) - 1;
	auto const& it_Y(std::lower_bound(_supp_pnts_y.begin(), _supp_pnts_y.end(),
		pnt_y_to_interpolate));
	std::size_t const y_interval_idx =
		std::distance(_supp_pnts_y.begin(), it_Y) - 1;
	const double x_supp_pnt_size = _supp_pnts_x.size();
	const double a00 = _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * y_interval_idx];
	const double a10= _values_at_supp_pnts[x_interval_idx + 1 + x_supp_pnt_size * y_interval_idx] 
		- _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * y_interval_idx];
	const double a01 = _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * (y_interval_idx + 1)]
		- _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * y_interval_idx];
	const double a11 = _values_at_supp_pnts[x_interval_idx + 1 + x_supp_pnt_size * (y_interval_idx + 1)]
		+ _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * y_interval_idx]
		- (_values_at_supp_pnts[x_interval_idx + 1 + x_supp_pnt_size * y_interval_idx]
			+ _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * (y_interval_idx + 1)]);
	return a10 + a11 *pnt_y_to_interpolate;

}
double PiecewiseBiLinearInterpolation::getBiDerivativeDy(double pnt_x_to_interpolate, double pnt_y_to_interpolate) const
{
	auto const& it_X(std::lower_bound(_supp_pnts_x.begin(), _supp_pnts_x.end(),
		pnt_x_to_interpolate));
	std::size_t const x_interval_idx =
		std::distance(_supp_pnts_x.begin(), it_X) - 1;
	auto const& it_Y(std::lower_bound(_supp_pnts_y.begin(), _supp_pnts_y.end(),
		pnt_y_to_interpolate));
	std::size_t const y_interval_idx =
		std::distance(_supp_pnts_y.begin(), it_Y) - 1;
	const double x_supp_pnt_size = _supp_pnts_x.size();
	const double a00 = _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * y_interval_idx];
	const double a10 = _values_at_supp_pnts[x_interval_idx + 1 + x_supp_pnt_size * y_interval_idx]
		- _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * y_interval_idx];
	const double a01 = _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * (y_interval_idx + 1)]
		- _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * y_interval_idx];
	const double a11 = _values_at_supp_pnts[x_interval_idx + 1 + x_supp_pnt_size * (y_interval_idx + 1)]
		+ _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * y_interval_idx]
		- (_values_at_supp_pnts[x_interval_idx + 1 + x_supp_pnt_size * y_interval_idx]
			+ _values_at_supp_pnts[x_interval_idx + x_supp_pnt_size * (y_interval_idx + 1)]);
	return a01 + a11 *pnt_x_to_interpolate;

}


double PiecewiseBiLinearInterpolation::getSupportMaxX() const
{
    assert(!_supp_pnts_x.empty());
    return _supp_pnts_x.front();
}
double PiecewiseBiLinearInterpolation::getSupportMinX() const
{
    assert(!_supp_pnts_x.empty());
    return _supp_pnts_x.back();
}
double PiecewiseBiLinearInterpolation::getSupportMaxY() const
{
    assert(!_supp_pnts_y.empty());
    return _supp_pnts_y.front();
}
double PiecewiseBiLinearInterpolation::getSupportMinY() const
{
    assert(!_supp_pnts_y.empty());
    return _supp_pnts_y.back();
}
}  // end MathLib
