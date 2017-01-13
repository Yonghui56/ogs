/**
 * @file TestPiecewiseLinearInterpolation.cpp
 * @author Thomas Fischer
 * @date Feb 12, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// stl
#include <limits>

// google test
#include "gtest/gtest.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseBiLinearInterpolation.h"

TEST(MathLibInterpolationAlgorithms, PiecewiseBiLinearInterpolation)
{
    std::vector<double> x_supp_pnts = {2, 3, 4, 5, 6};
    std::vector<double> y_supp_pnts = {20, 30, 40, 50, 60};
    std::vector<double> value_supp_pnts = {1,  2,  3,  4,  5,  6,  7,  8,  9,
                                           10, 11, 12, 13, 14, 15, 16, 17, 18,
                                           19, 20, 21, 22, 23, 24, 25};
    MathLib::PiecewiseBiLinearInterpolation interpolation{
        std::move(x_supp_pnts), std::move(y_supp_pnts),
        std::move(value_supp_pnts)};
    // Interpolation

    ASSERT_NEAR(19.25, interpolation.getBiValue(2.5, 55.5),
                std::numeric_limits<double>::epsilon());
}

TEST(MathLibInterpolationAlgorithms, PiecewiseBiLinearInterpolationDerivative)
{
    std::vector<double> x_supp_pnts = {2, 3, 4, 5, 6};
    std::vector<double> y_supp_pnts = {20, 30, 40, 50, 60};
    std::vector<double> value_supp_pnts = {1,  2,  3,  4,  5,  6,  7,  8,  9,
                                           10, 11, 12, 13, 14, 15, 16, 17, 18,
                                           19, 20, 21, 22, 23, 24, 25};

    MathLib::PiecewiseBiLinearInterpolation interpolation{
        std::move(x_supp_pnts), std::move(y_supp_pnts),
        std::move(value_supp_pnts)};
    // Interpolation

    ASSERT_NEAR(1,
                interpolation.getBiDerivativeDx(2.5, 55.5),
                std::numeric_limits<double>::epsilon());
}
