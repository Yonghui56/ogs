/**
 * @file TestPiecewiseBilinearInterpolation.cpp
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

#include "MathLib/InterpolationAlgorithms/PiecewiseBilinearInterpolation.h"

TEST(MathLibInterpolationAlgorithms, PiecewiseBilinearInterpolation)
{
    std::vector<double> x_supp_pnts = {2, 3, 4, 5, 6};
    std::vector<double> y_supp_pnts = {20, 30, 40, 50, 60};
    std::vector<double> value_supp_pnts = {1,  2,  3,  4,  5,  6,  7,  8,  9,
                                           10, 11, 12, 13, 14, 15, 16, 17, 18,
                                           19, 20, 21, 22, 23, 24, 25};
    MathLib::PiecewiseBilinearInterpolation interpolation{
        std::move(x_supp_pnts), std::move(y_supp_pnts),
        std::move(value_supp_pnts)};
    // Interpolation

    ASSERT_NEAR(19.25, interpolation.getValue(2.5, 55.5),
                std::numeric_limits<double>::epsilon());
}

TEST(MathLibInterpolationAlgorithms, PiecewiseBilinearInterpolationDerivative)
{
    std::vector<double> x_supp_pnts = {2, 3, 4, 5, 6};
    std::vector<double> y_supp_pnts = {20, 30, 40, 50, 60};
    std::vector<double> value_supp_pnts = {1,  2,  3,  4,  5,  6,  7,  8,  9,
                                           10, 11, 12, 13, 14, 15, 16, 17, 18,
                                           19, 20, 21, 22, 23, 24, 25};

    MathLib::PiecewiseBilinearInterpolation interpolation{
        std::move(x_supp_pnts), std::move(y_supp_pnts),
        std::move(value_supp_pnts)};
    // Interpolation

    ASSERT_NEAR(1,
                interpolation.getDerivativeDx(2.5, 55.5),
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(5,
                interpolation.getDerivativeDy(2.5, 55.5),
                std::numeric_limits<double>::epsilon());
}
