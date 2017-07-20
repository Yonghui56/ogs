/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateIterationNumberBasedAdaptiveTimeStepping.cpp
 *  Created on July 26, 2017, 4:43 PM
 */

#include "CreateIterationNumberBasedAdaptiveTimeStepping.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/makeVectorUnique.h"

#include "IterationNumberBasedAdaptiveTimeStepping.h"
#include "TimeStepAlgorithm.h"

namespace NumLib
{
class TimeStepAlgorithm;
std::unique_ptr<TimeStepAlgorithm> createIterationNumberBasedAdaptiveTimeStepping(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__time_stepping__type}
    config.checkConfigParameter("type", "IterationNumberBasedAdaptiveTimeStepping");

    //! \ogs_file_param{prj__time_loop__time_stepping__IterationNumberBasedAdaptiveTimeStepping__t_initial}
    auto const t0 = config.getConfigParameter<double>("t_initial");
    //! \ogs_file_param{prj__time_loop__time_stepping__IterationNumberBasedAdaptiveTimeStepping__t_end}
    auto const t_end = config.getConfigParameter<double>("t_end");

    //! \ogs_file_param{prj__time_loop__time_stepping__IterationNumberBasedAdaptiveTimeStepping__dt_min}
    auto const h_min = config.getConfigParameter<double>("dt_min");
    //! \ogs_file_param{prj__time_loop__time_stepping__IterationNumberBasedAdaptiveTimeStepping__dt_max}
    auto const h_max = config.getConfigParameter<double>("dt_max");
    //! \ogs_file_param{prj__time_loop__time_stepping__IterationNumberBasedAdaptiveTimeStepping__dt_initial_ts}
    auto const dt_initial_ts = config.getConfigParameter<double>("dt_initial_ts");

    auto iter_times_vector =
        //! \ogs_file_param{prj__time_loop__time_stepping__IterationNumberBasedAdaptiveTimeStepping__iter_times_vector}
        config.getConfigParameter<std::vector<std::size_t>>("iteration_number_vector",
                                                       std::vector<std::size_t>{});

    auto multiplier_vector =
        //! \ogs_file_param{prj__time_loop__time_stepping__IterationNumberBasedAdaptiveTimeStepping__multiplier_vector}
        config.getConfigParameter<std::vector<double>>("multiplier_coefficients",
            std::vector<double>{});

    return std::make_unique<IterationNumberBasedAdaptiveTimeStepping>(
        t0, t_end, h_min, h_max, dt_initial_ts, std::move(iter_times_vector), std::move(multiplier_vector));
}
}  // end of namespace NumLib
