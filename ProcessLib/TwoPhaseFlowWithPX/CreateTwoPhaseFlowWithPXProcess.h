/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_CREATETWOPHASEFLOWWITHPXPROCESS_H
#define OGS_CREATETWOPHASEFLOWWITHPXPROCESS_H

#include <memory>
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPX
{
std::unique_ptr<Process> CreateTwoPhaseFlowWithPXProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);
}  // end of namespace
}  // end of namespace

#endif /* CREATETWOPHASEFLOWWITHPXPROCESS_H */
