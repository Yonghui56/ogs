/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_CREATE_TWOPHASECOMPONENTIALPROCESS_H_
#define PROCESS_LIB_CREATE_TWOPHASECOMPONENTIALPROCESS_H_

#include <memory>
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace TwoPhaseComponential
{
std::unique_ptr<Process> createTwoPhaseComponentialProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config,
	std::map<std::string,
	std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
	curves);

}  // namespace TwoPhaseComponential
}  // namespace ProcessLib

#endif  // PROCESS_LIB_CREATE_TWOPHASECOMPONENTIALPROCESS_H_
