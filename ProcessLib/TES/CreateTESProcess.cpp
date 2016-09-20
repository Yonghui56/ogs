/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTESProcess.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "TESProcess.h"

namespace ProcessLib
{
namespace TES
{
std::unique_ptr<Process> createTESProcess(
    MeshLib::Mesh& mesh,
    Process::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<Process::TimeDiscretization>&& time_discretization,
    std::unique_ptr<NumLib::ConvergenceCriterion>&& convergence_criterion,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "TES");

    DBUG("Create TESProcess.");

    auto process_variables = findProcessVariables(
        variables, config,
        {"fluid_pressure", "temperature", "vapour_mass_fraction"});

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"TES_pressure", "TES_temperature", "TES_vapour_mass_fraction"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    //! \ogs_file_param{process__output}
    ProcessOutput process_output{config.getConfigSubtree("output")};

    return std::unique_ptr<Process>{new TESProcess{
        mesh, nonlinear_solver, std::move(time_discretization),
        std::move(convergence_criterion), parameters,
        std::move(process_variables), std::move(secondary_variables),
        std::move(process_output), std::move(named_function_caller), config}};
}

}  // namespace TES
}  // namespace ProcessLib
