/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateRichardsFlowProcess.h"

#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "RichardsFlowProcess.h"
#include "RichardsFlowProcessData.h"

namespace ProcessLib
{
namespace RichardsFlow
{
std::unique_ptr<Process> createRichardsFlowProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config,
	std::map<std::string,
	std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
	curves)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "RICHARDS_FLOW");

    DBUG("Create RichardsFlowProcess.");

    // Process variable.
    auto process_variables = findProcessVariables(
        variables, config,
        {//! \ogs_file_param_special{process__RICHARDS_FLOW__process_variables__process_variable}
         "process_variable"});

    // Hydraulic conductivity parameter.
    auto& intrinsic_permeability = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__RICHARDS_FLOW__hydraulic_conductivity}
        "intrinsic_permeability",
        parameters, 1);

    DBUG("Use \'%s\' as intrinsic permeability parameter.",
		intrinsic_permeability.name.c_str());

	// Porosity parameter.
	auto& porosity = findParameter<double>(
		config,
		//! \ogs_file_param_special{process__RICHARDS_FLOW__hydraulic_conductivity}
		"porosity",
		parameters, 1);

	DBUG("Use \'%s\' as porosity parameter.",
		porosity.name.c_str());

	// Viscosity parameter.
	auto& viscosity = findParameter<double>(
		config,
		//! \ogs_file_param_special{process__RICHARDS_FLOW__hydraulic_conductivity}
		"viscosity",
		parameters, 1);

	DBUG("Use \'%s\' as porosity parameter.",
		viscosity.name.c_str());

	//has  gravity
	auto grav = config.getConfigParameter<bool>("g");
	auto mass_lump = config.getConfigParameter<bool>("mass_lumping");

	RichardsFlowProcessData process_data{
		intrinsic_permeability,
		porosity,
		viscosity,
		grav,
		mass_lump,
		curves
	};
    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"RichardsFlow_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::unique_ptr<Process>{new RichardsFlowProcess{
        mesh, std::move(jacobian_assembler), parameters,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller)}};
}

}  // namespace RichardsFlow
}  // namespace ProcessLib
