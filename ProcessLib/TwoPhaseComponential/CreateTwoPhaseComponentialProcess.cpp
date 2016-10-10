/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTwoPhaseComponentialProcess.h"

#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "TwoPhaseComponentialProcess.h"
#include "TwoPhaseComponentialProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseComponential
{
std::unique_ptr<Process> createTwoPhaseComponentialProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
	unsigned const integration_order,
    BaseLib::ConfigTree const& config,
	std::map<std::string,
	std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
	curves)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "TWOPHASE_COMPONENTIAL");

    DBUG("Create TwoPhaseComponentialProcess.");

    // Process variable.
    auto process_variables = findProcessVariables(
        variables, config,
        {"gas_pressure", "molar_fraction_h", "molar_fraction_ch4", "molar_fraction_co2", "capillary_pressure"});

    // thermal conductivity parameter.
    auto& porosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HEAT_CONDUCTION__thermal_conductivity}
        "porosity", parameters, 1);

    DBUG("Use \'%s\' as Porosity parameter.",
		porosity.name.c_str());

    // heat capacity parameter.
    auto& intrinsic_permeability = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HEAT_CONDUCTION__heat_capacity}
        "intrinsic_permeability", parameters, 1);

    DBUG("Use \'%s\' as Intrinsic Permeability parameter.", 
		intrinsic_permeability.name.c_str());

    // density parameter.
    auto& water_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HEAT_CONDUCTION__density}
        "water_density", parameters, 1);

    DBUG("Use \'%s\' as density parameter.", 
		water_density.name.c_str());

	auto& mat_id = findParameter<double>(
		config,
		//! \ogs_file_param_special{process__HEAT_CONDUCTION__density}
		"material_id", parameters, 1);

	DBUG("Use \'%s\' as material id.",
		mat_id.name.c_str());

	auto& viscosity_gas = findParameter<double>(
		config,
		//! \ogs_file_param_special{process__HEAT_CONDUCTION__density}
		"viscosity_gas", parameters, 1);

	DBUG("Use \'%s\' as the viscosity of gas.",
		mat_id.name.c_str());

	auto& viscosity_liquid = findParameter<double>(
		config,
		//! \ogs_file_param_special{process__HEAT_CONDUCTION__density}
		"viscosity_liquid", parameters, 1);

	DBUG("Use \'%s\' as the viscosity of gas.",
		mat_id.name.c_str());

	//has  gravity
	auto grav = config.getConfigParameter<bool>("g");

	TwoPhaseComponentialProcessData process_data{ porosity, intrinsic_permeability,
		water_density
		, mat_id
		, viscosity_gas
		, viscosity_liquid
	    , grav
	    ,curves};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"Liquid_Pressure","Gas_Saturation"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::unique_ptr<Process>{new TwoPhaseComponentialProcess{
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller)}};
}

}  // namespace TWOPHASECOMPONENTIAL
}  // namespace ProcessLib
