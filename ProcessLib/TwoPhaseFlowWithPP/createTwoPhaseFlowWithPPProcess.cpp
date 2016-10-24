/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   createTwoPhaseFlowWithPPProcess.cpp
 *
 * Created on August 19, 2016, 1:30 PM
 */
#include <cassert>
#include "createTwoPhaseFlowWithPPProcess.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "TwoPhaseFlowWithPPProcess.h"
#include "TwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
std::unique_ptr<Process> createTwoPhaseFlowWithPPProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
	std::map<std::string,
	std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
	curves)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "TWOPHASE_FLOW_PP");

    DBUG("Create TwoPhaseFlowProcess with PP model.");

    // Process variable.
    auto process_variables = findProcessVariables(
        variables, config,
        {//! \ogs_file_param_special{process__TWOPHASE_FLOW__WITHPP__process_variables__process_variable}
			"capillary_pressure", "gas_pressure" });
         //"process_variable"});

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller({"TwoPhaseFlow_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    auto const gravity_param = config.getConfigParameter("gravitational_term");
    auto const axis =
        //! \ogs_file_attr{process__TWOPHASE_FLOW__gravitational_term__axis}
        gravity_param.getConfigAttributeOptional<std::string>("axis");
    // Gravitational acceleration
    auto const g =
        //! \ogs_file_attr{process__TWOPHASE_FLOW__gravitational_term__g}
        gravity_param.getConfigAttributeOptional<double>("g");

    int gravity_axis_id = -1;
    if (*axis == "x")
        gravity_axis_id = 0;
    else if (*axis == "y")
        gravity_axis_id = 1;
    else if (*axis == "z")
        gravity_axis_id = 2;

	// Specific body force parameter.
	auto& specific_body_force = findParameter<double>(
		config,
		//! \ogs_file_param_special{process__TWOPHASE_FLOW_specific_body_force}
		"specific_body_force", parameters, mesh.getDimension());
	DBUG("Use \'%s\' as specific body force parameter.",
		specific_body_force.name.c_str());

	// Assume constant parameter, then check the norm at arbitrary
	// SpatialPosition and time.
	assert(dynamic_cast<ConstantParameter<double>*>(&specific_body_force));
	bool const has_gravity =
		MathLib::toVector(specific_body_force(0, SpatialPosition{})).norm() > 0;

	// has mass lumping
	auto mass_lump = config.getConfigParameter<bool>("mass_lumping");

	TwoPhaseFlowWithPPProcessData process_data{ 
		specific_body_force,
		has_gravity,
		mass_lump};
    //! \ogs_file_param{process__TWOPHASE_FLOW__material_property}

    auto const& mat_config = config.getConfigSubtree("material_property");

    auto const& mat_ids =
        mesh.getProperties().getPropertyVector<int>("MaterialIDs");
	auto is_heterogeneous = config.getConfigParameter<bool>("heterogeneous");
    if (is_heterogeneous)//mat_ids
    {
        INFO("The twophase flow is in heterogeneous porous media.");
        return std::unique_ptr<Process>{new TwoPhaseFlowWithPPProcess{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(process_data), std::move(secondary_variables),
            std::move(named_function_caller), *mat_ids, gravity_axis_id, *g,
            mat_config,curves}};
    }
    else
    {
        INFO("The twophase flow is in homogeneous porous media.");

        MeshLib::Properties dummy_property;
        auto const& dummy_property_vector =
            dummy_property.createNewPropertyVector<int>(
                "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
        return std::unique_ptr<Process>{new TwoPhaseFlowWithPPProcess{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(process_data), std::move(secondary_variables),
            std::move(named_function_caller), *dummy_property_vector,
            gravity_axis_id, *g, mat_config, curves}};
    }
}

}  // end of namespace
}  // end of namespace
