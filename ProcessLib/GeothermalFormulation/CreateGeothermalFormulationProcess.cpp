/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "CreateGeothermalFormulationProcess.h"
#include <cassert>

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "CreateGeothermalFormulationMaterialProperties.h"
#include "GeothermalFormulationMaterialProperties.h"
#include "GeothermalFormulationProcess.h"
#include "GeothermalFormulationProcessData.h"
namespace ProcessLib
{
namespace GeothermalFormulation
{
std::unique_ptr<Process> createGeothermalFormulationProcess(
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
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "Geothermal_Formulation");

    DBUG("Create GeothermalFormulation with pressure-enthalpy model.");
    //! \ogs_file_param{prj__processes__process__Geothermal_Formulation__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__Geothermal_Formulation__process_variables__gas_pressure}
         "gas_pressure",
         //! \ogs_file_param_special{prj__processes__process__Geothermal_Formulation__process_variables__overall_enthalpy}
         "overall_enthalpy"});

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"GeothermalFormulation_pressureenthalpy"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);
    // Parameter for the density of the solid.
    auto& density_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__Geothermal_Formulation__density_solid}
        "density_solid", parameters, 1);
    DBUG("Use \'%s\' as density_solid parameter.", density_solid.name.c_str());

    // Parameter for the specific heat capacity of the solid.
    auto& specific_heat_capacity_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__Geothermal_Formulation__specific_heat_capacity_solid}
        "specific_heat_capacity_solid", parameters, 1);
    DBUG("Use \'%s\' as specific_heat_capacity_solid parameter.",
        specific_heat_capacity_solid.name.c_str());

    // Parameter for the thermal conductivity of the solid (only one scalar per
    // element, i.e., the isotropic case is handled at the moment)
    auto& thermal_conductivity_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__Geothermal_Formulation__thermal_conductivity_solid}
        "thermal_conductivity_solid", parameters, 1);
    DBUG("Use \'%s\' as thermal_conductivity_solid parameter.",
        thermal_conductivity_solid.name.c_str());

    // Parameter for the thermal conductivity of the fluid.
    auto& thermal_conductivity_fluid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__Geothermal_Formulation__thermal_conductivity_fluid}
        "thermal_conductivity_fluid", parameters, 1);
    DBUG("Use \'%s\' as thermal_conductivity_fluid parameter.",
        thermal_conductivity_fluid.name.c_str());

    // Specific body force
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__Geothermal_Formulation__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    Eigen::VectorXd specific_body_force(b.size());
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
        std::copy_n(b.data(), b.size(), specific_body_force.data());

    //! \ogs_file_param{prj__processes__process__Geothermal_Formulation__mass_lumping}
    auto const mass_lumping = config.getConfigParameter<bool>("mass_lumping");

    auto& temperature = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__Geothermal_Formulation__temperature}
        "temperature", parameters, 1);

    //! \ogs_file_param{prj__processes__process__Geothermal_Formulation__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    boost::optional<MeshLib::PropertyVector<int> const&> material_ids;
    if (mesh.getProperties().existsPropertyVector<int>("MaterialIDs"))
    {
        INFO("The twophase flow is in heterogeneous porous media.");
        material_ids =
            *mesh.getProperties().getPropertyVector<int>("MaterialIDs");
    }
    else
    {
        INFO("The twophase flow is in homogeneous porous media.");
    }
    std::unique_ptr<GeothermalFormulationMaterialProperties>
    material = createGeothermalFormulationMaterialProperties(mat_config, material_ids);

    GeothermalFormulationProcessData process_data{
        specific_body_force, has_gravity, mass_lumping, temperature, std::move(material),
    density_solid, specific_heat_capacity_solid,thermal_conductivity_solid,thermal_conductivity_fluid};

    return std::make_unique<GeothermalFormulationProcess>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller),
        mat_config, curves);
}

}  // end of namespace
}  // end of namespace
