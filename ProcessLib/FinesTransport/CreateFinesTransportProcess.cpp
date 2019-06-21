/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateFinesTransportProcess.h"

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/PorousMedium/CreatePorousMediaProperties.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "FinesTransportProcess.h"
#include "FinesTransportMaterialProperties.h"
#include "FinesTransportLocalAssemblerInterface.h"

namespace ProcessLib
{
namespace FinesTransport
{
std::unique_ptr<Process> createFinesTransportProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::string const& output_directory,
    std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "FINESTRANSPORT");

    DBUG("Create FinesTransportProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__FinesTransport__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__FinesTransport__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__FinesTransport__process_variables__temperature}
             "nonwetting_pressure",
             //! \ogs_file_param_special{prj__processes__process__FinesTransport__process_variables__pressure}
             "saturation",
             //! \ogs_file_param_special{prj__processes__process__FinesTransport__process_variables__pressure}
             "salt_concentration" });
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        OGS_FATAL(
            "NOT IMPLEMENTED YET ");
    }
    // Process IDs, which are set according to the appearance order of the
    // process variables.
    const int _heat_transport_process_id = 0;
    const int _hydraulic_process_id = 1;
    MaterialLib::PorousMedium::PorousMediaProperties porous_media_properties{
        MaterialLib::PorousMedium::createPorousMediaProperties(
            mesh, config, parameters) };

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__FinesTransport__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    if (b.size() < mesh.getDimension())
    {
        OGS_FATAL(
            "specific body force (gravity vector) has %d components, mesh "
            "dimension is %d",
            b.size(), mesh.getDimension());
    }
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        specific_body_force.resize(b.size());
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    std::unique_ptr<FinesTransportMaterialProperties> material_properties =
        std::make_unique<FinesTransportMaterialProperties>(
            std::move(media_map),
            std::move(porous_media_properties),
            specific_body_force,
            has_gravity);

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"FinesTransport_twophase_salt"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<FinesTransportProcess>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(material_properties),
        std::move(secondary_variables), std::move(named_function_caller),
        use_monolithic_scheme,
        _heat_transport_process_id, _hydraulic_process_id);
}

}  // namespace FinesTransport
}  // namespace ProcessLib
