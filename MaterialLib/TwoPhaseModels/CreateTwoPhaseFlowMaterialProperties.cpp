/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   createPorosityModel.cpp
 *
 * Created on August 16, 2016, 1:16 PM
 */

#include <logog/include/logog.hpp>

#include "BaseLib/reorderVector.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"

#include "CreateTwoPhaseFlowMaterialProperties.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"
#include "TwoPhaseFlowWithPPMaterialProperties.h"

namespace MaterialLib
{
namespace TwoPhaseFlowWithPP
{
std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties>
CreateTwoPhaseFlowMaterialProperties(
    BaseLib::ConfigTree const& config,
    bool const has_material_ids,
    MeshLib::PropertyVector<int> const& material_ids,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves_)
{
    DBUG("Reading material properties of two-phase flow process.");

    //! \ogs_file_param{prj__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    // Get fluid properties
    //! \ogs_file_param{prj__material_property__fluid__density}
    auto const& rho_conf = fluid_config.getConfigSubtree("liquiddensity");
    auto _liquid_density =
        MaterialLib::Fluid::createFluidDensityModel(rho_conf);
    auto const& rho_gas_conf = fluid_config.getConfigSubtree("gasdensity");
    auto _gas_density =
        MaterialLib::Fluid::createFluidDensityModel(rho_gas_conf);
    //! \ogs_file_param{prj__material_property__fluid__viscosity}
    auto const& mu_conf = fluid_config.getConfigSubtree("liquidviscosity");
    auto _viscosity = MaterialLib::Fluid::createViscosityModel(mu_conf);
    //! \ogs_file_param{prj__material_property__fluid__gas__viscosity}
    auto const& mu_gas_conf = fluid_config.getConfigSubtree("gasviscosity");
    auto _gas_viscosity = MaterialLib::Fluid::createViscosityModel(mu_gas_conf);
    // Get porous properties
    int cap_pressure_model;
    int rel_wet_perm_model;
    int rel_nonwet_perm_model;
    std::vector<int> mat_ids;
    std::vector<Eigen::MatrixXd> _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _storage_models;
    std::array<double, 4> cap_pressure_value;
    std::array<double, 4> rel_wet_perm_value;
    std::array<double, 4> rel_nonwet_perm_value;
    //! \ogs_file_param{prj__material_property__porous_medium}
    auto const& poro_config = config.getConfigSubtree("porous_medium");
    //! \ogs_file_param{prj__material_property__porous_medium__porous_medium}
    for (auto const& conf : poro_config.getConfigSubtreeList("porous_medium"))
    {
        //! \ogs_file_attr{prj__material_property__porous_medium__porous_medium__id}
        auto const id = conf.getConfigAttributeOptional<int>("id");
        mat_ids.push_back(*id);

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__permeability}
        auto const& perm_conf = conf.getConfigSubtree("permeability");
        _intrinsic_permeability_models.emplace_back(
            MaterialLib::PorousMedium::createPermeabilityModel(perm_conf));

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__porosity}
        auto const& poro_conf = conf.getConfigSubtree("porosity");
        auto n = MaterialLib::PorousMedium::createPorosityModel(poro_conf);
        _porosity_models.emplace_back(std::move(n));

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__storage}
        auto const& stora_conf = conf.getConfigSubtree("storage");
        auto beta = MaterialLib::PorousMedium::createStorageModel(stora_conf);
        _storage_models.emplace_back(std::move(beta));

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__cap_pressure}
        auto const& cap_pressure_conf =
            conf.getConfigSubtree("capillary_pressure");
        auto const cap_pressure_type =
            cap_pressure_conf.getConfigParameter<std::string>("type");
        if (cap_pressure_type == "Curve")
            cap_pressure_model = 0;
        else if (cap_pressure_type == "van_Genuchten")
            cap_pressure_model = 1;
        else if (cap_pressure_type == "Brooks_Corey")
        {
            cap_pressure_model = 2;
            cap_pressure_value = {
                {//! \ogs_file_param{material__fluid__density__linear_temperature__rho0}
                 cap_pressure_conf.getConfigParameter<double>("entry_pressure"),
                 //! \ogs_file_param{material__fluid__density__linear_temperature__temperature0}
                 cap_pressure_conf.getConfigParameter<double>(
                     "res_saturation_wet"),
                 //! \ogs_file_param{material__fluid__density__linear_temperature__beta}
                 cap_pressure_conf.getConfigParameter<double>(
                     "res_saturation_nonwet"),
                 //! \ogs_file_param{material__fluid__density__linear_temperature__beta}
                 cap_pressure_conf.getConfigParameter<double>("lambda")}};
        }
        else if (cap_pressure_type == "Liakopoulos")
            cap_pressure_model = 10;
        else
            OGS_FATAL("This model has not been implemented yet");
        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__rel_permeability_liquid}
        auto const& rel_wet_perm_conf =
            conf.getConfigSubtree("relative_wet_permeability");
        auto const rel_wet_perm_type =
            rel_wet_perm_conf.getConfigParameter<std::string>("type");
        if (rel_wet_perm_type == "Curve")
            rel_wet_perm_model = 0;
        else if (rel_wet_perm_type == "van_Genuchten")
            rel_wet_perm_model = 1;
        else if (rel_wet_perm_type == "Brooks_Corey")
        {
            rel_wet_perm_model = 2;
            rel_wet_perm_value = {
                {//! \ogs_file_param{material__fluid__density__linear_temperature__rho0}
                 rel_wet_perm_conf.getConfigParameter<double>(
                     "res_saturation_wet"),
                 //! \ogs_file_param{material__fluid__density__linear_temperature__temperature0}
                 rel_wet_perm_conf.getConfigParameter<double>(
                     "res_saturation_nonwet"),
                 //! \ogs_file_param{material__fluid__density__linear_temperature__beta}
                 rel_wet_perm_conf.getConfigParameter<double>("lambda"),
                 //! \ogs_file_param{material__fluid__density__linear_temperature__beta}
                 rel_wet_perm_conf.getConfigParameter<double>(
                     "minimum_value")}};
        }
        else if (rel_wet_perm_type == "Liakopoulos")
            rel_wet_perm_model = 10;
        else
            OGS_FATAL("This model has not been implemented yet");

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__rel_permeability_gas}
        auto const& rel_nonwet_perm_conf =
            conf.getConfigSubtree("relative_nonwet_permeability");
        auto const rel_nonwet_perm_type =
            rel_nonwet_perm_conf.getConfigParameter<std::string>("type");
        if (rel_nonwet_perm_type == "Curve")
            rel_nonwet_perm_model = 0;
        else if (rel_nonwet_perm_type == "van_Genuchten")
            rel_nonwet_perm_model = 1;
        else if (rel_nonwet_perm_type == "Brooks_Corey")
        {
            rel_nonwet_perm_model = 2;
            rel_nonwet_perm_value = {
                {//! \ogs_file_param{material__fluid__density__linear_temperature__rho0}
                 rel_nonwet_perm_conf.getConfigParameter<double>(
                     "res_saturation_wet"),
                 //! \ogs_file_param{material__fluid__density__linear_temperature__temperature0}
                 rel_nonwet_perm_conf.getConfigParameter<double>(
                     "res_saturation_nonwet"),
                 //! \ogs_file_param{material__fluid__density__linear_temperature__beta}
                 rel_nonwet_perm_conf.getConfigParameter<double>("lambda"),
                 //! \ogs_file_param{material__fluid__density__linear_temperature__beta}
                 rel_nonwet_perm_conf.getConfigParameter<double>(
                     "minimum_value")}};
        }
        else
            OGS_FATAL("This model has not been implemented yet");
    }

    BaseLib::reorderVector(_intrinsic_permeability_models, mat_ids);
    BaseLib::reorderVector(_porosity_models, mat_ids);
    BaseLib::reorderVector(_storage_models, mat_ids);

    return std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties>{
        new TwoPhaseFlowWithPPMaterialProperties{
            has_material_ids, material_ids, std::move(_liquid_density),
            std::move(_viscosity), std::move(_gas_density),
            std::move(_gas_viscosity), _intrinsic_permeability_models,
            std::move(_porosity_models), std::move(_storage_models),
            cap_pressure_model, rel_wet_perm_model, rel_nonwet_perm_model,
            cap_pressure_value, rel_wet_perm_value, rel_nonwet_perm_value,
            curves_}};
}

}  // end namespace
}  // end namespace
