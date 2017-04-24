/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermalTwoPhaseFlowWithPPMaterialProperties.h"
#include <logog/include/logog.hpp>
#include "BaseLib/reorderVector.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/SpecificHeatCapacity/CreateSpecificFluidHeatCapacityModel.h"
#include "MaterialLib/Fluid/ThermalConductivity/CreateFluidThermalConductivityModel.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MaterialLib/TwoPhaseModels/CreateTwoPhaseFlowMaterialProperties.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"
#include "ThermalTwoPhaseFlowWithPPMaterialProperties.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
std::unique_ptr<ThermalTwoPhaseFlowWithPPMaterialProperties>
createThermalTwoPhaseFlowWithPPMaterialProperties(
    BaseLib::ConfigTree const& config,
    boost::optional<MeshLib::PropertyVector<int> const&> material_ids)
{
    DBUG("Reading material properties of nonisothermal two-phase flow process.");

    auto two_phase_material_model
        = MaterialLib::TwoPhaseFlowWithPP::createTwoPhaseFlowMaterialProperties(config, material_ids);
    //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    // Get fluid properties
    //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__liquid_density}
    auto const& rho_conf = fluid_config.getConfigSubtree("liquid_density");
    auto liquid_density =
        MaterialLib::Fluid::createFluidDensityModel(rho_conf);
    //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__gas_density}
    auto const& rho_gas_conf = fluid_config.getConfigSubtree("gas_density");
    auto gas_density =
        MaterialLib::Fluid::createFluidDensityModel(rho_gas_conf);

    auto const& spec_heat_capacity_solid_conf =
        //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__specific_heat_capacity_solid}
        fluid_config.getConfigSubtree("specific_heat_capacity_solid");
    auto specific_heat_capacity_solid =
        MaterialLib::Fluid::createSpecificFluidHeatCapacityModel(spec_heat_capacity_solid_conf);
    auto const& spec_heat_capacity_water_conf =
        //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__specific_heat_capacity_water}
        fluid_config.getConfigSubtree("specific_heat_capacity_water");
    auto specific_heat_capacity_water =
        MaterialLib::Fluid::createSpecificFluidHeatCapacityModel(spec_heat_capacity_water_conf);
    auto const& spec_heat_capacity_air_conf =
        //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__specific_heat_capacity_air}
        fluid_config.getConfigSubtree("specific_heat_capacity_air");
    auto specific_heat_capacity_air =
        MaterialLib::Fluid::createSpecificFluidHeatCapacityModel(spec_heat_capacity_air_conf);
    auto const& spec_heat_capacity_vapor_conf =
        //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__specific_heat_capacity_vapor}
        fluid_config.getConfigSubtree("specific_heat_capacity_vapor");
    auto specific_heat_capacity_vapor =
        MaterialLib::Fluid::createSpecificFluidHeatCapacityModel(spec_heat_capacity_vapor_conf);

    auto const& thermal_conductivity_dry_solid_conf =
        //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__thermal_conductivity_dry_solid}
        fluid_config.getConfigSubtree("thermal_conductivity_dry_solid");
    auto thermal_conductivity_dry_solid =
        MaterialLib::Fluid::createFluidThermalConductivityModel(thermal_conductivity_dry_solid_conf);
    auto const& thermal_conductivity_wet_solid_conf =
        //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__thermal_conductivity_wet_solid}
        fluid_config.getConfigSubtree("thermal_conductivity_wet_solid");
    auto thermal_conductivity_wet_solid =
        MaterialLib::Fluid::createFluidThermalConductivityModel(thermal_conductivity_wet_solid_conf);

    //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__liquid_viscosity}
    auto const& mu_conf = fluid_config.getConfigSubtree("liquid_viscosity");
    auto viscosity = MaterialLib::Fluid::createViscosityModel(mu_conf);
    //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__gas_viscosity}
    auto const& mu_gas_conf = fluid_config.getConfigSubtree("gas_viscosity");
    auto gas_viscosity = MaterialLib::Fluid::createViscosityModel(mu_gas_conf);

    // Get porous properties
    std::vector<int> mat_ids;
    std::vector<int> mat_krel_ids;
    std::vector<Eigen::MatrixXd> intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        storage_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
        capillary_pressure_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>
        relative_permeability_models;

    //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__porous_medium}
    auto const& poro_config = config.getConfigSubtree("porous_medium");
    //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium}
    for (auto const& conf : poro_config.getConfigSubtreeList("porous_medium"))
    {
        //! \ogs_file_attr{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__id}
        auto const id = conf.getConfigAttributeOptional<int>("id");
        mat_ids.push_back(*id);

        //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__permeability}
        auto const& permeability_conf = conf.getConfigSubtree("permeability");
        intrinsic_permeability_models.emplace_back(
            MaterialLib::PorousMedium::createPermeabilityModel(
                permeability_conf));

        //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__porosity}
        auto const& porosity_conf = conf.getConfigSubtree("porosity");
        auto n = MaterialLib::PorousMedium::createPorosityModel(porosity_conf);
        porosity_models.emplace_back(std::move(n));

        //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__storage}
        auto const& storage_conf = conf.getConfigSubtree("storage");
        auto beta = MaterialLib::PorousMedium::createStorageModel(storage_conf);
        storage_models.emplace_back(std::move(beta));

        auto const& capillary_pressure_conf =
            //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__capillary_pressure}
            conf.getConfigSubtree("capillary_pressure");
        auto pc = MaterialLib::PorousMedium::createCapillaryPressureModel(
            capillary_pressure_conf);
        capillary_pressure_models.emplace_back(std::move(pc));

        auto const& krel_config =
            //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__relative_permeability}
            conf.getConfigSubtree("relative_permeability");
        for (
            auto const& krel_conf :
            //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__relative_permeability__relative_permeability}
            krel_config.getConfigSubtreeList("relative_permeability"))
        {
            auto const krel_id =
                //! \ogs_file_param{prj__processes__process__THERMAL_TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__relative_permeability__relative_permeability__id}
                krel_conf.getConfigAttributeOptional<int>("id");
            mat_krel_ids.push_back(*krel_id);
            auto krel_n =
                MaterialLib::PorousMedium::createRelativePermeabilityModel(
                    krel_conf);
            relative_permeability_models.emplace_back(std::move(krel_n));
        }
        BaseLib::reorderVector(relative_permeability_models, mat_krel_ids);
    }

    BaseLib::reorderVector(intrinsic_permeability_models, mat_ids);
    BaseLib::reorderVector(porosity_models, mat_ids);
    BaseLib::reorderVector(storage_models, mat_ids);

    return std::unique_ptr<ThermalTwoPhaseFlowWithPPMaterialProperties>{
        new ThermalTwoPhaseFlowWithPPMaterialProperties{
            std::move(two_phase_material_model),
            std::move(specific_heat_capacity_solid),
            std::move(specific_heat_capacity_water),std::move(specific_heat_capacity_air),
            std::move(specific_heat_capacity_vapor),std::move(thermal_conductivity_dry_solid),
            std::move(thermal_conductivity_wet_solid),
            std::move(relative_permeability_models)}};
}

}  // end namespace
}  // end namespace
