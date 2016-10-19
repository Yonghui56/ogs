/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TwoPhaseFlowWithPPMaterialProperties.cpp
 *
 * Created on August 18, 2016, 11:49 AM
 */

#include "TwoPhaseFlowWithPPMaterialProperties.h"

#include <logog/include/logog.hpp>

//#include "BaseLib/reorderVector.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
TwoPhaseFlowWithPPMaterialProperties::TwoPhaseFlowWithPPMaterialProperties(
    BaseLib::ConfigTree const& config,
    MeshLib::PropertyVector<int> const& material_ids)
    : _material_ids(material_ids)
{
    DBUG("Reading material properties of liquid flow process.");

    //! \ogs_file_param{prj__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    // Get fluid properties
    //! \ogs_file_param{prj__material_property__fluid__density}
    auto const& rho_conf = fluid_config.getConfigSubtree("density");
    _liquid_density = MaterialLib::Fluid::createFluidDensityModel(rho_conf);
	auto const& rho_gas_conf = fluid_config.getConfigSubtree("gasdensity");
	_gas_density = MaterialLib::Fluid::createFluidDensityModel(rho_gas_conf);
    //! \ogs_file_param{prj__material_property__fluid__viscosity}
    auto const& mu_conf = fluid_config.getConfigSubtree("viscosity");
    _viscosity = MaterialLib::Fluid::createViscosityModel(mu_conf);
	//! \ogs_file_param{prj__material_property__fluid__gas__viscosity}
	auto const& mu_gas_conf = fluid_config.getConfigSubtree("gasviscosity");
	_gas_viscosity = MaterialLib::Fluid::createViscosityModel(mu_gas_conf);
    // Get porous properties
    std::vector<int> mat_ids;
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
    }

    //BaseLib::reorderVector(_intrinsic_permeability_models, mat_ids);
    //BaseLib::reorderVector(_porosity_models, mat_ids);
    //BaseLib::reorderVector(_storage_models, mat_ids);
}

void TwoPhaseFlowWithPPMaterialProperties::setMaterialID(const SpatialPosition& pos)
{
    if (_material_ids.empty())
    {
        assert(pos.getElementID().get() < _material_ids.size());
        _current_material_id = _material_ids[pos.getElementID().get()];
    }
}

double TwoPhaseFlowWithPPMaterialProperties::getLiquidDensity(const double p,
                                                      const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pl)] = p;
    return _liquid_density->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getGasDensity(const double p,
	const double T) const
{
	ArrayType vars;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pl)] = p;
	return _gas_density->getValue(vars);
}
double TwoPhaseFlowWithPPMaterialProperties::getViscosity(const double p,
                                                  const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pl)] = p;
    return _viscosity->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getGasViscosity(const double p,
	const double T) const
{
	ArrayType vars;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pl)] = p;
	return _viscosity->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getMassCoefficient(
    const double /*t*/, const SpatialPosition& /*pos*/, const double p,
    const double T, const double porosity_variable,
    const double storage_variable) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pl)] = p;
    const double drho_dp = _liquid_density->getdValue(
        vars, MaterialLib::Fluid::PropertyVariableType::pl);
    const double rho = _liquid_density->getValue(vars);

    const double porosity =
        _porosity_models[_current_material_id]->getValue(porosity_variable, T);
    const double storage =
        _storage_models[_current_material_id]->getValue(storage_variable);
    return porosity * drho_dp / rho + storage;
}

Eigen::MatrixXd const& TwoPhaseFlowWithPPMaterialProperties::getPermeability(
    const double /*t*/, const SpatialPosition& /*pos*/, const int /*dim*/) const
{
    return _intrinsic_permeability_models[_current_material_id];
}

}  // end of namespace
}  // end of namespace
