/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#include "TwoPhaseComponentialFlowMaterialProperties.h"
#include <logog/include/logog.hpp>
#include "BaseLib/reorderVector.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "NumLib/NewtonRaphson.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"

using MaterialLib::PhysicalConstant::MolarMass::H2;
using MaterialLib::PhysicalConstant::IdealGasConstant;
using MaterialLib::PhysicalConstant::HenryConstant::HenryConstantH2;
namespace ProcessLib
{
namespace TwoPhaseComponentialFlow
{
    TwoPhaseComponentialFlowMaterialProperties::TwoPhaseComponentialFlowMaterialProperties(
    boost::optional<MeshLib::PropertyVector<int> const&> const material_ids,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>
        liquid_density,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>
        viscosity,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>
        gas_density,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>
        gas_viscosity,
    std::vector<Eigen::MatrixXd>
        intrinsic_permeability_models,
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>&&
        porosity_models,
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>&&
        storage_models,
    std::vector<std::unique_ptr<
        MaterialLib::PorousMedium::CapillaryPressureSaturation>>&&
        capillary_pressure_models,
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>&&
        relative_permeability_models)
    : _liquid_density(std::move(liquid_density)),
      _viscosity(std::move(viscosity)),
      _gas_density(std::move(gas_density)),
      _gas_viscosity(std::move(gas_viscosity)),
      _material_ids(material_ids),
      _intrinsic_permeability_models(intrinsic_permeability_models),
      _porosity_models(std::move(porosity_models)),
      _storage_models(std::move(storage_models)),
      _capillary_pressure_models(std::move(capillary_pressure_models)),
      _relative_permeability_models(std::move(relative_permeability_models))
{
    DBUG("Create material properties for Two-Phase flow with multi-component model.");
}

int TwoPhaseComponentialFlowMaterialProperties::getMaterialID(
    const std::size_t element_id)
{
    if (!_material_ids)
    {
        return 0;
    }

    assert(element_id < _material_ids->size());
    return (*_material_ids)[element_id];
}

double TwoPhaseComponentialFlowMaterialProperties::getLiquidDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _liquid_density->getValue(vars);
}

double TwoPhaseComponentialFlowMaterialProperties::getGasDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _gas_density->getValue(vars);
}

double TwoPhaseComponentialFlowMaterialProperties::getDerivGasDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;

    return _gas_density->getdValue(vars,
                                   MaterialLib::Fluid::PropertyVariableType::p);
}
double TwoPhaseComponentialFlowMaterialProperties::getLiquidViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _viscosity->getValue(vars);
}

double TwoPhaseComponentialFlowMaterialProperties::getGasViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _gas_viscosity->getValue(vars);
}

Eigen::MatrixXd const& TwoPhaseComponentialFlowMaterialProperties::getPermeability(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const int /*dim*/) const
{
    return _intrinsic_permeability_models[material_id];
}

double TwoPhaseComponentialFlowMaterialProperties::getPorosity(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double T, const double porosity_variable) const
{
    return _porosity_models[material_id]->getValue(porosity_variable, T);
}

double TwoPhaseComponentialFlowMaterialProperties::getNonwetRelativePermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    return _relative_permeability_models[0]->getValue(saturation);
}

double TwoPhaseComponentialFlowMaterialProperties::getWetRelativePermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    return _relative_permeability_models[1]->getValue(saturation);
}

double TwoPhaseComponentialFlowMaterialProperties::getCapillaryPressure(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double /*T*/, const double saturation) const
{
    return _capillary_pressure_models[material_id]->getCapillaryPressure(
        saturation);
}

double TwoPhaseComponentialFlowMaterialProperties::getCapillaryPressureDerivative(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double /*T*/, const double saturation) const
{
    return _capillary_pressure_models[material_id]->getdPcdS(saturation);
}

double TwoPhaseComponentialFlowMaterialProperties::getSaturation(const int material_id,
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double pc) const
{
    const double saturation =
        _capillary_pressure_models[material_id]->getSaturation(pc);
    return saturation;
}

double TwoPhaseComponentialFlowMaterialProperties::getDerivSaturation(const int material_id,
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    const double dpcdsw =
        _capillary_pressure_models[material_id]->getdPcdS(saturation);//
    const double dswdpc = 1 / dpcdsw;
    return dswdpc;
}

}  // end of namespace
}  // end of namespace
