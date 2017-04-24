/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#include "ThermalTwoPhaseFlowWithPPMaterialProperties.h"
#include <logog/include/logog.hpp>
#include "BaseLib/reorderVector.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MaterialLib/TwoPhaseModels/TwoPhaseFlowWithPPMaterialProperties.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "NumLib/NewtonRaphson.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"

/**
* Regarding the details of the material propertiy functions
* please refer to: Class H, Helmig R, Bastian P. Numerical simulation of
* non- isothermal multiphase multicomponent processes in Porous Media¨C¨C1. An
* efficient solution technique. Adv Water Resour, 2002.
*/

using MaterialLib::PhysicalConstant::MolarMass::Water;
using MaterialLib::PhysicalConstant::MolarMass::Air;
using MaterialLib::PhysicalConstant::IdealGasConstant;
namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
ThermalTwoPhaseFlowWithPPMaterialProperties::
    ThermalTwoPhaseFlowWithPPMaterialProperties(
        std::unique_ptr<MaterialLib::TwoPhaseFlowWithPP::TwoPhaseFlowWithPPMaterialProperties>
        two_phase_material_model,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            specific_heat_capacity_solid,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            specific_heat_capacity_water,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            specific_heat_capacity_air,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            specific_heat_capacity_vapor,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            thermal_conductivity_dry_solid,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            thermal_conductivity_wet_solid,
        std::vector<std::unique_ptr<
            MaterialLib::PorousMedium::CapillaryPressureSaturation>>&&
            capillary_pressure_models,
        std::vector<
            std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>&&
            relative_permeability_models)
    : _two_phase_material_model(std::move(two_phase_material_model)),
      _specific_heat_capacity_solid(std::move(specific_heat_capacity_solid)),
      _specific_heat_capacity_water(std::move(specific_heat_capacity_water)),
      _specific_heat_capacity_air(std::move(specific_heat_capacity_air)),
      _specific_heat_capacity_vapor(std::move(specific_heat_capacity_vapor)),
      _thermal_conductivity_dry_solid(std::move(thermal_conductivity_dry_solid)),
      _thermal_conductivity_wet_solid(std::move(thermal_conductivity_wet_solid)),
      _capillary_pressure_models(std::move(capillary_pressure_models)),
      _relative_permeability_models(std::move(relative_permeability_models))
{
    DBUG("Create material properties for non-isothermal two-phase flow model.");
}

int ThermalTwoPhaseFlowWithPPMaterialProperties::getMaterialID(
    const std::size_t element_id)
{
    if (!_material_ids)
    {
        return 0;
    }

    assert(element_id < _material_ids->size());
    return (*_material_ids)[element_id];
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getSpecificHeatCapacitySolid(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _specific_heat_capacity_solid->getValue(vars);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getSpecificHeatCapacityWater(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _specific_heat_capacity_water->getValue(vars);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getSpecificHeatCapacityAir(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _specific_heat_capacity_air->getValue(vars);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getSpecificHeatCapacityVapor(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _specific_heat_capacity_vapor->getValue(vars);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getThermalConductivityDrySolid(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _thermal_conductivity_dry_solid->getValue(vars);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getThermalConductivityWetSolid(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _thermal_conductivity_wet_solid->getValue(vars);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getNonwetRelativePermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    const double Se = (saturation - 0.15) / (1 - 0.15);
    if (saturation < 0.15)
        return 1.0;
    return pow(1 - Se, 3);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getWetRelativePermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    const double Se = (saturation - 0.15) / (1 - 0.15);
    if (Se < 0)
        return 0.0;
    return pow(Se, 3);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getSaturation(const int material_id,
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double pc) const
{
    return _capillary_pressure_models[material_id]->getSaturation(pc);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getCapillaryPressure(const int material_id,
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    return _capillary_pressure_models[material_id]
        ->getCapillaryPressure(saturation);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getSaturationDerivative(const int material_id,
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    const double dpcdsw =
        _capillary_pressure_models[material_id]->getdPcdS(saturation);
    const double dswdpc = 1 / dpcdsw;
    return dswdpc;
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::calculateUnsatHeatConductivity(
    double const /*t*/, ProcessLib::SpatialPosition const& /*x*/,
    double const Sw, double const lambda_pm_dry,
    double const lambda_pm_wet) const
{
    double lambda_pm =
        lambda_pm_dry + std::pow(Sw, 0.5) * (lambda_pm_wet - lambda_pm_dry);
    if (Sw > 1)
        lambda_pm = lambda_pm_wet;
    else if (Sw < 0)
        lambda_pm = lambda_pm_dry;
    return lambda_pm;
}
double
ThermalTwoPhaseFlowWithPPMaterialProperties::calculateSaturatedVaporPressure(
    const double T) const
{
    const double T_0 = 373.15;
    const double p_0 = 101325.0;
    const double h_wg = 2258000.0;
    return p_0 * exp(((1 / T_0) - (1 / T)) * Water * h_wg / IdealGasConstant);
}
double
ThermalTwoPhaseFlowWithPPMaterialProperties::calculateVaporPressureNonwet(
    const double pc, const double T, const double rho_mass_h2o) const
{
    const double p_sat = calculateSaturatedVaporPressure(T);
    const double c_w = Water / IdealGasConstant / T;
    return p_sat * exp(-pc * c_w / rho_mass_h2o);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculateDerivativedPsatdT(
    const double T) const
{
    const double T_0 = 373.15;
    const double p_0 = 101325;
    const double h_wg = 2258000.0;
    return p_0 * (Water * h_wg / IdealGasConstant) * (1. / T / T) *
           exp(((1. / T_0) - (1. / T)) * Water * h_wg / IdealGasConstant);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculateDerivativedPgwdT(
    const double pc, const double T, const double rho_mass_h2o) const
{
    const double c_w = Water / IdealGasConstant / T;
    const double p_sat = calculateSaturatedVaporPressure(T);
    const double dPsatdT = calculateDerivativedPsatdT(T);
    return dPsatdT * exp(-pc * c_w / rho_mass_h2o) +
           p_sat * exp(-pc * c_w / rho_mass_h2o) *
               (pc * Water / rho_mass_h2o / IdealGasConstant / T / T);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculateDerivativedPgwdPC(
    const double pc, const double T, const double rho_mass_h2o) const
{
    const double c_w = Water / IdealGasConstant / T;
    const double p_sat = calculateSaturatedVaporPressure(T);
    return p_sat * exp(-pc * c_w / rho_mass_h2o) * (-c_w / rho_mass_h2o);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculatedRhoNonwetdT(
    const double p_air_nonwet, const double p_vapor_nonwet, const double PC,
    const double T, const double rho_mass_h2o) const
{
    const double dPgwdT = calculateDerivativedPgwdT(PC, T, rho_mass_h2o);
    return -((p_air_nonwet * Air + p_vapor_nonwet * Water) / IdealGasConstant /
             T / T) +
           (Water - Air) * dPgwdT / IdealGasConstant / T;
}

}  // end of namespace
}  // end of namespace
