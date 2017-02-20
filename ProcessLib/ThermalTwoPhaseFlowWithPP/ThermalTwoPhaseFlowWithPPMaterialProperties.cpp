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
    DBUG("Create material properties for Two-Phase flow with Prho model.");
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

double ThermalTwoPhaseFlowWithPPMaterialProperties::getLiquidDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _liquid_density->getValue(vars);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getGasDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _gas_density->getValue(vars);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getDerivativeGasDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;

    return _gas_density->getdValue(vars,
                                   MaterialLib::Fluid::PropertyVariableType::p);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::getLiquidViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _viscosity->getValue(vars);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getGasViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _gas_viscosity->getValue(vars);
}

Eigen::MatrixXd const&
ThermalTwoPhaseFlowWithPPMaterialProperties::getPermeability(const int material_id,
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const int /*dim*/) const
{
    return _intrinsic_permeability_models[material_id];
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getPorosity(const int material_id,
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double T, const double porosity_variable) const
{
    return
        _porosity_models[material_id]->getValue(porosity_variable, T);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getNonwetRelativePermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    /*const double nonwet_krel =
        _relative_permeability_models[0]->getValue(saturation);
    return nonwet_krel;*/
    const double Se = (saturation - 0.15) / (1 - 0.15);
    if (saturation < 0.15)
        return 1.0;
    return pow(1 - Se, 3);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getWetRelativePermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    /*const double wet_krel =
        _relative_permeability_models[1]->getValue(saturation);
    return wet_krel;*/
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
    const double C_v = Air / IdealGasConstant / T;
    const double dPgwdT = calculateDerivativedPgwdT(PC, T, rho_mass_h2o);
    return -((p_air_nonwet * Air + p_vapor_nonwet * Water) / IdealGasConstant /
             T / T) +
           (Water - Air) * dPgwdT / IdealGasConstant / T;
}

}  // end of namespace
}  // end of namespace
