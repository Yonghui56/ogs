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

#include "BaseLib/reorderVector.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

namespace MaterialLib
{
namespace TwoPhaseFlowWithPP
{
TwoPhaseFlowWithPPMaterialProperties::TwoPhaseFlowWithPPMaterialProperties(
    bool const has_material_ids,
    MeshLib::PropertyVector<int> const& material_ids,
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
    int cap_pressure_model,
    int rel_wet_perm_model,
    int rel_nonwet_perm_model,
    std::array<double, 4>
        cap_pressure_value,
    std::array<double, 4>
        rel_wet_perm_value,
    std::array<double, 4>
        rel_nonwet_perm_value,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves_)
    : _has_material_ids(has_material_ids),
      _material_ids(material_ids),
      _liquid_density(std::move(liquid_density)),
      _viscosity(std::move(viscosity)),
      _gas_density(std::move(gas_density)),
      _gas_viscosity(std::move(gas_viscosity)),
      _intrinsic_permeability_models(intrinsic_permeability_models),
      _porosity_models(std::move(porosity_models)),
      _storage_models(std::move(storage_models)),
      _cap_pressure_model(cap_pressure_model),
      _rel_wet_perm_model(rel_wet_perm_model),
      _rel_nonwet_perm_model(rel_nonwet_perm_model),
      _cap_pressure_value(cap_pressure_value),
      _rel_wet_perm_value(rel_wet_perm_value),
      _rel_nonwet_perm_value(rel_nonwet_perm_value),
      _curves(curves_)
{
    DBUG("Create material properties for Two-Phase flow with PP model.");
}

void TwoPhaseFlowWithPPMaterialProperties::setMaterialID(
    const ProcessLib::SpatialPosition& pos)
{
    if (!_has_material_ids)
    {
        _current_material_id = 0;
        return;
    }

    assert(pos.getElementID().get() < _material_ids.size());
    _current_material_id = _material_ids[pos.getElementID().get()];
}

double TwoPhaseFlowWithPPMaterialProperties::getLiquidDensity(
    const double p, const double T) const
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
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pg)] = p;
    return _gas_density->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getDerivGasDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pg)] = p;

    return _gas_density->getdValue(
        vars, MaterialLib::Fluid::PropertyVariableType::pg);
}
double TwoPhaseFlowWithPPMaterialProperties::getLiquidViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pl)] = p;
    return _viscosity->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getGasViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pg)] = p;
    return _gas_viscosity->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getSaturation(double pc) const
{
    double Sw;
    /// TODO waiting for a better way to implemente the PC-S curve
    switch (_cap_pressure_model)
    {
        case 0:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        /// van Genuchten
        case 1:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        /// Brooks-Corey
        case 2:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        /// Liakopoulos
        case 10:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        default:
            OGS_FATAL("This model has not been implemented yet");
            break;
    }
    return Sw;
}

double TwoPhaseFlowWithPPMaterialProperties::getDerivSaturation(
    double const pc) const
{
    double dSwdPc;
    switch (_cap_pressure_model)
    {
        case 0:
            OGS_FATAL("This model has been implemented elsewhere");
            break;

        case 1:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        case 2:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        case 10:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        default:
            OGS_FATAL("This model has not been implemented yet");
            break;
    }
    return dSwdPc;
}

double TwoPhaseFlowWithPPMaterialProperties::getrelativePermeability_liquid(
    double const sw) const
{
    /// TODO waiting for a better way to implemente the Kr-S curve
    /*
    MathLib::PiecewiseLinearInterpolation const& interpolated_Kr =
        *curves.at("curveB");
    return interpolated_Kr.getValue(sw);
    */
    double rel_wet_perm;
    switch (_rel_wet_perm_model)
    {
        case 0:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        case 1:
            break;
        case 2:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        case 10:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        default:
            OGS_FATAL("This model has not been implemented yet");
            break;
    }
    return rel_wet_perm;
}

double TwoPhaseFlowWithPPMaterialProperties::getrelativePermeability_gas(
    double sw) const
{
    double rel_nonwet_perm;
    switch (_rel_nonwet_perm_model)
    {
        case 0:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        case 1:
            break;
        case 2:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        case 10:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        default:
            OGS_FATAL("This model has not been implemented yet");
            break;
    }
    return rel_nonwet_perm;
}

Eigen::MatrixXd const& TwoPhaseFlowWithPPMaterialProperties::getPermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const int /*dim*/) const
{
    return _intrinsic_permeability_models[_current_material_id];
}

double TwoPhaseFlowWithPPMaterialProperties::getPorosity(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double p, const double T, const double porosity_variable) const
{
    const double porosity =
        _porosity_models[_current_material_id]->getValue(porosity_variable, T);

    return porosity;
}

double TwoPhaseFlowWithPPMaterialProperties::getDissolvedGas(
    double const pg) const
{
    double const hen = 2e-6;     //
    double const M_air = 0.029;  // unit kg/mol
    return pg * hen * M_air;
}
/*!
* \brief Binary diffusion coefficent [m^2/s] for molecular water and CO2.
*
* To calculate the values, the \ref fullerMethod is used.
*/
double TwoPhaseFlowWithPPMaterialProperties::getgasDiffCoeff(
    double const pressure, double const temperature) const
{
    double const k = 1.3806504e-23;  // Boltzmann constant
    double const c = 4;  // slip parameter, can vary between 4 (slip condition)
                         // and 6 (stick condition)
    double const R_h = 1.72e-10;  // hydrodynamic radius of the solute
    const double mu = gasViscosity(temperature, pressure);  // CO2 viscosity
    return k / (c * M_PI * R_h) * (temperature / mu);
}

void TwoPhaseFlowWithPPMaterialProperties::calculateMoleFractions(
    double const pressure, double const temperature, double const salinity,
    double& xlCO2, double& ygH2O, double const rho_co2) const
{
	const double A = computeA_(pressure, temperature, rho_co2);
	/* salinity: conversion from mass fraction to mol fraction */
	const double x_NaCl = salinityToMolFrac_(salinity);
	const double p_sat = saturationPressure(temperature);
	if (pressure < 5e+4) {
		xlCO2 = pressure / henry(temperature);
		//virtual equilibrium mole fraction of water in the non-existing gas phase
		ygH2O = 0.0;// A * (1 - xlCO2 - x_NaCl);//virtual
		if (p_sat > pressure)
			ygH2O = 0.0;// if the pressure is smaller than the saturation pressure of pure water
	}
	else {
		double const molalityNaCl =
			moleFracToMolality_(x_NaCl);  // molality of NaCl //CHANGED
		double const m0_CO2 = molalityCO2inPureWater_(
			pressure, temperature, rho_co2);  // molality of CO2 in pure water
		double const gammaStar = activityCoefficient_(
			pressure, temperature,
			molalityNaCl);                  // activity coefficient of CO2 in brine
		double m_CO2 = m0_CO2 / gammaStar;  // molality of CO2 in brine
		xlCO2 = m_CO2 /
			(molalityNaCl + 55.508 + m_CO2);  // mole fraction of CO2 in brine
		ygH2O = 0.0;
			//A * (1 - xlCO2 - x_NaCl);  // mole fraction of water in the gas phase
	}
}

double TwoPhaseFlowWithPPMaterialProperties::moleFracCO2InBrine_duan(
    double temperature, double pressure, double const salinity,
    double const rhoCO2) const
{
    // regularisations:
    if (pressure > 2.5e8)
    {
        pressure = 2.5e8;
    }
    if (pressure < 2.e5)
    {
        pressure = 2.e5;
    }
    if (temperature < 275.)
    {
        temperature = 275;
    }
    if (temperature > 600.)
    {
        temperature = 600;
    }
	
    const double Ms = 58.8e-3; /* molecular weight of NaCl  [kg/mol] */
    /* salinity: conversion from mass fraction to mole fraction */
    const double x_NaCl = salinityToMolFrac_(salinity);
    // salinity: conversion from mole fraction to molality
    const double mol_NaCl = -55.56 * x_NaCl / (x_NaCl - 1);
    const double A =
        computeA_duan(temperature, pressure); /* mu_{CO2}^{l(0)}/RT */
    const double B =
        computeB_duan(temperature, pressure); /* lambda_{CO2-Na+} */
    const double C =
        computeC_duan(temperature, pressure); /* Xi_{CO2-Na+-Cl-} */
    const double pgCO2 = partialPressureCO2_(temperature, pressure);
    const double phiCO2 = fugacityCoeffCO2_(temperature, pressure, rhoCO2);
	//const double phiCO2 = fugacityCoefficientCO2(pressure, temperature, rhoCO2);
    const double exponent =
        A - log(phiCO2) + 2 * B * mol_NaCl + C * pow(mol_NaCl, 2);
    const double mol_CO2w =
        pgCO2 / (1e5 * exp(exponent)); /* paper: equation (6) */
    const double x_CO2w =
        mol_CO2w /
        (mol_CO2w + 55.56); /* conversion: molality to mole fraction */
    // Scalar X_CO2w = x_CO2w*MCO2/(x_CO2w*MCO2 + (1-x_CO2w)*Mw);   /*
    // conversion: mole fraction to mass fraction */
    return x_CO2w;
}

/*!
    * \brief The dynamic viscosity [Pa s] of CO2.
    *
    * Equations given in: - Vesovic et al., 1990
    *                        - Fenhour etl al., 1998
    */
double TwoPhaseFlowWithPPMaterialProperties::gasViscosity(
    double const pressure, const double temperature) const
{
    const double a0 = 0.235156;
    const double a1 = -0.491266;
    const double a2 = 5.211155e-2;
    const double a3 = 5.347906e-2;
    const double a4 = -1.537102e-2;

    const double d11 = 0.4071119e-2;
    const double d21 = 0.7198037e-4;
    const double d64 = 0.2411697e-16;
    const double d81 = 0.2971072e-22;
    const double d82 = -0.1627888e-22;

    const double ESP = 251.196;

    const double TStar = temperature / ESP;

    // mu0: viscosity in zero-density limit
    const double& logTStar = std::log(TStar);
    double SigmaStar = std::exp(
        a0 +
        logTStar * (a1 + logTStar * (a2 + logTStar * (a3 + logTStar * a4))));

    double mu0 = 1.00697 * std::sqrt(temperature) / SigmaStar;

    // const double rho = gasDensity(temperature, pressure); // CO2 mass density
    // [kg/m^3]
    const double rho = 1000;
    // dmu : excess viscosity at elevated density
    double dmu = d11 * rho + d21 * rho * rho +
                 d64 * std::pow(rho, 6.0) / (TStar * TStar * TStar) +
                 d81 * std::pow(rho, 8.0) + d82 * std::pow(rho, 8.0) / TStar;

    return (mu0 + dmu) / 1.0e6;  // conversion to [Pa s]
}
/*!
* \brief Henry coefficent \f$[N/m^2]\f$  for molecular CO2 in liquid water.
*  unit \f$Pa\f$
* See:
*
* IAPWS: "Guideline on the Henry's Constant and Vapor-Liquid
* Distribution Constant for Gases in H2O and D2O at High
* Temperatures"
* http://www.iapws.org/relguide/HenGuide.pdf
*/
double TwoPhaseFlowWithPPMaterialProperties::henry(const double temperature) const
{
	const double E = 1672.9376;
	const double F = 28.1751;
	const double G = -112.4619;
	const double H = 85.3807;

	return henryIAPWS(E, F, G, H, temperature);
}
double TwoPhaseFlowWithPPMaterialProperties::saturationvaporpressure(double const temperature) const
{
	const std::array<double, 10> n = {
		0.11670521452767e4,  -0.72421316703206e6, -0.17073846940092e2,
		0.12020824702470e5,  -0.32325550322333e7, 0.14915108613530e2,
		-0.48232657361591e4, 0.40511340542057e6,  -0.23855557567849,
		0.65017534844798e3 };
	double sigma = temperature + n[8] / (temperature - n[9]);
	double A = (sigma + n[0]) * sigma + n[1];
	double B = (n[2] * sigma + n[3]) * sigma + n[4];
	double C = (n[5] * sigma + n[6]) * sigma + n[7];
	double tmp = 2 * C / (std::sqrt(B * B - 4 * A * C) - B);
	tmp *= tmp;
	tmp *= tmp;
	return 1e6 * tmp;
}

}  // end of namespace
}  // end of namespace
