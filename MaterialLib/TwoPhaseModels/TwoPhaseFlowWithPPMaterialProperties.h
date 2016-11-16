/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TWOPHASEFLOWMaterialProperties.h
 *
 * Created on August 18, 2016, 11:03 AM
 */

#ifndef OGS_TWOPHASEFLOWWITHPPMATERIALPROPERTIES_H
#define OGS_TWOPHASEFLOWWITHPPMATERIALPROPERTIES_H

#include <iostream>
#include <memory>
#include <vector>
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace ProcessLib
{
class SpatialPosition;
}
namespace MaterialLib
{
namespace PorousMedium
{
class Porosity;
class Storage;
}
}

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace MaterialLib
{
namespace TwoPhaseFlowWithPP
{
class TwoPhaseFlowWithPPMaterialProperties
{
public:
    typedef MaterialLib::Fluid::FluidProperty::ArrayType ArrayType;

    TwoPhaseFlowWithPPMaterialProperties(
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
            curves_);

    void setMaterialID(const ProcessLib::SpatialPosition& pos);

    /**
     * \brief Compute the coefficient of the mass term by
     *      \f[
     *           n \frac{partial \rho_l}{\partial p} + \beta_s
     *      \f]
     *     where \f$n\f$ is the porosity, \f$rho_l\f$ is the liquid density,
     *     \f$bata_s\f$ is the storage.
     * \param t                  Time.
     * \param pos                Position of element.
     * \param porosity_variable  The first variable for porosity model, and it
     *                           passes a double type value that could be
     *                           saturation, and invariant of stress or strain.
     * \param storage_variable   Variable for storage model.
     * \param p                  Pressure value
     * \param T                  Temperature value
     */

    Eigen::MatrixXd const& getPermeability(
        const double t,
        const ProcessLib::SpatialPosition& pos,
        const int dim) const;

    double getPorosity(const double t, const ProcessLib::SpatialPosition& pos,
                       const double p, const double T,
                       const double porosity_variable) const;

    double getLiquidDensity(const double p, const double T) const;
    double getGasDensity(const double p, const double T) const;
    double getGasViscosity(const double p, const double T) const;
    double getLiquidViscosity(const double p, const double T) const;
    double getDerivGasDensity(double const p, double const T) const;
    double getDissolvedGas(double const pg) const;

    virtual double getSaturation(double pc) const;
    virtual double getrelativePermeability_liquid(double const sw) const;
    virtual double getrelativePermeability_gas(double const sw) const;
    virtual double getDerivSaturation(double const pc) const;

    // EOS function
    /*!
    * \brief Binary diffusion coefficent [m^2/s] for molecular water and CO2.
    *
    * To calculate the values, the \ref fullerMethod is used. Spycher and Pruess
    */
    double getgasDiffCoeff(double const pressure,
                           double const temperature, double const rho_co2) const;
    void calculateMoleFractions(double const pressure, double const temperature,
                                double const salinity, double& xlCO2,
                                double& ygH2O, double const rho_co2) const;
    double gasViscosity(double const pressure, const double temperature, const double rho_co2/*co2 mass density*/) const;
    /*!
    * \brief Binary diffusion coefficent [m^2/s] for molecular water and CO2.
    *
    * To calculate the values, the \ref fullerMethod is used.
    * Duan 2003
    */
    double moleFracCO2InBrine_duan(double temperature, double pressure,
                                   double const salinity,
                                   double const rhoCO2) const;
	double henry(double const temperature) const;
	double saturationvaporpressure(double const temperature) const;
protected:
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _liquid_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _viscosity;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _gas_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _gas_viscosity;

    int _cap_pressure_model;
    std::array<double, 4> _cap_pressure_value;

    int _rel_wet_perm_model;
    std::array<double, 4> _rel_wet_perm_value;

    int _rel_nonwet_perm_model;
    std::array<double, 4> _rel_nonwet_perm_value;

    /// A flag to indicate whether the reference member, _material_ids,
    /// is not assigned.
    const bool _has_material_ids;

    /** Use porous medium models for different material zones.
     *  Material IDs must be given as mesh element properties.
     */
    MeshLib::PropertyVector<int> const& _material_ids;

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        _curves;

    int _current_material_id = 0;
    std::vector<Eigen::MatrixXd> _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _storage_models;

    // Note: For the statistical data of porous media, they could be read from
    // vtu files directly. This can be done by using property vectors directly.
    // Such property vectors will be added here if they are needed.
private:
    /*
    * molar mass of each component
    */
    const double molar_mass_co2 = MaterialLib::PhysicalConstant::MolarMass::CO2;
    const double molar_mass_h2o =
        MaterialLib::PhysicalConstant::MolarMass::Water;
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
    double criticalPressure_co2() const { return 73.8e5; /* [Pa] */ }
    double criticalTemperature_co2() const { return 273.15 + 30.95; /* [K] */ }
    /*!
    * \brief Returns the molality of NaCl (mol NaCl / kg water) for a given mole
    * fraction
    *
    * \param salinity the salinity [kg NaCl / kg solution]
    */
    double salinityToMolFrac_(double salinity) const
    {
        const double Ms = 58.8e-3; /* molecular weight of NaCl  [kg/mol] */

        const double X_NaCl = salinity;
        /* salinity: conversion from mass fraction to mol fraction */
        const double x_NaCl =
            -molar_mass_h2o * X_NaCl / ((Ms - molar_mass_h2o) * X_NaCl - Ms);
        return x_NaCl;
    }
    /*!
    * \brief Returns the molality of NaCl (mol NaCl / kg water) for a given mole
    * fraction (mol NaCl / mol solution)
    *
    * \param x_NaCl mole fraction of NaCL in brine [mol/mol]
    */
    double moleFracToMolality_(const double x_NaCl) const
    {
        // conversion from mol fraction to molality (dissolved CO2 neglected)
        return 55.508 * x_NaCl / (1 - x_NaCl);
    }

    /*!
    * \brief Returns the equilibrium molality of CO2 (mol CO2 / kg water) for a
    * CO2-water mixture at a given pressure and temperature
    *
    * \param temperature The temperature [K]
    * \param pg The gas phase pressure [Pa]
    */
    double molalityCO2inPureWater_(const double pressure,
                                   const double temperature,
                                   const double rho_co2) const
    {
        const double A = computeA_(
            pressure, temperature,
            rho_co2);  // according to Spycher, Pruess and Ennis-King (2003)
        const double B = computeB_(
            pressure, temperature,
            rho_co2);  // according to Spycher, Pruess and Ennis-King (2003)
        const double yH2OinGas =
            (1 - B) /
            (1. / A - B);  // equilibrium mol fraction of H2O in the gas phase
        const double xCO2inWater =
            B *
            (1 -
             yH2OinGas);  // equilibrium mol fraction of CO2 in the water phase
        return (xCO2inWater * 55.508) / (1 - xCO2inWater);  // CO2 molality
    }
    /*!
    * \brief Returns the activity coefficient of CO2 in brine for a
    *           molal description. According to "Duan and Sun 2003"
    *           given in "Spycher and Pruess 2005"
    *
    * \param temperature the temperature [K]
    * \param pg the gas phase pressure [Pa]
    * \param molalityNaCl molality of NaCl (mol NaCl / kg water)
    */
    double activityCoefficient_(const double pressure,
                                const double temperature,
                                double molalityNaCl) const
    {
        const double lambda =
            computeLambda_(pressure, temperature);  // lambda_{CO2-Na+}
        const double xi =
            computeXi_(pressure, temperature);  // Xi_{CO2-Na+-Cl-}
        const double lnGammaStar =
            2 * molalityNaCl * lambda + xi * molalityNaCl * molalityNaCl;
        return std::exp(lnGammaStar);
    }
    /*!
    * \brief Returns the paramater A for the calculation of
    * them mutual solubility in the water-CO2 system.
    * Given in Spycher, Pruess and Ennis-King (2003)
    *
    * \param T the temperature [K]
    * \param pg the gas phase pressure [Pa]
    */
    double computeA_(double const pressure, double const temperature,
                     double const rho_co2) const
    {
        const double& deltaP =
            pressure / 1e5 -
            1;  // pressure range [bar] from p0 = 1bar to pg[bar]
        double const v_av_H2O =
            18.1;  // average partial molar volume of H2O [cm^3/mol]
        double R = 8.314 * 10;
        const double k0_H2O = equilibriumConstantH2O_(
            temperature);  // equilibrium constant for H2O at 1 bar
        const double& phi_H2O = fugacityCoefficientH2O(
            pressure, temperature,
            rho_co2);  // fugacity coefficient of H2O for the water-CO2 system
        const double& pg_bar = pressure / 1.e5;
        return k0_H2O / (phi_H2O * pg_bar) *
               std::exp(deltaP * v_av_H2O / (R * temperature));
    }
    /*!
    * \brief Returns the paramater B for the calculation of
    * the mutual solubility in the water-CO2 system.
    * Given in Spycher, Pruess and Ennis-King (2003)
    *
    * \param temperature the temperature [K]
    * \param pg the gas phase pressure [Pa]
    */
    double computeB_(const double pressure, const double temperature,
                     const double rho_co2) const
    {
        const double deltaP =
            pressure / 1e5 -
            1;  // pressure range [bar] from p0 = 1bar to pg[bar]
        const double v_av_CO2 =
            32.6;  // average partial molar volume of CO2 [cm^3/mol]
        const double R = MaterialLib::PhysicalConstant::IdealGasConstant * 10;
        const double k0_CO2 = equilibriumConstantCO2_(
            temperature);  // equilibrium constant for CO2 at 1 bar
        const double phi_CO2 = fugacityCoefficientCO2(
            pressure, temperature,
            rho_co2);  // fugacity coefficient of CO2 for the water-CO2 system
        const double pg_bar = pressure / 1.e5;
        return phi_CO2 * pg_bar / (55.508 * k0_CO2) *
               std::exp(-(deltaP * v_av_CO2) / (R * temperature));
    }

    /*!
    * \brief Returns the parameter lambda, which is needed for the
    * calculation of the CO2 activity coefficient in the brine-CO2 system.
    * Given in Spycher and Pruess (2005)
    * \param temperature the temperature [K]
    * \param pg the gas phase pressure [Pa]
    */
    double computeLambda_(const double pressure, const double temperature) const
    {
        static const std::array<double, 6> c = {-0.411370585, 6.07632013E-4,
                                                97.5347708,   -0.0237622469,
                                                0.0170656236, 1.41335834E-5};

        double pg_bar = pressure / 1.0E5; /* conversion from Pa to bar */
        return c[0] + c[1] * temperature + c[2] / temperature +
               c[3] * pg_bar / temperature +
               c[4] * pg_bar / (630.0 - temperature) +
               c[5] * temperature * std::log(pg_bar);
    }
    /*!
    * \brief Returns the parameter xi, which is needed for the
    * calculation of the CO2 activity coefficient in the brine-CO2 system.
    * Given in Spycher and Pruess (2005)
    * \param temperature the temperature [K]
    * \param pg the gas phase pressure [Pa]
    */
    double computeXi_(const double pressure, const double temperature) const
    {
        static const std::array<double, 6> c = {3.36389723E-4, -1.98298980E-5,
                                                2.12220830E-3, -5.24873303E-3};

        double pg_bar = pressure / 1.0E5;
        return c[0] + c[1] * temperature + c[2] * pg_bar / temperature +
               c[3] * pg_bar / (630.0 - temperature);
    }
    /*!
    * \brief Returns the equilibrium constant for CO2, which is needed for the
    * calculation of the mutual solubility in the water-CO2 system
    * Given in Spycher, Pruess and Ennis-King (2003)
    * \param temperature the temperature [K]
    */
    double equilibriumConstantCO2_(double const temperature) const
    {
        double temperatureCelcius = temperature - 273.15;
        std::array<double, 3> c = {1.189, 1.304e-2, -5.446e-5};
        double logk0_CO2 =
            c[0] + temperatureCelcius * (c[1] + temperatureCelcius * c[2]);
        double k0_CO2 = std::pow(10.0, logk0_CO2);
        return k0_CO2;
    }
    /*!
    * \brief Returns the equilibrium constant for H2O, which is needed for the
    * calculation of the mutual solubility in the water-CO2 system
    * Given in Spycher, Pruess and Ennis-King (2003)
    * \param temperature the temperature [K]
    */
    static double equilibriumConstantH2O_(double const temperature)
    {
        double temperatureCelcius = temperature - 273.15;
        std::array<double, 4> c = {-2.209, 3.097e-2, -1.098e-4, 2.048e-7};
        double logk0_H2O =
            c[0] +
            temperatureCelcius *
                (c[1] +
                 temperatureCelcius * (c[2] + temperatureCelcius * c[3]));
        return std::pow(10.0, logk0_H2O);
    }
    /*!
    * \brief Returns the fugacity coefficient of the CO2 component in a
    * water-CO2 mixture
    *
    * (given in Spycher, Pruess and Ennis-King (2003))
    *
    * \param T the temperature [K]
    * \param pg the gas phase pressure [Pa]
    */
    double fugacityCoefficientCO2(const double pressure,
                                  const double temperature,
                                  const double rho_co2) const
    {
        const double molar_mass_co2 =
            MaterialLib::PhysicalConstant::MolarMass::CO2;
        const double V =
            1 / (rho_co2 / molar_mass_co2) * 1.e6;  // molar volume in cm^3/mol
        const double pg_bar = pressure / 1.e5;      // gas phase pressure in bar
        const double a_CO2 =
            (7.54e7 -
             4.13e4 *
                 temperature);  // mixture parameter of  Redlich-Kwong equation
        const double b_CO2 =
            27.8;  // mixture parameter of Redlich-Kwong equation
        const double R = MaterialLib::PhysicalConstant::IdealGasConstant *
                         10.;  // ideal gas constant with unit bar cm^3 /(K mol)
        double lnPhiCO2 = std::log(V / (V - b_CO2));
        lnPhiCO2 += b_CO2 / (V - b_CO2);
        lnPhiCO2 -= 2 * a_CO2 / (R * std::pow(temperature, 1.5) * b_CO2) *
                    log((V + b_CO2) / V);
        lnPhiCO2 += a_CO2 * b_CO2 /
                    (R * std::pow(temperature, 1.5) * b_CO2 * b_CO2) *
                    (std::log((V + b_CO2) / V) - b_CO2 / (V + b_CO2));
        lnPhiCO2 -= std::log(pg_bar * V / (R * temperature));

        return std::exp(lnPhiCO2);  // fugacity coefficient of CO2
    }
    double fugacityCoefficientH2O(const double pressure,
                                  const double temperature,
                                  const double rho_co2) const
    {
        const double molar_co2 = MaterialLib::PhysicalConstant::MolarMass::CO2;
        const double V =
            1 / (rho_co2 / molar_co2) * 1.e6;   // molar volume in cm^3/mol
        const double pg_bar = pressure / 1.e5;  // gas phase pressure in bar
        const double a_CO2 =
            (7.54e7 -
             4.13e4 *
                 temperature);  // mixture parameter of  Redlich-Kwong equation
        const double a_CO2_H2O =
            7.89e7;  // mixture parameter of Redlich-Kwong equation
        const double b_CO2 =
            27.8;  // mixture parameter of Redlich-Kwong equation
        const double b_H2O =
            18.18;  // mixture parameter of Redlich-Kwong equation
        const double R = MaterialLib::PhysicalConstant::IdealGasConstant *
                         10.;  // ideal gas constant with unit bar cm^3 /(K mol)
        const double lnPhiH2O =
            std::log(V / (V - b_CO2)) + b_H2O / (V - b_CO2) -
            2 * a_CO2_H2O / (R * std::pow(temperature, 1.5) * b_CO2) *
                std::log((V + b_CO2) / V) +
            a_CO2 * b_H2O / (R * std::pow(temperature, 1.5) * b_CO2 * b_CO2) *
                (std::log((V + b_CO2) / V) - b_CO2 / (V + b_CO2)) -
            std::log(pg_bar * V / (R * temperature));
        return std::exp(lnPhiH2O);  // fugacity coefficient of H2O
    }
    /* EoS for Duan 2003
    * calculate the solubility of co2 in brine
    */
    double computeA_duan(double T, double pg) const
    {
        static const std::array<double, 10> c = {
            28.9447706,    -0.0354581768, -4770.67077,    1.02782768E-5,
            33.8126098,    9.04037140E-3, -1.14934031E-3, -0.307405726,
            -0.0907301486, 9.32713393E-4,
        };
        const double pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        const double Tr = 630.0 - T;
        return c[0] + c[1] * T + c[2] / T + c[3] * T * T + c[4] / Tr +
               c[5] * pg_bar + c[6] * pg_bar * log(T) + c[7] * pg_bar / T +
               c[8] * pg_bar / Tr + c[9] * pg_bar * pg_bar / (Tr * Tr);
    }
    /*!
    * \brief computation of B
    *
    * \param T the temperature \f$\mathrm{[K]}\f$
    * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
    */
    double computeB_duan(double T, double pg) const
    {
        const double c1 = -0.411370585;
        const double c2 = 6.07632013E-4;
        const double c3 = 97.5347708;
        const double c8 = -0.0237622469;
        const double c9 = 0.0170656236;
        const double c11 = 1.41335834E-5;
        const double pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        return c1 + c2 * T + c3 / T + c8 * pg_bar / T +
               c9 * pg_bar / (630.0 - T) + c11 * T * std::log(pg_bar);
    }
    /*!
    * \brief computation of C
    *
    * \param T the temperature \f$\mathrm{[K]}\f$
    * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
    */
    double computeC_duan(double T, double pg) const
    {
        const double c1 = 3.36389723E-4;
        const double c2 = -1.98298980E-5;
        const double c8 = 2.12220830E-3;
        const double c9 = -5.24873303E-3;
        const double pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        return c1 + c2 * T + c8 * pg_bar / T + c9 * pg_bar / (630.0 - T);
    }
    /*!
    * \brief computation of partial pressure CO2
    *
    * We assume that the partial pressure of brine is its vapor pressure.
    * \warning: Strictly this is assumption is invalid for CO2 because the
    *           mole fraction of CO2 in brine can be considerable
    *
    * \param temperature the temperature \f$\mathrm{[K]}\f$
    * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
    */
    double partialPressureCO2_(double temperature, double pg) const
    {
        double T = std::min(T, criticalTemperature_co2());
        T = std::max(T, 273.15);
        return pg - saturationPressure(T);
    }
	double saturationPressure(double temperature) const
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
 
    /*!
    * \brief The fugacity coefficent of CO2 for a CO2-H2O mixture.
    *
    * \param temperature the temperature \f$\mathrm{[K]}\f$
    * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
    * \param rhoCO2 the density of CO2 for the critical volume
    * \f$\mathrm{[kg/m^3]}\f$
    */
    double fugacityCoeffCO2_(double temperature, double pg, double rhoCO2) const
    {
        const std::array<double, 15> a = {
            8.99288497E-2,  -4.94783127E-1, 4.77922245E-2, 1.03808883E-2,
            -2.82516861E-2, 9.49887563E-2,  5.20600880E-4, -2.93540971E-4,
            -1.77265112E-3, -2.51101973E-5, 8.93353441E-5, 7.88998563E-5,
            -1.66727022E-2, 1.3980,         2.96000000E-2};
        // reduced temperature
        const double Tr = temperature / criticalTemperature_co2();
        // reduced pressure
        const double pr = pg / criticalPressure_co2();
        // reduced molar volume. ATTENTION: Vc is _NOT_ the critical
        // molar volume of CO2. See the reference!
        const double Vc =
            R * criticalTemperature_co2() / criticalPressure_co2();
        const double Vr =
            // molar volume of CO2 at (temperature, pg)
            molar_mass_co2 / rhoCO2 *
            // "pseudo-critical" molar volume
            1.0 / Vc;
        // the Z coefficient
        const double Z = pr * Vr / Tr;
        const double A = a[0] + a[1] / (Tr * Tr) + a[2] / (Tr * Tr * Tr);
        const double B = a[3] + a[4] / (Tr * Tr) + a[5] / (Tr * Tr * Tr);
        const double C = a[6] + a[7] / (Tr * Tr) + a[8] / (Tr * Tr * Tr);
        const double D = a[9] + a[10] / (Tr * Tr) + a[11] / (Tr * Tr * Tr);
        const double lnphiCO2 = Z - 1 - log(Z) + A / Vr + B / (2 * Vr * Vr) +
                                C / (4 * Vr * Vr * Vr * Vr) +
                                D / (5 * Vr * Vr * Vr * Vr * Vr) +
                                a[12] / (2 * Tr * Tr * Tr * a[14]) *
                                    (a[13] + 1 -
                                     (a[13] + 1 + a[14] / (Vr * Vr)) *
                                         std::exp(-a[14] / (Vr * Vr)));
        return std::exp(lnphiCO2);
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
	double henryIAPWS(double E,
		double F,
		double G,
		double H,
		const double temperature) const
	{

		double Tr = temperature / 647.096;// H2O::criticalTemperature();
		double tau = 1 - Tr;

		const std::array<double,6> c = {
			1.99274064, 1.09965342, -0.510839303,
			-1.75493479,-45.5170352, -6.7469445e5
		};
		const std::array<double, 6> d = {
			1 / 3.0, 2 / 3.0, 5 / 3.0,
			16 / 3.0, 43 / 3.0, 110 / 3.0
		};
		const double q = -0.023767;

		double f = 0;
		for (int i = 0; i < 6; ++i) {
			f += c[i] * std::pow(tau, d[i]);
		}

		const double& exponent =
			q*F +
			E / temperature*f +
			(F +
				G*std::pow(tau, 2.0 / 3) +
				H*tau)*
			std::exp((273.15 - temperature) / 100);
		// CAUTION: K_D is formulated in mole fractions. We have to
		// multiply it with the vapor pressure of water in order to get
		// derivative of the partial pressure.
		return std::exp(exponent)*saturationPressure(temperature);
	}
};

}  // end of namespace
}  // end of namespace
#endif /* TWOPHASEFLOWWITHPPMATERIALPROPERTIES_H */
