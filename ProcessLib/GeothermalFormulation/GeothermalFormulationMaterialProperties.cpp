/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#include "GeothermalFormulationMaterialProperties.h"
#include <logog/include/logog.hpp>
#include <utility>
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
#include "Numlib/NewtonRaphson.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"
namespace ProcessLib
{
namespace GeothermalFormulation
{
GeothermalFormulationMaterialProperties::
    GeothermalFormulationMaterialProperties(
        boost::optional<MeshLib::PropertyVector<int> const&> const material_ids,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            liquid_density,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            liquid_viscosity,
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
      _liquid_viscosity(std::move(liquid_viscosity)),
      _gas_density(std::move(gas_density)),
      _gas_viscosity(std::move(gas_viscosity)),
      _material_ids(material_ids),
      _intrinsic_permeability_models(std::move(intrinsic_permeability_models)),
      _porosity_models(std::move(porosity_models)),
      _storage_models(std::move(storage_models)),
      _capillary_pressure_models(std::move(capillary_pressure_models)),
      _relative_permeability_models(std::move(relative_permeability_models))
{
    DBUG("Create material properties for Two-Phase flow with PP model.");
}

int GeothermalFormulationMaterialProperties::getMaterialID(
    const std::size_t element_id)
{
    if (!_material_ids)
    {
        return 0;
    }

    assert(element_id < _material_ids->size());
    return (*_material_ids)[element_id];
}

double GeothermalFormulationMaterialProperties::getLiquidDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _liquid_density->getValue(vars);
}

double GeothermalFormulationMaterialProperties::getGasDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _gas_density->getValue(vars);
}

double GeothermalFormulationMaterialProperties::getGasDensityDerivative(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;

    return _gas_density->getdValue(vars,
                                   MaterialLib::Fluid::PropertyVariableType::p);
}
double GeothermalFormulationMaterialProperties::getLiquidViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _liquid_viscosity->getValue(vars);
}

double GeothermalFormulationMaterialProperties::getGasViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _gas_viscosity->getValue(vars);
}

Eigen::MatrixXd const& GeothermalFormulationMaterialProperties::getPermeability(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const int /*dim*/) const
{
    return _intrinsic_permeability_models[material_id];
}

double GeothermalFormulationMaterialProperties::getPorosity(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double T, const double porosity_variable) const
{
    return _porosity_models[material_id]->getValue(porosity_variable, T);
}

double GeothermalFormulationMaterialProperties::getNonwetRelativePermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    return _relative_permeability_models[0]->getValue(saturation);
}

double GeothermalFormulationMaterialProperties::getWetRelativePermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    return _relative_permeability_models[1]->getValue(saturation);
}

double GeothermalFormulationMaterialProperties::getSaturation(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double /*T*/, const double pc) const
{
    return _capillary_pressure_models[material_id]->getSaturation(pc);
}
double GeothermalFormulationMaterialProperties::getSaturationDerivative(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double /*T*/, const double saturation) const
{
    const double dpcdsw =
        _capillary_pressure_models[material_id]->getdPcdS(saturation);
    return 1 / dpcdsw;
}

double GeothermalFormulationMaterialProperties::getViscositySteam(
    const double /*p*/, const double T) const
{
    // here T is in [K]
    const double mu_steam = 1e-4 * (0.407 * (T - 273.15) + 80.4);
    return mu_steam;
}

double GeothermalFormulationMaterialProperties::getViscosityWater(
    const double /*p*/, const double T) const
{
    // HERE T is in [K]
    const double mu_water =
        1e-4 * (241.4 * std::pow(10, 247.8 / ((T - 273.15) + 133.15)));
    return mu_water;
}

double GeothermalFormulationMaterialProperties::getDensityWater(
    const double p, const double h) const
{
    // The density is in kg/m3
    const double rho_water =
        1000 * (1.00207 + 4.42607 * 1e-10 * p - 5.47456e-8 * h +
                5.02875e-16 * h * p - 1.24791e-13 * h * h);
    return rho_water;
}

double GeothermalFormulationMaterialProperties::getDerivWaterDensity_dPressure(
    const double p, const double h) const
{
    // The density is in kg/m3
    const double d_rhowater_d_p =
        1000 * (4.42607 * 1e-10 +
            5.02875e-16 * h );
    return d_rhowater_d_p;
}

double GeothermalFormulationMaterialProperties::getDerivWaterDensity_dEnthalpy(
    const double p, const double h) const
{
    // The density is in kg/m3
    const double d_rhowater_d_h =
        1000 * (-5.47456e-8 +
            5.02875e-16 * p - 2 * 1.24791e-13 * h);
    return d_rhowater_d_h;
}


double GeothermalFormulationMaterialProperties::getDensitySteam(
    const double p, const double h) const
{
    // The density is in kg/m3
    const double rho_steam =
        1000 * (-2.26162e-5 + 4.38411e-8 * p - 1.79088e-14 * p * h +
                3.69276e-32 * std::pow(p, 4) + 5.17644e-28 * p * pow(h, 3));
    return rho_steam;
}

double GeothermalFormulationMaterialProperties::getDerivSteamDensity_dPressure(
    const double p, const double h) const
{
    // The density is in kg/m3
    const double d_rhosteam_d_p =
        1000 * (4.38411e-8 - 1.79088e-14 * h +
            4* 3.69276e-32 * std::pow(p, 3) + 5.17644e-28 * pow(h, 3));
    return d_rhosteam_d_p;
}

double GeothermalFormulationMaterialProperties::getDerivSteamDensity_dEnthalpy(
    const double p, const double h) const
{
    // The density is in kg/m3
    const double d_rhosteam_d_h =
        1000 * (- 1.79088e-14 * h +
             3*5.17644e-28 *p* pow(h, 2));
    return d_rhosteam_d_h;
}


double GeothermalFormulationMaterialProperties::getTemperatureSteamRegion(
    const double p, const double h) const
{
    // temperature is in K
    double const temperature =
        273.15 - 374.669 + 4.79921e-5* p - 6.33606e-13 * p * p +
        7.39386e-11 * h * h - 3.3372e+24 * std::pow(h, -2) * std::pow(p, -2) +
        3.57154e+16* std::pow(p, -3) - 1.1725e-24 * std::pow(h, 3) * p -
        2.26861e+27 * std::pow(h, -4);
    return temperature;
}

double GeothermalFormulationMaterialProperties::getTemperatureWaterRegion(
    const double p, const double h) const
{
    // temperature is in K
    // SI unit
    double const temperature = 273.15 - 2.41231 + 2.56222e-4 * h -
                               9.31415e-15 * std::pow(p, 2) -
                               2.2568e-11 * std::pow(h, 2);
    return temperature;
}

double GeothermalFormulationMaterialProperties::getSatSteamEnthalpy(
    const double p, const double /*h*/
) const
{
    double const steamEnthalpy = 1e-7*(2.82282e+10 - 3.91952e+5*std::pow(p, -1)
        + 2.54342e+21*std::pow(p, -1) - 9.38879e-8*std::pow(p, 2));
    return steamEnthalpy;
}

double GeothermalFormulationMaterialProperties::getSatWaterEnthalpy(
    const double p, const double /*h*/
) const
{
    double const waterEnthalpy = 1e-7*(7.30984e+9 + 1.29239e+2*p - 1.00333e-6*std::pow(p, 2)+
        3.9881e-15*std::pow(p, 3) - 9.90697e+15*std::pow(p, -1) + 1.29267e+22*std::pow(p, -2) 
        - 6.28359e+27*std::pow(p, -3));
    return waterEnthalpy;
}


/*
* calculate the secondary variable based on a series of complementary condition
*/
bool GeothermalFormulationMaterialProperties::computeConstitutiveRelation(
    double const t,
    ProcessLib::SpatialPosition const& x_position,
    const int material_id,
    double const p,
    double const h,
    double& Sw,
    double& h_s, /*enthalpy of the steam*/
    double& h_w, /*enthalpy of the water*/
    double& dsw_dp,
    double& dsw_dh,
    double& dhs_dp,
    double& dhs_dh,
    double& dhw_dp,
    double& dhw_dh)
{
    {  // Local Newton solver
        using LocalJacobianMatrix =
            Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
        using LocalResidualVector = Eigen::Matrix<double, 3, 1>;
        using LocalUnknownVector = Eigen::Matrix<double, 3, 1>;
        LocalJacobianMatrix J_loc;

        Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(3);
        auto const update_residual = [&](LocalResidualVector& residual) {
            calculateResidual(material_id, p, h, Sw, h_s, h_w, residual);
        };

        auto const update_jacobian = [&](LocalJacobianMatrix& jacobian) {
            calculateJacobian(material_id, t, x_position, p, h, jacobian, Sw,
                              h_s,h_w);  // for solution dependent Jacobians
        };

        auto const update_solution = [&](LocalUnknownVector const& increment) {
            // increment solution vectors
            Sw += increment[0];
            h_s += increment[1];
            h_w += increment[2];
        };

        // TODO Make the following choice of maximum iterations and convergence
        // criteria available from the input file configuration. See Ehlers
        // material model implementation for the example.
        const int maximum_iterations(20);
        const double tolerance(1.e-14);

        auto newton_solver = NumLib::NewtonRaphson<
            decltype(linear_solver), LocalJacobianMatrix,
            decltype(update_jacobian), LocalResidualVector,
            decltype(update_residual), decltype(update_solution)>(
            linear_solver, update_jacobian, update_residual, update_solution,
            {maximum_iterations, tolerance});

        auto const success_iterations = newton_solver.solve(J_loc);

        if (!success_iterations)
            return false;
    }
    return true;
}
void GeothermalFormulationMaterialProperties::calculateResidual(
    const int material_id, double const p, double const h,
    double Sw, double enthalpy_steam, double enthalpy_water,
    ResidualVector& res)
{
    // const double pg =
    // pl + _capillary_pressure_models[material_id]->getCapillaryPressure(Sw);
    const double molar_density_water =
        GeothermalFormulationMaterialProperties::getDensityWater(p, h);
    const double molar_density_steam =
        GeothermalFormulationMaterialProperties::getDensitySteam(p, h);
    const double overall_molar_density = Sw*molar_density_water + (1 - Sw)*molar_density_steam;
    res(0) = h*overall_molar_density - 
        (Sw*molar_density_water*enthalpy_water + (1 - Sw)*molar_density_steam*enthalpy_steam);
        // calculating residual
    res(1) = std::min(Sw, enthalpy_steam - getSatSteamEnthalpy(p,h));
    res(2) = std::min(1-Sw, getSatWaterEnthalpy(p,h)-enthalpy_water);
}

void GeothermalFormulationMaterialProperties::calculateJacobian(
    const int material_id, double const /*t*/,
    ProcessLib::SpatialPosition const& /*x*/, double const p,
    double const h, JacobianMatrix& Jac, double Sw,
    double enthalpy_steam, double enthalpy_water)
{
    const double molar_density_water =
        GeothermalFormulationMaterialProperties::getDensityWater(p, h);
    const double molar_density_steam =
        GeothermalFormulationMaterialProperties::getDensitySteam(p, h);
    Jac.setZero();
    double const dF1dSw = h*(molar_density_water - molar_density_steam) -
        (molar_density_water*enthalpy_water - molar_density_steam*enthalpy_steam);
    double const dF1dhw = -molar_density_water*Sw;
    double const dF1dhs = -molar_density_steam*(1 - Sw);
    if (Sw > (enthalpy_steam - getSatSteamEnthalpy(p, h))) {
        double const dF2dSw = 0.0;
        double const dF2dhw = 0.0;
        double const dF2dhs = 1.0;
    }
    else {
        double const dF2dSw = 1.0;
        double const dF2dhw = 0.0;
        double const dF2dhs = 0.0;
    }
    if (1-Sw > (getSatWaterEnthalpy(p, h) - enthalpy_water)) {
        double const dF3dSw = 0.0;
        double const dF3dhw = -1.0;
        double const dF3dhs = 0.0;
    }
    else {
        double const dF3dSw = -1.0;
        double const dF3dhw = 0.0;
        double const dF3dhs = 0.0;
    }


}

}  // end of namespace
}  // end of namespace
