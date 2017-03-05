/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include <memory>
#include <vector>
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
class SpatialPosition;
namespace ThermalTwoPhaseFlowWithPP
{
class ThermalTwoPhaseFlowWithPPMaterialProperties
{
public:
    static int const JacobianResidualSize = 2;
    using ResidualVector = Eigen::Matrix<double, JacobianResidualSize, 1>;
    using JacobianMatrix = Eigen::Matrix<double, JacobianResidualSize,
                                         JacobianResidualSize, Eigen::RowMajor>;
    using UnknownVector = Eigen::Matrix<double, JacobianResidualSize, 1>;

public:
    using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

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
            relative_permeability_models);

    int getMaterialID(const std::size_t element_id);

    Eigen::MatrixXd const& getPermeability(const int material_id,
        const double t,
        const ProcessLib::SpatialPosition& pos,
        const int dim) const;

    double getPorosity(const int material_id,
        const double t, const ProcessLib::SpatialPosition& pos,
                       const double p, const double T,
                       const double porosity_variable) const;

    double getNonwetRelativePermeability(const double t,
                                         const ProcessLib::SpatialPosition& pos,
                                         const double p, const double T,
                                         const double saturation) const;
    double getWetRelativePermeability(const double t,
                                      const ProcessLib::SpatialPosition& pos,
                                      const double p, const double T,
                                      const double saturation) const;
    double getSaturation(const int material_id,
        const double t, const ProcessLib::SpatialPosition& pos,
                         const double p, const double T, const double pc) const;
    double getSaturationDerivative(const int material_id, const double t,
                                   const ProcessLib::SpatialPosition& pos,
                                   const double p, const double T,
                                   const double saturation) const;
    double getCapillaryPressure(const int material_id, const double t,
                                const ProcessLib::SpatialPosition& pos,
                                const double p, const double T,
                                const double saturation) const;
    double getLiquidDensity(const double p, const double T) const;
    double getGasDensity(const double p, const double T) const;
    double getGasViscosity(const double p, const double T) const;
    double getLiquidViscosity(const double p, const double T) const;
    double getSpecificHeatCapacitySolid(const double p, const double T) const;
    double getSpecificHeatCapacityWater(const double p, const double T) const;
    double getSpecificHeatCapacityAir(const double p, const double T) const;
    double getSpecificHeatCapacityVapor(const double p, const double T) const;
    double getThermalConductivityDrySolid(const double p, const double T) const;
    double getThermalConductivityWetSolid(const double p, const double T) const;
    double getDerivativeGasDensity(double const p, double const T) const;
    /// Calculates the unsaturated heat conductivity
    double calculateUnsatHeatConductivity(double const t,
                                          ProcessLib::SpatialPosition const& x,
                                          double const Sw,
                                          double const lambda_pm_dry,
                                          double const lambda_pm_wet) const;
    /// water vapor saturation pressure
    double calculateSaturatedVaporPressure(const double T) const;
    /// partial water vapor pressure in nonwetting phase
    /// Kelvin equation
    double calculateVaporPressureNonwet(const double pc, const double T, const double rho_mass_h2o) const;
    /// Derivative of SaturatedVaporPressure in terms of T
    double calculateDerivativedPsatdT(const double T) const;
    /// Derivative of partial vapor pressure in terms of T
    double calculateDerivativedPgwdT(const double pc, const double T, const double rho_mass_h2o) const;
    /// Derivative of partial vapor pressure in terms of PC
    double calculateDerivativedPgwdPC(const double pc, const double T, const double rho_mass_h2o) const;
    ///
    double calculatedRhoNonwetdT(const double p_air_nonwet,
                                 const double p_vapor_nonwetconst, double pc,
                                 const double T, const double rho_mass_h2o) const;

protected:

    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _liquid_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _viscosity;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _gas_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _gas_viscosity;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _specific_heat_capacity_solid;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _specific_heat_capacity_water;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _specific_heat_capacity_air;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _specific_heat_capacity_vapor;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _thermal_conductivity_dry_solid;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _thermal_conductivity_wet_solid;

    /** Use two phase models for different material zones.
    *  Material IDs must be given as mesh element properties.
    */
    boost::optional<MeshLib::PropertyVector<int> const&> const _material_ids;

    std::vector<Eigen::MatrixXd> _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _storage_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
        _capillary_pressure_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>
        _relative_permeability_models;

};

}  // end of namespace
}  // end of namespace