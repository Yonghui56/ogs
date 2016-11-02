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
        curves;

    int _current_material_id = 0;
    std::vector<Eigen::MatrixXd> _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _storage_models;

    // Note: For the statistical data of porous media, they could be read from
    // vtu files directly. This can be done by using property vectors directly.
    // Such property vectors will be added here if they are needed.
};

}  // end of namespace
}  // end of namespace
#endif /* TWOPHASEFLOWWITHPPMATERIALPROPERTIES_H */
