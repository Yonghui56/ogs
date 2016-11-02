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

#ifndef OGS_TWOPHASEFLOWWITHPPMATERIALPROPERTIES_H_L
#define OGS_TWOPHASEFLOWWITHPPMATERIALPROPERTIES_H_L

#include <iostream>
#include <memory>
#include <vector>
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

#include "TwoPhaseFlowWithPPMaterialProperties.h"

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
class Liakopoulos : public TwoPhaseFlowWithPPMaterialProperties
{
public:
    typedef MaterialLib::Fluid::FluidProperty::ArrayType ArrayType;

    Liakopoulos(
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

    double getSaturation(double pc) const override;
    double getrelativePermeability_liquid(double sw) const override;
    double getrelativePermeability_gas(double const sw) const override;
};

}  // end of namespace
}  // end of namespace
#endif /* TWOPHASEFLOWWITHPPMATERIALPROPERTIES_H_L */
