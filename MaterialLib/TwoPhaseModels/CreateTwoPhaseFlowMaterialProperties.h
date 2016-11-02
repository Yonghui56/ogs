/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   createPorosityModel.h
 *
 * Created on August 16, 2016, 1:16 PM
 */

#ifndef OGS_CREATETWOPHASEFLOWMATERIALPROPERTIES_H
#define OGS_CREATETWOPHASEFLOWMATERIALPROPERTIES_H

#include <memory>
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/TwoPhaseModels/TwoPhaseFlowWithPPMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace TwoPhaseFlowWithPP
{
/** Create a porosity model
 *  @param config  ConfigTree object has a tag of <porosity>
 */
std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties>
CreateTwoPhaseFlowMaterialProperties(
    BaseLib::ConfigTree const& config,
    bool const has_material_ids,
    MeshLib::PropertyVector<int> const& material_ids,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves_);

}  // end namespace
}  // end namespace

#endif /* CREATEPOROSITYMODEL_H */