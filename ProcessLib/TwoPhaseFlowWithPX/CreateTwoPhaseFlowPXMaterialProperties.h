/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_CREATETWOPHASEFLOWPXMATERIALPROPERTIES_H
#define OGS_CREATETWOPHASEFLOWPXMATERIALPROPERTIES_H

#include <memory>
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "ProcessLib/TwoPhaseFlowWithPX/TwoPhaseFlowWithPXMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
namespace TwoPhaseFlowWithPX
{
std::unique_ptr<TwoPhaseFlowWithPXMaterialProperties>
CreateTwoPhaseFlowPXMaterialProperties(
    BaseLib::ConfigTree const& config,
    bool const has_material_ids,
    MeshLib::PropertyVector<int> const& material_ids);

}  // end namespace
}  // end namespace

#endif /* CREATETWOPHASEFLOWPXMATERIALPROPERTIES_H */
