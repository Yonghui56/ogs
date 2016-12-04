/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_CREATETWOPHASEFLOWPRHOMATERIALPROPERTIES_H
#define OGS_CREATETWOPHASEFLOWPRHOMATERIALPROPERTIES_H

#include <memory>
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "ProcessLib/TwoPhaseFlowWithPrho/TwoPhaseFlowWithPrhoMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
std::unique_ptr<TwoPhaseFlowWithPrhoMaterialProperties>
CreateTwoPhaseFlowPrhoMaterialProperties(
    BaseLib::ConfigTree const& config,
    bool const has_material_ids,
    MeshLib::PropertyVector<int> const& material_ids);

}  // end namespace
}  // end namespace

#endif /* CREATETWOPHASEFLOWPRHOMATERIALPROPERTIES_H */
