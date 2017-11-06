/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once


#include "GeothermalFormulationMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
namespace GeothermalFormulation
{
std::unique_ptr<GeothermalFormulationMaterialProperties>
createGeothermalFormulationMaterialProperties(
    BaseLib::ConfigTree const& config,
    boost::optional<MeshLib::PropertyVector<int> const&>
        material_ids);

}  // end namespace
}  // end namespace
