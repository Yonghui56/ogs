/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once


#include "TwoPhaseCarbonationMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
namespace TwoPhaseCarbonation
{
std::unique_ptr<TwoPhaseCarbonationMaterialProperties>
createTwoPhaseCarbonationMaterialProperties(
    BaseLib::ConfigTree const& config,
    boost::optional<MeshLib::PropertyVector<int> const&>
        material_ids);

}  // end namespace
}  // end namespace
