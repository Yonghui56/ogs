/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/PorousMediaProperties.h"

namespace ProcessLib
{
template <typename ReturnType>
struct Parameter;

namespace FinesTransport
{
struct FinesTransportMaterialProperties final
{
    FinesTransportMaterialProperties(
        std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>&&
            media_map_,
        MaterialLib::PorousMedium::PorousMediaProperties&&
        porous_media_properties_,
        bool const has_fluid_thermal_expansion_,
        ParameterLib::Parameter<double> const& solid_intrinsic_permeability_,
        ParameterLib::Parameter<double> const& porosity_constant_,
        Eigen::VectorXd specific_body_force_,
        bool const has_gravity_)
        : media_map(std::move(media_map_)),
          porous_media_properties(std::move(porous_media_properties_)),
          has_fluid_thermal_expansion(has_fluid_thermal_expansion_),
          solid_intrinsic_permeability(solid_intrinsic_permeability_),
          porosity_constant(porosity_constant_),
          specific_body_force(std::move(specific_body_force_)),
          has_gravity(has_gravity_)
    {
    }

    FinesTransportMaterialProperties(FinesTransportMaterialProperties&&) = delete;
    FinesTransportMaterialProperties(FinesTransportMaterialProperties const&) = delete;
    void operator=(FinesTransportMaterialProperties&&) = delete;
    void operator=(FinesTransportMaterialProperties const&) = delete;

    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    MaterialLib::PorousMedium::PorousMediaProperties porous_media_properties;
    bool const has_fluid_thermal_expansion;
    ParameterLib::Parameter<double> const& solid_intrinsic_permeability;//intrinsic_permeability
    ParameterLib::Parameter<double> const& porosity_constant;//porosity

    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
    double dt = 0.0;
};

}  // namespace FinesTransport
}  // namespace ProcessLib
