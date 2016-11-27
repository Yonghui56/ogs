/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_TWOPHASEFLOWWITHPX_CREATEEOSODEALMIX_H_
#define PROCESSLIB_TWOPHASEFLOWWITHPX_CREATEEOSODEALMIX_H_

#include <logog/include/logog.hpp>

#include "EoS_IdealMix.h"
#include "EoSBase.h"
#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter

namespace ProcessLib
{
namespace TwoPhaseFlowWithPX
{
std::unique_ptr<EoSBase> CreateEoS_IdealMix(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
	//! \ogs_file_param{process__SMALL_DEFORMATION__constitutive_relation__type}
	config.checkConfigParameter("type", "EoS_IdealMix");
	DBUG("Create EoS_IdealMix model");

	// Kelvin shear modulus.
	auto& henry_constant = ProcessLib::findParameter<double>(
		//! \ogs_file_param_special{process__SMALL_DEFORMATION__constitutive_relation__Lubby2__kelvin_shear_modulus}
		config, "henry_constant", parameters, 1);

	DBUG("Use '%s' as kelvin shear modulus parameter.",
		henry_constant.name.c_str());
    typename EoS_IdealMix::EoS_IdealMix_Properties eos_mp{
		henry_constant};

    return std::unique_ptr<EoSBase>{
        new EoS_IdealMix{ eos_mp }};
}

}  // namespace Solids
}  // namespace MaterialLib

#endif  // MATERIALLIB_SOLIDMODELS_CREATELUBBY2_H_
