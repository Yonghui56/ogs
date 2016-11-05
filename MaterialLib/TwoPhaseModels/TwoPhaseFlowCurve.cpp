/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TwoPhaseFlowWithPPMaterialProperties.cpp
 *
 * Created on August 18, 2016, 11:49 AM
 */

#include "TwoPhaseFlowCurve.h"

#include <logog/include/logog.hpp>

#include "BaseLib/reorderVector.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

namespace MaterialLib
{
namespace TwoPhaseFlowWithPP
{
TwoPhaseFlowCurve::TwoPhaseFlowCurve(
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
        curves_)
    : TwoPhaseFlowWithPPMaterialProperties(has_material_ids,
                                           material_ids,
                                           std::move(liquid_density),
                                           std::move(viscosity),
                                           std::move(gas_density),
                                           std::move(gas_viscosity),
                                           intrinsic_permeability_models,
                                           std::move(porosity_models),
                                           std::move(storage_models),
                                           cap_pressure_model,
                                           rel_wet_perm_model,
                                           rel_nonwet_perm_model,
                                           cap_pressure_value,
                                           rel_wet_perm_value,
                                           rel_nonwet_perm_value,
                                           curves_)
{
    DBUG("TwoPhaseFlowCurve.");
}

double TwoPhaseFlowCurve::getSaturation(double pc) const
{
    double Sw;
    /// TODO waiting for a better way to implemente the PC-S curve
    assert(_cap_pressure_model == 0);
	MathLib::PiecewiseLinearInterpolation const& interpolated_Pc =
		*_curves.at("curve_PC_S");
	Sw = interpolated_Pc.getValue(pc);
	return Sw;
}

double TwoPhaseFlowCurve::getDerivSaturation(double const pc) const
{
	double dSwdPc;
	MathLib::PiecewiseLinearInterpolation const& interpolated_Pc =
		*_curves.at("curve_PC_S");
	dSwdPc = interpolated_Pc.getDerivative(pc);
	if (pc > interpolated_Pc.getSupportMax())
		dSwdPc = interpolated_Pc.getDerivative(
			interpolated_Pc.getSupportMax());
	else if (pc < interpolated_Pc.getSupportMin())
		dSwdPc = interpolated_Pc.getDerivative(
			interpolated_Pc.getSupportMin());
	return dSwdPc;
}

double TwoPhaseFlowCurve::getrelativePermeability_liquid(
    double const sw) const
{
	double rel_wet_perm;
	MathLib::PiecewiseLinearInterpolation const& interpolated_Kr =
		*_curves.at("curve_S_Krel_wet");
	rel_wet_perm = interpolated_Kr.getValue(sw);
	return rel_wet_perm;
}

double TwoPhaseFlowCurve::getrelativePermeability_gas(double const sw) const
{
    double rel_nonwet_perm;
    double const slr = _rel_nonwet_perm_value[0];
    double const sgr = _rel_nonwet_perm_value[1];
    double const slm = 1.0;
    double const m = _rel_nonwet_perm_value[2];
    double const kr_min = _rel_nonwet_perm_value[3];
    double S_le = (sw - slr) / (slm - slr - sgr);
    rel_nonwet_perm =
        std::pow(1.0 - S_le, 2) * (1.0 - std::pow(S_le, 1.0 + 2.0 / m));
    if (rel_nonwet_perm < kr_min)
        rel_nonwet_perm = kr_min;

    return rel_nonwet_perm;
}
}  // end of namespace
}  // end of namespace
