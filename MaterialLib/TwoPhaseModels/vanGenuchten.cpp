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

#include "vanGenuchten.h"

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
vanGenuchten::vanGenuchten(
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
    DBUG("vanGenuchten.");
}

double vanGenuchten::getSaturation(double pc) const
{
	double Sg;
	double const Pb_van_Genuchten = _cap_pressure_value[0];
	double const S_lr = _cap_pressure_value[1];
	double const S_gr = _cap_pressure_value[2];
	double const n_van_Genuchten= _cap_pressure_value[3];
	double const m = 1.0 - 1.0 / n_van_Genuchten;
	if (pc > 0) {
		Sg = 1 - (((1 - S_gr - S_lr) / pow((pow((pc / Pb_van_Genuchten), n_van_Genuchten) + 1), m)) + S_lr);
	}
	else
		Sg = 0.0;//just for test
	return 1 - Sg;
}

double vanGenuchten::getDerivSaturation(double const pc) const
{
	double dPC;
	double const Pb_van_Genuchten = _cap_pressure_value[0];
	double const S_lr = _cap_pressure_value[1];
	double const S_gr = _cap_pressure_value[2];
	double const n_van_Genuchten = _cap_pressure_value[3];
	double const m = 1.0 - 1.0 / n_van_Genuchten;
	dPC = m*n_van_Genuchten*
		(1 - S_gr - S_lr)*(1 / Pb_van_Genuchten)*(pow((pc / Pb_van_Genuchten), (n_van_Genuchten - 1)))*pow(((pow((pc / Pb_van_Genuchten), n_van_Genuchten)) + 1), (-m - 1));
	if (pc <= 0)
	{
		dPC = 0.0;//just for test
	}
	return -dPC;
}

double vanGenuchten::getrelativePermeability_liquid(
    double const sw) const
{
	double Kr_L;
	double const S_lr = _rel_wet_perm_value[0];
	double const S_gr = _rel_wet_perm_value[1];
	double const n_van_Genuchten = _rel_wet_perm_value[2];
	double const m = 1.0 - 1.0 / n_van_Genuchten;
	double const kr_min = _rel_wet_perm_value[3];
	double EffectSat_l = (sw - S_lr) / (1 - S_gr - S_lr);
	Kr_L = sqrt(EffectSat_l)*pow(1 - pow(1 - pow(EffectSat_l, 1 / m), m), 2);
	if (sw < 0)
		Kr_L = kr_min;
	else if (sw > 1)
		Kr_L = 1;
	return Kr_L;
}


double vanGenuchten::getrelativePermeability_gas(double const sw) const
{
	double Kr_G;
	double const S_lr = _rel_nonwet_perm_value[0];
	double const S_gr = _rel_nonwet_perm_value[1];
	double const n_van_Genuchten = _rel_nonwet_perm_value[2];
	double const m = 1.0 - 1.0 / n_van_Genuchten;
	double const kr_min = _rel_wet_perm_value[3];
	double EffectSat_g = (1 - sw - S_gr) / (1 - S_gr - S_lr);
	Kr_G = sqrt(EffectSat_g)*pow(1 - pow(1 - EffectSat_g, 1 / m), 2 * m);
	if (sw < 0)
		Kr_G = 1;
	else if (sw > 1)
		Kr_G = kr_min;
	return Kr_G;;
}
}  // end of namespace
}  // end of namespace
