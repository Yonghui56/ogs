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

#include "TwoPhaseFlowWithPPMaterialProperties.h"

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

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
TwoPhaseFlowWithPPMaterialProperties::TwoPhaseFlowWithPPMaterialProperties(
    BaseLib::ConfigTree const& config,
    MeshLib::PropertyVector<int> const& material_ids,
	std::map<std::string,
	std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
	curves_)
    : _material_ids(material_ids),
	curves(curves_)
{
    DBUG("Reading material properties of two-phase flow process.");

    //! \ogs_file_param{prj__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    // Get fluid properties
    //! \ogs_file_param{prj__material_property__fluid__density}
    auto const& rho_conf = fluid_config.getConfigSubtree("density");
    _liquid_density = MaterialLib::Fluid::createFluidDensityModel(rho_conf);
	auto const& rho_gas_conf = fluid_config.getConfigSubtree("gasdensity");
	_gas_density = MaterialLib::Fluid::createFluidDensityModel(rho_gas_conf);
	auto const& rho_dissolve_gas_conf = fluid_config.getConfigSubtree("dissolvegasdensity");
	_dissolve_gas_rho = MaterialLib::Fluid::createFluidDensityModel(rho_dissolve_gas_conf);
    //! \ogs_file_param{prj__material_property__fluid__viscosity}
    auto const& mu_conf = fluid_config.getConfigSubtree("viscosity");
    _viscosity = MaterialLib::Fluid::createViscosityModel(mu_conf);
	//! \ogs_file_param{prj__material_property__fluid__gas__viscosity}
	auto const& mu_gas_conf = fluid_config.getConfigSubtree("gasviscosity");
	_gas_viscosity = MaterialLib::Fluid::createViscosityModel(mu_gas_conf);
    // Get porous properties
    std::vector<int> mat_ids;
    //! \ogs_file_param{prj__material_property__porous_medium}
    auto const& poro_config = config.getConfigSubtree("porous_medium");
    //! \ogs_file_param{prj__material_property__porous_medium__porous_medium}
    for (auto const& conf : poro_config.getConfigSubtreeList("porous_medium"))
    {
        //! \ogs_file_attr{prj__material_property__porous_medium__porous_medium__id}
        auto const id = conf.getConfigAttributeOptional<int>("id");
        mat_ids.push_back(*id);

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__permeability}
        auto const& perm_conf = conf.getConfigSubtree("permeability");
        _intrinsic_permeability_models.emplace_back(
            MaterialLib::PorousMedium::createPermeabilityModel(perm_conf));

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__porosity}
        auto const& poro_conf = conf.getConfigSubtree("porosity");
        auto n = MaterialLib::PorousMedium::createPorosityModel(poro_conf);
        _porosity_models.emplace_back(std::move(n));

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__storage}
        auto const& stora_conf = conf.getConfigSubtree("storage");
        auto beta = MaterialLib::PorousMedium::createStorageModel(stora_conf);
        _storage_models.emplace_back(std::move(beta));
    }

    BaseLib::reorderVector(_intrinsic_permeability_models, mat_ids);
    BaseLib::reorderVector(_porosity_models, mat_ids);
    BaseLib::reorderVector(_storage_models, mat_ids);
}

void TwoPhaseFlowWithPPMaterialProperties::setMaterialID(const SpatialPosition& pos)
{
    if (_material_ids.empty())
    {
        assert(pos.getElementID().get() < _material_ids.size());
        _current_material_id = _material_ids[pos.getElementID().get()];
    }
}

double TwoPhaseFlowWithPPMaterialProperties::getLiquidDensity(const double p,
                                                      const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pl)] = p;
    return _liquid_density->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getGasDensity(const double p,
	const double T) const
{
	ArrayType vars;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pg)] = p;
	return _gas_density->getValue(vars);

}

double TwoPhaseFlowWithPPMaterialProperties::getDerivGasDensity(const double p,
	const double T) const
{
	ArrayType vars;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pg)] = p;
	
	return _gas_density->getdValue(vars, MaterialLib::Fluid::PropertyVariableType::pg);

}
double TwoPhaseFlowWithPPMaterialProperties::getLiquidViscosity(const double p,
                                                  const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pl)] = p;
    return _viscosity->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getGasViscosity(const double p,
	const double T) const
{
	ArrayType vars;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pg)] = p;
	return _gas_viscosity->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getSaturation(double pc) const
{
	/*// read from curve
	MathLib::PiecewiseLinearInterpolation const& interpolated_Pc =
		*curves.at("curveA");
	return interpolated_Pc.getValue(pc);
	*/
	/* for liakopoulos test
	if (pc<0)
		return 1 - (1.9722e-11)*pow(0.0, 2.4279);
	return 1 - (1.9722e-11)*pow(pc, 2.4279);
	*/
	/*
	 * Brooks Corey model
	 */
	double const pb = 5000;
	double const slr = 0.0;
	double const slm = 1.0;
	double const m = 2.0; // always >= 1.0
	if (pc < pb)
		pc = pb;
	double const se = pow(pc / pb, -m);
	double sl = se * (slm - slr) + slr;
	sl =MRange(slr + DBL_EPSILON, sl, slm - DBL_EPSILON);
	return sl;
}

double TwoPhaseFlowWithPPMaterialProperties::getrelativePermeability_liquid(double const sw) const
{   
	/*
	 * read from curve
	 */
	/*MathLib::PiecewiseLinearInterpolation const& interpolated_Kr =
		*curves.at("curveB");
	return interpolated_Kr.getValue(sw);*/
	/*
	* Brooks Corey model
	*/
	double const slr = 0.0;
	double const slm = 1.0;
	double const m = 2.0;
	double sl = sw;
	sl= this->MRange(slr, sl, slm);
	double const se = (sl - slr) / (slm - slr);
	//
	double kr = pow(se, 3.0 + 2.0 / m);
	if (kr < 1e-9)
		kr = 1e-9;
	return kr;
}

double TwoPhaseFlowWithPPMaterialProperties::getDerivSaturation(double pc) const
{
	/*// read from curve
	MathLib::PiecewiseLinearInterpolation const& interpolated_Pc =
		*curves.at("curveA");
	double dSwdPc = interpolated_Pc.getDerivative(pc);
	if (pc > interpolated_Pc.getSupportMax())
		dSwdPc = interpolated_Pc.getDerivative(
			interpolated_Pc.getSupportMax());
	else if (pc < interpolated_Pc.getSupportMin())
		dSwdPc = interpolated_Pc.getDerivative(
			interpolated_Pc.getSupportMin());
	return dSwdPc;
	*/
	/* // for Liakopoulos
	if (pc < 0)
		return -(1.9722e-11)*2.4279*pow(0.0, 1.4279);
	return -(1.9722e-11)*2.4279*pow(pc, 1.4279);
	*/
	// for Brooks Corey
	double const pb = 5000;
	double const slr = 0.0;
	double const slm = 1.0;
	double const m = 2.0; // always >= 1.0
	if (pc < pb)
		pc = pb;
	double const se = pow(pc / pb, -m);
	double sl = se * (slm - slr) + slr;
	sl = MRange(slr + DBL_EPSILON, sl, slm - DBL_EPSILON);
	
	double v1 = pow(((sl - slr) / (slm - slr)), (-1.0 / m));
	double const dpds = (pb * v1) / (m * (slr - sl));
	return 1/dpds;
}

double TwoPhaseFlowWithPPMaterialProperties::getrelativePermeability_gas(double sw) const
{
	/*for liakopoulos
	double k_rG = 0.0;
	double k_min = 1e-5;
	double Lambda_Brook = 3.;
	double S_gr = 0.;
	double S_lr = 0.2;
	double S_le = (sw - S_lr) / (1 - S_lr - S_gr);
	k_rG = pow(1.0 - S_le, 2)*(1.0 - pow(S_le, 1.0 + 2.0 / Lambda_Brook));
	if (k_rG < k_min)
		return k_min;
	return k_rG;
	*/
	double const slr = 0.0; // slr = 1.0 - sgm
	double const slm = 1.0; // slm = 1.0 - sgr
	double const m = 2;
	double sl = sw;
	sl = MRange(slr, sl, slm);
	double se = (sl - slr) / (slm - slr);
	//
	double kr = pow(1.0 - se, 2) * (1.0 - pow(se, 1.0 + 2.0 / m));
	if (kr < 1e-9)
		kr = 1e-9;
	return kr;
}


Eigen::MatrixXd const& TwoPhaseFlowWithPPMaterialProperties::getPermeability(
    const double /*t*/, const SpatialPosition& /*pos*/, const int /*dim*/) const
{
    return _intrinsic_permeability_models[_current_material_id];
}

double TwoPhaseFlowWithPPMaterialProperties::getPorosity(
	const double /*t*/, const SpatialPosition& /*pos*/, const double p,
	const double T, const double porosity_variable) const
{

	const double porosity =
		_porosity_models[_current_material_id]->getValue(porosity_variable, T);
	
	return porosity;
}

double TwoPhaseFlowWithPPMaterialProperties::getDissolvedGas(double const pg, double const T, double const molg) const
{
	ArrayType vars;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pg)] = pg;
	vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::molarg)] = molg;
	return _dissolve_gas_rho->getValue(vars);
}







}  // end of namespace
}  // end of namespace
