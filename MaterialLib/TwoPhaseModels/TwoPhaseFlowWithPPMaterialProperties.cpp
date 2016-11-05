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

namespace MaterialLib
{
namespace TwoPhaseFlowWithPP
{
TwoPhaseFlowWithPPMaterialProperties::TwoPhaseFlowWithPPMaterialProperties(
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
    : _has_material_ids(has_material_ids),
      _material_ids(material_ids),
      _liquid_density(std::move(liquid_density)),
      _viscosity(std::move(viscosity)),
      _gas_density(std::move(gas_density)),
      _gas_viscosity(std::move(gas_viscosity)),
      _intrinsic_permeability_models(intrinsic_permeability_models),
      _porosity_models(std::move(porosity_models)),
      _storage_models(std::move(storage_models)),
      _cap_pressure_model(cap_pressure_model),
      _rel_wet_perm_model(rel_wet_perm_model),
      _rel_nonwet_perm_model(rel_nonwet_perm_model),
      _cap_pressure_value(cap_pressure_value),
      _rel_wet_perm_value(rel_wet_perm_value),
      _rel_nonwet_perm_value(rel_nonwet_perm_value),
	  _curves(curves_)
{
    DBUG("Create material properties for Two-Phase flow with PP model.");
}

void TwoPhaseFlowWithPPMaterialProperties::setMaterialID(
    const ProcessLib::SpatialPosition& pos)
{
    if (!_has_material_ids)
    {
        _current_material_id = 0;
        return;
    }

    assert(pos.getElementID().get() < _material_ids.size());
    _current_material_id = _material_ids[pos.getElementID().get()];
}

double TwoPhaseFlowWithPPMaterialProperties::getLiquidDensity(
    const double p, const double T) const
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

double TwoPhaseFlowWithPPMaterialProperties::getDerivGasDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pg)] = p;

    return _gas_density->getdValue(
        vars, MaterialLib::Fluid::PropertyVariableType::pg);
}
double TwoPhaseFlowWithPPMaterialProperties::getLiquidViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pl)] = p;
    return _viscosity->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getGasViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::pg)] = p;
    return _gas_viscosity->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getSaturation(double pc) const
{
    double Sw;
    /// TODO waiting for a better way to implemente the PC-S curve
    switch (_cap_pressure_model)
    {
        case 0:
			OGS_FATAL("This model has been implemented elsewhere");
			break;
        /// van Genuchten
        case 1:
			OGS_FATAL("This model has been implemented elsewhere");
			break;
        /// Brooks-Corey
        case 2:
			OGS_FATAL("This model has been implemented elsewhere");
			break;
        /// Liakopoulos
        case 10:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        default:
            OGS_FATAL("This model has not been implemented yet");
            break;
    }
    return Sw;
}

double TwoPhaseFlowWithPPMaterialProperties::getDerivSaturation(
    double const pc) const
{
    double dSwdPc;
    switch (_cap_pressure_model)
    {
        case 0:
			OGS_FATAL("This model has been implemented elsewhere");
			break;

        case 1:
			OGS_FATAL("This model has been implemented elsewhere");
			break;
        case 2:
			OGS_FATAL("This model has been implemented elsewhere");
			break;
        case 10:
			OGS_FATAL("This model has been implemented elsewhere");
			break;
        default:
            OGS_FATAL("This model has not been implemented yet");
            break;
    }
    return dSwdPc;
}

double TwoPhaseFlowWithPPMaterialProperties::getrelativePermeability_liquid(
    double const sw) const
{
    /// TODO waiting for a better way to implemente the Kr-S curve
    /*
    MathLib::PiecewiseLinearInterpolation const& interpolated_Kr =
        *curves.at("curveB");
    return interpolated_Kr.getValue(sw);
    */
    double rel_wet_perm;
    switch (_rel_wet_perm_model)
    {
        case 0:
			OGS_FATAL("This model has been implemented elsewhere");
			break;
        case 1:
            break;
        case 2:
			OGS_FATAL("This model has been implemented elsewhere");
			break;
        case 10:
            OGS_FATAL("This model has been implemented elsewhere");
            break;
        default:
            OGS_FATAL("This model has not been implemented yet");
            break;
    }
    return rel_wet_perm;
}

double TwoPhaseFlowWithPPMaterialProperties::getrelativePermeability_gas(
    double sw) const
{
    double rel_nonwet_perm;
    switch (_rel_nonwet_perm_model)
    {
        case 0:
			OGS_FATAL("This model has been implemented elsewhere");
			break;
        case 1:
            break;
        case 2:
			OGS_FATAL("This model has been implemented elsewhere");
			break;
		case 10:
			OGS_FATAL("This model has been implemented elsewhere");
			break;
        default:
			OGS_FATAL("This model has not been implemented yet");
            break;
    }
    return rel_nonwet_perm;
}

Eigen::MatrixXd const& TwoPhaseFlowWithPPMaterialProperties::getPermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const int /*dim*/) const
{
    return _intrinsic_permeability_models[_current_material_id];
}

double TwoPhaseFlowWithPPMaterialProperties::getPorosity(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double p, const double T, const double porosity_variable) const
{
    const double porosity =
        _porosity_models[_current_material_id]->getValue(porosity_variable, T);

    return porosity;
}

double TwoPhaseFlowWithPPMaterialProperties::getDissolvedGas(
    double const pg) const
{
    double const hen = 2e-6;     //
    double const M_air = 0.029;  // unit kg/mol
    return pg * hen * M_air;
}

}  // end of namespace
}  // end of namespace
