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
      curves(curves_),
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
      _rel_nonwet_perm_value(rel_nonwet_perm_value)
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
        {
            MathLib::PiecewiseLinearInterpolation const& interpolated_Pc =
                *curves.at("curve_PC_S");
            Sw = interpolated_Pc.getValue(pc);
            break;
        }
        /// van Genuchten
        case 1:
        {
            const double pb = _cap_pressure_value[0];
            const double slr = _cap_pressure_value[1];
            const double slm = 1.0;
            const double sgr = _cap_pressure_value[2];
            const double m =
                _cap_pressure_value[3];  // always <= 1.0.  Input is
                                         // exponent = 1 /
                                         // (1-lambda)
            assert(m <= 1);              //
            if (pc < 0.0)
                pc = 0.0;
            double effect_sw = std::pow(pc / pb, 1.0 / (1.0 - m)) + 1.0;
            effect_sw = std::pow(effect_sw, -m);
            Sw = effect_sw * (slm - slr) + slr;
            //
            break;
        }
        /// Brooks-Corey
        case 2:
        {
            const double pb = _cap_pressure_value[0];
            const double slr = _cap_pressure_value[1];
            const double slm = 1.0;
            const double sgr = _cap_pressure_value[2];
            const double lambda = _cap_pressure_value[3];  // always >= 1.0
            if (pc < pb)
                pc = pb;
            double const effect_sw = std::pow(pc / pb, -lambda);
            Sw = effect_sw * (slm - sgr - slr) + slr;
            break;
        }
        /// Liakopoulos
        case 10:
            Sw = 1 - (1.9722e-11) * std::pow(pc, 2.4279);
            if (pc < 0)
                // return 1 - (1.9722e-11)*std::pow(0.0, 2.4279);
                // extend
                Sw = 1 + pc * getDerivSaturation(0.0);
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
        {
            MathLib::PiecewiseLinearInterpolation const& interpolated_Pc =
                *curves.at("curve_PC_S");
            dSwdPc = interpolated_Pc.getDerivative(pc);
            if (pc > interpolated_Pc.getSupportMax())
                dSwdPc = interpolated_Pc.getDerivative(
                    interpolated_Pc.getSupportMax());
            else if (pc < interpolated_Pc.getSupportMin())
                dSwdPc = interpolated_Pc.getDerivative(
                    interpolated_Pc.getSupportMin());
            break;
        }

        case 1:
        {
            double const pb = _cap_pressure_value[0];
            double const slr = _cap_pressure_value[1];
            double const slm = 1.0;
            double const sgr = _cap_pressure_value[2];
            double const m = _cap_pressure_value[3];  // always <= 1.0.

            assert(m <= 1);
            // pc = MRange(FLT_EPSILON, pc, capillary_pressure_values[4]);

            double const v1 = std::pow((pc / pb), (1.0 / (1.0 - m)));
            double const v2 = std::pow((1.0 + v1), (-1.0 - m));
            dSwdPc = (m * v1 * v2 * (slm - slr)) / ((m - 1.0) * pc);
            break;
        }
        case 2:
        {
            double const pb = _cap_pressure_value[0];
            double const slr = _cap_pressure_value[1];
            double const slm = 1.0;
            double const sgr = _cap_pressure_value[2];
            const double lambda = _cap_pressure_value[3];  // always >= 1.0
                                                           //
            double const v1 = std::pow((pc / pb), -lambda);
            dSwdPc = (lambda * v1 * (slr - slm)) / pc;
            break;
        }
        case 10:
            dSwdPc = -(1.9722e-11) * 2.4279 * std::pow(pc, 1.4279);
            if (pc <= 0)
                dSwdPc = -3.7901e-7;
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
        {
            MathLib::PiecewiseLinearInterpolation const& interpolated_Kr =
                *curves.at("curve_S_Krel_wet");
            rel_wet_perm = interpolated_Kr.getValue(sw);
            break;
        }
        case 1:
            break;
        case 2:
        {
            double const slr = _rel_wet_perm_value[0];
            double const sgr = _rel_wet_perm_value[1];
            double const slm = 1.0;
            double const m = _rel_wet_perm_value[2];
            double const kr_min = _rel_wet_perm_value[3];
            double sl = sw;
            double const se = (sl - slr) / (slm - slr - sgr);
            //
            rel_wet_perm = std::pow(se, 3.0 + 2.0 / m);
            if (rel_wet_perm < kr_min)
                rel_wet_perm = kr_min;
            break;
        }
        case 10:
            rel_wet_perm = 1 - 2.207 * std::pow((1 - sw), 1.0121);  //
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
        {
            MathLib::PiecewiseLinearInterpolation const& interpolated_Kr =
                *curves.at("curve_S_Krel_nonwet");
            rel_nonwet_perm = interpolated_Kr.getValue(sw);
            break;
        }
        case 1:
            break;
        case 2:
        {
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
            break;
        }
        default:
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
