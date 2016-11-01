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
    BaseLib::ConfigTree const& config,
    bool const has_material_ids,
    MeshLib::PropertyVector<int> const& material_ids,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves_)
    : _has_material_ids(has_material_ids),
      _material_ids(material_ids),
      curves(curves_)
{
    DBUG("Reading material properties of two-phase flow process.");

    //! \ogs_file_param{prj__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    // Get fluid properties
    //! \ogs_file_param{prj__material_property__fluid__density}
    auto const& rho_conf = fluid_config.getConfigSubtree("liquiddensity");
    _liquid_density = MaterialLib::Fluid::createFluidDensityModel(rho_conf);
    auto const& rho_gas_conf = fluid_config.getConfigSubtree("gasdensity");
    _gas_density = MaterialLib::Fluid::createFluidDensityModel(rho_gas_conf);
    //! \ogs_file_param{prj__material_property__fluid__viscosity}
    auto const& mu_conf = fluid_config.getConfigSubtree("liquidviscosity");
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

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__cap_pressure}
        auto const& cap_pressure_conf =
            conf.getConfigSubtree("capillary_pressure");
        auto const cap_pressure_type =
            cap_pressure_conf.getConfigParameter<std::string>("type");
        if (cap_pressure_type == "Curve")
            cap_pressure_model = 0;
        else if (cap_pressure_type == "van_Genuchten")
            cap_pressure_model = 1;
        else if (cap_pressure_type == "Brooks_Corey")
        {
            cap_pressure_model = 2;
            cap_pressure_value[0] =
                cap_pressure_conf.getConfigParameter<double>("entry_pressure");
            cap_pressure_value[1] =
                cap_pressure_conf.getConfigParameter<double>(
                    "res_saturation_wet");
            cap_pressure_value[2] =
                cap_pressure_conf.getConfigParameter<double>(
                    "res_saturation_nonwet");
            cap_pressure_value[3] =
                cap_pressure_conf.getConfigParameter<double>("lambda");
        }
        else if (cap_pressure_type == "Liakopoulos")
            cap_pressure_model = 10;
        else
            OGS_FATAL("This model has not been implemented yet");

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__rel_permeability_liquid}
        auto const& rel_wet_perm_conf =
            conf.getConfigSubtree("relative_wet_permeability");
        auto const rel_wet_perm_type =
            rel_wet_perm_conf.getConfigParameter<std::string>("type");
        if (rel_wet_perm_type == "Curve")
            rel_wet_perm_model = 0;
        else if (rel_wet_perm_type == "van_Genuchten")
            rel_wet_perm_model = 1;
        else if (rel_wet_perm_type == "Brooks_Corey")
        {
            rel_wet_perm_model = 2;
            rel_wet_perm_value[0] =
                rel_wet_perm_conf.getConfigParameter<double>(
                    "res_saturation_wet");
            rel_wet_perm_value[1] =
                rel_wet_perm_conf.getConfigParameter<double>(
                    "res_saturation_nonwet");
            rel_wet_perm_value[2] =
                rel_wet_perm_conf.getConfigParameter<double>("lambda");
            rel_wet_perm_value[3] =
                rel_wet_perm_conf.getConfigParameter<double>("minimum_value");
        }
        else if (rel_wet_perm_type == "Liakopoulos")
            rel_wet_perm_model = 10;
        else
            OGS_FATAL("This model has not been implemented yet");

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__rel_permeability_gas}
        auto const& rel_nonwet_perm_conf =
            conf.getConfigSubtree("relative_nonwet_permeability");
        auto const rel_nonwet_perm_type =
            rel_nonwet_perm_conf.getConfigParameter<std::string>("type");
        if (rel_nonwet_perm_type == "Curve")
            rel_nonwet_perm_model = 0;
        else if (rel_nonwet_perm_type == "van_Genuchten")
            rel_nonwet_perm_model = 1;
        else if (rel_nonwet_perm_type == "Brooks_Corey")
        {
            rel_nonwet_perm_model = 2;
            rel_nonwet_perm_value[0] =
                rel_nonwet_perm_conf.getConfigParameter<double>(
                    "res_saturation_wet");
            rel_nonwet_perm_value[1] =
                rel_nonwet_perm_conf.getConfigParameter<double>(
                    "res_saturation_nonwet");
            rel_nonwet_perm_value[2] =
                rel_nonwet_perm_conf.getConfigParameter<double>("lambda");
            rel_nonwet_perm_value[3] =
                rel_nonwet_perm_conf.getConfigParameter<double>(
                    "minimum_value");
        }
        else
            OGS_FATAL("This model has not been implemented yet");
    }

    BaseLib::reorderVector(_intrinsic_permeability_models, mat_ids);
    BaseLib::reorderVector(_porosity_models, mat_ids);
    BaseLib::reorderVector(_storage_models, mat_ids);
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
    switch (cap_pressure_model)
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
            const double pb = cap_pressure_value[0];
            const double slr = cap_pressure_value[1];
            const double slm = 1.0;
            const double sgr = cap_pressure_value[2];
            const double m = cap_pressure_value[3];  // always <= 1.0.  Input is
                                                     // exponent = 1 /
                                                     // (1-lambda)
            assert(m <= 1);                          //
            if (pc < 0.0)
                pc = 0.0;
            double effect_sw = pow(pc / pb, 1.0 / (1.0 - m)) + 1.0;
            effect_sw = pow(effect_sw, -m);
            Sw = effect_sw * (slm - slr) + slr;
            //
            break;
        }
        /// Brooks-Corey
        case 2:
        {
            const double pb = cap_pressure_value[0];
            const double slr = cap_pressure_value[1];
            const double slm = 1.0;
            const double sgr = cap_pressure_value[2];
            const double lambda = cap_pressure_value[3];  // always >= 1.0
            if (pc < pb)
                pc = pb;
            double const effect_sw = pow(pc / pb, -lambda);
            Sw = effect_sw * (slm - sgr - slr) + slr;
            break;
        }
        /// Liakopoulos
        case 10:
            Sw = 1 - (1.9722e-11) * pow(pc, 2.4279);
            if (pc < 0)
                // return 1 - (1.9722e-11)*pow(0.0, 2.4279);
                // extend
                Sw = 1 + pc * getDerivSaturation(0.0);
            break;
        default:
            break;
    }
    return Sw;
}

double TwoPhaseFlowWithPPMaterialProperties::getDerivSaturation(
    double const pc) const
{
    double dSwdPc;
    switch (cap_pressure_model)
    {
        case 0:
        {
            MathLib::PiecewiseLinearInterpolation const& interpolated_Pc =
                *curves.at("curve_PC_S");
            double dSwdPc = interpolated_Pc.getDerivative(pc);
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
            double const pb = cap_pressure_value[0];
            double const slr = cap_pressure_value[1];
            double const slm = 1.0;
            double const sgr = cap_pressure_value[2];
            double const m = cap_pressure_value[3];  // always <= 1.0.  Input is
                                                     // exponent = 1 /
                                                     // (1-lambda)
            assert(m <= 1);
            // pc = MRange(FLT_EPSILON, pc, capillary_pressure_values[4]);
            //
            //
            double const v1 = pow((pc / pb), (1.0 / (1.0 - m)));
            double const v2 = pow((1.0 + v1), (-1.0 - m));
            dSwdPc = (m * v1 * v2 * (slm - slr)) / ((m - 1.0) * pc);
            break;
        }
        case 2:
        {
            double const pb = cap_pressure_value[0];
            double const slr = cap_pressure_value[1];
            double const slm = 1.0;
            double const sgr = cap_pressure_value[2];
            const double lambda = cap_pressure_value[3];  // always >= 1.0
                                                          //
            double const v1 = pow((pc / pb), -lambda);
            dSwdPc = (lambda * v1 * (slr - slm)) / pc;
            break;
        }
        case 10:
            dSwdPc = -(1.9722e-11) * 2.4279 * pow(pc, 1.4279);
            if (pc <= 0)
                dSwdPc = -3.7901e-7;  // -(1.9722e-11)*2.4279*pow(0.0, 1.4279);
            break;
        default:
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
    switch (rel_wet_perm_model)
    {
        case 0:
        {
            MathLib::PiecewiseLinearInterpolation const& interpolated_Kr =
                *curves.at("curve_S_Krel_wet");
            rel_wet_perm = interpolated_Kr.getValue(sw);
            break;
        }
        case 1:
        case 2:
        {
            double const slr = rel_wet_perm_value[0];
            double const sgr = rel_wet_perm_value[1];
            double const slm = 1.0;
            double const m = rel_wet_perm_value[2];
            double const kr_min = rel_wet_perm_value[3];
            double sl = sw;
            double const se = (sl - slr) / (slm - slr - sgr);
            //
            rel_wet_perm = pow(se, 3.0 + 2.0 / m);
            if (rel_wet_perm < kr_min)
                rel_wet_perm = kr_min;
            break;
        }
        case 10:
            rel_wet_perm = 1 - 2.207 * pow((1 - sw), 1.0121);  //
            break;
        default:
            break;
    }
    return rel_wet_perm;
}

double TwoPhaseFlowWithPPMaterialProperties::getrelativePermeability_gas(
    double sw) const
{
    double rel_nonwet_perm;
    switch (rel_nonwet_perm_model)
    {
        case 0:
        {
            MathLib::PiecewiseLinearInterpolation const& interpolated_Kr =
                *curves.at("curve_S_Krel_nonwet");
            rel_nonwet_perm = interpolated_Kr.getValue(sw);
            break;
        }
        case 1:
        case 2:
        {
            double const slr = rel_nonwet_perm_value[0];
            double const sgr = rel_nonwet_perm_value[1];
            double const slm = 1.0;
            double const m = rel_nonwet_perm_value[2];
            double const kr_min = rel_nonwet_perm_value[3];
            double S_le = (sw - slr) / (slm - slr - sgr);
            rel_nonwet_perm =
                pow(1.0 - S_le, 2) * (1.0 - pow(S_le, 1.0 + 2.0 / m));
            if (rel_nonwet_perm < kr_min)
                rel_nonwet_perm = kr_min;
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
