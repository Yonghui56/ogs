/**
* \copyright
* Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#include "TwoPhaseFlowWithPrhoMaterialProperties.h"
#include <logog/include/logog.hpp>
#include "BaseLib/reorderVector.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"

namespace ProcessLib
{
    namespace TwoPhaseFlowWithPrho
    {
        TwoPhaseFlowWithPrhoMaterialProperties::TwoPhaseFlowWithPrhoMaterialProperties(
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
            std::vector<std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>&&
            capillary_pressure_models,
            std::vector<std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>&&
            nonwet_relative_permeability_models,
            std::vector<std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>&&
            wet_relative_permeability_models)
            : _has_material_ids(has_material_ids),
            _liquid_density(std::move(liquid_density)),
            _viscosity(std::move(viscosity)),
            _gas_density(std::move(gas_density)),
            _gas_viscosity(std::move(gas_viscosity)),
            _material_ids(material_ids),
            _intrinsic_permeability_models(intrinsic_permeability_models),
            _porosity_models(std::move(porosity_models)),
            _storage_models(std::move(storage_models)),
            _capillary_pressure_models(std::move(capillary_pressure_models)),
            _nonwet_relative_permeability_models(std::move(nonwet_relative_permeability_models)),
            _wet_relative_permeability_models(std::move(wet_relative_permeability_models))

        {
            DBUG("Create material properties for Two-Phase flow with PP model.");
        }

        void TwoPhaseFlowWithPrhoMaterialProperties::setMaterialID(
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

        double TwoPhaseFlowWithPrhoMaterialProperties::getLiquidDensity(
            const double p, const double T) const
        {
            ArrayType vars;
            vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
            vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
            return _liquid_density->getValue(vars);
        }

        double TwoPhaseFlowWithPrhoMaterialProperties::getGasDensity(const double p,
            const double T) const
        {
            ArrayType vars;
            vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
            vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
            return _gas_density->getValue(vars);
        }

        double TwoPhaseFlowWithPrhoMaterialProperties::getDerivGasDensity(
            const double p, const double T) const
        {
            ArrayType vars;
            vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
            vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;

            return _gas_density->getdValue(
                vars, MaterialLib::Fluid::PropertyVariableType::p);
        }
        double TwoPhaseFlowWithPrhoMaterialProperties::getLiquidViscosity(
            const double p, const double T) const
        {
            ArrayType vars;
            vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
            vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
            return _viscosity->getValue(vars);
        }

        double TwoPhaseFlowWithPrhoMaterialProperties::getGasViscosity(
            const double p, const double T) const
        {
            ArrayType vars;
            vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
            vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
            return _gas_viscosity->getValue(vars);
        }

        Eigen::MatrixXd const& TwoPhaseFlowWithPrhoMaterialProperties::getPermeability(
            const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
            const int /*dim*/) const
        {
            return _intrinsic_permeability_models[_current_material_id];
        }

        double TwoPhaseFlowWithPrhoMaterialProperties::getPorosity(
            const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
            const double /*p*/, const double T, const double porosity_variable) const
        {
            const double porosity =
                _porosity_models[_current_material_id]->getValue(porosity_variable, T);

            return porosity;
        }

        double TwoPhaseFlowWithPrhoMaterialProperties::getNonwetRelativePermeability(
            const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
            const double /*p*/, const double /*T*/, const double saturation) const
        {
            const double nonwet_krel = _nonwet_relative_permeability_models[0]->getValue(saturation);
                 

            return nonwet_krel;
        }

        double TwoPhaseFlowWithPrhoMaterialProperties::getWetRelativePermeability(
            const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
            const double /*p*/, const double /*T*/, const double saturation) const
        {
            const double wet_krel = _nonwet_relative_permeability_models[1]->getValue(saturation);


            return wet_krel;
        }
        
        double TwoPhaseFlowWithPrhoMaterialProperties::getCapillaryPressure(
            const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
            const double /*p*/, const double /*T*/, const double saturation) const
        {
            const double pc = _capillary_pressure_models[_current_material_id]->getCapillaryPressure(saturation);
            return pc;
        }

		double TwoPhaseFlowWithPrhoMaterialProperties::getRegularizedCapillaryPressure(
			const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
			const double /*p*/, const double /*T*/, const double saturation) const
		{
			const double pc = _capillary_pressure_models[_current_material_id]->getRegularizedCapillaryPressure(saturation);
			return pc;
		}

        double TwoPhaseFlowWithPrhoMaterialProperties::getDerivCapillaryPressure(
            const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
            const double /*p*/, const double /*T*/, const double saturation) const
        {
            const double dpcdsw = _capillary_pressure_models[_current_material_id]->getdPcdS(saturation);
            return dpcdsw;
        }
		double TwoPhaseFlowWithPrhoMaterialProperties::getRegularizedDerivCapillaryPressure(
			const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
			const double /*p*/, const double /*T*/, const double saturation) const
		{
			const double dpcdsw = _capillary_pressure_models[_current_material_id]->getRegularizedPcdS(saturation);
			return dpcdsw;
		}
    }  // end of namespace
}  // end of namespace
