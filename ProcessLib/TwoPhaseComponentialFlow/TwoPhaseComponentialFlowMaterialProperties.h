/**
* \copyright
* Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MeshLib
{
	template <typename PROP_VAL_TYPE>
	class PropertyVector;
}

namespace ProcessLib
{
	class SpatialPosition;
	namespace TwoPhaseComponentialFlow
	{
		class TwoPhaseComponentialFlowMaterialProperties
		{
		public:
			static int const JacobianResidualSize = 2;
			using ResidualVector = Eigen::Matrix<double, JacobianResidualSize, 1>;
			using JacobianMatrix = Eigen::Matrix<double, JacobianResidualSize,
				JacobianResidualSize, Eigen::RowMajor>;
			using UnknownVector = Eigen::Matrix<double, JacobianResidualSize, 1>;

		public:
			using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

			TwoPhaseComponentialFlowMaterialProperties(
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
				std::vector<std::unique_ptr<
				MaterialLib::PorousMedium::CapillaryPressureSaturation>>&&
				capillary_pressure_models,
				std::vector<
				std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>&&
				relative_permeability_models);

			void setMaterialID(const ProcessLib::SpatialPosition& pos);
			double getMaterialID(const ProcessLib::SpatialPosition& pos) const;
			Eigen::MatrixXd const& getPermeability(
				const double t,
				const ProcessLib::SpatialPosition& pos,
				const int dim) const;

			double getPorosity(const double t, const ProcessLib::SpatialPosition& pos,
				const double p, const double T,
				const double porosity_variable) const;

			double getNonwetRelativePermeability(const double t,
				const ProcessLib::SpatialPosition& pos,
				const double p, const double T,
				const double saturation) const;
			double getWetRelativePermeability(const double t,
				const ProcessLib::SpatialPosition& pos,
				const double p, const double T,
				const double saturation) const;
			double getSaturation(const double t,
				const ProcessLib::SpatialPosition& pos,
				const double p, const double T,
				const double pc) const;
			double getDerivSaturation(const double t,
				const ProcessLib::SpatialPosition& pos,
				const double p, const double T,
				const double pc) const;
			double getLiquidDensity(const double p, const double T) const;
			double getGasDensity(const double p, const double T) const;
			double getGasViscosity(const double p, const double T) const;
			double getLiquidViscosity(const double p, const double T) const;
			double getDerivGasDensity(double const p, double const T) const;

		protected:
			/// A flag to indicate whether the reference member, _material_ids,
			/// is not assigned.
			const bool _has_material_ids;

			std::unique_ptr<MaterialLib::Fluid::FluidProperty> _liquid_density;
			std::unique_ptr<MaterialLib::Fluid::FluidProperty> _viscosity;
			std::unique_ptr<MaterialLib::Fluid::FluidProperty> _gas_density;
			std::unique_ptr<MaterialLib::Fluid::FluidProperty> _gas_viscosity;

			/** Use two phase models for different material zones.
			*  Material IDs must be given as mesh element properties.
			*/
			MeshLib::PropertyVector<int> const& _material_ids;

			int _current_material_id = 0;
			std::vector<Eigen::MatrixXd> _intrinsic_permeability_models;
			std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
				_porosity_models;
			std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
				_storage_models;
			std::vector<
				std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
				_capillary_pressure_models;
			std::vector<
				std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>
				_relative_permeability_models;

		private:
			double const Hen = 7.65e-6;  // mol/Pa./m3
			double const molar_mass_h2o =
				MaterialLib::PhysicalConstant::MolarMass::Water;
			double const molar_mass_h2 = MaterialLib::PhysicalConstant::MolarMass::H2;
			/**
			* mass density of water
			*/
			double const rho_mass_h20 = 1000;

			double const R = MaterialLib::PhysicalConstant::IdealGasConstant;
		};

	}  // end of namespace
}  // end of namespace
