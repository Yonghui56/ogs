/**
* \copyright
* Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#ifndef OGS_TWOPHASEFLOWWITHPRHOMATERIALPROPERTIES_H
#define OGS_TWOPHASEFLOWWITHPRHOMATERIALPROPERTIES_H

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
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MeshLib
{
	template <typename PROP_VAL_TYPE>
	class PropertyVector;
}

namespace ProcessLib
{
	class SpatialPosition;
	namespace TwoPhaseFlowWithPrho
	{
		class TwoPhaseFlowWithPrhoMaterialProperties
		{
		public:
			static int const JacobianResidualSize = 2;
			using ResidualVector = Eigen::Matrix<double, JacobianResidualSize, 1>;
			using JacobianMatrix = Eigen::Matrix<double,
				JacobianResidualSize,
				JacobianResidualSize,
				Eigen::RowMajor>;
			using UnknownVector = Eigen::Matrix<double, JacobianResidualSize, 1>;
		public:
			using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

			TwoPhaseFlowWithPrhoMaterialProperties(
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
				wet_relative_permeability_models);

			void setMaterialID(const ProcessLib::SpatialPosition& pos);

			Eigen::MatrixXd const& getPermeability(
				const double t,
				const ProcessLib::SpatialPosition& pos,
				const int dim) const;

			double getPorosity(const double t, const ProcessLib::SpatialPosition& pos,
				const double p, const double T,
				const double porosity_variable) const;

			double getNonwetRelativePermeability(const double t, const ProcessLib::SpatialPosition& pos,
				const double p, const double T,
				const double saturation) const;
			double getWetRelativePermeability(const double t, const ProcessLib::SpatialPosition& pos,
				const double p, const double T,
				const double saturation) const;
			double getCapillaryPressure(const double t, const ProcessLib::SpatialPosition& pos,
				const double p, const double T,
				const double saturation) const;
			double getDerivCapillaryPressure(const double t, const ProcessLib::SpatialPosition& pos,
				const double p, const double T,
				const double saturation) const;
			double getLiquidDensity(const double p, const double T) const;
			double getGasDensity(const double p, const double T) const;
			double getGasViscosity(const double p, const double T) const;
			double getLiquidViscosity(const double p, const double T) const;
			double getDerivGasDensity(double const p, double const T) const;
			bool computeConstitutiveRelation(
				double const t,
				ProcessLib::SpatialPosition const& x_position,
				double const PG,
				double const X,
				double& Sw,
				double& X_m,
				double& dsw_dpg,
				double& dsw_dX,
				double& dxm_dpg,
				double& dxm_dX
			);
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
			std::vector<std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
				_capillary_pressure_models;
			std::vector<std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>
				_nonwet_relative_permeability_models;
			std::vector<std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>
				_wet_relative_permeability_models;
		private:
			/// Calculates the 18x1 residual vector.
			void calculateResidual(double const PG, double const X, double Sw, double X_m,
				ResidualVector& res);

			/// Calculates the 18x18 Jacobian.
			void calculateJacobian(double const t,
				ProcessLib::SpatialPosition const& x,
				double const PG, double const X,
				JacobianMatrix& Jac, double Sw, double X_m);

			/**
			* Complementary condition 1
			* for calculating molar fraction of light component in the liquid phase
			*/
			const double Calc_equili_rho_m(double const PG, double const Sw, double const rho_h2_wet)
			{
			    double const rho_equili_h2_wet = PG * Hen * molar_mass_h2;
				return std::min(1 - Sw, rho_equili_h2_wet - rho_h2_wet);
			}
			/**
			* Complementary condition 2
			* for calculating the saturation
			*/
			const double Calc_Saturation(double PL, double X, double Sw, double rho_h2_wet,
				double rho_h2_nonwet, double T)
			{
				return X - (Sw*rho_h2_wet + (1 - Sw)*rho_h2_nonwet);
			}
			/**
			* Calculate the derivatives using the analytical way
			*/
			const double Calculate_dSwdP(double PL, double S, double rho_h2_wet,
				double T = 303.15)
			{
				const double PG = PL + _capillary_pressure_models[_current_material_id]->getCapillaryPressure(S);
				double const rho_equili_h2_wet = PG * Hen * molar_mass_h2;
				const double rho_h2_nonwet = PG*molar_mass_h2 / R / 303.15;
				if ((1 - S) < (rho_equili_h2_wet - rho_h2_wet))
				{
					return 0.0;
				}
				else {
					double const drhoh2wet_dpg = Hen*molar_mass_h2;
					double const drhoh2nonwet_dpg = molar_mass_h2 / R / T;
					double const alpha = ((drhoh2nonwet_dpg - drhoh2wet_dpg)*(1 - S) + drhoh2wet_dpg);
					double const beta= (drhoh2nonwet_dpg - drhoh2wet_dpg)*PG;//NOTE here should be PG^h, but we ignore vapor
					double const dPC_dSw = _capillary_pressure_models[_current_material_id]->getdPcdS(S);
					return -alpha / (beta - alpha*dPC_dSw);
				}
			}
			/*
			* Calculate the derivative using the analytical way
			*/
			const double Calculate_dSwdX(double PL, double X, double S, double rho_h2_wet,
				double T = 303.15)
			{
				const double PG = PL + _capillary_pressure_models[_current_material_id]->getCapillaryPressure(S);
				double const rho_equili_h2_wet = PG * Hen * molar_mass_h2;
				const double rho_h2_nonwet = PG*molar_mass_h2 / R / 303.15;
				if ((1 - S) < (rho_equili_h2_wet - rho_h2_wet))
				{
					return 0.0;
				}
				else {
					double const drhoh2wet_dpg = Hen*molar_mass_h2;
					double const drhoh2nonwet_dpg = molar_mass_h2 / R / T;
					double const alpha = ((drhoh2nonwet_dpg - drhoh2wet_dpg)*(1 - S) + drhoh2wet_dpg);
					double const beta = (drhoh2nonwet_dpg - drhoh2wet_dpg)*PG;//NOTE here should be PG^h, but we ignore vapor
					double const dPC_dSw = _capillary_pressure_models[_current_material_id]->getdPcdS(S);
					return -1 / (beta - alpha*dPC_dSw);
				}
			}
			/*
			* Calculate the derivative using the analytical way
			*/
			const double Calculate_dX_mdX(double PL, double Sw, double rho_h2_wet, double dSwdX)
			{
				const double PG = PL + _capillary_pressure_models[_current_material_id]->getCapillaryPressure(Sw);
				double const rho_equili_h2_wet = PG * Hen * molar_mass_h2;
				double const dPC_dSw = _capillary_pressure_models[_current_material_id]->getdPcdS(Sw);
				if ((1 - Sw) < (rho_equili_h2_wet - rho_h2_wet))
					return 1.0;
				return Hen * molar_mass_h2*dPC_dSw*dSwdX;
			}
			/*
			* Calculate the derivative using the analytical way
			*/
			const double Calculate_dX_mdP(double PL, double Sw, double rho_h2_wet, double dSwdP)
			{
				const double PG = PL + _capillary_pressure_models[_current_material_id]->getCapillaryPressure(Sw);
				double const rho_equili_h2_wet = PG * Hen * molar_mass_h2;
				double const dPC_dSw = _capillary_pressure_models[_current_material_id]->getdPcdS(Sw);
				if ((1 - Sw) < (rho_equili_h2_wet - rho_h2_wet))
					return 0.0;
				return Hen * molar_mass_h2*(1+ dPC_dSw*dSwdP);
			}

		private:
			double const Hen = 7.65e-6;  // mol/Pa./m3
										 /**
										 * Molar mass of water
										 */
			double const molar_mass_h2o =
				MaterialLib::PhysicalConstant::MolarMass::Water;
			double const molar_mass_h2 = MaterialLib::PhysicalConstant::MolarMass::H2;
			/**
			* mass density of water
			*/
			double const rho_mass_h20 = 1000;
			/**
			* molar density of water
			*/
			double const rho_mol_h2o = rho_mass_h20 / molar_mass_h2o;
			double const R = MaterialLib::PhysicalConstant::IdealGasConstant;
		};

	}  // end of namespace
}  // end of namespace
#endif /* TWOPHASEFLOWWITHPPMATERIALPROPERTIES_H */
