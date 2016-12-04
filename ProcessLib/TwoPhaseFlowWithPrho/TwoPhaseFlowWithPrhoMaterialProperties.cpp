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
#include "NewtonRaphson.h"
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
		

		double TwoPhaseFlowWithPrhoMaterialProperties::getDerivCapillaryPressure(
			const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
			const double /*p*/, const double /*T*/, const double saturation) const
		{
			const double dpcdsw = _capillary_pressure_models[_current_material_id]->getdPcdS(saturation);
			return dpcdsw;
		}

		bool TwoPhaseFlowWithPrhoMaterialProperties::computeConstitutiveRelation(
			double const t,
			ProcessLib::SpatialPosition const& x,
			double const PG,
			double const X,
			double& Sw,
			double& X_m,
			double& dsw_dpg,
			double& dsw_dX,
			double& dxm_dpg,
			double& dxm_dX)
		{
			using LocalJacobianMatrix =
				Eigen::Matrix<double, 2, 2,
				Eigen::RowMajor>;
			LocalJacobianMatrix J_loc;
			// Linear solver for the newton loop is required after the loop with the
			// same matrix. This saves one decomposition.
			Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(2);

			// Different solvers are available for the solution of the local system.
			// TODO Make the following choice of linear solvers available from the
			// input file configuration:
			//      K_loc.partialPivLu().solve(-res_loc);
			//      K_loc.fullPivLu().solve(-res_loc);
			//      K_loc.householderQr().solve(-res_loc);
			//      K_loc.colPivHouseholderQr().solve(res_loc);
			//      K_loc.fullPivHouseholderQr().solve(-res_loc);
			//      K_loc.llt().solve(-res_loc);
			//      K_loc.ldlt().solve(-res_loc);

			{  // Local Newton solver
				using LocalResidualVector =
					Eigen::Matrix<double, 2, 1>;
				using LocalUnknownVector = Eigen::Matrix<double, 2, 1>;
				LocalJacobianMatrix J_loc;
				//LocalUnknownVector sec_var_unknown;
				//sec_var_unknown(0) = Sw;
				//sec_var_unknown(1) = X_m;
				auto const update_residual = [&](LocalResidualVector& residual) {
					calculateResidual(PG, X, Sw, X_m, residual);
				};

				auto const update_jacobian = [&](LocalJacobianMatrix& jacobian) {
					calculateJacobian(
						t, x, PG, X, jacobian, Sw, X_m);  // for solution dependent Jacobians
				};

				auto const update_solution = [&](LocalResidualVector const& increment) {
					// increment solution vectors
					Sw += increment[0];
					X_m += increment[1];
				};

				// TODO Make the following choice of maximum iterations and convergence
				// criteria available from the input file configuration:
				const int maximum_iterations(20);
				const double tolerance(1.e-14);

				auto newton_solver =
					NewtonRaphson<decltype(linear_solver), LocalJacobianMatrix,
					decltype(update_jacobian), LocalResidualVector,
					decltype(update_residual), decltype(update_solution)>(
						linear_solver, update_jacobian, update_residual,
						update_solution, maximum_iterations, tolerance);

				auto const success_iterations = newton_solver.solve(J_loc);

				if (!success_iterations)
					return false;

				// If the Newton loop didn't run, the linear solver will not be
				// initialized.
				// This happens usually for the first iteration of the first timestep.
				if (*success_iterations == 0)
					linear_solver.compute(J_loc);
			}
			dsw_dpg =Calculate_dSwdP(PG, Sw, X_m);
			dsw_dX = Calculate_dSwdX(PG, X, Sw, X_m);
			dxm_dpg = Calculate_dX_mdP(PG, Sw, X_m, dsw_dpg);
			dxm_dX = Calculate_dX_mdX(PG, Sw, X_m, dsw_dX);
			return true;
		}
		void TwoPhaseFlowWithPrhoMaterialProperties::calculateResidual(double const PL, double const X, double Sw, double rho_h2_wet,
			ResidualVector& res)
		{
			// getting unknowns
			const double PG = PL + _capillary_pressure_models[_current_material_id]->getCapillaryPressure(Sw);
			const double rho_h2_nonwet = PG*molar_mass_h2 / R / 303.15;

			// calculating residual
			res(0) = Calc_equili_rho_m(PG, Sw, rho_h2_wet);
			res(1) = Calc_Saturation(PL, X, Sw, rho_h2_wet, rho_h2_nonwet, 303.15);
		}

		void TwoPhaseFlowWithPrhoMaterialProperties::calculateJacobian(double const t,
			ProcessLib::SpatialPosition const& x,
			double const PL, double const X,
			JacobianMatrix& Jac, double Sw, double rho_h2_wet)
		{
			// getting unknowns
			const double PG = PL + _capillary_pressure_models[_current_material_id]->getCapillaryPressure(Sw);
			const double rho_h2_nonwet = PG*molar_mass_h2 / R / 303.15;
			double const rho_equili_h2_wet = PG * Hen * molar_mass_h2;
			double const dPC_dSw = _capillary_pressure_models[_current_material_id]->getdPcdS(Sw);
			double const drhoh2wet_dpg = Hen*molar_mass_h2;
			double const drhoh2nonwet_dpg = molar_mass_h2 / R / 303.15;
			const double RT = R * 303.15;
			// evaluate J
			Jac.setZero();
			if ((1 - Sw) < (rho_equili_h2_wet - rho_h2_wet))
			{
				Jac(0, 0) = -1;
				Jac(0, 1) = 0.0;
			}
			else
			{
				Jac(0, 0) = drhoh2wet_dpg*dPC_dSw;
				Jac(0, 1) = -1;
			}

			Jac(1, 0) = rho_h2_nonwet - rho_h2_wet;// -(1 - Sw)*drhoh2nonwet_dpg*dPC_dSw;
			//Jac(1, 1) = -Sw*N_L + Sw*(X - X_m)*rho_mol_h2o / std::pow(1 - X_m, 2);
			Jac(1, 1) = -Sw;
			//std::cout << Jac << std::endl;
		}
	}  // end of namespace
}  // end of namespace
