/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TwoPhaseFlowWithPPLocalAssembler.cpp
 *
 * Created on October 19, 2016, 2:28 PM
 */

#ifndef OGS_TWOPHASEFLOWWITHPPLOCALASSEMBLER_IMPL_H
#define OGS_TWOPHASEFLOWWITHPPLOCALASSEMBLER_IMPL_H

#include <iostream>
#include "TwoPhaseFlowWithPPLocalAssembler.h"

#include "NumLib/Function/Interpolation.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "TwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseFlowWithPPLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assemble(double const t, std::vector<double> const& local_x,
             std::vector<double>& local_M_data,
             std::vector<double>& local_K_data,
             std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();
	
	assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
	int const n_nodes = ShapeFunction::NPOINTS;
	Eigen::MatrixXd mass_mat_coeff =
		Eigen::MatrixXd::Zero(NUM_NODAL_DOF, NUM_NODAL_DOF);
	Eigen::MatrixXd K_mat_coeff = Eigen::MatrixXd::Zero(NUM_NODAL_DOF, NUM_NODAL_DOF);
	Eigen::VectorXd H_vec_coeff = Eigen::VectorXd::Zero(NUM_NODAL_DOF);
	Eigen::VectorXd F_vec_coeff = Eigen::VectorXd::Zero(NUM_NODAL_DOF);
	Eigen::MatrixXd localMass_tmp = Eigen::MatrixXd::Zero(
		ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
	Eigen::MatrixXd localDispersion_tmp = Eigen::MatrixXd::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
	Eigen::VectorXd localGravity_tmp = Eigen::VectorXd::Zero(ShapeFunction::NPOINTS);

    auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
	auto local_b = MathLib::createZeroedVector<NodalVectorType>(
		local_b_data, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    _material_properties.setMaterialID(pos);

    const Eigen::MatrixXd& perm =
        _material_properties.getPermeability(t, pos, _element.getDimension());

    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to perm.rows() == _element->getDimension()
    assert(perm.rows() == GlobalDim || perm.rows() == 1);

	/*Eigen::MatrixXd mass_mat_coeff = Eigen::MatrixXd::Zero(N_size, N_size);
	Eigen::MatrixXd K_mat_coeff = Eigen::MatrixXd::Zero(N_size, N_size);
	Eigen::VectorXd H_vec_coeff = Eigen::VectorXd::Zero(N_size);
	Eigen::MatrixXd localMass_tmp = Eigen::MatrixXd::Zero(n_nodes, n_nodes);
	Eigen::MatrixXd localDispersion_tmp = Eigen::MatrixXd::Zero(n_nodes, n_nodes);
	Eigen::VectorXd localGravity_tmp = Eigen::VectorXd::Zero(n_nodes);*/
	//MathLib::PiecewiseLinearInterpolation const&  interpolated_Pc = *_process_data.curves.at("curveA");
	//MathLib::PiecewiseLinearInterpolation const&  interpolated_Kr_L = *_process_data.curves.at("curveB");
    // TODO: The following two variables should be calculated inside the
    //       the integration loop for non-constant porosity and storage models.
    double porosity_variable = 0.;
    double storage_variable = 0.;
	_temperature = 293.15;
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);

        double pc_int_pt = 0.;
		double pg_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pc_int_pt, pg_int_pt);

        // TODO : compute _temperature from the heat transport pcs
		double const pl = pg_int_pt - pc_int_pt;
        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();
		double const rho_gas = _material_properties.getGasDensity(pg_int_pt, _temperature);
		double const rho_w =
			_material_properties.getLiquidDensity(pl, _temperature);
		double const mu_gas = _material_properties.getGasViscosity(pg_int_pt, _temperature);
		double const mu_liquid = _material_properties.getLiquidViscosity(pl, _temperature);

		double Sw = _material_properties.getSaturation(pc_int_pt);//pc
		double const dSwdPc = _material_properties.getDerivSaturation(pc_int_pt);
		double const k_rel_L = _material_properties.getrelativePermeability_liquid(Sw);
		
		double const k_rel_G = _material_properties.getrelativePermeability_gas(Sw);//Sw
		double const drhogas_dpg = _material_properties.getDerivGasDensity(pg_int_pt, _temperature);
		double const poro = _material_properties.getPorosity(t, pos, pg_int_pt, _temperature, porosity_variable);
		double const rho_L_air = _material_properties.getDissolvedGas(pg_int_pt);
		double const drhoLairdPG = 0.029*2e-6;
        // Assemble mass matrix, M
		_saturation[ip] = Sw;
		//air
		mass_mat_coeff(0, 0) = -poro*rho_gas*dSwdPc+ poro*rho_L_air*dSwdPc;
		mass_mat_coeff(0, 1) = poro* (1 - Sw)*drhogas_dpg+poro*Sw*drhoLairdPG;
		mass_mat_coeff(1, 0) = poro*rho_w*dSwdPc;
		mass_mat_coeff(1, 1) = 0.0;

		//std::cout << mass_mat_coeff << std::endl;
		//assembly the mass matrix
		for (int ii = 0; ii < NUM_NODAL_DOF; ii++) {
			for (int jj = 0; jj < NUM_NODAL_DOF; jj++) {
				localMass_tmp.setZero();
				localMass_tmp.noalias() = sm.N.transpose() *
					mass_mat_coeff(ii, jj) *  sm.N *
					sm.detJ * sm.integralMeasure * wp.getWeight();
				local_M.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes).noalias() += localMass_tmp;
			}
		}
		//std::cout << local_M << std::endl;
		/*
		*construct the K matrix
		*/
		K_mat_coeff(0, 0) = -rho_L_air*perm(0, 0)*k_rel_L / mu_liquid;// _process_data.intrinsic_permeability(_element)*k_rel / _process_data.viscosity(_element);
		K_mat_coeff(0, 1) = rho_gas*perm(0, 0)*k_rel_G / mu_gas + rho_L_air*perm(0, 0)*k_rel_L / mu_liquid;
		K_mat_coeff(1, 0) = -rho_w*perm(0, 0)*k_rel_L / mu_liquid;
		K_mat_coeff(1, 1) = rho_w*perm(0, 0)*k_rel_L / mu_liquid;
		//std::cout << K_mat_coeff << std::endl;
		//assembly the mass matrix
		for (int ii = 0; ii < NUM_NODAL_DOF; ii++) {
			for (int jj = 0; jj < NUM_NODAL_DOF; jj++) {
				localDispersion_tmp.setZero();
				localDispersion_tmp.noalias() = sm.dNdx.transpose() *
					K_mat_coeff(ii, jj) * sm.dNdx *
					sm.detJ * sm.integralMeasure * wp.getWeight();
				local_K.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes).noalias() += localDispersion_tmp;

			}
		}
		//std::cout << local_K << std::endl;
		H_vec_coeff(0) = rho_gas*rho_gas*perm(0, 0)*k_rel_G / mu_gas
			+ rho_L_air *rho_w*perm(0, 0)*k_rel_L / mu_liquid;
		H_vec_coeff(1) = rho_w*rho_w*perm(0, 0)*k_rel_L / mu_liquid;
		//std::cout << H_vec_coeff << std::endl;
		if (_process_data.has_gravity)
		{
			auto const body_force =
				_process_data.specific_body_force(t, pos);
			assert(body_force.size() == GlobalDim);
			auto const b = MathLib::toVector<GlobalDimVectorType>(
				body_force, GlobalDim);
			localGravity_tmp.noalias() = sm.dNdx.transpose() * H_vec_coeff(0)*
					b * sm.detJ * sm.integralMeasure *
					wp.getWeight();
			
			local_b.block<n_nodes, 1>(0, 0).noalias() += localGravity_tmp;
			localGravity_tmp.setZero();
			localGravity_tmp.noalias() = sm.dNdx.transpose() * H_vec_coeff(1)*
				b * sm.detJ * sm.integralMeasure *
				wp.getWeight();
			
			local_b.block<n_nodes, 1>(n_nodes, 0).noalias() += localGravity_tmp;
			
		}// end of has gravity
		//std::cout << local_b << std::endl;
    }// end of GP
	if (_process_data.has_mass_lumping)
	{
		for (int idx_ml = 0; idx_ml < local_M.cols(); idx_ml++)
		{
			double const mass_lump_val = local_M.col(idx_ml).sum();
			local_M.col(idx_ml).setZero();
			local_M(idx_ml, idx_ml) = mass_lump_val;
		}
	}  // end of mass lumping
}

}  // end of namespace
}  // end of namespace

#endif
