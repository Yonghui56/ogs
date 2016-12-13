/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_TWOPHASEFLOWWITHPRHOLOCALASSEMBLER_IMPL_H
#define OGS_TWOPHASEFLOWWITHPRHOLOCALASSEMBLER_IMPL_H

#include "TwoPhaseFlowWithPrhoLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseFlowWithPrhoProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseFlowWithPrhoLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, std::vector<double> const& local_x,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();

    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    NodalMatrixType mass_operator;
    mass_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Mgp =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mgx = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);

    auto Mlp = local_M.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Mlx = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    NodalMatrixType laplace_operator;
    laplace_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kgp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Kgx = local_K.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);

    auto Klp = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Klx = local_K.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    auto Bg = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);

    auto Bl =
        local_b.template segment<cap_pressure_size>(cap_pressure_matrix_index);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    _process_data._material->setMaterialID(pos);

    const Eigen::MatrixXd& perm = _process_data._material->getPermeability(
        t, pos, _element.getDimension());
    assert(perm.rows() == GlobalDim || perm.rows() == 1);
    GlobalDimMatrixType permeability =
        GlobalDimMatrixType::Zero(GlobalDim, GlobalDim);
    if (perm.rows() == GlobalDim)
        permeability = perm;
    else if (perm.rows() == 1)
        permeability.diagonal().setConstant(perm(0, 0));
    // Note: currently only isothermal case is considered, so the temperature is
    // assumed to be const
    // the variation of temperature will be taken into account in future
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];

        double pl_int_pt = 0.;
        double totalrho_int_pt = 0.;  // total mass density of the light component
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pl_int_pt,
                                         totalrho_int_pt);

        auto const& wp = _integration_method.getWeightedPoint(ip);
        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();

        double const rho_gas =
            _process_data._material->getGasDensity(pl_int_pt, _temperature);
        double const rho_h2o = _process_data._material->getLiquidDensity(
            _pressure_wetting[ip], _temperature);

        double& Sw = _ip_data[ip]._sw;
        double const X_h2_nonwet = 1.0;        // TODO
        double& rho_h2_wet = _ip_data[ip]._rho_m;  // TODO
        double& dSwdP_gp = _ip_data[ip]._dsw_dpg;
		double& dSwdrho_gp = _ip_data[ip]._dsw_drho;
        double& drhoh2wet_dp = _ip_data[ip]._drhom_dpg;
        double& drhoh2wet_drho = _ip_data[ip]._drhom_drho;
		if (!_ip_data[ip]._mat_property.computeConstitutiveRelation(
			t,
			pos,
			pl_int_pt,
			totalrho_int_pt,
			Sw,
			rho_h2_wet,
			dSwdP_gp,
			dSwdrho_gp,
			drhoh2wet_dp,
			drhoh2wet_drho))
			OGS_FATAL("Computation of local constitutive relation failed.");
        double const pc =
            _process_data._material->getCapillaryPressure(t, pos, pl_int_pt,
           _temperature, Sw);
		double const rho_wet = rho_h2o + rho_h2_wet;
        
        _saturation[ip] = Sw;
        _pressure_wetting[ip] = pl_int_pt - pc;

        // Assemble M matrix
        // nonwetting
        double const drhogas_dpg = _process_data._material->getDerivGasDensity(
            pl_int_pt, _temperature);

        double dPC_dSw_gp =
            _process_data._material->getDerivCapillaryPressure(
                t, pos, pl_int_pt, _temperature, Sw);

        /*double dPC_dSw_gp =
            _process_data._material->getRegularizedDerivCapillaryPressure(
                t, pos, pl_int_pt, _temperature, Sw);*/

        double const porosity = _process_data._material->getPorosity(
            t, pos, pl_int_pt, _temperature, 0);
        mass_operator.noalias() = sm.N.transpose() * sm.N * integration_factor;
       
        Mgx.noalias() +=
            porosity *
            mass_operator;

		Mlp.noalias() += porosity*rho_h2o*dSwdP_gp*mass_operator;

		Mlx.noalias() += porosity*(1 + dSwdrho_gp*rho_h2o) * mass_operator;
        double const k_rel_G =
            _process_data._material->getNonwetRelativePermeability(
                t, pos, pl_int_pt, _temperature, Sw);
        double const mu_gas =
            _process_data._material->getGasViscosity(pl_int_pt, _temperature);
        double const lambda_G = k_rel_G / mu_gas;
        double const diffusion_coeff_componenth2 =
            _process_data._diffusion_coeff_componentb(t, pos)[0];
        double const diffusion_coeff_componentw =
            _process_data._diffusion_coeff_componenta(t, pos)[0];
        // wet
        double const k_rel_L =
            _process_data._material->getWetRelativePermeability(
                t, pos, pl_int_pt, _temperature,
                Sw); 
        double const mu_liquid = _process_data._material->getLiquidViscosity(
            _pressure_wetting[ip], _temperature);
        double const lambda_L = k_rel_L / mu_liquid;

        laplace_operator.noalias() =
            sm.dNdx.transpose() * permeability * sm.dNdx * integration_factor;

        Kgp.noalias() +=
            (rho_gas * X_h2_nonwet * lambda_G *(1+ dPC_dSw_gp * dSwdP_gp)+
             rho_h2_wet * lambda_L ) *
                laplace_operator +
            (Sw * porosity * diffusion_coeff_componenth2 * (rho_h2o/rho_wet) *
             drhoh2wet_dp )*
                sm.dNdx.transpose() * sm.dNdx * integration_factor;
        Kgx.noalias() +=
            (rho_gas * X_h2_nonwet * lambda_G * dPC_dSw_gp * dSwdrho_gp) *
                laplace_operator +
			(Sw * porosity * diffusion_coeff_componenth2 * (rho_h2o / rho_wet) * 
				drhoh2wet_drho )*
				sm.dNdx.transpose() * sm.dNdx * integration_factor;
		Klp.noalias() +=
			(rho_gas * lambda_G * (1+ dPC_dSw_gp * dSwdP_gp) +
				rho_wet * lambda_L) *
			laplace_operator;

		Klx.noalias() +=
			(rho_gas * lambda_G  * dPC_dSw_gp * dSwdrho_gp) *
			laplace_operator;

        // rho_mol_nonwet * ((1 - x_h2o_in_nonwet) * molar_mass_co2 +
        // x_h2o_in_nonwet * molar_mass_h2o);
        if (_process_data._has_gravity)
        {
            auto const& b = _process_data._specific_body_force;
            Bg.noalias() +=
                (rho_gas * rho_gas * lambda_G +
                 rho_h2_wet * rho_wet * lambda_L) *
                sm.dNdx.transpose() * permeability * b * integration_factor;
            Bl.noalias() +=
                (rho_wet * lambda_L * rho_wet +
					rho_gas * rho_gas  * lambda_G) *
                sm.dNdx.transpose() * permeability * b * integration_factor;

        }  // end of has gravity
    }      // end of GP
    if (_process_data._has_mass_lumping)
    {
        for (unsigned row = 0; row < Mgp.cols(); row++)
        {
            for (unsigned column = 0; column < Mgp.cols(); column++)
            {
                if (row != column)
                {
                    Mgx(row, row) += Mgx(row, column);
                    Mgx(row, column) = 0.0;
                    Mgp(row, row) += Mgp(row, column);
                    Mgp(row, column) = 0.0;
                    Mlx(row, row) += Mlx(row, column);
                    Mlx(row, column) = 0.0;
                    Mlp(row, row) += Mlp(row, column);
                    Mlp(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // end of namespace
}  // end of namespace

#endif
