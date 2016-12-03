/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_TWOPHASEFLOWWITHPXLOCALASSEMBLER_IMPL_H
#define OGS_TWOPHASEFLOWWITHPXLOCALASSEMBLER_IMPL_H

#include "TwoPhaseFlowWithPXLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseFlowWithPXProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPX
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseFlowWithPXLocalAssembler<
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

    // NodalMatrixType laplace_diffusion_operator;
    // laplace_diffusion_operator.setZero(ShapeFunction::NPOINTS,
    // ShapeFunction::NPOINTS);

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

        double pg_int_pt = 0.;
        double totalX_int_pt = 0.;  // total molar
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pg_int_pt,
                                         totalX_int_pt);

        auto const& wp = _integration_method.getWeightedPoint(ip);
        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();

        double const rho_gas =
            _process_data._material->getGasDensity(pg_int_pt, _temperature);
        double const rho_w = _process_data._material->getLiquidDensity(
            _pressure_wetting[ip], _temperature);
        double const rho_mol_nonwet = rho_gas / molar_mass_h2;
        double const rho_mol_w_std = rho_w / molar_mass_h2o;

        double& Sw = _ip_data[ip]._sw;
        double const X_h2_nonwet = 1.0;        // TODO
        double& X_h2_wet = _ip_data[ip]._x_m;  // TODO
        double& dSwdP_gp = _ip_data[ip]._dsw_dpg;
        double& dSwdX_gp = _ip_data[ip]._dsw_dX;
        double& dXh2wet_dpg = _ip_data[ip]._dxm_dpg;
        double& dXh2wet_dX = _ip_data[ip]._dxm_dX;
        if (!_ip_data[ip]._EoS_material.computeConstitutiveRelation(
                t,
                pos,
                pg_int_pt,
                totalX_int_pt,
                Sw,
                X_h2_wet,
                dSwdP_gp,
                dSwdX_gp,
                dXh2wet_dpg,
                dXh2wet_dX))
            OGS_FATAL("Computation of local constitutive relation failed.");
        /*double const pc =
            _process_data._material->getCapillaryPressure(t, pos, pg_int_pt,
           _temperature, Sw);*/
        double const pc =
            _process_data._material->getRegularizedCapillaryPressure(
                t, pos, pg_int_pt, _temperature, Sw);
        _saturation[ip] = Sw;  // there is no need
        _pressure_wetting[ip] = pg_int_pt - pc;
        double const rho_mol_wet = rho_mol_w_std / (1 - X_h2_wet);
        // Assemble M matrix
        // nonwetting
        double const drhogas_dpg = _process_data._material->getDerivGasDensity(
            pg_int_pt, _temperature);
        double const drhomolwet_dpg =
            rho_mol_w_std * dXh2wet_dpg / pow(-1 + X_h2_wet, 2);  // dN_LdPG
        double const drhomolwet_dX =
            rho_mol_w_std * dXh2wet_dX / pow(-1 + X_h2_wet, 2);        // TODO
        double const drhomolnonwet_dpg = drhogas_dpg / molar_mass_h2;  //
        double const drhomolnonwet_dX = 0.0;                           // TODO
        double const dXh2nonwet_dpg = 0.0;                             // TODO
        double const dXh2nonwet_dX = 0.0;                              // TODO
        double dPC_dSw_gp_test =
            _process_data._material->getDerivCapillaryPressure(
                t, pos, pg_int_pt, _temperature, Sw);

        double dPC_dSw_gp =
            _process_data._material->getRegularizedDerivCapillaryPressure(
                t, pos, pg_int_pt, _temperature, Sw);

        double const porosity = _process_data._material->getPorosity(
            t, pos, pg_int_pt, _temperature, 0);
        mass_operator.noalias() = sm.N.transpose() * sm.N * integration_factor;
        // X*(N_G*(1-Sw)+N_L*Sw)
        /*double const mgp_coeff = porosity *
            (totalX_int_pt *
            (-dSwdP_gp * rho_mol_nonwet + (1 - Sw) * drhomolnonwet_dpg +
                Sw * drhomolwet_dpg + dSwdP_gp * rho_mol_wet));
        double const mgx_coeff = porosity *
            ((1 - Sw) * rho_mol_nonwet + Sw * rho_mol_wet +
                totalX_int_pt *
                ((1 - Sw) * drhomolnonwet_dX - rho_mol_nonwet * dSwdX_gp +
                    Sw * drhomolwet_dX + dSwdX_gp * rho_mol_wet));
        double const mlp_coeff = porosity *
            ((1 - totalX_int_pt) *
            (-dSwdP_gp * rho_mol_nonwet + (1 - Sw) * drhomolnonwet_dpg +
                Sw * drhomolwet_dpg + dSwdP_gp * rho_mol_wet));
        double const mlx_coeff = porosity* (-((1 - Sw) * rho_mol_nonwet + Sw *
        rho_mol_wet) +
            (1 - totalX_int_pt) *
            ((1 - Sw) * drhomolnonwet_dX - rho_mol_nonwet * dSwdX_gp +
                Sw * drhomolwet_dX + dSwdX_gp * rho_mol_wet));
        std::cout << mgp_coeff << std::endl;
        std::cout << mgx_coeff << std::endl;
        std::cout << mlp_coeff << std::endl;
        std::cout << mlx_coeff << std::endl;*/
        Mgp.noalias() +=
            porosity *
            (totalX_int_pt *
             (-dSwdP_gp * rho_mol_nonwet + (1 - Sw) * drhomolnonwet_dpg +
              Sw * drhomolwet_dpg + dSwdP_gp * rho_mol_wet)) *
            mass_operator;
        Mgx.noalias() +=
            porosity *
            ((1 - Sw) * rho_mol_nonwet + Sw * rho_mol_wet +
             totalX_int_pt *
                 ((1 - Sw) * drhomolnonwet_dX - rho_mol_nonwet * dSwdX_gp +
                  Sw * drhomolwet_dX + dSwdX_gp * rho_mol_wet)) *
            mass_operator;
        //(1-X)*(N_G*(1-Sw)+N_L*Sw)
        Mlp.noalias() +=
            porosity *
            ((1 - totalX_int_pt) *
             (-dSwdP_gp * rho_mol_nonwet + (1 - Sw) * drhomolnonwet_dpg +
              Sw * drhomolwet_dpg + dSwdP_gp * rho_mol_wet)) *
            mass_operator;
        Mlx.noalias() +=
            porosity *
            (-((1 - Sw) * rho_mol_nonwet + Sw * rho_mol_wet) +
             (1 - totalX_int_pt) *
                 ((1 - Sw) * drhomolnonwet_dX - rho_mol_nonwet * dSwdX_gp +
                  Sw * drhomolwet_dX + dSwdX_gp * rho_mol_wet)) *
            mass_operator;
        /*Mgp.noalias() +=
            porosity *
            (totalX_int_pt *
            (-0.0 * rho_mol_nonwet + (1 - Sw) * drhomolnonwet_dpg +
                Sw * drhomolwet_dpg + 0.0 * rho_mol_wet)) *
            mass_operator;
        Mgx.noalias() +=
            porosity *
            ((1 - Sw) * rho_mol_nonwet + Sw * rho_mol_wet +
                totalX_int_pt *
                ((1 - Sw) * drhomolnonwet_dX - rho_mol_nonwet * 0.0 +
                    Sw * drhomolwet_dX + 0.0 * rho_mol_wet)) *
            mass_operator;
        //(1-X)*(N_G*(1-Sw)+N_L*Sw)
        Mlp.noalias() +=
            porosity *
            ((1 - totalX_int_pt) *
            (-0.0 * rho_mol_nonwet + (1 - Sw) * drhomolnonwet_dpg +
                Sw * drhomolwet_dpg + 0.0 * rho_mol_wet)) *
            mass_operator;
        Mlx.noalias() +=
            porosity *
            (-((1 - Sw) * rho_mol_nonwet + Sw * rho_mol_wet) +
            (1 - totalX_int_pt) *
                ((1 - Sw) * drhomolnonwet_dX - rho_mol_nonwet * 0.0 +
                    Sw * drhomolwet_dX + 0.0 * rho_mol_wet)) *
            mass_operator;*/
        double const k_rel_G =
            _process_data._material->getNonwetRelativePermeability(
                t, pos, pg_int_pt, _temperature, Sw);
        double const mu_gas =
            _process_data._material->getGasViscosity(pg_int_pt, _temperature);
        double const lambda_G = k_rel_G / mu_gas;
        double const diffusion_coeff_componenth2 =
            _process_data._diffusion_coeff_componentb(t, pos)[0];
        double const diffusion_coeff_componentw =
            _process_data._diffusion_coeff_componenta(t, pos)[0];
        // wet
        double const k_rel_L =
            _process_data._material->getWetRelativePermeability(
                t, pos, pg_int_pt, _temperature,
                Sw);  // interpolated_Kr_wet.getValue(Sw);
        double const mu_liquid = _process_data._material->getLiquidViscosity(
            _pressure_wetting[ip], _temperature);
        double const lambda_L = k_rel_L / mu_liquid;

        laplace_operator.noalias() =
            sm.dNdx.transpose() * permeability * sm.dNdx * integration_factor;
        // N_L*X_L^h*K*lambda_L*(PG-PC)+N_G*X_G^h*K*lambda_G+poro*Sw*D_L^h*d(X_L^h*N_L)
        /*double const kgp_coeff = (rho_mol_nonwet * X_h2_nonwet * lambda_G +
            rho_mol_wet * X_h2_wet * lambda_L * (1 - dPC_dSw_gp *
        dSwdP_gp))*5e-20
            + (rho_mol_wet * Sw * porosity * diffusion_coeff_componenth2 *
                dXh2wet_dpg);
        double const Kgx_coeff = (-rho_mol_wet * X_h2_wet * lambda_L *
        dPC_dSw_gp * dSwdX_gp) * 5e-20
            + (rho_mol_wet * Sw * porosity * diffusion_coeff_componenth2 *
                dXh2wet_dX);
        double const Klp_coeff = (rho_mol_nonwet * (1 - X_h2_nonwet) * lambda_G
        +
            rho_mol_w_std * lambda_L * (1 - dPC_dSw_gp * dSwdP_gp)) *5e-20
            + (-rho_mol_wet * Sw * porosity * diffusion_coeff_componenth2 *
                dXh2wet_dpg);
        double const Klx_coeff = (-rho_mol_w_std * lambda_L * dPC_dSw_gp *
        dSwdX_gp)*5e-20
            + (-
                rho_mol_wet * Sw * porosity *
                diffusion_coeff_componenth2 * dXh2wet_dX);
        std::cout << kgp_coeff << std::endl;
        std::cout << Kgx_coeff << std::endl;
        std::cout << Klp_coeff << std::endl;
        std::cout << Klx_coeff << std::endl;*/
        Kgp.noalias() +=
            (rho_mol_nonwet * X_h2_nonwet * lambda_G +
             rho_mol_wet * X_h2_wet * lambda_L * (1 - dPC_dSw_gp * dSwdP_gp)) *
                laplace_operator +
            (Sw * porosity * diffusion_coeff_componenth2 *
             (rho_mol_wet * dXh2wet_dpg + X_h2_wet * drhomolwet_dpg)) *
                sm.dNdx.transpose() * sm.dNdx * integration_factor;
        Kgx.noalias() +=
            (-rho_mol_wet * X_h2_wet * lambda_L * dPC_dSw_gp * dSwdX_gp) *
                laplace_operator +
            (Sw * porosity * diffusion_coeff_componenth2 *
             (rho_mol_wet * dXh2wet_dX + X_h2_wet * drhomolwet_dX)) *
                sm.dNdx.transpose() * sm.dNdx * integration_factor;
        Klp.noalias() +=
            (rho_mol_nonwet * (1 - X_h2_nonwet) * lambda_G +
             rho_mol_w_std * lambda_L * (1 - dPC_dSw_gp * dSwdP_gp)) *
                laplace_operator -
            (Sw * porosity * diffusion_coeff_componenth2 *
             (rho_mol_wet * dXh2wet_dpg + X_h2_wet * drhomolwet_dpg)) *
                sm.dNdx.transpose() * sm.dNdx * integration_factor;

        Klx.noalias() +=
            (-rho_mol_w_std * lambda_L * dPC_dSw_gp * dSwdX_gp) *
                laplace_operator -
            (Sw * porosity * diffusion_coeff_componenth2 *
             (rho_mol_wet * dXh2wet_dX + X_h2_wet * drhomolwet_dX)) *
                sm.dNdx.transpose() * sm.dNdx * integration_factor;

        double const rho_mass_wet =
            rho_mol_wet *
            (X_h2_wet * molar_mass_h2 + (1 - X_h2_wet) * molar_mass_h2o);
        double const rho_mass_nonwet = rho_gas;
        // rho_mol_nonwet * ((1 - x_h2o_in_nonwet) * molar_mass_co2 +
        // x_h2o_in_nonwet * molar_mass_h2o);
        if (_process_data._has_gravity)
        {
            auto const& b = _process_data._specific_body_force;
            Bg.noalias() +=
                (rho_mol_nonwet * X_h2_nonwet * rho_mass_nonwet * lambda_G +
                 rho_mol_wet * X_h2_wet * lambda_L * rho_mass_wet) *
                sm.dNdx.transpose() * permeability * b * integration_factor;
            Bl.noalias() +=
                ((rho_w / molar_mass_h2o) * lambda_L * rho_mass_wet +
                 rho_mol_nonwet * (1 - X_h2_nonwet) * lambda_G *
                     rho_mass_nonwet) *
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
