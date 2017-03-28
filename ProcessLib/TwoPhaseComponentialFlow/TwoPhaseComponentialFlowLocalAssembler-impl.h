/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TwoPhaseComponentialFlowLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseComponentialFlowProcessData.h"

using MaterialLib::PhysicalConstant::MolarMass::Water;
namespace ProcessLib
{
namespace TwoPhaseComponentialFlow
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseComponentialFlowLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, std::vector<double> const& local_x,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();
    auto const n_nodes = ShapeFunction::NPOINTS;
    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    Eigen::MatrixXd mass_mat_coeff =
        Eigen::MatrixXd::Zero(NUM_NODAL_DOF, NUM_NODAL_DOF);
    Eigen::MatrixXd K_mat_coeff = Eigen::MatrixXd::Zero(NUM_NODAL_DOF, NUM_NODAL_DOF);
    Eigen::VectorXd H_vec_coeff = Eigen::VectorXd::Zero(NUM_NODAL_DOF);
    Eigen::VectorXd F_vec_coeff = Eigen::VectorXd::Zero(NUM_NODAL_DOF);

    NodalMatrixType localMass_tmp; 
    localMass_tmp.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
    NodalMatrixType localDispersion_tmp;
    localDispersion_tmp.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
    Eigen::VectorXd localGravity_tmp =
        Eigen::VectorXd::Zero(ShapeFunction::NPOINTS);
    Eigen::VectorXd localSource_tmp =
        Eigen::VectorXd::Zero(ShapeFunction::NPOINTS);

    NodalMatrixType mass_operator;
    mass_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    NodalMatrixType laplace_operator;
    laplace_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
    NodalMatrixType diffusive_operator;
    diffusive_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id =
        _process_data._material->getMaterialID(pos.getElementID().get());

    const Eigen::MatrixXd& perm = _process_data._material->getPermeability(
        material_id, t, pos, _element.getDimension());
    assert(perm.rows() == _element.getDimension() || perm.rows() == 1);
    GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
        _element.getDimension(), _element.getDimension());
    if (perm.rows() == _element.getDimension())
        permeability = perm;
    else if (perm.rows() == 1)
        permeability.diagonal().setConstant(perm(0, 0));
    MathLib::PiecewiseLinearInterpolation const& interpolated_Q_slow =
        _process_data._interpolated_Q_slow;
    MathLib::PiecewiseLinearInterpolation const& interpolated_Q_fast =
        _process_data._interpolated_Q_fast;
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];

        double pg_int_pt = 0.0;
        double X1_int_pt = 0.0;
        double X2_int_pt = 0.0;
        double X3_int_pt = 0.0;
        double PC_int_pt = 0.0;

        NumLib::shapeFunctionInterpolate(local_x, sm.N, pg_int_pt, X1_int_pt,
            X2_int_pt, X3_int_pt, PC_int_pt);
        //PC_int_pt = 3.1245e+7;
        _pressure_wetting[ip] = pg_int_pt - PC_int_pt;
        const double dt = _process_data._dt;
        auto const& wp = _integration_method.getWeightedPoint(ip);
        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();

        const double temperature = _process_data._temperature(t, pos)[0];
        const double kelvin_term = exp(PC_int_pt*0.018 / rho_l_std / R / temperature);
        double X_L_h_gp = pg_int_pt * X1_int_pt / Hen_L_h;  // Henry law
        double X_L_c_gp = pg_int_pt * X2_int_pt / Hen_L_c;
        double X_L_co2_gp = pg_int_pt * X3_int_pt / Hen_L_co2;

        double P_sat_gp = get_P_sat(temperature);

        double K_G_w = pg_int_pt / P_sat_gp;  // henry law ratio
        double K_G_air = pg_int_pt / Hen_L_air;
        double L = 1 - (X_L_h_gp + X_L_c_gp + X_L_co2_gp);
        double G = 1 - X1_int_pt - X2_int_pt - X3_int_pt;
        double x_nonwet_air = get_x_nonwet_air_gp(pg_int_pt, X1_int_pt, X2_int_pt,
            X3_int_pt, P_sat_gp);

        /*double X_G_h2o_gp = get_X_G_h2o_gp(pg_int_pt, X1_int_pt, X2_int_pt,
        X3_int_pt, P_sat_gp);*/
        double const x_nonwet_h2o = get_x_nonwet_h2o(
            pg_int_pt, X1_int_pt, X2_int_pt, X3_int_pt, P_sat_gp,kelvin_term);

        double const x_wet_h2o = pg_int_pt * x_nonwet_h2o*kelvin_term / P_sat_gp;
        double const x_wet_air = pg_int_pt * x_nonwet_air / Hen_L_air;
        /*double const rho_gas =
            _process_data._material->getGasDensity(pg_int_pt, temperature);*/
        double const rho_w = _process_data._material->getLiquidDensity(
            _pressure_wetting[ip], temperature);

        double const Sw = _process_data._material->getSaturation(material_id,
            t, pos, pg_int_pt, temperature, PC_int_pt);
        double const S_G_gp = 1 - Sw;

        _saturation[ip] = Sw;
        double dSwdPc = _process_data._material->getDerivSaturation(material_id,
            t, pos, pg_int_pt, temperature, Sw);
        const double dSgdPC = -dSwdPc;

        const double rho_mol_nonwet = pg_int_pt / R / temperature;
        const double rho_mol_water = rho_l_std / M_L;
        const double rho_mass_G_gp =
            rho_mol_nonwet *
            (X1_int_pt * M_H + X2_int_pt * M_C + X3_int_pt * M_CO2 +
                x_nonwet_air * M_AIR + x_nonwet_h2o * M_L);

        double dLdPG =
            -X1_int_pt / Hen_L_h - X2_int_pt / Hen_L_c - X3_int_pt / Hen_L_co2;
        double d_x_nonwet_air_d_pg =
            ((G / P_sat_gp - dLdPG) * (K_G_w - K_G_air) -
            (K_G_w * G - L) * (1 / P_sat_gp - 1 / Hen_L_air)) /
                (K_G_w - K_G_air) / (K_G_w - K_G_air);

        double d_x_nonwet_h2o_d_pg =
            -((G / Hen_L_air - dLdPG) * (K_G_w - K_G_air) +
            (L - K_G_air * G) * (1 / P_sat_gp - 1 / Hen_L_air)) /
                (K_G_w - K_G_air) / (K_G_w - K_G_air);

        
        double const d_x_nonwet_air_d_x1 =
            ((-K_G_w + pg_int_pt / Hen_L_h) / (K_G_w - K_G_air));

        double d_x_nonwet_air_d_x2 =
            ((-K_G_w + pg_int_pt / Hen_L_c) / (K_G_w - K_G_air));

        double d_x_nonwet_air_d_x3 =
            ((-K_G_w + pg_int_pt / Hen_L_co2) / (K_G_w - K_G_air));

        double d_x_nonwet_h2o_d_x1 =
            ((K_G_air - pg_int_pt / Hen_L_h) / (K_G_w - K_G_air));
        /// double const d_x_nonwet_h2o_d_x1 =
        /// get_derivative_x_nonwet_h2o_d_x1(pg_int_pt, X1_int_pt, X2_int_pt,
        /// X3_int_pt, P_sat_gp);
        double d_x_nonwet_h2o_d_x2 =
            ((K_G_air - pg_int_pt / Hen_L_c) / (K_G_w - K_G_air));
        /// double const d_x_nonwet_h2o_d_x2 =
        /// get_derivative_x_nonwet_h2o_d_x2(pg_int_pt, X1_int_pt, X2_int_pt,
        /// X3_int_pt, P_sat_gp);
        double d_x_nonwet_h2o_d_x3 =
            ((K_G_air - pg_int_pt / Hen_L_co2) / (K_G_w - K_G_air));
        /// double const d_x_nonwet_h2o_d_x3 =
        /// get_derivative_x_nonwet_h2o_d_x3(pg_int_pt, X1_int_pt, X2_int_pt,
        /// X3_int_pt, P_sat_gp);
        double porosity = _process_data._material->getPorosity(material_id,
            t, pos, pg_int_pt, temperature, 0);
        //should be only valid for material 0
        double& porosity2 = _ip_data[ip].porosity;
        // Assemble M matrix
        // nonwetting
        double const d_rho_mol_nonwet_d_pg = 1 / R / temperature;

        double const rho_mol_wet = rho_mol_water/ x_wet_h2o;
        const double rho_mass_wet =
            rho_mol_wet *
            (X_L_h_gp * M_H + X_L_c_gp * M_C + X_L_co2_gp * M_CO2 +
                x_wet_air * M_AIR + x_wet_h2o * M_L);
        double const d_rho_mol_wet_d_pg = 
            -rho_mol_water * (x_nonwet_h2o / P_sat_gp +
                pg_int_pt * d_x_nonwet_h2o_d_pg / P_sat_gp) /
            x_wet_h2o / x_wet_h2o;
        double const d_rho_mol_wet_d_x1 = 
            -rho_mol_water * (pg_int_pt * d_x_nonwet_h2o_d_x1 / P_sat_gp) /
            x_wet_h2o / x_wet_h2o;
        double const d_rho_mol_wet_d_x2 = 
            -rho_mol_water * (pg_int_pt * d_x_nonwet_h2o_d_x2 / P_sat_gp) /
            x_wet_h2o / x_wet_h2o;
        double const d_rho_mol_wet_d_x3 = 
            -rho_mol_water * (pg_int_pt * d_x_nonwet_h2o_d_x3 / P_sat_gp) /
            x_wet_h2o / x_wet_h2o;

        //calculate the carbonation and ASR source/sink term
        double& rho_mol_sio2_wet = _ip_data[ip].rho_mol_sio2;
        if (_process_data._material->getMaterialID(pos.getElementID().get()) == 0)
        {
            porosity = porosity2;//use look-up table value  
        }
        double& rho_mol_co2_cumul_total =
            _ip_data[ip].rho_mol_co2_cumul_total;  // get cumulative co2
                                                   // calculate the current ammount of co2
        double rho_mol_total_co2 =
            _ip_data[ip].porosity_prev *
            (rho_mol_nonwet * X3_int_pt * (1 - Sw) +
                rho_mol_wet * X_L_co2_gp * Sw);

        //_process_data._material->getPorosity(material_id, t, pos, pn_int_pt, temperature, 0);
        double const pH = bi_interpolation(_ip_data[ip].rho_mol_sio2_prev,
            _ip_data[ip].rho_mol_co2_cumul_total_prev, _pH_at_supp_pnt);
        _pH_value[ip] = pH;
        //
        // double const d_rho_mol_wet_d_pg=
        mass_operator.noalias() = sm.N.transpose() * sm.N * integration_factor;
        // H2
        mass_mat_coeff(0, 0) =
            porosity * ((1 - Sw) * X1_int_pt * d_rho_mol_nonwet_d_pg +
                Sw * rho_mol_wet * X1_int_pt / Hen_L_h +
                Sw * X_L_h_gp * d_rho_mol_wet_d_pg);
        mass_mat_coeff(0, 1) =
            porosity *
            (rho_mol_nonwet * (1 - Sw) + Sw * rho_mol_wet * pg_int_pt / Hen_L_h +
                Sw * X_L_h_gp * d_rho_mol_wet_d_x1);
        mass_mat_coeff(0, 2) =
            porosity * (Sw * X_L_h_gp * d_rho_mol_wet_d_x2);
        mass_mat_coeff(0, 3) =
            porosity * (Sw * X_L_h_gp * d_rho_mol_wet_d_x3);
        mass_mat_coeff(0, 4) = porosity * (rho_mol_nonwet * X1_int_pt * dSgdPC -
            rho_mol_wet * X_L_h_gp * dSgdPC);
        //CH4
        mass_mat_coeff(1, 0) = porosity *
            ((1 - Sw) * X2_int_pt * d_rho_mol_nonwet_d_pg +
                Sw * rho_mol_wet * X2_int_pt / Hen_L_c +
                Sw * X_L_c_gp * d_rho_mol_wet_d_pg);
        mass_mat_coeff(1, 1) =
            porosity * (Sw * X_L_c_gp * d_rho_mol_wet_d_x1);
        mass_mat_coeff(1, 2) =
            porosity *
            (rho_mol_nonwet * (1 - Sw) + Sw * rho_mol_wet * pg_int_pt / Hen_L_c +
                Sw * X_L_c_gp * d_rho_mol_wet_d_x2);
        mass_mat_coeff(1, 3) =
            porosity * (Sw * X_L_c_gp * d_rho_mol_wet_d_x3);

        mass_mat_coeff(1, 4) = porosity * (rho_mol_nonwet * X2_int_pt * dSgdPC -
            rho_mol_wet * X_L_c_gp * dSgdPC);
        // co2
        mass_mat_coeff(2, 0) = porosity *
            ((1 - Sw) * X3_int_pt * d_rho_mol_nonwet_d_pg +
                Sw * rho_mol_wet * X3_int_pt / Hen_L_co2 +
                Sw * X_L_co2_gp * d_rho_mol_wet_d_pg);
        mass_mat_coeff(2, 1) = porosity *
            (Sw * X_L_co2_gp * d_rho_mol_wet_d_x1);
        mass_mat_coeff(2, 2) =
            porosity * (Sw * X_L_co2_gp * d_rho_mol_wet_d_x2);
        mass_mat_coeff(2, 3) = porosity *
            (rho_mol_nonwet * (1 - Sw) +
                Sw * rho_mol_wet * pg_int_pt / Hen_L_co2 +
                Sw * X_L_co2_gp * d_rho_mol_wet_d_x3);
        mass_mat_coeff(2, 4) = porosity * (rho_mol_nonwet * X3_int_pt * dSgdPC -
            rho_mol_wet * X_L_co2_gp * dSgdPC);
        // air
        mass_mat_coeff(3, 0) = porosity *
            ((1 - Sw) * x_nonwet_air * d_rho_mol_nonwet_d_pg +
            (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_pg +
                Sw * rho_mol_wet * d_x_nonwet_air_d_pg * K_G_air +
                Sw * rho_mol_wet * x_nonwet_air / Hen_L_air +
                Sw * x_wet_air * d_rho_mol_wet_d_pg);
        mass_mat_coeff(3, 1) =
            porosity * ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x1 +
                Sw * rho_mol_wet * K_G_air * d_x_nonwet_air_d_x1 +
                Sw * x_wet_air * d_rho_mol_wet_d_x1);
        mass_mat_coeff(3, 2) =
            porosity * ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x2 +
                Sw * rho_mol_wet * K_G_air * d_x_nonwet_air_d_x2 +
                Sw * x_wet_air * d_rho_mol_wet_d_x2);
        mass_mat_coeff(3, 3) =
            porosity * ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x3 +
                Sw * rho_mol_wet * K_G_air * d_x_nonwet_air_d_x3 +
                Sw * x_wet_air * d_rho_mol_wet_d_x3);
        mass_mat_coeff(3, 4) = porosity *
            (rho_mol_nonwet * x_nonwet_air * dSgdPC -
                rho_mol_wet * K_G_air * x_nonwet_air * dSgdPC);

                            // h2o
        mass_mat_coeff(4, 0) = porosity *
            ((1 - Sw) * d_rho_mol_nonwet_d_pg
                + Sw * d_rho_mol_wet_d_pg);
        mass_mat_coeff(4, 1) = porosity * Sw * d_rho_mol_wet_d_x1;
        mass_mat_coeff(4, 2) = porosity * Sw * d_rho_mol_wet_d_x2;
        mass_mat_coeff(4, 3) = porosity * Sw * d_rho_mol_wet_d_x3;
        mass_mat_coeff(4, 4) = porosity *
            (rho_mol_nonwet * dSgdPC -
                rho_mol_wet * dSgdPC);
        //-------------debugging------------------------
        // std::cout << "mass_mat_coeff=" << std::endl;
        // std::cout << mass_mat_coeff << std::endl;
        //--------------end debugging-------------------
        for (int ii = 0; ii < NUM_NODAL_DOF; ii++)
        {
            for (int jj = 0; jj < NUM_NODAL_DOF; jj++)
            {
                localMass_tmp.setZero();
                localMass_tmp.noalias() = mass_mat_coeff(ii, jj) * mass_operator;
                local_M.block(n_nodes * ii, n_nodes * jj, n_nodes, n_nodes)
                    .noalias() += localMass_tmp;
            }
        }

        double const k_rel_G =
            _process_data._material->getNonwetRelativePermeability(
                t, pos, pg_int_pt, temperature, Sw);
        double const mu_gas =
            _process_data._material->getGasViscosity(pg_int_pt, temperature);
        double const lambda_G = k_rel_G / mu_gas;
        // diffusion coefficient in water phase
        double const D_L = _process_data._diffusion_coeff_component_b(t, pos)[0];

        // wet
        double const k_rel_L =
            _process_data._material->getWetRelativePermeability(
                t, pos, _pressure_wetting[ip], temperature, Sw);
        double const mu_liquid = _process_data._material->getLiquidViscosity(
            _pressure_wetting[ip], temperature);
        double const lambda_L = k_rel_L / mu_liquid;
        // diffusion coefficient in gas phase
        double const D_G = _process_data._diffusion_coeff_component_a(t, pos)[0];

        laplace_operator.noalias() =
            sm.dNdx.transpose() * permeability * sm.dNdx * integration_factor;
        diffusive_operator.noalias() = sm.dNdx.transpose() * sm.dNdx * integration_factor;

        K_mat_coeff(0, 0) =
            (lambda_G * rho_mol_nonwet * X1_int_pt +
                lambda_L * rho_mol_wet * X_L_h_gp) * permeability(0, 0) +
                (porosity * D_L * Sw * rho_mol_wet * X1_int_pt / Hen_L_h);
        K_mat_coeff(0, 1) =
            (porosity * D_G * (1 - Sw) * rho_mol_nonwet +
                porosity * D_L * Sw * rho_mol_wet * pg_int_pt / Hen_L_h);
        K_mat_coeff(0, 2) = 0.0 ;
        K_mat_coeff(0, 3) = 0.0 ;
        K_mat_coeff(0, 4) =
            (-lambda_L * rho_mol_wet * X_L_h_gp)*permeability(0, 0);
        // ch4
        K_mat_coeff(1, 0) =
            (lambda_G * rho_mol_nonwet * X2_int_pt +
                lambda_L * rho_mol_wet * X_L_c_gp) * permeability(0, 0) +
                (porosity * D_L * Sw * rho_mol_wet * X2_int_pt / Hen_L_c);
        K_mat_coeff(1, 1) = 0.0;
        K_mat_coeff(1, 2) =
            (porosity * D_G * (1 - Sw) * rho_mol_nonwet +
                porosity * D_L * Sw * rho_mol_wet * pg_int_pt / Hen_L_c);
        K_mat_coeff(1, 3) = 0.0;
        K_mat_coeff(1, 4) =
            (-lambda_L * rho_mol_wet * X_L_c_gp) * permeability(0, 0);
        // co2
        K_mat_coeff(2, 0) =
            (lambda_G * rho_mol_nonwet * X3_int_pt +
                lambda_L * rho_mol_wet * X_L_co2_gp) * permeability(0,0)+
            (porosity * D_L * Sw * rho_mol_wet * X3_int_pt / Hen_L_co2);
        K_mat_coeff(2, 1) = 0.0;
        K_mat_coeff(2, 2) = 0.0;
        K_mat_coeff(2, 3) =
            (porosity * D_G * (1 - Sw) * rho_mol_nonwet +
                porosity * D_L * Sw * rho_mol_wet * pg_int_pt / Hen_L_co2);
        K_mat_coeff(2, 4) =
            (-lambda_L * rho_mol_wet * X_L_co2_gp) * permeability(0,0);
        // air
        K_mat_coeff(3, 0) =
            (lambda_G * rho_mol_nonwet * x_nonwet_air +
                lambda_L * rho_mol_wet * x_wet_air) * permeability(0, 0) +
                (porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_pg +
                    porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_pg * K_G_air +
                    porosity * D_L * Sw * rho_mol_wet * x_nonwet_air / Hen_L_air);
        K_mat_coeff(3, 1) =
            (porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x1 +
                porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_x1 *
                K_G_air);
        K_mat_coeff(3, 2) =
            (porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x2 +
                porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_x2 *
                K_G_air);
        K_mat_coeff(3, 3) =
            (porosity * D_G * S_G_gp * rho_mol_nonwet * d_x_nonwet_air_d_x3 +
                porosity * D_L * (1 - S_G_gp) * rho_mol_wet * d_x_nonwet_air_d_x3 *
                K_G_air);
        K_mat_coeff(3, 4) =
            (-lambda_L * rho_mol_wet * K_G_air * x_nonwet_air) * permeability(0,0);
        // h2o
        K_mat_coeff(4, 0) =
            (lambda_G * rho_mol_nonwet +
                lambda_L * rho_mol_wet) *
            permeability(0, 0);

        K_mat_coeff(4, 1) = 0.0;
        K_mat_coeff(4, 2) = 0.0;
        K_mat_coeff(4, 3) = 0.0;
        K_mat_coeff(4, 4) =
            (-lambda_L * rho_mol_wet) * permeability(0,0);

        //-------------debugging------------------------
        // std::cout << "K_mat_coeff=" << std::endl;
        // std::cout << K_mat_coeff << std::endl;
        //--------------end debugging-------------------

        for (int ii = 0; ii < NUM_NODAL_DOF; ii++)
        {
            for (int jj = 0; jj < NUM_NODAL_DOF; jj++)
            {
                localDispersion_tmp.setZero();
                localDispersion_tmp.noalias() =
                    K_mat_coeff(ii, jj) * diffusive_operator;
                local_K.block(n_nodes * ii, n_nodes * jj, n_nodes, n_nodes)
                    .noalias() += localDispersion_tmp;
            }
        }

        if (_process_data._has_gravity)
        {
            auto const& b = _process_data._specific_body_force;
            NodalVectorType gravity_operator = sm.dNdx.transpose() *
                permeability * b *
                integration_factor;
            H_vec_coeff(0) =
                (-lambda_G * rho_mol_nonwet * X1_int_pt * rho_mass_G_gp -
                    lambda_L * rho_mol_wet * X1_int_pt * pg_int_pt * rho_mass_wet /
                    Hen_L_h);
            H_vec_coeff(1) =
                (-lambda_G * rho_mol_nonwet * X2_int_pt * rho_mass_G_gp -
                    lambda_L * rho_mol_wet * X2_int_pt * pg_int_pt * rho_mass_wet /
                    Hen_L_c);
            H_vec_coeff(2) =
                (-lambda_G * rho_mol_nonwet * X3_int_pt * rho_mass_G_gp -
                    lambda_L * rho_mol_wet * X3_int_pt * pg_int_pt * rho_mass_wet /
                    Hen_L_co2);
            H_vec_coeff(3) =
                (-lambda_G * rho_mol_nonwet * x_nonwet_air * rho_mass_G_gp -
                    lambda_L * rho_mol_wet * x_nonwet_air* pg_int_pt * rho_mass_wet /
                    Hen_L_air);
            H_vec_coeff(4) =
                (-lambda_G * rho_mol_nonwet  * rho_mass_G_gp
                    - lambda_L * rho_mol_wet * rho_mass_wet);
            for (int idx = 0; idx < NUM_NODAL_DOF; idx++)
            {
                // since no primary vairable involved
                // directly assemble to the Right-Hand-Side
                // F += dNp^T * K * gz
                localGravity_tmp.setZero();
                localGravity_tmp.noalias() = H_vec_coeff(idx) * gravity_operator;
                local_b.block(n_nodes * idx, 0, n_nodes, 1).noalias() +=
                    localGravity_tmp;
            }
        }  // end of hasGravityEffect
           // load the source term
        double const flag_carbon = bi_interpolation(
            _ip_data[ip].rho_mol_sio2_prev,
            _ip_data[ip].rho_mol_co2_cumul_total_prev, _flag_carbon_suppt_pnt);
        double& fluid_volume = _ip_data[ip].fluid_volume;
        fluid_volume = bi_interpolation(
            _ip_data[ip].rho_mol_sio2_prev,
            _ip_data[ip].rho_mol_co2_cumul_total_prev, _fluid_volume_suppt_pnt);
        double quartz_dissolute_rate = bi_interpolation(
            _ip_data[ip].rho_mol_sio2_prev,
            _ip_data[ip].rho_mol_co2_cumul_total_prev, _quartz_rate_suppt_pnt);
        // quartz_dissolute_rate is always nonpositive.
        if (quartz_dissolute_rate > 0)
            quartz_dissolute_rate = 0;
        if (Sw > 0.3 && dt > 0)
        {
            Eigen::VectorXd F_vec_coeff = Eigen::VectorXd::Zero(NUM_NODAL_DOF);
            double Q_organic_slow_co2_ini =
                interpolated_Q_slow.getValue(t);  // read from curves
            double Q_organic_fast_co2_ini =
                interpolated_Q_fast.getValue(t);  // read from curves
            if (_process_data._material->getMaterialID(pos.getElementID().get()) == 1)
            {

                F_vec_coeff(0) = Q_steel;

                const double Q_organic_slow_co2 =
                    Q_organic_slow_co2_ini * para_slow;

                const double Q_organic_fast_co2 =
                    Q_organic_fast_co2_ini * para_fast;

                F_vec_coeff(1) += (Q_organic_slow_co2 * 31 / 19);
                F_vec_coeff(2) += Q_organic_slow_co2;

                F_vec_coeff(1) += Q_organic_fast_co2;

                F_vec_coeff(2) += Q_organic_fast_co2;

                F_vec_coeff(4) = 
                    (Q_organic_slow_co2 * 12 / 19) + (Q_organic_fast_co2 * 5 / 3);
                _porosity_value[ip] = porosity;
            }
            else if (_process_data._material->getMaterialID(pos.getElementID().get()) == 0)
            {
                const double rho_co2_ele =
                    porosity * ((1 - Sw) * rho_mol_nonwet * X3_int_pt +
                        Sw * rho_mol_wet * X_L_co2_gp);

                double const fluid_volume_rate =
                    (fluid_volume - _ip_data[ip].fluid_volume_prev) / dt;
                if (_ip_data[ip].rho_mol_co2_cumul_total_prev >= 3800)  // means carbonation stops
                    rho_mol_total_co2 = 0.0;
                // update the current cumulated co2 consumption
                rho_mol_co2_cumul_total =
                    _ip_data[ip].rho_mol_co2_cumul_total_prev + rho_mol_total_co2;
                // co2 consumption
                F_vec_coeff(2) -= (rho_mol_total_co2 / dt);
                    // water source/sink term
                F_vec_coeff(4) += (fluid_volume_rate)-(rho_mol_total_co2 / dt);
                // update the amount of dissolved sio2
                rho_mol_sio2_wet =
                    _ip_data[ip].rho_mol_sio2_prev -
                    quartz_dissolute_rate * dt;  // cumulative dissolved sio2
                //porosity update
                porosity2 =
                    bi_interpolation(_ip_data[ip].rho_mol_sio2_prev,
                        _ip_data[ip].rho_mol_co2_cumul_total_prev, _porosity_at_supp_pnts);//porosity update
                _porosity_value[ip] = porosity2;
            }
            for (int idx = 0; idx < NUM_NODAL_DOF; idx++)
            {
                // since no primary vairable involved
                // directly assemble to the Right-Hand-Side
                // F += dNp^T * K * gz
                localSource_tmp.setZero();
                localSource_tmp.noalias() = sm.N.transpose() * F_vec_coeff(idx) * integration_factor;
                local_b.block(n_nodes * idx, 0, n_nodes, 1).noalias() +=
                    localSource_tmp;
            }
        }
    }// end of GP asm
    
    if (_process_data._has_mass_lumping)
    {
        auto Mhpg =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);

        auto Mhmolh2 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                nonwet_pressure_matrix_index,
                nonwet_pressure_size * mol_fraction_h_coeff_index);

        auto Mhmolch4 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                nonwet_pressure_matrix_index,
                nonwet_pressure_size * mol_fraction_ch4_coeff_index);
        auto Mhmolco2 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                nonwet_pressure_matrix_index,
                nonwet_pressure_size * mol_fraction_co2_coeff_index);

        auto Mhpc =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                nonwet_pressure_matrix_index,
                nonwet_pressure_size * cap_pressure_coeff_index);

        auto Mch4pg =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                nonwet_pressure_size, nonwet_pressure_matrix_index);

        auto Mch4molh2 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_h_coeff_index);

        auto Mch4molch4 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_ch4_coeff_index);
        auto Mch4molco2 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_co2_coeff_index);

        auto Mch4pc =
            local_M.block<nonwet_pressure_size, cap_pressure_size>(
                nonwet_pressure_size,
                nonwet_pressure_size * cap_pressure_coeff_index);

        auto Mco2pg =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                2 * nonwet_pressure_size, nonwet_pressure_matrix_index);

        auto Mco2molh2 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                2 * nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_h_coeff_index);

        auto Mco2molch4 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                2 * nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_ch4_coeff_index);
        auto Mco2molco2 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                2 * nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_co2_coeff_index);

        auto Mco2pc =
            local_M.block<nonwet_pressure_size, cap_pressure_size>(
                2 * nonwet_pressure_size,
                nonwet_pressure_size * cap_pressure_coeff_index);

        auto Mairpg =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                3 * nonwet_pressure_size, nonwet_pressure_matrix_index);

        auto Mairmolh2 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                3 * nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_h_coeff_index);

        auto Mairmolch4 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                3 * nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_ch4_coeff_index);
        auto Mairmolco2 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                3 * nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_co2_coeff_index);

        auto Mairpc =
            local_M.block<nonwet_pressure_size, cap_pressure_size>(
                3 * nonwet_pressure_size,
                nonwet_pressure_size * cap_pressure_coeff_index);

        auto Mh2opg =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                4 * nonwet_pressure_size, nonwet_pressure_matrix_index);

        auto Mh2omolh2 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                4 * nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_h_coeff_index);

        auto Mh2omolch4 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                4 * nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_ch4_coeff_index);
        auto Mh2omolco2 =
            local_M.block<nonwet_pressure_size, nonwet_pressure_size>(
                4 * nonwet_pressure_size,
                nonwet_pressure_size * mol_fraction_co2_coeff_index);

        auto Mh2opc =
            local_M.block<nonwet_pressure_size, cap_pressure_size>(
                4 * nonwet_pressure_size,
                nonwet_pressure_size * cap_pressure_coeff_index);
        for (unsigned row = 0; row < Mhpg.cols(); row++)
        {
            for (unsigned column = 0; column < Mhpg.cols(); column++)
            {
                if (row != column)
                {
                    Mhpg(row, row) += Mhpg(row, column);
                    Mhpg(row, column) = 0.0;
                    Mhmolh2(row, row) += Mhmolh2(row, column);
                    Mhmolh2(row, column) = 0.0;
                    Mhmolch4(row, row) += Mhmolch4(row, column);
                    Mhmolch4(row, column) = 0.0;
                    Mhmolco2(row, row) += Mhmolco2(row, column);
                    Mhmolco2(row, column) = 0.0;
                    Mhpc(row, row) += Mhpc(row, column);
                    Mhpc(row, column) = 0.0;

                    Mch4pg(row, row) += Mch4pg(row, column);
                    Mch4pg(row, column) = 0.0;
                    Mch4molh2(row, row) += Mch4molh2(row, column);
                    Mch4molh2(row, column) = 0.0;
                    Mch4molch4(row, row) += Mch4molch4(row, column);
                    Mch4molch4(row, column) = 0.0;
                    Mch4molco2(row, row) += Mch4molco2(row, column);
                    Mch4molco2(row, column) = 0.0;
                    Mch4pc(row, row) += Mch4pc(row, column);
                    Mch4pc(row, column) = 0.0;

                    Mco2pg(row, row) += Mco2pg(row, column);
                    Mco2pg(row, column) = 0.0;
                    Mco2molh2(row, row) += Mco2molh2(row, column);
                    Mco2molh2(row, column) = 0.0;
                    Mco2molch4(row, row) += Mco2molch4(row, column);
                    Mco2molch4(row, column) = 0.0;
                    Mco2molco2(row, row) += Mco2molco2(row, column);
                    Mco2molco2(row, column) = 0.0;
                    Mco2pc(row, row) += Mco2pc(row, column);
                    Mco2pc(row, column) = 0.0;

                    Mairpg(row, row) += Mairpg(row, column);
                    Mairpg(row, column) = 0.0;
                    Mairmolh2(row, row) += Mairmolh2(row, column);
                    Mairmolh2(row, column) = 0.0;
                    Mairmolch4(row, row) += Mairmolch4(row, column);
                    Mairmolch4(row, column) = 0.0;
                    Mairmolco2(row, row) += Mairmolco2(row, column);
                    Mairmolco2(row, column) = 0.0;
                    Mairpc(row, row) += Mairpc(row, column);
                    Mairpc(row, column) = 0.0;

                    Mh2opg(row, row) += Mh2opg(row, column);
                    Mh2opg(row, column) = 0.0;
                    Mh2omolh2(row, row) += Mh2omolh2(row, column);
                    Mh2omolh2(row, column) = 0.0;
                    Mh2omolch4(row, row) += Mh2omolch4(row, column);
                    Mh2omolch4(row, column) = 0.0;
                    Mh2omolco2(row, row) += Mh2omolco2(row, column);
                    Mh2omolco2(row, column) = 0.0;
                    Mh2opc(row, row) += Mh2opc(row, column);
                    Mh2opc(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // end of namespace
}  // end of namespace
