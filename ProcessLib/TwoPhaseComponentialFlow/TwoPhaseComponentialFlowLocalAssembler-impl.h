/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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

    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    NodalMatrixType mass_operator;
    mass_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Mhpg =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Mhmolh2 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index,
            nonwet_pressure_size * mol_fraction_h_coeff_index);

    auto Mhmolch4 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index,
            nonwet_pressure_size * mol_fraction_ch4_coeff_index);
    auto Mhmolco2 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index,
            nonwet_pressure_size * mol_fraction_co2_coeff_index);

    auto Mhpc =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index,
            nonwet_pressure_size * cap_pressure_coeff_index);

    auto Mch4pg =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_size, nonwet_pressure_matrix_index);

    auto Mch4molh2 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_h_coeff_index);

    auto Mch4molch4 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_ch4_coeff_index);
    auto Mch4molco2 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_co2_coeff_index);

    auto Mch4pc =
        local_M.template block<nonwet_pressure_size, cap_pressure_size>(
            nonwet_pressure_size,
            nonwet_pressure_size * cap_pressure_coeff_index);

    auto Mco2pg =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            2 * nonwet_pressure_size, nonwet_pressure_matrix_index);

    auto Mco2molh2 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            2 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_h_coeff_index);

    auto Mco2molch4 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            2 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_ch4_coeff_index);
    auto Mco2molco2 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            2 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_co2_coeff_index);

    auto Mco2pc =
        local_M.template block<nonwet_pressure_size, cap_pressure_size>(
            2 * nonwet_pressure_size,
            nonwet_pressure_size * cap_pressure_coeff_index);

    auto Mairpg =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            3 * nonwet_pressure_size, nonwet_pressure_matrix_index);

    auto Mairmolh2 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            3 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_h_coeff_index);

    auto Mairmolch4 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            3 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_ch4_coeff_index);
    auto Mairmolco2 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            3 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_co2_coeff_index);

    auto Mairpc =
        local_M.template block<nonwet_pressure_size, cap_pressure_size>(
            3 * nonwet_pressure_size,
            nonwet_pressure_size * cap_pressure_coeff_index);

    auto Mh2opg =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            4 * nonwet_pressure_size, nonwet_pressure_matrix_index);

    auto Mh2omolh2 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            4 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_h_coeff_index);

    auto Mh2omolch4 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            4 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_ch4_coeff_index);
    auto Mh2omolco2 =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            4 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_co2_coeff_index);

    auto Mh2opc =
        local_M.template block<nonwet_pressure_size, cap_pressure_size>(
            4 * nonwet_pressure_size,
            nonwet_pressure_size * cap_pressure_coeff_index);

    NodalMatrixType laplace_operator;
    laplace_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
    NodalMatrixType diffusive_operator;
    diffusive_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
    auto Khpg =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Khmolh2 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index,
            nonwet_pressure_size * mol_fraction_h_coeff_index);

    auto Khmolch4 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index,
            nonwet_pressure_size * mol_fraction_ch4_coeff_index);
    auto Khmolco2 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index,
            nonwet_pressure_size * mol_fraction_co2_coeff_index);

    auto Khpc = local_K.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index,
        nonwet_pressure_size * cap_pressure_coeff_index);

    auto Kch4pg =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_size, nonwet_pressure_matrix_index);

    auto Kch4molh2 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_h_coeff_index);

    auto Kch4molch4 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_ch4_coeff_index);
    auto Kch4molco2 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_co2_coeff_index);

    auto Kch4pc =
        local_K.template block<nonwet_pressure_size, cap_pressure_size>(
            nonwet_pressure_size,
            nonwet_pressure_size * cap_pressure_coeff_index);

    auto Kco2pg =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            2 * nonwet_pressure_size, nonwet_pressure_matrix_index);

    auto Kco2molh2 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            2 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_h_coeff_index);

    auto Kco2molch4 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            2 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_ch4_coeff_index);
    auto Kco2molco2 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            2 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_co2_coeff_index);

    auto Kco2pc =
        local_K.template block<nonwet_pressure_size, cap_pressure_size>(
            2 * nonwet_pressure_size,
            nonwet_pressure_size * cap_pressure_coeff_index);

    auto Kairpg =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            3 * nonwet_pressure_size, nonwet_pressure_matrix_index);

    auto Kairmolh2 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            3 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_h_coeff_index);

    auto Kairmolch4 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            3 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_ch4_coeff_index);
    auto Kairmolco2 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            3 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_co2_coeff_index);

    auto Kairpc =
        local_K.template block<nonwet_pressure_size, cap_pressure_size>(
            3 * nonwet_pressure_size,
            nonwet_pressure_size * cap_pressure_coeff_index);

    auto Kh2opg =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            4 * nonwet_pressure_size, nonwet_pressure_matrix_index);

    auto Kh2omolh2 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            4 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_h_coeff_index);

    auto Kh2omolch4 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            4 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_ch4_coeff_index);
    auto Kh2omolco2 =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            4 * nonwet_pressure_size,
            nonwet_pressure_size * mol_fraction_co2_coeff_index);

    auto Kh2opc =
        local_K.template block<nonwet_pressure_size, cap_pressure_size>(
            4 * nonwet_pressure_size,
            nonwet_pressure_size * cap_pressure_coeff_index);

    auto Bh2 = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);

    auto Bch4 = local_b.template segment<nonwet_pressure_size>(
        mol_fraction_h2_matrix_index);

    auto Bco2 = local_b.template segment<nonwet_pressure_size>(
        mol_fraction_ch4_matrix_index);
    auto Bair = local_b.template segment<nonwet_pressure_size>(
        mol_fraction_co2_matrix_index);
    auto Bh2o = local_b.template segment<nonwet_pressure_size>(
        cap_pressure_matrix_index);

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
    MathLib::PiecewiseLinearInterpolation const& interpolated_Q_slow =
        _process_data._interpolated_Q_slow;
    MathLib::PiecewiseLinearInterpolation const& interpolated_Q_fast =
        _process_data._interpolated_Q_fast;

    // Note: currently only isothermal case is considered, so the temperature is
    // assumed to be const
    // the variation of temperature will be taken into account in future
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

        _pressure_wetting[ip] = pg_int_pt - PC_int_pt;
        const double dt = _process_data._dt;
        auto const& wp = _integration_method.getWeightedPoint(ip);
        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();

        double X_L_h_gp = pg_int_pt * X1_int_pt / Hen_L_h;  // Henry law
        double X_L_c_gp = pg_int_pt * X2_int_pt / Hen_L_c;
        double X_L_co2_gp = pg_int_pt * X3_int_pt / Hen_L_co2;

        double P_sat_gp = get_P_sat(T_0);

        double K_G_w = pg_int_pt / P_sat_gp;  // henry law ratio
        double K_G_air = pg_int_pt / Hen_L_air;
        double L = 1 - (X_L_h_gp + X_L_c_gp + X_L_co2_gp);
        double G = 1 - X1_int_pt - X2_int_pt - X3_int_pt;
        double x_nonwet_air = get_x_nonwet_air_gp(pg_int_pt, X1_int_pt, X2_int_pt,
                                           X3_int_pt, P_sat_gp);

        /*double X_G_h2o_gp = get_X_G_h2o_gp(pg_int_pt, X1_int_pt, X2_int_pt,
                                           X3_int_pt, P_sat_gp);*/
        double const x_nonwet_h2o = get_x_nonwet_h2o(
            pg_int_pt, X1_int_pt, X2_int_pt, X3_int_pt, P_sat_gp);

        double const x_wet_h2o = pg_int_pt * x_nonwet_h2o / P_sat_gp;
        double const x_wet_air = pg_int_pt * x_nonwet_air / Hen_L_air;
        double const rho_gas =
            _process_data._material->getGasDensity(pg_int_pt, _temperature);
        double const rho_w = _process_data._material->getLiquidDensity(
            _pressure_wetting[ip], _temperature);

        double const Sw = _process_data._material->getSaturation(
            t, pos, pg_int_pt, _temperature, PC_int_pt);
        double const S_G_gp = 1 - Sw;

        _saturation[ip] = Sw;
        double dSwdPc = _process_data._material->getDerivSaturation(
            t, pos, pg_int_pt, _temperature, Sw);
        const double dSgdPC = -dSwdPc;
        /*if (pc_int_pt > interpolated_Pc.getSupportMax())
            dSwdPc =
                interpolated_Pc.getDerivative(interpolated_Pc.getSupportMax());
        else if (pc_int_pt < interpolated_Pc.getSupportMin())
            dSwdPc =
                interpolated_Pc.getDerivative(interpolated_Pc.getSupportMin());*/
        const double rho_mol_nonwet = pg_int_pt / R / T_0;
        const double rho_mol_water = rho_l_std / M_L;
        const double rho_mass_G_gp =
            rho_mol_nonwet *
            (X1_int_pt * M_H + X2_int_pt * M_C + X3_int_pt * M_CO2 +
             x_nonwet_air * M_AIR + x_nonwet_h2o * M_L);
        const double rho_mass_L_gp = rho_l_std;

        double Q_organic_slow_co2_ini =
            interpolated_Q_slow.getValue(t);  // read from curves
        double Q_organic_fast_co2_ini =
            interpolated_Q_fast.getValue(t);  // read from curves

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
        double const porosity = _process_data._material->getPorosity(
            t, pos, pg_int_pt, _temperature, 0);

        // Assemble M matrix
        // nonwetting
        double const drhogas_dpg = _process_data._material->getDerivGasDensity(
            pg_int_pt, _temperature);
        double const d_rho_mol_nonwet_d_pg = 1 / R / T_0;
        
        double const rho_mol_wet = rho_mol_water / x_wet_h2o;
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

        // double const d_rho_mol_wet_d_pg=
        mass_operator.noalias() = sm.N.transpose() * sm.N * integration_factor;
        // H2
        Mhpg.noalias() +=
            porosity * ((1 - Sw) * X1_int_pt * d_rho_mol_nonwet_d_pg +
                         Sw * rho_mol_wet * X1_int_pt / Hen_L_h +
                         Sw * X_L_h_gp * d_rho_mol_wet_d_pg) *
            mass_operator;  // dPG
        Mhmolh2.noalias() +=
            porosity *
            (rho_mol_nonwet * (1 - Sw) + Sw * rho_mol_wet * pg_int_pt / Hen_L_h +
             Sw * X_L_h_gp * d_rho_mol_wet_d_x1) *
            mass_operator;  // dX1 h2
        Mhmolch4.noalias() +=
            porosity * (Sw * X_L_h_gp * d_rho_mol_wet_d_x2) * mass_operator;
        Mhmolco2.noalias() +=
            porosity * (Sw * X_L_h_gp * d_rho_mol_wet_d_x3) * mass_operator;
        Mhpc.noalias() += porosity * (rho_mol_nonwet * X1_int_pt * dSgdPC -
                                      rho_mol_wet * X_L_h_gp * dSgdPC) *
                          mass_operator;  // dPC
        // CH4
        Mch4pg.noalias() += porosity *
                            ((1 - Sw) * X2_int_pt * d_rho_mol_nonwet_d_pg +
                             Sw * rho_mol_wet * X2_int_pt / Hen_L_c +
                             Sw * X_L_c_gp * d_rho_mol_wet_d_pg) *
                            mass_operator;  // dPG
        Mch4molh2.noalias() +=
            porosity * (Sw * X_L_c_gp * d_rho_mol_wet_d_x1) * mass_operator;
        Mch4molch4.noalias() +=
            porosity *
            (rho_mol_nonwet * (1 - Sw) + Sw * rho_mol_wet * pg_int_pt / Hen_L_c +
             Sw * X_L_c_gp * d_rho_mol_wet_d_x2) *
            mass_operator;  // Dx2 ch4
        Mch4molco2.noalias() +=
            porosity * (Sw * X_L_c_gp * d_rho_mol_wet_d_x3) * mass_operator;
        
        Mch4pc.noalias() += porosity * (rho_mol_nonwet * X2_int_pt * dSgdPC -
                                        rho_mol_wet * X_L_c_gp * dSgdPC) *
                            mass_operator;  // dPC
        // co2
        Mco2pg.noalias() += porosity *
                            ((1 - Sw) * X3_int_pt * d_rho_mol_nonwet_d_pg +
                             Sw * rho_mol_wet * X3_int_pt / Hen_L_co2 +
                             Sw * X_L_co2_gp * d_rho_mol_wet_d_pg) *
                            mass_operator;  // dPG
        Mco2molh2.noalias() += porosity *
                               (Sw * X_L_co2_gp * d_rho_mol_wet_d_x1) *
                               mass_operator;  // dX1 h2
        Mco2molch4.noalias() +=
            porosity * (Sw * X_L_co2_gp * d_rho_mol_wet_d_x2) * mass_operator;
        Mco2molco2.noalias() += porosity *
                                (rho_mol_nonwet * (1 - Sw) +
                                 Sw * rho_mol_wet * pg_int_pt / Hen_L_co2 +
                                 Sw * X_L_co2_gp * d_rho_mol_wet_d_x3) *
                                mass_operator;  // dx3 CO2
        Mco2pc.noalias() += porosity * (rho_mol_nonwet * X3_int_pt * dSgdPC -
                                        rho_mol_wet * X_L_co2_gp * dSgdPC) *
                            mass_operator;  // dPC
        // air
        Mairpg.noalias() += porosity *
                            ((1 - Sw) * x_nonwet_air * d_rho_mol_nonwet_d_pg +
                             (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_pg +
                             Sw * rho_mol_wet * d_x_nonwet_air_d_pg * K_G_air +
                             Sw * rho_mol_wet * x_nonwet_air / Hen_L_air +
                             Sw * x_wet_air * d_rho_mol_wet_d_pg) *
                            mass_operator;  // dPG
        Mairmolh2.noalias() +=
            porosity * ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x1 +
                        Sw * rho_mol_wet * K_G_air * d_x_nonwet_air_d_x1 +
                        Sw * x_wet_air * d_rho_mol_wet_d_x1) *
            mass_operator;  // dX1 h2
        Mairmolch4.noalias() +=
            porosity * ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x2 +
                        Sw * rho_mol_wet * K_G_air * d_x_nonwet_air_d_x2 +
                        Sw * x_wet_air * d_rho_mol_wet_d_x2) *
            mass_operator;  // dX2 CH4
        Mairmolco2.noalias() +=
            porosity * ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x3 +
                        Sw * rho_mol_wet * K_G_air * d_x_nonwet_air_d_x3 +
                        Sw * x_wet_air * d_rho_mol_wet_d_x3) *
            mass_operator;  // dX3 co2
        Mairpc.noalias() += porosity *
                            (rho_mol_nonwet * x_nonwet_air * dSgdPC -
                             rho_mol_wet * K_G_air * x_nonwet_air * dSgdPC) *
                            mass_operator;  // dPC

        // h2o
        Mh2opg.noalias() += porosity *
                            ((1 - Sw) * x_nonwet_h2o * d_rho_mol_nonwet_d_pg +
                             (1 - Sw) * rho_mol_nonwet * d_x_nonwet_h2o_d_pg +
                             Sw * rho_mol_wet * d_x_nonwet_h2o_d_pg * K_G_w +
                             Sw * rho_mol_wet * x_nonwet_h2o / P_sat_gp +
                             Sw * x_wet_h2o * d_rho_mol_wet_d_pg) *
                            mass_operator;  // dPG
        Mh2omolh2.noalias() += porosity *
                               ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_h2o_d_x1 +
                                Sw * rho_mol_wet * K_G_w * d_x_nonwet_h2o_d_x1 +
                                Sw * x_wet_h2o * d_rho_mol_wet_d_x1) *
                               mass_operator;  // dX1 h2
        Mh2omolch4.noalias() +=
            porosity * ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_h2o_d_x2 +
                        Sw * rho_mol_wet * K_G_w * d_x_nonwet_h2o_d_x2 +
                        Sw * x_wet_h2o * d_rho_mol_wet_d_x2) *
            mass_operator;  // dX2 CH4
        Mh2omolco2.noalias() +=
            porosity * ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_h2o_d_x3 +
                        Sw * rho_mol_wet * K_G_w * d_x_nonwet_h2o_d_x3 +
                        Sw * x_wet_h2o * d_rho_mol_wet_d_x3) *
            mass_operator;  // dX3 co2
        Mh2opc.noalias() += porosity *
                            (rho_mol_nonwet * x_nonwet_h2o * dSgdPC -
                             rho_mol_wet * K_G_w * x_nonwet_h2o * dSgdPC) *
                            mass_operator;  // dPC

        // Assemble M matrix
        // nonwet
        double const k_rel_G =
            _process_data._material->getNonwetRelativePermeability(
                t, pos, pg_int_pt, _temperature, Sw);
        double const mu_gas =
            _process_data._material->getGasViscosity(pg_int_pt, _temperature);
        double const lambda_G = k_rel_G / mu_gas;
        // diffusion coefficient in water phase
        double const D_L = _process_data._diffusion_coeff_componentb(t, pos)[0];

        // wet
        double const k_rel_L =
            _process_data._material->getWetRelativePermeability(
                t, pos, _pressure_wetting[ip], _temperature, Sw);
        double const mu_liquid = _process_data._material->getLiquidViscosity(
            _pressure_wetting[ip], _temperature);
        double const lambda_L = k_rel_L / mu_liquid;
        // diffusion coefficient in gas phase
        double const D_G = _process_data._diffusion_coeff_componenta(t, pos)[0];

        laplace_operator.noalias() =
            sm.dNdx.transpose() * permeability * sm.dNdx * integration_factor;
        diffusive_operator.noalias() = sm.dNdx.transpose() * sm.dNdx * integration_factor;
        Khpg.noalias() +=
            (lambda_G * rho_mol_nonwet * X1_int_pt +
             lambda_L * rho_mol_wet * X_L_h_gp) *
                laplace_operator +
            (porosity * D_L * Sw * rho_mol_wet * X1_int_pt / Hen_L_h) *
            diffusive_operator;
        Khmolh2.noalias() +=
            (porosity * D_G * (1 - Sw) * rho_mol_nonwet +
             porosity * D_L * Sw * rho_mol_wet * pg_int_pt / Hen_L_h) *
            diffusive_operator;
        Khmolch4.noalias() += 0.0 * laplace_operator;
        Khmolco2.noalias() += 0.0 * laplace_operator;
        Khpc.noalias() +=
            (-lambda_L * rho_mol_wet * X_L_h_gp) * laplace_operator;
        // ch4
        Kch4pg.noalias() +=
            (lambda_G * rho_mol_nonwet * X2_int_pt +
             lambda_L * rho_mol_wet * X_L_c_gp) *
                laplace_operator +
            (porosity * D_L * Sw * rho_mol_wet * X2_int_pt / Hen_L_c) *
            diffusive_operator;
        Kch4molh2.noalias() += 0.0 * laplace_operator;
        Kch4molch4.noalias() +=
            (porosity * D_G * (1 - Sw) * rho_mol_nonwet +
             porosity * D_L * Sw * rho_mol_wet * pg_int_pt / Hen_L_c) *
            diffusive_operator;
        Kch4molco2.noalias() += 0.0 * laplace_operator;
        Kch4pc.noalias() +=
            (-lambda_L * rho_mol_wet * X_L_c_gp) * laplace_operator;
        // co2
        Kco2pg.noalias() +=
            (lambda_G * rho_mol_nonwet * X3_int_pt +
             lambda_L * rho_mol_wet * X_L_co2_gp) *
                laplace_operator +
            (porosity * D_L * Sw * rho_mol_wet * X3_int_pt / Hen_L_co2) *
            diffusive_operator;
        Kco2molh2.noalias() += 0.0 * laplace_operator;
        Kco2molch4.noalias() += 0.0 * laplace_operator;
        Kco2molco2.noalias() +=
            (porosity * D_G * (1 - Sw) * rho_mol_nonwet +
             porosity * D_L * Sw * rho_mol_wet * pg_int_pt / Hen_L_co2) *
            diffusive_operator;
        Kco2pc.noalias() +=
            (-lambda_L * rho_mol_wet * X_L_co2_gp) * laplace_operator;
        // air
        Kairpg.noalias() +=
            (lambda_G * rho_mol_nonwet * x_nonwet_air +
             lambda_L * rho_mol_wet * K_G_air * x_nonwet_air) *
                laplace_operator +
            (porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_pg +
             porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_pg * K_G_air +
             porosity * D_L * Sw * rho_mol_wet * x_nonwet_air / Hen_L_air) *
            diffusive_operator;
        Kairmolh2.noalias() +=
            (porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x1 +
             porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_x1 *
                 K_G_air) *
            diffusive_operator;
        Kairmolch4.noalias() +=
            (porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x2 +
             porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_x2 *
                 K_G_air) *
            diffusive_operator;
        Kairmolco2.noalias() +=
            (porosity * D_G * S_G_gp * rho_mol_nonwet * d_x_nonwet_air_d_x3 +
             porosity * D_L * (1 - S_G_gp) * rho_mol_wet * d_x_nonwet_air_d_x3 *
                 K_G_air) *
            diffusive_operator;
        Kairpc.noalias() +=
            (-lambda_L * rho_mol_wet * K_G_air * x_nonwet_air) * laplace_operator;
        // h2o
        Kh2opg.noalias() +=
            (lambda_G * rho_mol_nonwet * x_nonwet_h2o +
             lambda_L * rho_mol_water) *
                laplace_operator -
            ((porosity * D_L * Sw * rho_mol_wet *
              (X1_int_pt / Hen_L_h + X2_int_pt / Hen_L_c +
               X3_int_pt / Hen_L_co2)) +
             porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_pg * K_G_air +
             porosity * D_L * Sw * rho_mol_wet * x_nonwet_air / Hen_L_air +
             porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_pg) *
            diffusive_operator;

        Kh2omolh2.noalias() +=
            (-porosity * D_G * (1 - Sw) * rho_mol_nonwet -
             porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x1 -
             porosity * D_L * Sw * rho_mol_wet * pg_int_pt / Hen_L_h -
             porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_x1 *
                 K_G_air) *
            diffusive_operator;
        Kh2omolch4.noalias() +=
            (-porosity * D_G * (1 - Sw) * rho_mol_nonwet -
             porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x2 -
             porosity * D_L * Sw * rho_mol_wet * pg_int_pt / Hen_L_c -
             porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_x2 *
                 K_G_air) *
            diffusive_operator;
        Kh2omolco2.noalias() +=
            (-porosity * D_G * (1 - Sw) * rho_mol_nonwet -
             porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x3 -
             porosity * D_L * Sw * rho_mol_wet * pg_int_pt /
                 Hen_L_co2 -
             porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_x3 *
                 K_G_air) *
            diffusive_operator;
        Kh2opc.noalias() +=
            (-lambda_L * rho_mol_wet * K_G_w * x_nonwet_h2o) * laplace_operator;

        if (_process_data._has_gravity)
        {
            auto const& b = _process_data._specific_body_force;
            NodalVectorType gravity_operator = sm.dNdx.transpose() *
                permeability * b *
                integration_factor;
            Bh2.noalias() +=
                (-lambda_G * rho_mol_nonwet * X1_int_pt * rho_mass_G_gp -
                 lambda_L * rho_mol_wet * X1_int_pt * pg_int_pt * rho_mass_wet /
                     Hen_L_h) *
                gravity_operator;
            Bch4.noalias() +=
                (-lambda_G * rho_mol_nonwet * X2_int_pt * rho_mass_G_gp -
                 lambda_L * rho_mol_wet * X2_int_pt * pg_int_pt * rho_mass_wet /
                     Hen_L_c) *
                gravity_operator;
            Bco2.noalias() +=
                (-lambda_G * rho_mol_nonwet * X3_int_pt * rho_mass_G_gp -
                 lambda_L * rho_mol_wet * X3_int_pt * pg_int_pt * rho_mass_wet /
                     Hen_L_co2) *
                gravity_operator;
            Bair.noalias() +=
                (-lambda_G * rho_mol_nonwet * x_nonwet_air * rho_mass_G_gp -
                 lambda_L * rho_mol_wet * x_nonwet_air * pg_int_pt *
                     rho_mass_wet / Hen_L_air) *
                gravity_operator;
            Bh2o.noalias() +=
                (-lambda_G * rho_mol_nonwet * x_nonwet_h2o * rho_mass_G_gp -
                 lambda_L * rho_mol_wet * x_nonwet_h2o * pg_int_pt *
                     rho_mass_wet / P_sat_gp) *
                gravity_operator;

        }  // end of has gravity
        // apply source and sink term
        if (_process_data._material->getMaterialID(pos) == 1 && dt>0)
        {
            Bh2.noalias() += sm.N.transpose() * Q_steel * integration_factor;
            Bh2o.noalias() +=
                sm.N.transpose() * (-Q_steel) * integration_factor;

            const double Q_organic_slow_co2 =
                Q_organic_slow_ch4_ini * para_slow;

            const double Q_organic_fast_co2 =
                Q_organic_fast_co2_ini * para_fast;
            Bch4.noalias() +=
                sm.N.transpose() * (Q_organic_slow_co2 * 31 / 19) *
                integration_factor;  // ch4 organic degradation slow
            Bco2.noalias() +=
                sm.N.transpose() * Q_organic_slow_co2 *
                integration_factor;  // co2 organic degradation slow
            Bh2o.noalias() += sm.N.transpose() * 2 * (-Q_organic_slow_co2) *
                              integration_factor;  // h2o consume slow
            Bch4.noalias() +=
                sm.N.transpose() * Q_organic_fast_co2 *
                integration_factor;  // ch4 organic degradation fast
            Bco2.noalias() +=
                sm.N.transpose() * Q_organic_fast_co2 *
                integration_factor;  // co2 organic degradation fast
            Bh2o.noalias() += sm.N.transpose() * (-Q_organic_fast_co2 / 3) *
                              integration_factor;  // h2o consume fast
        }
        else if (_process_data._material->getMaterialID(pos) == 0 && dt>0)
        {
            //co2 consumption
            // curve_slow->eval(time.getTime(), Q_organic_slow_co2);
            const double Q_organic_slow_co2 =
                Q_organic_slow_co2_ini * para_slow;
            // curve_fast->eval(time.getTime(), Q_organic_fast_co2);
            const double Q_organic_fast_co2 =
                Q_organic_fast_co2_ini * para_fast;
            const double rho_co2_ele =
                porosity * ((1 - Sw) * rho_mol_nonwet * X3_int_pt +
                            Sw * rho_mol_wet * X_L_co2_gp);

            Bco2.noalias() += sm.N.transpose() * (-rho_co2_ele) *
                              integration_factor;  // carbonation consume co2
            Bh2o.noalias() +=
                sm.N.transpose() * (rho_co2_ele) * integration_factor;

            Bh2o.noalias() -= sm.N.transpose() * 2.57635 *
                              integration_factor;  // concrete degradation
        }
    }  // end of GP
    if (_process_data._has_mass_lumping)
    {
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
