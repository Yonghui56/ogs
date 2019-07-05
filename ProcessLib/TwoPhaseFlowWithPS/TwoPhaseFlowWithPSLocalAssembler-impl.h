/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
* common nomenclature
* --------------primary variable----------------------
* pn_int_pt    pressure for nonwetting phase at each integration point
* pc_int_pt    capillary pressure at each integration point
* --------------secondary variable--------------------
* temperature              capillary pressure
* Sw wetting               phase saturation
* dSw_dpc                  derivative of wetting phase saturation with respect
* to capillary pressure
* rho_nonwet               density of nonwetting phase
* drhononwet_dpn           derivative of nonwetting phase density with respect
*to nonwetting phase pressure
* rho_wet                  density of wetting phase
* k_rel_nonwet             relative permeability of nonwetting phase
* mu_nonwet                viscosity of nonwetting phase
* lambda_nonwet            mobility of nonwetting phase
* k_rel_wet                relative permeability of wetting phase
* mu_wet                   viscosity of wetting phase
* lambda_wet               mobility of wetting phase
*/
#pragma once

#include "TwoPhaseFlowWithPSLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseFlowWithPSProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPS
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseFlowWithPSLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, std::vector<double> const& local_x,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();

    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
    auto const num_nodes = ShapeFunction::NPOINTS;
    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    auto Mpp =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mps = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);

    auto Msp = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mss = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    NodalMatrixType laplace_operator =
        NodalMatrixType::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kpp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
     
    auto Ksp = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Kss = local_K.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    auto Bp = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);

    auto Bs =
        local_b.template segment<cap_pressure_size>(cap_pressure_matrix_index);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto p_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[0], num_nodes);
    auto s_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[num_nodes], num_nodes);
    auto const& b = _process_data.specific_body_force;

    std::vector<double> testdata = {0, 0, 0, 0};
    double min_length =
        CalcEffectiveElementLength(GlobalDim, _element, num_nodes);
    const double& dt = _process_data.dt;
    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id =
        _process_data.material->getMaterialID(pos.getElementID().get());

    const Eigen::MatrixXd& perm = _process_data.material->getPermeability(
        material_id, t, pos, _element.getDimension());
    assert(perm.rows() == _element.getDimension() || perm.rows() == 1);
    GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
        _element.getDimension(), _element.getDimension());
    if (perm.rows() == _element.getDimension())
    {
        permeability = perm;
    }
    else if (perm.rows() == 1)
    {
        permeability.diagonal().setConstant(perm(0, 0));
    }

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& ip_data = this->_ip_data[ip];
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        auto const& w = ip_data.integration_weight;
        auto v_dN = testdata;

        double sn_int_pt = 0.;
        double pn_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, pn_int_pt,
                                         sn_int_pt);
        bool sign = false;
        if (sn_int_pt > 1)
        {
            sn_int_pt = 0.99999;
            sign = true;
        }




        const double temperature = _process_data.temperature(t, pos)[0];
        double const rho_nonwet =
            _process_data.material->getGasDensity(pn_int_pt, temperature);
        double const rho_wet = _process_data.material->getLiquidDensity(
            _pressure_wet[ip], temperature);

        double const porosity = _process_data.material->getPorosity(
            material_id, t, pos, pn_int_pt, temperature, 0);

        // Assemble M matrix
        // nonwetting
        double const drhononwet_dpn =
            _process_data.material->getGasDensityDerivative(pn_int_pt,
                                                            temperature);

        double const compressibility = 4.4e-10;// pre-assume
        double const sn_r = 0.3;
        double const sw_r = 0.01;
        double Sn_eff = (sn_int_pt - sn_r) / (1 - sn_r - sw_r);
        // nonwet
        double const k_rel_nonwet = Sn_eff * Sn_eff;
        double const mu_nonwet =
            _process_data.material->getGasViscosity(pn_int_pt, temperature);
        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;

        // wet
        double const k_rel_wet = (1 - Sn_eff) * (1 - Sn_eff);
        double const mu_wet = _process_data.material->getLiquidViscosity(
            _pressure_wet[ip], temperature);
        double const lambda_wet = k_rel_wet / mu_wet;

        double const lambda_total = lambda_nonwet + lambda_wet;

        //the fractional flow function of phase water/nonwetting
        double f_l = lambda_wet / lambda_total;
        double f_n = 1 - f_l;

        //calculate the derivatives of


        double  dLambda_n_dS = 2 * Sn_eff / mu_nonwet / (1 - sn_r - sw_r);
        double  dLambda_w_dS =
            -2 * (1 - Sn_eff) / mu_wet / (1 - sn_r - sw_r);
        //fractional flow function vesus nonwet phase sn
        double  df_n_dsn = dLambda_n_dS / (lambda_nonwet + lambda_wet) -
                       lambda_nonwet/ std::pow((lambda_nonwet + lambda_wet),2) * (dLambda_n_dS + dLambda_w_dS);
        if (sign)
        {
            dLambda_w_dS = 0;
            dLambda_n_dS = 0;
            df_n_dsn = 0;
        }
        laplace_operator.noalias() = _ip_data[ip].dNdx.transpose() *
                                     permeability * _ip_data[ip].dNdx *
                                     _ip_data[ip].integration_weight;
        //calculate the phase velocity
        GlobalDimVectorType const velocity_total =
            _process_data.has_gravity
                ? GlobalDimVectorType(
                      -permeability * lambda_total *
                                      (dNdx * p_nodal_values - (rho_nonwet*f_n+rho_wet*f_l)*b))
                : GlobalDimVectorType(-permeability * lambda_total * dNdx *
                                      p_nodal_values);
        GlobalDimVectorType const velocity_nonwet =
            _process_data.has_gravity
                ? GlobalDimVectorType(-permeability * lambda_nonwet *
                                      (dNdx * p_nodal_values -
                                       (rho_nonwet * f_n + rho_wet * f_l) * b))
                : GlobalDimVectorType(-permeability * lambda_nonwet * dNdx *
                                      p_nodal_values);
        Mpp.noalias() +=
            porosity * compressibility * _ip_data[ip].massOperator;
        Kpp.noalias() += lambda_total * laplace_operator;

        Mss.noalias() +=
            porosity *_ip_data[ip].massOperator;
        Kss.noalias() += 
        _ip_data[ip].N.transpose() * velocity_total.transpose() *
            df_n_dsn * _ip_data[ip].dNdx *
                         _ip_data[ip].integration_weight;
        Ksp.noalias() += lambda_nonwet * laplace_operator;
        // calculate the weighting function
        double u_norm = GlobalDim > 2 ? (std::pow(velocity_nonwet(0), 2) +
                                         std::pow(velocity_nonwet(1), 2) +
                                         std::pow(velocity_nonwet(2), 2))
                                      : (std::pow(velocity_nonwet(0), 2) +
                                         std::pow(velocity_nonwet(1), 2));
        u_norm = std::sqrt(u_norm);
        for (int i = 0; i < num_nodes; i++)
            for (int k = 0; k < GlobalDim; k++)
                v_dN[i] += dNdx(k * num_nodes + i) * velocity_nonwet(k);
        auto tau2 = dt == 0
                       ? 0
                       : std::pow(1 / (0.5 * dt) + 2.0 * u_norm / min_length +
                                      4 * 1e-15 / pow(min_length, 2.0),
                                  -1);

        auto tau = std::pow(1 / (dt * dt) + std::pow(2.0 * u_norm / min_length,2),
                                  -1/2);
        /*CalcSUPGCoefficient(u_norm, ip, min_length,
            diffusion_coeff_component_salt,
            dt);*/
        for (int i = 0; i < num_nodes; i++)
            for (int j = 0; j < num_nodes; j++)
                Kss(i, j) += w * tau * v_dN[i] * v_dN[j] * df_n_dsn
                             ;
        for (int i = 0; i < num_nodes; i++)
            for (int j = 0; j < num_nodes; j++)
                Mss(i, j) +=
                    w * tau * (porosity) * v_dN[i] * N(j);

        if (_process_data.has_gravity)
        {
            NodalVectorType gravity_operator = _ip_data[ip].dNdx.transpose() *
                                               permeability * b *
                                               _ip_data[ip].integration_weight;
            Bp.noalias() +=
                rho_nonwet * rho_nonwet * lambda_nonwet * gravity_operator;
            Bs.noalias() += rho_wet * rho_wet * lambda_wet * gravity_operator;
        }  // end of has gravity
    }
    if (true)
    {
        for (unsigned row = 0; row < Mpp.cols(); row++)
        {
            for (unsigned column = 0; column < Mpp.cols(); column++)
            {
                if (row != column)
                {
                    Mpp(row, row) += Mpp(row, column);
                    Mpp(row, column) = 0.0;
                    Mss(row, row) += Mss(row, column);
                    Mss(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // namespace TwoPhaseFlowWithPS
}  // namespace ProcessLib
