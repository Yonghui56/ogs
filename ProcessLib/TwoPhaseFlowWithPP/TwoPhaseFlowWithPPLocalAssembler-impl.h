/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowLocalAssembler.cpp
 *
 * Created on October 19, 2016, 2:28 PM
 */

#ifndef OGS_TWOPHASEFLOWWITHPPLOCALASSEMBLER_IMPL_H
#define OGS_TWOPHASEFLOWWITHPPLOCALASSEMBLER_IMPL_H

#include "TwoPhaseFlowWithPPLocalAssembler.h"

#include "NumLib/Function/Interpolation.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

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
    assert(local_matrix_size == ShapeFunction::NPOINTS);

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

	Eigen::MatrixXd mass_mat_coeff = Eigen::MatrixXd::Zero(N_size, N_size);
	Eigen::MatrixXd K_mat_coeff = Eigen::MatrixXd::Zero(N_size, N_size);
	Eigen::VectorXd H_vec_coeff = Eigen::VectorXd::Zero(N_size);
	Eigen::MatrixXd localMass_tmp = Eigen::MatrixXd::Zero(n_nodes, n_nodes);
	Eigen::MatrixXd localDispersion_tmp = Eigen::MatrixXd::Zero(n_nodes, n_nodes);
	Eigen::VectorXd localGravity_tmp = Eigen::VectorXd::Zero(n_nodes);
	MathLib::PiecewiseLinearInterpolation const&  interpolated_Pc = *_process_data.curves.at("curveA");
	MathLib::PiecewiseLinearInterpolation const&  interpolated_Kr_L = *_process_data.curves.at("curveB");
    // TODO: The following two variables should be calculated inside the
    //       the integration loop for non-constant porosity and storage models.
    double porosity_variable = 0.;
    double storage_variable = 0.;
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);

        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, p);
        // TODO : compute _temperature from the heat transport pcs

        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();

        // Assemble mass matrix, M
        local_M.noalias() +=
            _material_properties.getMassCoefficient(
                t, pos, porosity_variable, storage_variable, p, _temperature) *
            sm.N.transpose() * sm.N * integration_factor;

        // Compute gas density:
        const double rho_g =
            _material_properties.getLiquidDensity(p, _temperature) *
            _gravitational_acceleration;
        // Compute viscosity:
        const double mu = _material_properties.getViscosity(p, _temperature);

        // Assemble Laplacian, K, and RHS by the gravitational term
        if (perm.size() == 1)
        {
            //  Use scalar number for isotropic permeability
            //  to save the computation time.
            const double K = perm(0, 0) / mu;
            const double fac = K * integration_factor;
            local_K.noalias() += fac * sm.dNdx.transpose() * sm.dNdx;
            if (_gravitational_axis_id >= 0)
            {
                local_b.noalias() -=
                    fac * sm.dNdx.transpose().col(_gravitational_axis_id) *
                    rho_g;
            }
        }
        else
        {
            const double fac = integration_factor / mu;
            local_K.noalias() += fac * sm.dNdx.transpose() * perm * sm.dNdx;
            if (_gravitational_axis_id >= 0)
            {
                local_b.noalias() -= fac * rho_g * sm.dNdx.transpose() *
                                     perm.col(_gravitational_axis_id);
            }
        }
    }
}

}  // end of namespace
}  // end of namespace

#endif
