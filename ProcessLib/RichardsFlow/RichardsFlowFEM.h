/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "RichardsFlowProcessData.h"

namespace ProcessLib
{
namespace RichardsFlow
{
template <typename NodalRowVectorType, typename GlobalDimNodalMatrixType,
          typename NodalMatrixType>
struct IntegrationPointData final
{
    IntegrationPointData(NodalRowVectorType const& N_,
                         GlobalDimNodalMatrixType const& dNdx_,
                         double const& integration_weight_,
                         NodalMatrixType const mass_operator_)
        : N(N_),
          dNdx(dNdx_),
          integration_weight(integration_weight_),
          mass_operator(mass_operator_)
    {
    }

    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;
    NodalMatrixType const mass_operator;
    double saturation;
    double tauSUPG;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
const unsigned NUM_NODAL_DOF = 1;

class RichardsFlowLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtTauSUPG(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public RichardsFlowLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       RichardsFlowProcessData const& process_data)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _saturation(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _tauSUPG(
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
        (void)local_matrix_size;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                element, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& sm = shape_matrices[ip];
            const double integration_factor = sm.integralMeasure * sm.detJ;
            _ip_data.emplace_back(
                sm.N, sm.dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    integration_factor,
                sm.N.transpose() * sm.N * integration_factor *
                    _integration_method.getWeightedPoint(ip).getWeight());
        }
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<NodalVectorType>(
            local_b_data, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());
        const int material_id =
            _process_data.material->getMaterialID(_element.getID());
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
        double const pressure_size = ShapeFunction::NPOINTS;
        auto p_L =
            Eigen::Map<const NodalVectorType>(local_x.data(),
                pressure_size);
        double const& dt = _process_data.dt;
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;
            double p_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, p_int_pt);
            double p_cap_ip;
            NumLib::shapeFunctionInterpolate(-p_L, N, p_cap_ip);//PC on each itegrate point
            const double& temperature = _process_data.temperature(t, pos)[0];
            auto const porosity = _process_data.material->getPorosity(
                material_id, t, pos, p_int_pt, temperature, 0);

            double const pc_int_pt = -p_int_pt;

            double const Sw = _process_data.material->getSaturation(
                material_id, t, pos, p_int_pt, temperature, pc_int_pt);
            double& S_L = _ip_data[ip].saturation;
            S_L = _process_data.material->getSaturation(
                material_id, t, pos, -p_cap_ip, temperature, p_cap_ip);

            double const dSw_dpc =
                _process_data.material->getSaturationDerivative(
                    material_id, t, pos, p_int_pt, temperature, Sw);

            auto const& body_force = _process_data.specific_body_force;
            auto const rho_LR = _process_data.material->getFluidDensity(
                -p_cap_ip, temperature);
            // \TODO Extend to pressure dependent density.
            double const drhow_dp(0.0);
            auto const storage = _process_data.material->getStorage(
                material_id, t, pos, p_int_pt, temperature, 0);
            double const mass_mat_coeff =
                storage * Sw + porosity * Sw * drhow_dp - porosity * dSw_dpc;


            local_M.noalias() += mass_mat_coeff * _ip_data[ip].mass_operator;

            double const k_rel =
                _process_data.material->getRelativePermeability(
                    t, pos, p_int_pt, temperature, Sw);
            auto const mu = _process_data.material->getFluidViscosity(
                p_int_pt, temperature);
            local_K.noalias() += _ip_data[ip].dNdx.transpose() * permeability *
                _ip_data[ip].dNdx *
                _ip_data[ip].integration_weight * (k_rel / mu);

            double const dk_rel_dS_l =
                _process_data.material->getRelativePermeabilityDerivative(
                    t, pos
                    , -p_cap_ip, temperature, Sw);
            typename ShapeMatricesType::GlobalDimVectorType const
                grad_p_cap = -dNdx * p_L;

            //velocity
            GlobalDimVectorType const velocity =
                _process_data.has_gravity
                ? GlobalDimVectorType(-permeability * (1 )*
                (dNdx * p_L -
                    rho_LR * body_force))
                : GlobalDimVectorType(-permeability * (1 ) * dNdx *
                    p_L);
            //velocity norm
            double u_norm = GlobalDim > 2 ? (std::pow(velocity(0), 2) +
                std::pow(velocity(1), 2) +
                std::pow(velocity(2), 2))
                : (std::pow(velocity(0), 2) +
                    std::pow(velocity(1), 2));
            u_norm = std::sqrt(u_norm);
            //double const min_length = 0.02;
            //tau for SUPG
            double& tau2 = _ip_data[ip].tauSUPG;
            /*tau2 = dt == 0
                ? 0
                : std::pow(1 / (0.5 * dt) + 2.0 * u_norm / min_length +
                    4 * 1e-15 / pow(min_length, 2.0),
                    -1);*/
            auto min_length = (this->_element.getContent());
            /*tau2 = dt == 0
                ? 0
                : 1/sqrt( pow(2.0 * u_norm / min_length, 2.0));*/
            auto vdN = dNdx.transpose() * velocity;//4*1 pow(0 / (0.5*dt),2.0)+
            /*auto norm_b = (2.0 * u_norm) / min_length;
            double alpha = 0.5 * u_norm * min_length / (2*permeability(0,0) * 0.1); // this is the Peclet number
            const double xi_tilde = cosh_relation(alpha);
            tau2 = xi_tilde / norm_b;*/
            //tau2 = 0;
            double alpha = 0.5 * u_norm * min_length / (2 * permeability(0, 0) *100); // this is the Peclet number
            const double xi_tilde = cosh_relation(alpha);
            tau2 = 0.5 * min_length / u_norm * xi_tilde;
            local_K.noalias() += dNdx.transpose() * (permeability/mu)
                * (-1)*(-grad_p_cap-body_force*rho_LR) *
                dk_rel_dS_l * dSw_dpc *vdN.transpose() *tau2* w;
            local_M.noalias() += tau2 * mass_mat_coeff*vdN*N*w;
            if (_process_data.has_gravity)
            {
                assert(body_force.size() == GlobalDim);
                NodalVectorType gravity_operator =
                    _ip_data[ip].dNdx.transpose() * permeability * body_force *
                    _ip_data[ip].integration_weight;
                local_b.noalias() += (k_rel / mu) * rho_LR * gravity_operator;
                //local_b.noalias() +=  dk_rel_dS_l * dSw_dpc* rho_w*(permeability/mu)* body_force*vdN.transpose()*tau2*w;
            }
        }
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

    void assembleWithJacobian(double const t, std::vector<double> const& local_x,
        std::vector<double> const& local_xdot,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_rhs_data,
        std::vector<double>& local_Jac_data)
    {
        double const pressure_size = ShapeFunction::NPOINTS;
        double const pressure_index = 0;
        assert(local_x.size() == pressure_size);
        auto p_L =
            Eigen::Map<const NodalVectorType>(local_x.data(),
                pressure_size);
        
        auto p_L_dot =
            Eigen::Map<const NodalVectorType>(local_xdot.data(),
                pressure_size);

        auto local_Jac = MathLib::createZeroedMatrix<
            NodalMatrixType>(
                local_Jac_data, pressure_size,
                pressure_size);

        auto local_rhs = MathLib::createZeroedVector<
            NodalVectorType >(
                local_rhs_data, pressure_size);

        typename ShapeMatricesType::NodalMatrixType laplace_p =
            ShapeMatricesType::NodalMatrixType::Zero(pressure_size,
                pressure_size);

        typename ShapeMatricesType::NodalMatrixType storage_p =
            ShapeMatricesType::NodalMatrixType::Zero(pressure_size,
                pressure_size);

        double const& dt = _process_data.dt;
        auto const material_id =
            _process_data.material->getMaterialID(_element.getID());;

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        const Eigen::MatrixXd& perm = _process_data.material->getPermeability(
            material_id, t, x_position, _element.getDimension());
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

        //loop for each gauss point
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;
            auto const x_coord =
                interpolateXCoordinate<ShapeFunction,
                ShapeMatricesType>(_element,
                    N);
            double p_cap_ip;
            NumLib::shapeFunctionInterpolate(-p_L, N, p_cap_ip);//PC on each itegrate point

            double p_cap_dot_ip;
            NumLib::shapeFunctionInterpolate(-p_L_dot, N, p_cap_dot_ip);
            double temperature = 293.15;
            auto const porosity = _process_data.material->getPorosity(
                material_id, t, x_position, -p_cap_ip, temperature, 0);

            auto const rho_LR = _process_data.material->getFluidDensity(
                -p_cap_ip, temperature);

            double& S_L = _ip_data[ip].saturation;
            S_L = _process_data.material->getSaturation(
                material_id, t, x_position, -p_cap_ip, temperature, p_cap_ip);

            double const dS_L_dp_cap =
                _process_data.material->getSaturationDerivative(
                    material_id, t, x_position, -p_cap_ip, temperature, S_L);

            double const d2S_L_dp_cap_2 =
                _process_data.material->getSaturationDerivative2(
                    material_id, t, x_position, -p_cap_ip, temperature, S_L);

            double const k_rel =
                _process_data.material->getRelativePermeability(
                    t, x_position, -p_cap_ip, temperature, S_L);
            auto const mu = _process_data.material->getFluidViscosity(
                -p_cap_ip, temperature);
            auto const& body_force = _process_data.specific_body_force;

            GlobalDimMatrixType const rho_Ki_over_mu =
                permeability *
                    (rho_LR / mu);
            //velocity
            GlobalDimVectorType const velocity =
                _process_data.has_gravity
                ? GlobalDimVectorType(-permeability * (k_rel / mu)*
                (dNdx * p_L -
                    rho_LR * body_force))
                : GlobalDimVectorType(-permeability * (k_rel / mu) * dNdx *
                    p_L);
            //velocity norm
            double u_norm = GlobalDim > 2 ? (std::pow(velocity(0), 2) +
                std::pow(velocity(1), 2) +
                std::pow(velocity(2), 2))
                : (std::pow(velocity(0), 2) +
                    std::pow(velocity(1), 2));
            u_norm = std::sqrt(u_norm);
            double const min_length = 0.02;
            //tau for SUPG
            auto tau2 = dt == 0
                ? 0
                : std::pow(1 / (0.5 * dt) + 2.0 * u_norm / min_length +
                    4 * 1e-15 / pow(min_length, 2.0),
                    -1);
            auto vdN = dNdx.transpose() * velocity;


            laplace_p.noalias() +=
                dNdx.transpose() * k_rel * rho_Ki_over_mu * dNdx * w;

            double const specific_storage =
                dS_L_dp_cap * (-1)*porosity;

            double const dspecific_storage_dp_cap =
                d2S_L_dp_cap_2 * (-1)*porosity;

            storage_p.noalias() +=
                N.transpose() * rho_LR * specific_storage * N * w;
            //mass change part
            local_Jac
                .noalias() += N.transpose() * rho_LR * p_cap_dot_ip *
                dspecific_storage_dp_cap * N * w;
            //flux change
            double const dk_rel_dS_l =
                _process_data.material->getRelativePermeabilityDerivative(
                    t, x_position, -p_cap_ip, temperature, S_L);
            typename ShapeMatricesType::GlobalDimVectorType const
                grad_p_cap = -dNdx * p_L;

            GlobalDimVectorType supg_kernel =
                dNdx.transpose() * rho_Ki_over_mu * grad_p_cap *
                dk_rel_dS_l * dS_L_dp_cap *vdN.transpose() *p_L * w;//CHECK FOR NEGATIVE OR POSITYVE
            local_rhs.noalias() += supg_kernel;

            local_Jac
                .noalias() += dNdx.transpose() * rho_Ki_over_mu * grad_p_cap *
                dk_rel_dS_l * dS_L_dp_cap * N * w;
           //gravitional part
            local_Jac
                .noalias() += dNdx.transpose() * rho_LR * rho_Ki_over_mu * body_force *
                dk_rel_dS_l * dS_L_dp_cap * N * w;
            //right hand side for the gravitional part
            local_rhs.noalias() +=
                dNdx.transpose() * rho_LR * k_rel * rho_Ki_over_mu * body_force * w;

        }
        // pressure equation, pressure part.
        local_Jac
            .noalias() += laplace_p + storage_p / dt;

        // pressure equation
        local_rhs.noalias() -=
            laplace_p * p_L + storage_p * p_L_dot;


    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

        std::vector<double> const& 
        getIntPtSaturation(const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
            std::vector<double>& cache) const
    {
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<
            Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(cache, 1,
                num_intpts);

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            cache_mat[ip] = _ip_data[ip].saturation;
        }

        return cache;
    }

        std::vector<double> const&
            getIntPtTauSUPG(const double /*t*/,
                GlobalVector const& /*current_solution*/,
                NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                std::vector<double>& cache) const
        {
            auto const num_intpts = _ip_data.size();

            cache.clear();
            auto cache_mat = MathLib::createZeroedMatrix<
                Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(cache, 1,
                    num_intpts);

            for (unsigned ip = 0; ip < num_intpts; ++ip)
            {
                cache_mat[ip] = _ip_data[ip].tauSUPG;
            }

            return cache;
        }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& cache) const override
    {
        auto const num_intpts = _shape_matrices.size();

        auto const indices = NumLib::getIndices(
            _element.getID(), dof_table);
        assert(!indices.empty());
        auto const local_x = current_solution.get(indices);

        cache.clear();
        auto cache_vec = MathLib::createZeroedMatrix<
            Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, GlobalDim, num_intpts);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());
        const int material_id =
            _process_data.material->getMaterialID(_element.getID());

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

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        auto const p_nodal_values = Eigen::Map<const NodalVectorType>(
            &local_x[0], ShapeFunction::NPOINTS);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            double p_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, p_int_pt);
            double const pc_int_pt = -p_int_pt;
            const double& temperature = _process_data.temperature(t, pos)[0];
            double const Sw = _process_data.material->getSaturation(
                material_id, t, pos, p_int_pt, temperature, pc_int_pt);
            double const k_rel =
                _process_data.material->getRelativePermeability(
                    t, pos, p_int_pt, temperature, Sw);
            auto const mu = _process_data.material->getFluidViscosity(
                p_int_pt, temperature);
            auto const K_mat_coeff = permeability * (k_rel / mu);
            cache_vec.col(ip).noalias() =
                -K_mat_coeff * _ip_data[ip].dNdx * p_nodal_values;
            if (_process_data.has_gravity)
            {
                auto const rho_w = _process_data.material->getFluidDensity(
                    p_int_pt, temperature);
                auto const& body_force = _process_data.specific_body_force;
                assert(body_force.size() == GlobalDim);
                // here it is assumed that the vector body_force is directed
                // 'downwards'
                cache_vec.col(ip).noalias() += K_mat_coeff * rho_w * body_force;
            }
        }

        return cache;
    }
    double cosh_relation(double alpha)
    {
        if (alpha >= 5.0 || alpha <= -5.0)
            return ((alpha > 0.0) ? 1.0 : -1.0) - 1.0 / alpha; // prevents overflows
        else if (alpha == 0)
            return 0.0;
        return 1.0 / std::tanh(alpha) - 1.0 / alpha;
    }
private:
    MeshLib::Element const& _element;
    RichardsFlowProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;
    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType,
                             NodalMatrixType>,
        Eigen::aligned_allocator<IntegrationPointData<
            NodalRowVectorType, GlobalDimNodalMatrixType, NodalMatrixType>>>
        _ip_data;
    std::vector<double> _saturation;
    std::vector<double> _tauSUPG;
};

}  // namespace RichardsFlow
}  // namespace ProcessLib
