/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   MonolithicHTFEM.h
 *  Created on October 11, 2017, 2:33 PM
 */

#pragma once

#include <Eigen/Dense>
#include <vector>

#include "HTMaterialProperties.h"
#include "MaterialLib/MPL/Medium.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "HTFEM.h"

namespace ProcessLib
{
namespace HT
{
const unsigned NUM_NODAL_DOF = 2;

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class MonolithicHTFEM
    : public HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalMatrixType = typename ShapeMatricesType::template MatrixType<
        NUM_NODAL_DOF * ShapeFunction::NPOINTS,
        NUM_NODAL_DOF * ShapeFunction::NPOINTS>;
    using LocalVectorType =
        typename ShapeMatricesType::template VectorType<NUM_NODAL_DOF *
                                                        ShapeFunction::NPOINTS>;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    MonolithicHTFEM(MeshLib::Element const& element,
                    std::size_t const local_matrix_size,
                    bool is_axially_symmetric,
                    unsigned const integration_order,
                    HTMaterialProperties const& material_properties)
        : HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>(
              element, local_matrix_size, is_axially_symmetric,
              integration_order, material_properties, NUM_NODAL_DOF)
    {
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

        auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<LocalVectorType>(
            local_b_data, local_matrix_size);

        auto const num_nodes = ShapeFunction::NPOINTS;

        auto Ktt = local_K.template block<num_nodes, num_nodes>(0, 0);
        auto Ktp = local_K.template block<num_nodes, num_nodes>(0, num_nodes);
        auto Mtt = local_M.template block<num_nodes, num_nodes>(0, 0);
        auto Kpp =
            local_K.template block<num_nodes, num_nodes>(num_nodes, num_nodes);
        auto Mpp =
            local_M.template block<num_nodes, num_nodes>(num_nodes, num_nodes);
        auto Bp = local_b.template block<num_nodes, 1>(num_nodes, 0);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(this->_element.getID());

        auto p_nodal_values =
            Eigen::Map<const NodalVectorType>(&local_x[num_nodes], num_nodes);

        auto const& process_data = this->_material_properties;
        auto const& medium =
            *process_data.media_map->getMedium(this->_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        auto const& solid_phase = medium.phase("Solid");

        auto const& b = process_data.specific_body_force;
        //Eigen::Vector2d b = Eigen::Vector2d::Zero();
        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

        MaterialPropertyLib::VariableArray vars;

        unsigned const n_integration_points =
            this->_integration_method.getNumberOfPoints();

        for (unsigned ip(0); ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);

            // \todo the argument to getValue() has to be changed for non
            // constant storage model
            auto const specific_storage =
                solid_phase.property(MaterialPropertyLib::PropertyType::storage)
                    .template value<double>(vars);

            auto const& ip_data = this->_ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;

            double T_int_pt = 0.0;
            double p_int_pt = 0.0;
            // Order matters: First T, then P!
            NumLib::shapeFunctionInterpolate(local_x, N, T_int_pt, p_int_pt);

            auto const porosity =
                process_data.porous_media_properties.getPorosity(t, pos)
                .getValue(t, pos, 0.0, p_int_pt);

            auto const& intrinsic_permeability =
                process_data.porous_media_properties.getIntrinsicPermeability(
                    t, pos).getValue(t, pos, 0.0, 0.0);
            const double K = std::pow(10, intrinsic_permeability(0, 0));
            /*auto const porosity =
                solid_phase
                    .property(MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(vars);

            auto const intrinsic_permeability =
                intrinsicPermeability<GlobalDim>(
                    solid_phase
                        .property(
                            MaterialPropertyLib::PropertyType::permeability)
                        .value(vars));*/

            vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
                T_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;

            auto const specific_heat_capacity_fluid =
                liquid_phase
                    .property(MaterialPropertyLib::specific_heat_capacity)
                    .template value<double>(vars);

            // Use the fluid density model to compute the density
            auto const fluid_density =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars);

            // Use the viscosity model to compute the viscosity
            auto const viscosity = liquid_phase
                    .property(MaterialPropertyLib::PropertyType::viscosity)
                    .template value<double>(vars);
            double const K_over_mu = K / viscosity;

            GlobalDimVectorType const velocity = 
                process_data.has_gravity
                    ? GlobalDimVectorType(-K_over_mu * (dNdx * p_nodal_values -
                                                        fluid_density * b))
                    : GlobalDimVectorType(-K_over_mu * dNdx * p_nodal_values);
            GlobalDimVectorType const velocity_supg =
                process_data.has_gravity
                ? GlobalDimVectorType(-K_over_mu * (dNdx * p_nodal_values -
                    fluid_density * b))
                : GlobalDimVectorType(-K_over_mu * dNdx * p_nodal_values);
            //velocity norm
            double u_norm = GlobalDim > 2 ? (std::pow(velocity_supg(0), 2) +
                std::pow(velocity_supg(1), 2) +
                std::pow(velocity_supg(2), 2))
                : (std::pow(velocity_supg(0), 2) +
                    std::pow(velocity_supg(1), 2));
            u_norm = std::sqrt(u_norm);
            double& tau2 = _ip_data[ip].tauSUPG;

            auto min_length = sqrt(this->_element.getContent());
            auto vdN = dNdx.transpose() * velocity_supg;//4*1 pow(0 / (0.5*dt),2.0)+
            double alpha = 0.5 * u_norm * min_length / (2 * K * 0.1); // this is the Peclet number
            const double xi_tilde = cosh_relation(alpha);
            tau2 = xi_tilde == 0
                ? 0
                : 0;
            /*: 0.5 * min_length / u_norm * xi_tilde;*/
            
            GlobalDimMatrixType const thermal_conductivity_dispersivity =
                this->getThermalConductivityDispersivity(
                    vars, porosity, fluid_density, specific_heat_capacity_fluid,
                    velocity, I);
            //tau
            /*tau2 = 1 / sqrt(pow(2.0 * u_norm / min_length, 2.0)
                + 9 * pow(4 * thermal_conductivity_dispersivity(0,0) / min_length / min_length, 2.0));*/
            // matrix assembly
            Ktt.noalias() +=
                (dNdx.transpose() * thermal_conductivity_dispersivity * dNdx +
                 N.transpose() * velocity.transpose() * dNdx * fluid_density *
                     specific_heat_capacity_fluid) *
                w;
            /*Ktp.noalias() +=
                dNdx.transpose()*fluid_density*specific_heat_capacity_fluid*K_over_mu*(T_int_pt-273.15)
                *dNdx*w;*/
            Kpp.noalias() += w * dNdx.transpose() * K_over_mu * dNdx;
            Mtt.noalias() +=
                w *
                this->getHeatEnergyCoefficient(vars, porosity, fluid_density,
                                               specific_heat_capacity_fluid) *
                N.transpose() * N;
            Mpp.noalias() += w * N.transpose() * specific_storage * N;
            //SUPG
            Ktt.noalias() += dNdx.transpose() *fluid_density*specific_heat_capacity_fluid
                * velocity
                * vdN.transpose() *tau2* w;
            Mtt.noalias() += tau2 * this->getHeatEnergyCoefficient(vars, porosity, fluid_density,
                specific_heat_capacity_fluid)*vdN*N*w;
            if (process_data.has_gravity)
            {
                Bp += w * fluid_density * dNdx.transpose() * K_over_mu * b;
            }
            /* with Oberbeck-Boussing assumption density difference only exists
             * in buoyancy effects */
        }
        if (true)
        {
            for (unsigned row = 0; row < Mtt.cols(); row++)
            {
                for (unsigned column = 0; column < Mtt.cols(); column++)
                {
                    if (row != column)
                    {
                        Mtt(row, row) += Mtt(row, column);
                        Mtt(row, column) = 0.0;
                        Mpp(row, row) += Mpp(row, column);
                        Mpp(row, column) = 0.0;
                    }
                }
            }
        }  // end of mass-lumping
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& cache) const override
    {
        auto const indices =
            NumLib::getIndices(this->_element.getID(), dof_table);
        assert(!indices.empty());
        auto local_x = current_solution.get(indices);

        std::vector<double> local_p(
            std::make_move_iterator(local_x.begin() + local_x.size() / 2),
            std::make_move_iterator(local_x.end()));
        // only T is kept in local_x
        local_x.erase(local_x.begin() + local_x.size() / 2, local_x.end());

        return this->getIntPtDarcyVelocityLocal(t, local_p, local_x, cache);
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
};

}  // namespace HT
}  // namespace ProcessLib
