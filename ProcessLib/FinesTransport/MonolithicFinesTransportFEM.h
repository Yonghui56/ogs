/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   MonolithicFinesTransportFEM.h
 *  Created on October 11, 2017, 2:33 PM
 */

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <iostream>

#include "FinesTransportMaterialProperties.h"
#include "MaterialLib/MPL/Medium.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "FinesTransportFEM.h"
#include <boost/array.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/numeric/odeint.hpp>


using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix=boost::phoenix;
typedef boost::numeric::ublas::vector<double> vector_type;
typedef runge_kutta_dopri5<double> stepper_type;
namespace ProcessLib
{
namespace FinesTransport
{
const unsigned NUM_NODAL_DOF = 3;

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class MonolithicFinesTransportFEM
    : public FinesTransportFEM<ShapeFunction, IntegrationMethod, GlobalDim>
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
    MonolithicFinesTransportFEM(MeshLib::Element const& element,
                    std::size_t const local_matrix_size,
                    bool is_axially_symmetric,
                    unsigned const integration_order,
                    FinesTransportMaterialProperties const& material_properties)
        : FinesTransportFEM<ShapeFunction, IntegrationMethod, GlobalDim>(
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

        auto Kpp = local_K.template block<num_nodes, num_nodes>(0, 0);
        auto Kss = local_K.template block<num_nodes, num_nodes>(num_nodes, num_nodes);
        auto Ksp = local_K.template block<num_nodes, num_nodes>(num_nodes, 0);

        auto Kcc = local_K.template block<num_nodes, num_nodes>(2 * num_nodes, 2 * num_nodes);
        auto Kcp = local_K.template block<num_nodes, num_nodes>(2 * num_nodes, 0);

        auto Mpp = local_M.template block<num_nodes, num_nodes>(0, 0);
        auto Mps = local_M.template block<num_nodes, num_nodes>(0, num_nodes);
        auto Mss =
            local_M.template block<num_nodes, num_nodes>(num_nodes, num_nodes);

        auto Mcc =
            local_M.template block<num_nodes, num_nodes>(2 * num_nodes, 2 * num_nodes);
        auto Mcs =
            local_M.template block<num_nodes, num_nodes>(2 * num_nodes, num_nodes);

        auto Bp = local_b.template block<num_nodes, 1>(0, 0);

        auto Bs = local_b.template block<num_nodes, 1>(num_nodes, 0);

        auto Bc = local_b.template block<num_nodes, 1>(2*num_nodes, 0);
        //auto Bf = local_b.template block<num_nodes, 1>(3*num_nodes, 0);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(this->_element.getID());

        auto p_nodal_values =
            Eigen::Map<const NodalVectorType>(&local_x[0], num_nodes);
        auto s_nodal_values =
            Eigen::Map<const NodalVectorType>(&local_x[num_nodes], num_nodes);
        auto c_nodal_values =
            Eigen::Map<const NodalVectorType>(&local_x[2*num_nodes], num_nodes);


        auto const& process_data = this->_material_properties;
        auto const& medium =
            *process_data.media_map->getMedium(this->_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        auto const& solid_phase = medium.phase("Solid");
        auto const& nonwetting_phase = medium.phase("NonAqueousLiquid");

        auto const& b = process_data.specific_body_force;

        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

        MaterialPropertyLib::VariableArray vars;

        unsigned const n_integration_points =
            this->_integration_method.getNumberOfPoints();
        const double& dt = process_data.dt;

        /*vector_type x(2, 1.0);

        size_t num_of_steps = integrate_const(
            make_dense_output<rosenbrock4<double>>(1.0e-6, 1.0e-6),
            make_pair(stiff_system(), stiff_system_jacobi()), x, 0.0, 50.0,
            0.01,
            std::cout << phoenix::arg_names::arg2 << " "
                 << phoenix::arg_names::arg1[0] << "\n" 
            );*/
        double x_pt(0.0);
        double x_d(0.0);
        //this->u_norm = 0.0037;
        //this->concentration = 1.4488e-10;0
        size_t num_of_steps_pt =
            integrate_adaptive(make_controlled(1E-12, 1E-12, stepper_type()),
            Rate_pt(6.2e-6,0.0037, 1.4488e-10), x_pt, 0.0, 86400.0, 0.1,
                               std::cout << phoenix::arg_names::arg2 << " "
                                         << phoenix::arg_names::arg1 << "\n");
        double const sigma_pt = x_pt;
        double ah = 3.8e-4;
        double acl = 0.0;
        double ad = 1.2e-4;
        double u_norm = 0.0037, net_u_norm = 0.0037, net_conc = 0.1,
               concentration = 1.4488e-10
            ;
        size_t num_of_steps_d = integrate_adaptive(
            make_controlled(1E-12, 1E-12, stepper_type()),
            Rate_d(ah, acl, ad, net_u_norm, net_conc, u_norm, concentration), x_d,
            0.0, 86400.0, 0.1,
            std::cout << phoenix::arg_names::arg2 << " "
                      << phoenix::arg_names::arg1 << "\n");
        double const sigma_d = x_d;

        for (unsigned ip(0); ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);

            // \todo the argument to getValue() has to be changed for non
            // constant storage model
            auto const specific_storage =
                solid_phase.property(MaterialPropertyLib::PropertyType::storage)
                    .template value<double>(vars);

            auto & ip_data = this->_ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;

            double p_int_pt = 0.0;//co2 pressure oil pressure
            double sw_int_pt = 0.0;//water/brine saturation
            double c_int_pt = 0.0;//salt concentration
            double c_porethoat_int_pt = 0.0;
            double c_porebody_int_pt = 0.0;
            double const particle_density = 2500.0;
            // Order matters: First p, then s,last c
            NumLib::shapeFunctionInterpolate(local_x, N, p_int_pt, sw_int_pt,c_int_pt);

            auto const porosity =
                process_data.porous_media_properties.getPorosity(t, pos)
                .getValue(t, pos, 0.0, p_int_pt);

            auto const& intrinsic_permeability =
                process_data.porous_media_properties.getIntrinsicPermeability(
                    t, pos).getValue(t, pos, 0.0, 0.0);
            //std::cout << intrinsic_permeability << std::endl;

            double& porosity0 = ip_data.porosity_curr;
            double const porosity_old = porosity0;
            double& permeability0 = ip_data.permeability_curr;
            const double K = std::pow(10, intrinsic_permeability(0, 0));
            /*auto const intrinsic_permeability =
                intrinsicPermeability<GlobalDim>
                (process_data.solid_intrinsic_permeability(t, pos));

            auto const porosity =
                intrinsicPermeability<GlobalDim>(process_data.porosity_constant(t, pos));
             */   

            /*vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
                T_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;*/

            // Use the fluid density model to compute the density
            auto const fluid_density =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars);
            auto const nonwet_density =
                nonwetting_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars);

            // Use the viscosity model to compute the viscosity
            auto const liquid_viscosity = liquid_phase
                    .property(MaterialPropertyLib::PropertyType::viscosity)
                    .template value<double>(vars);
            auto const nonwet_viscosity = nonwetting_phase
                .property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars);

            auto const liquid_residual_saturation = liquid_phase
                .property(MaterialPropertyLib::PropertyType::residual_liquid_saturation)
                .template value<double>(vars);
            auto const nonwet_residual_saturation = nonwetting_phase
                .property(MaterialPropertyLib::PropertyType::residual_gas_saturation)
                .template value<double>(vars);
            //relative permeability
            //nonwet phase
            double const swe= (sw_int_pt- liquid_residual_saturation)
                / (1 - nonwet_residual_saturation - liquid_residual_saturation);
            double K_rel_G = std::pow((1 - swe), 2);
            if (K_rel_G > 1)
                K_rel_G = 1;
            double K_rel_L = std::pow((swe), 2);
            if (K_rel_L > 1)
                K_rel_L = 1;
            const double K_over_mu_wet = K / liquid_viscosity;

            const double K_over_mu_nonwet = K / nonwet_viscosity;

            double const diffusion_coeff_component_salt = 1e-9;
            const double R = 8.314;
            const double temperature = 393.15;
            const double drho_gas_d_pg = 1e-8;// M_G_comp / R / temperature;
            // matrix assembly
            Mpp.noalias() += w * (porosity * (1 - sw_int_pt)*drho_gas_d_pg)*N.transpose() * N;
            Mps.noalias() -= w*(porosity * nonwet_density)*
                N.transpose() * N;

            Mss.noalias() += w * (porosity * fluid_density)*N.transpose() * N;
            Mcs.noalias() += w * (porosity*c_int_pt)*N.transpose() * N;
            Mcc.noalias() += w * (porosity*sw_int_pt)*N.transpose() * N;


            Kpp.noalias() +=
                (dNdx.transpose() * K_over_mu_nonwet*K_rel_G* nonwet_density * dNdx
                 ) *
                w;

            Ksp.noalias()+= (dNdx.transpose() * K_over_mu_wet*K_rel_L* fluid_density* dNdx
                ) * 
                w;
            /*Kss.noalias()+= (dNdx.transpose() * K_over_mu_wet*K_rel_L * fluid_density * dNdx
                ) *
                w;*/

            Kcp.noalias()+= (dNdx.transpose() * K_over_mu_wet*K_rel_L * c_int_pt * dNdx
                ) *
                w;
            Kcc.noalias()+= (dNdx.transpose() * diffusion_coeff_component_salt*porosity * sw_int_pt * dNdx
                ) *
                w;
            //phase velocity
            GlobalDimVectorType const velocity_wet =
                process_data.has_gravity
                    ? GlobalDimVectorType(
                          -K_over_mu_wet *
                          K_rel_L*
                          (dNdx * p_nodal_values - fluid_density * b))
                    : GlobalDimVectorType(-K_over_mu_wet * K_rel_L * dNdx *
                                          p_nodal_values);
            GlobalDimVectorType const velocity_nonwetting =
                process_data.has_gravity
                    ? GlobalDimVectorType(
                          -K_over_mu_nonwet * K_rel_G *
                          (dNdx * p_nodal_values - nonwet_density * b))
                    : GlobalDimVectorType(-K_over_mu_nonwet * K_rel_G * dNdx *
                                          p_nodal_values);
            double u_norm = GlobalDim>2 ? (std::pow(velocity_nonwetting(0), 2) +
                            std::pow(velocity_nonwetting(1), 2) + 
                                std::pow(velocity_nonwetting(2), 2))
                :(std::pow(velocity_nonwetting(0), 2) +
                            std::pow(velocity_nonwetting(1), 2));
            Kcc.noalias() +=
                N.transpose() * velocity_wet.transpose() * dNdx * w;

            if (process_data.has_gravity)
            {
                Bp.noalias() +=
                    w * nonwet_density * dNdx.transpose() * K_over_mu_nonwet *K_rel_G* nonwet_density*b;
                Bs.noalias() +=
                    w * fluid_density * dNdx.transpose() * K_over_mu_wet* K_rel_L* fluid_density*b;
                Bc.noalias()+=
                    w * c_int_pt * dNdx.transpose() * K_over_mu_wet* K_rel_L* fluid_density*b;

            }
            double rate_h = -ah * y.*net_u_norm;
            double rate_cl = -acl * y.*net_conc;
                            //Colloidal rel.rate 
            double rate_d = ad * u_norm.*c_int_pt;
            double rate_p =
                apt * u_norm * c_int_pt;
            /* with Oberbeck-Boussing assumption density difference only exists
             * in buoyancy effects */
            //update the porosity
            porosity0 = update_porosity(porosity0, c_porethoat_int_pt,
                                        c_porebody_int_pt, particle_density);
            //calculate the flow efficiency
            double const fe = calculate_fe(0.6, c_porethoat_int_pt);
            //update the permeability
            permeability0 = permeability0 * update_permeability_ratio(porosity0,porosity_old,0,fe,3);
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
                        Mps(row, row) += Mps(row, column);
                        Mps(row, column) = 0.0;
                        Mss(row, row) += Mss(row, column);
                        Mss(row, column) = 0.0;
                        Mcs(row, row) += Mcs(row, column);
                        Mcs(row, column) = 0.0;
                        Mcc(row, row) += Mcc(row, column);
                        Mcc(row, column) = 0.0;
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
            std::make_move_iterator(local_x.begin()),
            std::make_move_iterator(local_x.begin() + local_x.size() / 3));
        std::vector<double> local_s(
            std::make_move_iterator(local_x.begin() + local_x.size() / 3),
            std::make_move_iterator(local_x.begin() + 2*local_x.size() / 3));
        std::vector<double> local_c(
            std::make_move_iterator(local_x.begin() + 2 * local_x.size() / 3),
            std::make_move_iterator(local_x.end()));

        return this->getIntPtDarcyVelocityLocal(t, local_p, local_s, local_c, cache);
    }
};

}  // namespace FinesTransport
}  // namespace ProcessLib
