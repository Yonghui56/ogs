/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>
#include <vector>

#include "FinesTransportMaterialProperties.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/PropertyType.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include <boost/array.hpp>
#include "FinesTransportLocalAssemblerInterface.h"

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

namespace ProcessLib
{
namespace FinesTransport
{
template <int GlobalDim>
Eigen::Matrix<double, GlobalDim, GlobalDim> intrinsicPermeability(
    std::vector<double> const& values)
{
    if (boost::get<double>(&values))
    {
        return Eigen::Matrix<double, GlobalDim, GlobalDim>::Identity() *
               boost::get<double>(values);
    }
    if (boost::get<MaterialPropertyLib::Vector>(&values))
    {
        return Eigen::Map<Eigen::Matrix<double, GlobalDim, 1> const>(
                   boost::get<MaterialPropertyLib::Vector>(values).data(),
                   GlobalDim, 1)
            .asDiagonal();
    }
    if (boost::get<MaterialPropertyLib::Tensor2d>(&values))
    {
        return Eigen::Map<Eigen::Matrix<double, GlobalDim, GlobalDim> const>(
            boost::get<MaterialPropertyLib::Tensor2d>(values).data(), GlobalDim,
            GlobalDim);
    }
    if (boost::get<MaterialPropertyLib::Tensor>(&values))
    {
        return Eigen::Map<Eigen::Matrix<double, GlobalDim, GlobalDim> const>(
            boost::get<MaterialPropertyLib::Tensor>(values).data(), GlobalDim,
            GlobalDim);
    }
    OGS_FATAL(
        "Intrinsic permeability parameter values size is neither one nor %d "
        "nor %d squared.",
        GlobalDim, GlobalDim);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class FinesTransportFEM : public FinesTransportLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    FinesTransportFEM(
        MeshLib::Element const& element,
        std::size_t const local_matrix_size,
        bool is_axially_symmetric,
        unsigned const integration_order,
        FinesTransportMaterialProperties const& material_properties,
        const unsigned dof_per_node)
        : FinesTransportLocalAssemblerInterface(),
          _element(element),
          _material_properties(material_properties),
          _integration_method(integration_order)
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * dof_per_node);
        (void)local_matrix_size;
        (void)dof_per_node;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                element, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices[ip].N, shape_matrices[ip].dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ);
        }
    }

    typedef boost::numeric::ublas::vector<double> vector_type;
    typedef boost::numeric::ublas::matrix<double> matrix_type;

    struct stiff_system
    {
        void operator()(const vector_type& x, vector_type& dxdt, double /* t */)
        {
            dxdt[0] = -101.0 * x[0] - 100.0 * x[1];
            dxdt[1] = x[0];
        }
    };

    struct stiff_system_jacobi
    {
        void operator()(const vector_type& /* x */, matrix_type& J,
                        const double& /* t */, vector_type& dfdt)
        {
            J(0, 0) = -101.0;
            J(0, 1) = -100.0;
            J(1, 0) = 1.0;
            J(1, 1) = 0.0;
            dfdt[0] = 0.0;
            dfdt[1] = 0.0;
        }
    };

    struct Rate_pt  // pore throat
    {
        Rate_pt(double apt_ = 0.0, double u_norm_ = 0.0,
                double concentration_ = 0.0)
            : apt(apt_), u_norm(u_norm_), concentration(concentration_)
        {
        }
        void operator()(const double /*x*/, double& dxdt, double /* t */)
        {
            dxdt = apt * u_norm * concentration;
        }
        double apt, u_norm, concentration;
    };

    struct Rate_d  // pore body
    {
        Rate_d(double ah_ = 0.0, double acl_ = 0.0, double ad_ = 0.0,
               double net_u_norm_ = 0.0, double net_conc_ = 0.0,
               double u_norm_ = 0.0, double concentration_ = 0.0)
            : ah(ah_),
              acl(acl_),
              ad(ad_),
              net_u_norm(net_u_norm_),
              net_conc(net_conc_),
              u_norm(u_norm_),
              concentration(concentration_)
        {
        }
        void operator()(const double x, double& dxdt, double /* t */)
        {
            dxdt = -ah * x * net_u_norm            // Hydrodynamic rel.rate
                   - acl * x * net_conc            // Colloidal rel.rate
                   + ad * u_norm * concentration;  // Surface dep. rate
        }
        double ah, acl, ad, net_u_norm, u_norm, net_conc, concentration;
    };

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    /// Computes the flux in the point \c pnt_local_coords that is given in
    /// local coordinates using the values from \c local_x.
    Eigen::Vector3d getFlux(MathLib::Point3d const& pnt_local_coords,
                            double const t,
                            std::vector<double> const& local_x) const override
    {
        // eval dNdx and invJ at given point
        using FemType =
            NumLib::TemplateIsoparametric<ShapeFunction, ShapeMatricesType>;

        FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(
            &this->_element));

        typename ShapeMatricesType::ShapeMatrices shape_matrices(
            ShapeFunction::DIM, GlobalDim, ShapeFunction::NPOINTS);

        // Note: Axial symmetry is set to false here, because we only need dNdx
        // here, which is not affected by axial symmetry.
        fe.computeShapeFunctions(pnt_local_coords.getCoords(), shape_matrices,
                                 GlobalDim, false);

        // fetch permeability, viscosity, density
        ParameterLib::SpatialPosition pos;
        pos.setElementID(this->_element.getID());

        MaterialPropertyLib::VariableArray vars;
        auto const& process_data = this->_material_properties;
        // local_x contains the local pressure values

        double p_int_pt = 0.0;
        double s_int_pt = 0.0;
        double c_int_pt = 0.0;
        NumLib::shapeFunctionInterpolate(local_x, shape_matrices.N, p_int_pt,
                                         s_int_pt, c_int_pt);

        auto const& medium =
            *_material_properties.media_map->getMedium(_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        auto const& solid_phase = medium.phase("Solid");
        auto const& nonwetting_phase = medium.phase("NonAqueousLiquid");
        auto const porosity =
            process_data.porous_media_properties.getPorosity(t, pos).getValue(
                t, pos, 0.0, p_int_pt);
        auto const& intrinsic_permeability =
            process_data.porous_media_properties
                .getIntrinsicPermeability(t, pos)
                .getValue(t, pos, 0.0, 0.0);
        // Use the fluid density model to compute the density
        auto const fluid_density =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars);
        auto const nonwet_density =
            nonwetting_phase
                .property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars);

        auto const liquid_viscosity =
            liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars);
        auto const nonwet_viscosity =
            nonwetting_phase
                .property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars);
        GlobalDimMatrixType K_over_mu_wet =
            intrinsic_permeability / liquid_viscosity;
        GlobalDimMatrixType K_over_mu_nonwet =
            intrinsic_permeability / nonwet_viscosity;
        auto const p_nodal_values =
            Eigen::Map<const NodalVectorType>(&local_x[local_x.size() / 3], 0);
        // relative permeability
        // nonwet phase
        auto const liquid_residual_saturation =
            liquid_phase
                .property(MaterialPropertyLib::PropertyType::
                              residual_liquid_saturation)
                .template value<double>(vars);
        auto const nonwet_residual_saturation =
            nonwetting_phase
                .property(
                    MaterialPropertyLib::PropertyType::residual_gas_saturation)
                .template value<double>(vars);
        double const swe =
            (s_int_pt - liquid_residual_saturation) /
            (1 - nonwet_residual_saturation - liquid_residual_saturation);
        double const K_rel_G = std::pow((1 - swe), 2);
        double const K_rel_L = std::pow((swe), 2);
        GlobalDimVectorType q_wet =
            -K_over_mu_wet * K_rel_L * shape_matrices.dNdx * p_nodal_values;
        GlobalDimVectorType q_nonwet =
            -K_over_mu_nonwet * K_rel_G * shape_matrices.dNdx * p_nodal_values;
        if (this->_material_properties.has_gravity)
        {
            auto const b = this->_material_properties.specific_body_force;
            q_wet += K_over_mu_wet * K_rel_L * fluid_density * b;
            q_nonwet += K_over_mu_nonwet * K_rel_G * nonwet_density * b;
        }

        Eigen::Vector3d flux;
        flux.head<GlobalDim>() = q_wet;
        return flux;
    }
    // fe: flow efficiency factor
    // expressing the fraction of unplugged pores available  for flows
    // sigma_fe: coefficient of flow efficiency
                            //concentration_pt: concentration of particle deposit on the pore throat
    double const calculate_fe(double const sigma_fe,
                              double const concentration_pt)
    {
        return 1 - (sigma_fe * concentration_pt);
    }
    double const update_porosity(double const phi0,
                              double const concentration_porebody,
        double const concentration_porethroat,
        double const density_particle_suspension)
    {
        return phi0 - ((concentration_porebody + concentration_porethroat) /
                       density_particle_suspension);
    }
    double const update_permeability_ratio(double porosity0,
                                           double const porosity_new,
                                           double const kappa, double const fe,
                                           double const exponential_n)
    {
        return std::pow((1 - fe) * kappa + fe * (porosity_new / porosity0),
                        exponential_n);
    }
    /*double const diffusion_coef(, double const temperature, double const viscosity)
    {
        double const k = 1.3806504e-23;
        double const 
        return (k * temperature) / (3 * pi * viscosity * this.diameter);
    }*/

protected:
    MeshLib::Element const& _element;
    FinesTransportMaterialProperties const& _material_properties;

    double apt;
    double u_norm;
    double concentration;

    IntegrationMethod const _integration_method;
    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>,
        Eigen::aligned_allocator<
            IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>>>
        _ip_data;

    std::vector<double> const& getIntPtDarcyVelocityLocal(
        const double t, std::vector<double> const& local_p,
        std::vector<double> const& local_s, std::vector<double> const& local_c,
        std::vector<double>& cache) const
    {
        auto const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<
            Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, GlobalDim, n_integration_points);

        /*ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        MaterialPropertyLib::VariableArray vars;
        auto const& process_data = this->_material_properties;
        auto const p_nodal_values = Eigen::Map<const NodalVectorType>(
            &local_p[0], ShapeFunction::NPOINTS);

        auto const& medium =
            *_material_properties.media_map->getMedium(_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        auto const& solid_phase = medium.phase("Solid");
        auto const& nonwetting_phase = medium.phase("NonAqueousLiquid");

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;

            pos.setIntegrationPoint(ip);

            double p_int_pt = 0.0;
            double s_int_pt = 0.0;
            double c_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_p, N, p_int_pt);
            NumLib::shapeFunctionInterpolate(local_s, N, s_int_pt);
            NumLib::shapeFunctionInterpolate(local_c, N, c_int_pt);

            auto const porosity =
                process_data.porous_media_properties.getPorosity(t, pos)
                .getValue(t, pos, 0.0, p_int_pt);
            auto const& intrinsic_permeability =
                process_data.porous_media_properties.getIntrinsicPermeability(
                    t, pos).getValue(t, pos, 0.0, 0.0);
            // Use the fluid density model to compute the density
            auto const fluid_density =
                liquid_phase
                .property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars);
            auto const nonwet_density =
                nonwetting_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars);

            auto const liquid_viscosity = liquid_phase
                .property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars);
            auto const nonwet_viscosity = nonwetting_phase
                .property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars);
            GlobalDimMatrixType K_over_mu_wet = intrinsic_permeability /
        liquid_viscosity; GlobalDimMatrixType K_over_mu_nonwet =
        intrinsic_permeability / nonwet_viscosity; auto const
        liquid_residual_saturation = liquid_phase
                .property(MaterialPropertyLib::PropertyType::residual_liquid_saturation)
                .template value<double>(vars);
            auto const nonwet_residual_saturation = nonwetting_phase
                .property(MaterialPropertyLib::PropertyType::residual_gas_saturation)
                .template value<double>(vars);
            double const swe = (s_int_pt - liquid_residual_saturation)
                / (1 - nonwet_residual_saturation - liquid_residual_saturation);
            double const K_rel_G = std::pow((1 - swe), 2);
            double const K_rel_L = std::pow((swe), 2);
            cache_mat.col(ip).noalias() = -K_over_mu_wet * K_rel_L* dNdx *
        p_nodal_values;

            if (_material_properties.has_gravity)
            {
                auto const b = _material_properties.specific_body_force;
                // here it is assumed that the vector b is directed 'downwards'
                cache_mat.col(ip).noalias() += K_over_mu_wet * K_rel_L *
        nonwet_density * b;
            }
        }*/

        return cache;
    }
};

}  // namespace FinesTransport
}  // namespace ProcessLib
