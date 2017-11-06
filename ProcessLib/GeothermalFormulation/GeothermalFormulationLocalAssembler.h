/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "GeothermalFormulationMaterialProperties.h"
#include "GeothermalFormulationProcessData.h"

namespace ProcessLib
{
namespace GeothermalFormulation
{
template <typename NodalRowVectorType, typename GlobalDimNodalMatrixType,
          typename NodalMatrixType>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        NodalRowVectorType N_, GlobalDimNodalMatrixType dNdx_,
        GeothermalFormulationMaterialProperties& material_property_,
        double const& integration_weight_, NodalMatrixType const massOperator_)
        : N(std::move(N_)),
          dNdx(std::move(dNdx_)),
          mat_property(material_property_),
          integration_weight(integration_weight_),
          massOperator(massOperator_),
          sw(1.0),
          h_steam(0.0),
          h_water(0.0),
          dsw_dp(0.0),
          dsw_dh(0.0),
          d_h_steam_d_p(0.0),
          d_h_steam_d_h(0.0),
          d_h_water_d_p(0.0),
          d_h_water_d_h(0.0)

    {
    }
    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    GeothermalFormulationMaterialProperties& mat_property;
    const double integration_weight;
    NodalMatrixType const massOperator;
    double sw;
    double h_steam;
    double h_water;
    double dsw_dp;
    double dsw_dh;
    double d_h_steam_d_p;
    double d_h_steam_d_h;
    double d_h_water_d_p;
    double d_h_water_d_h;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
const unsigned NUM_NODAL_DOF = 2;

class GeothermalFormulationLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtWetPressure(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class GeothermalFormulationLocalAssembler
    : public GeothermalFormulationLocalAssemblerInterface
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
    using LocalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using LocalVectorType = typename LocalAssemblerTraits::LocalVector;

public:
    GeothermalFormulationLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        GeothermalFormulationProcessData const& process_data)
        : _element(element),
          _integration_method(integration_order),
          _process_data(process_data),
          _saturation(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _pressure_wet(
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
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
            _ip_data.emplace_back(
                sm.N, sm.dNdx, *_process_data.material,
                sm.integralMeasure * sm.detJ *
                    _integration_method.getWeightedPoint(ip).getWeight(),
                sm.N.transpose() * sm.N * sm.integralMeasure * sm.detJ *
                    _integration_method.getWeightedPoint(ip).getWeight());
        }
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSaturation(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_saturation.empty());
        return _saturation;
    }

    std::vector<double> const& getIntPtWetPressure(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_pressure_wet.empty());
        return _pressure_wet;
    }

private:
    MeshLib::Element const& _element;

    IntegrationMethod const _integration_method;

    GeothermalFormulationProcessData const& _process_data;
    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType,
                             NodalMatrixType>,
        Eigen::aligned_allocator<IntegrationPointData<
            NodalRowVectorType, GlobalDimNodalMatrixType, NodalMatrixType>>>
        _ip_data;

    // output vector for wetting phase saturation with
    // respect to each integration point
    std::vector<double> _saturation;
    // output vector for wetting phase pressure with respect
    // to each integration point
    std::vector<double> _pressure_wet;
    static const int nonwet_pressure_coeff_index = 0;
    static const int cap_pressure_coeff_index = 1;

    static const int nonwet_pressure_matrix_index = 0;
    static const int cap_pressure_matrix_index = ShapeFunction::NPOINTS;

    static const int nonwet_pressure_size = ShapeFunction::NPOINTS;
    static const int cap_pressure_size = ShapeFunction::NPOINTS;
private:
    const double get_P_sat(double T)
    {
        double P_sat(0.0);
        double T_0 = 373.15;
        double P_0 = 101325.0;
        double h_wg = 2258000.0;
        double R = 8.314;
        double m_w = 0.018;
        return P_0 * exp(((1 / T_0) - (1 / T)) * m_w * h_wg / R);
    }
};

}  // end of namespace
}  // end of namespace

#include "GeothermalFormulationLocalAssembler-impl.h"
