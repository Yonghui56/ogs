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

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "TwoPhaseFlowWithPSMaterialProperties.h"
#include "TwoPhaseFlowWithPSProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPS
{
template <typename NodalRowVectorType, typename GlobalDimNodalMatrixType,
          typename NodalMatrixType>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        NodalRowVectorType N_, GlobalDimNodalMatrixType dNdx_,
        TwoPhaseFlowWithPSMaterialProperties& material_property_,
        double const& integration_weight_, NodalMatrixType const massOperator_)
        : N(std::move(N_)),
          dNdx(std::move(dNdx_)),
          mat_property(material_property_),
          integration_weight(integration_weight_),
          massOperator(massOperator_)

    {
    }
    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    TwoPhaseFlowWithPSMaterialProperties const& mat_property;
    const double integration_weight;
    NodalMatrixType const massOperator;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
const unsigned NUM_NODAL_DOF = 2;

class TwoPhaseFlowWithPSLocalAssemblerInterface
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
class TwoPhaseFlowWithPSLocalAssembler
    : public TwoPhaseFlowWithPSLocalAssemblerInterface
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
    TwoPhaseFlowWithPSLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        TwoPhaseFlowWithPSProcessData const& process_data)
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
    double CalcEffectiveElementLength(unsigned GlobalDim,
                                      MeshLib::Element const& element,
                                      double const num_nodes)
    {
        double L = 0;
        if (GlobalDim <= 1.1)
            L = element.getContent();
        if (GlobalDim == 2)
        {
            auto* edge = element.getEdge(0);
            double max = edge->getContent();
            delete edge;
            for (int k = 0; k < num_nodes; k++)
            {
                auto* new_edge = element.getEdge(k);
                L = new_edge->getContent();
                if (L > max)
                    max = L;
                delete new_edge;
            }
            L = max;
            // delete edge;
        }
        return L;
    }

     double MLangevin(double v)
    {
        double s = 0.0;
        if (v < 0.01)
            s = v * (1.0 / 3.0 + v * v * (-1.0 / 45.0 + 18.0 / 8505.0 * v * v));
        else if (0.01 <= v && v < 20)
            s = (exp(v) + exp(-v)) / (exp(v) - exp(-v)) - 1 / v;
        //    s = coth(v)-1/v;
        else if (20 <= v)
            s = 1.0;

        return s;
    }
    double CalcSUPGCoefficient(double v_mag, int ip, const double ele_length,
                               const double diff_tensor, double const dt)
    {
        //--------------------------------------------------------------------
        // Collect following information to determine SUPG coefficient
        // + flow velocity
        // + diffusivity coefficient (scalar)
        // + characteristic element length
        // + (Peclet number)
        // vel is a double array

        // Characteristic element length
        // Diffusivity = (effective heat conductivity) / (fluid heat capacity)
        const double dispersion_tensor = diff_tensor;

        double diff = dispersion_tensor;

        //--------------------------------------------------------------------
        // Here calculates SUPG coefficient (tau)
        double tau = 0.0;
        switch (2)
        {
            case 1:
            {
                // this coefficient matches with the analytical solution in 1D
                // steady state case
                double alpha = 0.5 * v_mag * ele_length / diff;  // 0.5*Pe
                double func = MLangevin(alpha);
                // tau = 0.5 * ele_length / v_mag * func;
                tau = func / v_mag;
            }
            break;
            case 2:
            {
                // taking into account time step
                //          tau = 1.0 / sqrt(pow(2.0/dt
                //          ,2.0)+pow(2.0*v_mag/ele_len,2.0));
                tau = 1.0 /
                      std::sqrt((2.0 / dt) * (2.0 / dt) +
                                (2.0 * v_mag / ele_length) *
                                    (2.0 * v_mag / ele_length) +
                                (4.0 * diff / (ele_length * ele_length)) *
                                    (4.0 * diff / (ele_length * ele_length)));
            }
            break;
        }

        return tau;
    }

private:
    MeshLib::Element const& _element;

    IntegrationMethod const _integration_method;

    TwoPhaseFlowWithPSProcessData const& _process_data;
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
};

}  // namespace TwoPhaseFlowWithPS
}  // namespace ProcessLib

#include "TwoPhaseFlowWithPSLocalAssembler-impl.h"
