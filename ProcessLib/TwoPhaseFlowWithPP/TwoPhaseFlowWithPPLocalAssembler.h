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

#include "TwoPhaseFlowWithPPMaterialProperties.h"
#include "TwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
template <typename NodalMatrixType>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        TwoPhaseFlowWithPPMaterialProperties& material_property_,
        double const& integration_weight_, NodalMatrixType const massOperator_, NodalMatrixType const diffusion_operator_)
        : mat_property(material_property_),
          integration_weight(integration_weight_),
          massOperator(massOperator_),
        diffusion_operator(diffusion_operator_)

    {
    }
    TwoPhaseFlowWithPPMaterialProperties const& mat_property;

    const double integration_weight;
    NodalMatrixType const massOperator;
    NodalMatrixType const diffusion_operator;
};
const unsigned NUM_NODAL_DOF = 2;

class TwoPhaseFlowWithPPLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtWetPressure(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class TwoPhaseFlowWithPPLocalAssembler
    : public TwoPhaseFlowWithPPLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using LocalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using LocalVectorType = typename LocalAssemblerTraits::LocalVector;

public:
    TwoPhaseFlowWithPPLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        TwoPhaseFlowWithPPProcessData const& process_data)
        : _element(element),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _process_data(process_data),
          _saturation(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _pressure_wet(
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& sm = _shape_matrices[ip];
            _ip_data.emplace_back(
                *_process_data.material,
                sm.integralMeasure * sm.detJ *
                    _integration_method.getWeightedPoint(ip).getWeight(),
                sm.N.transpose() * sm.N  * sm.integralMeasure * sm.detJ *
                _integration_method.getWeightedPoint(ip).getWeight(),
                sm.dNdx.transpose() * sm.dNdx * sm.integralMeasure * sm.detJ *
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
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSaturation(
        std::vector<double>& /*cache*/) const override
    {
        assert(_saturation.size() > 0);
        return _saturation;
    }

    std::vector<double> const& getIntPtWetPressure(
        std::vector<double>& /*cache*/) const override
    {
        assert(_pressure_wet.size() > 0);
        return _pressure_wet;
    }

private:
    MeshLib::Element const& _element;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;

    TwoPhaseFlowWithPPProcessData const& _process_data;
    std::vector<IntegrationPointData<NodalMatrixType>> _ip_data;

    // output vector for wetting phase saturation with
    // respect to each integration point
    std::vector<double> _saturation;
    // output vector for wetting phase pressure with respect
    // to each integration point
    std::vector<double> _pressure_wet;
    const double _hen_L_air = 9.077e+9; // Henry constant in [Pa]
    const double _rho_l_std = 1000;
    const double& _molar_mass_water = MaterialLib::PhysicalConstant::MolarMass::Water;
    const double& _ideal_gas_const = MaterialLib::PhysicalConstant::IdealGasConstant;
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
        return P_0 * exp(((1 / T_0) - (1 / T)) * _molar_mass_water * h_wg / _ideal_gas_const);
    }

    const double get_x_nonwet_vapor(double pg, double p_sat, double kelvin_term)
    {

        return (1 / pg - 1/_hen_L_air)/(kelvin_term/p_sat -1/_hen_L_air);
    }

    const double get_derivative_x_nonwet_h2o_d_pg(double const pg,
        double const p_sat,
        double const kelvin_term)
    {
        return (1 / (kelvin_term / p_sat - 1 / _hen_L_air)) * (-1 / pg / pg);
    }

    const double get_derivative_x_nonwet_h2o_d_kelvin(double const pg,
        double const p_sat,
        double const kelvin_term)
    {
        return -(1 / pg - 1 / _hen_L_air) / (kelvin_term / p_sat - 1 / _hen_L_air) / (kelvin_term / p_sat - 1 / _hen_L_air) / p_sat;
    }

};

}  // end of namespace
}  // end of namespace

#include "TwoPhaseFlowWithPPLocalAssembler-impl.h"
