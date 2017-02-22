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
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "TwoPhaseComponentialFlowProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseComponentialFlow
{
template <typename NodalMatrixType>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        TwoPhaseComponentialFlowMaterialProperties& material_property_)
        : mat_property(material_property_),
        sw(1.0),
        rho_m(0.0),
        dsw_dpg(0.0),
        dsw_drho(0.0),
        drhom_dpg(0.0),
        drhom_drho(0.0)
    {
    }
    TwoPhaseComponentialFlowMaterialProperties& mat_property;
    double sw;
    double rho_m;
    double dsw_dpg;
    double dsw_drho;
    double drhom_dpg;
    double drhom_drho;
    double pressure_nonwetting;
    double integration_weight;
    NodalMatrixType massOperator;
    NodalMatrixType diffusionOperator;
};
const unsigned NUM_NODAL_DOF = 5;

class TwoPhaseComponentialFlowLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtWettingPressure(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class TwoPhaseComponentialFlowLocalAssembler
    : public TwoPhaseComponentialFlowLocalAssemblerInterface
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
    TwoPhaseComponentialFlowLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        TwoPhaseComponentialFlowProcessData const& process_data)
        : _element(element),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _process_data(process_data),
          _saturation(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _pressure_wetting(
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(*_process_data._material);
            auto const& sm = _shape_matrices[ip];
            _ip_data[ip].integration_weight =
                sm.integralMeasure * sm.detJ *
                _integration_method.getWeightedPoint(ip).getWeight();
            _ip_data[ip].massOperator.setZero(ShapeFunction::NPOINTS,
                                              ShapeFunction::NPOINTS);
            _ip_data[ip].diffusionOperator.setZero(ShapeFunction::NPOINTS,
                                                   ShapeFunction::NPOINTS);
            _ip_data[ip].massOperator.noalias() =
                sm.N.transpose() * sm.N * _ip_data[ip].integration_weight;
            _ip_data[ip].diffusionOperator.noalias() =
                sm.dNdx.transpose() * sm.dNdx * _ip_data[ip].integration_weight;
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

    std::vector<double> const& getIntPtWettingPressure(
        std::vector<double>& /*cache*/) const override
    {
        assert(_pressure_wetting.size() > 0);
        return _pressure_wetting;
    }

private:
    MeshLib::Element const& _element;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;

    TwoPhaseComponentialFlowProcessData const& _process_data;
    std::vector<IntegrationPointData<NodalMatrixType>,
                Eigen::aligned_allocator<IntegrationPointData<NodalMatrixType>>>
        _ip_data;

    std::vector<double> _saturation;  /// used for secondary variable output
    std::vector<double> _pressure_wetting;
    static const int nonwet_pressure_coeff_index = 0;
    static const int mol_fraction_h_coeff_index = 1;
    static const int mol_fraction_ch4_coeff_index = 2;
    static const int mol_fraction_co2_coeff_index = 3;
    static const int cap_pressure_coeff_index = 4;

    static const int nonwet_pressure_matrix_index = 0;
    static const int mol_fraction_h2_matrix_index = ShapeFunction::NPOINTS;
    static const int mol_fraction_ch4_matrix_index = 2 * ShapeFunction::NPOINTS;
    static const int mol_fraction_co2_matrix_index = 3 * ShapeFunction::NPOINTS;
    static const int cap_pressure_matrix_index = 4 * ShapeFunction::NPOINTS;

    static const int nonwet_pressure_size = ShapeFunction::NPOINTS;
    static const int cap_pressure_size = ShapeFunction::NPOINTS;
private:
    const double Hen_L_h = 7.26e+9;     // Henry constant in [Pa]
    const double Hen_L_c = 4.13e+9;     // Henry constant in [Pa]
    const double Hen_L_air = 9.077e+9;  // Henry constant in [Pa]
    const double Hen_L_co2 = 0.163e+9;  // Henry constant in [Pa]
    const double rho_l_std = 1000.0;
    const double M_H = 0.002;// MaterialLib::PhysicalConstant::MolarMass::H2;
    const double M_L = MaterialLib::PhysicalConstant::MolarMass::Water;
    const double M_C = 0.016;// MaterialLib::PhysicalConstant::MolarMass::CH4;
    const double M_AIR = MaterialLib::PhysicalConstant::MolarMass::Air;
    const double M_CO2 = 0.044;// MaterialLib::PhysicalConstant::MolarMass::CO2;
    const double& R = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double Q_steel = 7.8591;           // generate H2
    const double para_slow = 401.55;
    const double para_fast = 191.47286;
private:
    const double get_P_sat(double T)
    {
        double P_sat(0.0);
        double T_0 = 373.15;
        double P_0 = 101325.0;
        double h_wg = 2258000.0;
        return P_0 * exp(((1 / T_0) - (1 / T)) * M_L * h_wg / R);
    }
    const double get_x_nonwet_air_gp(double PG, double X1, double X2, double X3,
        double P_sat)
    {
        double K_G_w = PG / P_sat;
        double K_G_air = PG / Hen_L_air;
        double L =
            1 - (X1 * PG / Hen_L_h + X2 * PG / Hen_L_c + X3 * PG / Hen_L_co2);
        double G = 1 - X1 - X2 - X3;
        double X_G_air = (K_G_w * G - L) / (K_G_w - K_G_air);
        return X_G_air;
    }
    const double get_X_G_h2o_gp(double PG, double X1, double X2, double X3,
        double P_sat)
    {
        double K_G_w = PG / P_sat;
        double K_G_air = PG / Hen_L_air;
        double L =
            1 - (X1 * PG / Hen_L_h + X2 * PG / Hen_L_c + X3 * PG / Hen_L_co2);
        double G = 1 - X1 - X2 - X3;
        double X_G_h2o = (L - K_G_air * G) / (K_G_w - K_G_air);
        return X_G_h2o;
    }
    const double get_x_nonwet_h2o(double const pg, double const x1,
        double const x2, double const x3,
        double const p_sat)
    {
        double const A = 1 - x1 - x2 - x3;
        double const B =
            1 - (x1 * pg / Hen_L_h + x2 * pg / Hen_L_c + x3 * pg / Hen_L_co2);
        return (B / pg - A / Hen_L_air) / (1 / p_sat - 1 / Hen_L_air);
    }
    const double get_derivative_x_nonwet_h2o_d_pg(double const pg,
        double const /*x1*/,
        double const /*x2*/,
        double const /*x3*/,
        double const p_sat)
    {
        return (1 / (1 / p_sat - 1 / Hen_L_air)) * (-1 / pg / pg);
    }

    const double get_derivative_x_nonwet_h2o_d_x1(double const pg,
        double const /*x1*/,
        double const /*x2*/,
        double const /*x3*/,
        double const p_sat)
    {
        return (1 / (1 / p_sat - 1 / Hen_L_air)) *
            (1 / Hen_L_air - 1 / Hen_L_h);
    }
    const double get_derivative_x_nonwet_h2o_d_x2(double const pg,
        double const /*x1*/,
        double const /*x2*/,
        double const /*x3*/,
        double const p_sat)
    {
        return (1 / (1 / p_sat - 1 / Hen_L_air)) *
            (1 / Hen_L_air - 1 / Hen_L_c);
    }
    const double get_derivative_x_nonwet_h2o_d_x3(double const pg,
        double const /*x1*/,
        double const /*x2*/,
        double const /*x3*/,
        double const p_sat)
    {
        return (1 / (1 / p_sat - 1 / Hen_L_air)) *
            (1 / Hen_L_air - 1 / Hen_L_co2);
    }


};

}  // end of namespace
}  // end of namespace

#include "TwoPhaseComponentialFlowLocalAssembler-impl.h"
