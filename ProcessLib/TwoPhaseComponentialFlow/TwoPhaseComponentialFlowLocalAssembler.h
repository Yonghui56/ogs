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
        rho_mol_sio2(0.0),
        rho_mol_sio2_prev(0.0),
        porosity(0.0763),
        porosity_prev(0.0763),
        rho_mol_co2_cumul_total(0.0),
        rho_mol_co2_cumul_total_prev(0.0),
        fluid_volume(0.0),
        fluid_volume_prev(0.0)
    {
    }
    TwoPhaseComponentialFlowMaterialProperties& mat_property;
    double rho_mol_sio2;
    double rho_mol_sio2_prev;
    double porosity;
    double porosity_prev;
    double rho_mol_co2_cumul_total;
    double rho_mol_co2_cumul_total_prev;
    double fluid_volume;
    double fluid_volume_prev;
    double integration_weight;
    NodalMatrixType massOperator;
    NodalMatrixType diffusionOperator;

    void pushBackState()
    {
        rho_mol_sio2_prev = rho_mol_sio2;
        porosity_prev = porosity;
        rho_mol_co2_cumul_total_prev = rho_mol_co2_cumul_total;
        fluid_volume_prev = fluid_volume;
    }
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
    virtual std::vector<double> const& getIntPtpHValue(
        std::vector<double>& /*cache*/) const = 0;
    virtual std::vector<double> const& getIntPtPorosityValue(
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
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _pH_value(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _porosity_value(
            std::vector<double>(_integration_method.getNumberOfPoints())),
          _mol_fraction_nonwet_vapor(
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
            _ip_data[ip].massOperator.noalias() =
                sm.N.transpose() * sm.N * _ip_data[ip].integration_weight;
            _ip_data[ip].diffusionOperator.setZero(ShapeFunction::NPOINTS,
                ShapeFunction::NPOINTS);
            _ip_data[ip].diffusionOperator.noalias() =
                sm.dNdx.transpose() * sm.dNdx * _ip_data[ip].integration_weight;
        }
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    void postTimestepConcrete(std::vector<double> const& /*local_x*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

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

    std::vector<double> const& getIntPtpHValue(
        std::vector<double>& /*cache*/) const override
    {
        assert(_pH_value.size() > 0);
        return _pH_value;
    }
    std::vector<double> const& getIntPtPorosityValue(
        std::vector<double>& /*cache*/) const override
    {
        assert(_porosity_value.size() > 0);
        return _porosity_value;
    }

    std::vector<double> const& getIntPtMolFracNonwetVapor(
        std::vector<double>& /*cache*/) const override
    {
        assert(_mol_fraction_nonwet_vapor.size() > 0);
        return _mol_fraction_nonwet_vapor;
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
    std::vector<double> _pH_value;
    std::vector<double> _porosity_value;
    std::vector<double> _mol_fraction_nonwet_vapor;
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
    const double M_H = 0.002;  // MaterialLib::PhysicalConstant::MolarMass::H2;
    const double M_L = MaterialLib::PhysicalConstant::MolarMass::Water;
    const double M_C = 0.016;  // MaterialLib::PhysicalConstant::MolarMass::CH4;
    const double M_AIR = MaterialLib::PhysicalConstant::MolarMass::Air;
    const double M_CO2 =
        0.044;  // MaterialLib::PhysicalConstant::MolarMass::CO2;
    const double& R = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double Q_steel = 7.8591 * 4/3;  // generate H2
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
                                  double const p_sat, double const kelvin_term)
    {
        double const A = 1 - x1 - x2 - x3;
        double const B =
            1 - (x1 * pg / Hen_L_h + x2 * pg / Hen_L_c + x3 * pg / Hen_L_co2);
        return (B / pg - A / Hen_L_air) / (kelvin_term / p_sat - 1 / Hen_L_air);
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

private:
    const double bi_interpolation(double pnt_x_to_interpolate,
                                  double pnt_y_to_interpolate,
                                  std::vector<double>
                                      _values_at_supp_pnts)
    {
        if (pnt_x_to_interpolate <= _supp_pnts_x.front() &&
            pnt_y_to_interpolate <= _supp_pnts_y.front())
        {
            return _values_at_supp_pnts.front();
        }
        if (pnt_x_to_interpolate >= _supp_pnts_x.back() &&
            pnt_y_to_interpolate >= _supp_pnts_y.back())
        {
            return _values_at_supp_pnts.back();
        }

        const double x_supp_pnt_size = _supp_pnts_x.size();
        if (pnt_x_to_interpolate <= _supp_pnts_x.front() &&
            pnt_y_to_interpolate >= _supp_pnts_y.back())
        {
            return _values_at_supp_pnts[x_supp_pnt_size *
                                        (x_supp_pnt_size - 1)];
        }
        if (pnt_x_to_interpolate >= _supp_pnts_x.back() &&
            pnt_y_to_interpolate <= _supp_pnts_y.front())
        {
            return _values_at_supp_pnts[x_supp_pnt_size - 1];
        }

        auto const& it_X(std::lower_bound(
            _supp_pnts_x.begin(), _supp_pnts_x.end(), pnt_x_to_interpolate));
        std::size_t x_interval_idx =
            std::distance(_supp_pnts_x.begin(), it_X) - 1;
        if (_supp_pnts_x.begin() == it_X)
        {
            x_interval_idx = 0;
            pnt_x_to_interpolate = _supp_pnts_x.front();
        }
        auto const& it_Y(std::lower_bound(
            _supp_pnts_y.begin(), _supp_pnts_y.end(), pnt_y_to_interpolate));
        std::size_t y_interval_idx =
            std::distance(_supp_pnts_y.begin(), it_Y) - 1;
        if (_supp_pnts_y.begin() == it_Y)
        {
            y_interval_idx = 0;
            pnt_y_to_interpolate = _supp_pnts_y.front();
        }

        const double f_r1 =
            _values_at_supp_pnts[x_interval_idx +
                                 x_supp_pnt_size * y_interval_idx] *
                (_supp_pnts_x[x_interval_idx + 1] - pnt_x_to_interpolate) /
                (_supp_pnts_x[x_interval_idx + 1] -
                 _supp_pnts_x[x_interval_idx]) +
            _values_at_supp_pnts[x_interval_idx + 1 +
                                 x_supp_pnt_size * y_interval_idx] *
                (pnt_x_to_interpolate - _supp_pnts_x[x_interval_idx]) /
                (_supp_pnts_x[x_interval_idx + 1] -
                 _supp_pnts_x[x_interval_idx]);
        const double f_r2 =
            _values_at_supp_pnts[x_interval_idx +
                                 x_supp_pnt_size * (y_interval_idx + 1)] *
                (_supp_pnts_x[x_interval_idx + 1] - pnt_x_to_interpolate) /
                (_supp_pnts_x[x_interval_idx + 1] -
                 _supp_pnts_x[x_interval_idx]) +
            _values_at_supp_pnts[x_interval_idx + 1 +
                                 x_supp_pnt_size * (y_interval_idx + 1)] *
                (pnt_x_to_interpolate - _supp_pnts_x[x_interval_idx]) /
                (_supp_pnts_x[x_interval_idx + 1] -
                 _supp_pnts_x[x_interval_idx]);
        const double f_p =
            f_r1 * (_supp_pnts_y[y_interval_idx + 1] - pnt_y_to_interpolate) /
                (_supp_pnts_y[y_interval_idx + 1] -
                 _supp_pnts_y[y_interval_idx]) +
            f_r2 * (pnt_y_to_interpolate - _supp_pnts_y[y_interval_idx]) /
                (_supp_pnts_y[y_interval_idx + 1] -
                 _supp_pnts_y[y_interval_idx]);
        const double f_val_11 =
            _values_at_supp_pnts[x_interval_idx +
                                 x_supp_pnt_size * y_interval_idx] *
            (_supp_pnts_x[x_interval_idx + 1] - pnt_x_to_interpolate) *
            (_supp_pnts_y[y_interval_idx + 1] - pnt_y_to_interpolate) /
            (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]) /
            (_supp_pnts_y[y_interval_idx + 1] - _supp_pnts_y[y_interval_idx]);
        const double f_val_21 =
            _values_at_supp_pnts[x_interval_idx + 1 +
                                 x_supp_pnt_size * y_interval_idx] *
            (pnt_x_to_interpolate - _supp_pnts_x[x_interval_idx]) *
            (_supp_pnts_y[y_interval_idx + 1] - pnt_y_to_interpolate) /
            (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]) /
            (_supp_pnts_y[y_interval_idx + 1] - _supp_pnts_y[y_interval_idx]);
        const double f_val_12 =
            _values_at_supp_pnts[x_interval_idx +
                                 x_supp_pnt_size * (y_interval_idx + 1)] *
            (_supp_pnts_x[x_interval_idx + 1] - pnt_x_to_interpolate) *
            (pnt_y_to_interpolate - _supp_pnts_y[y_interval_idx]) /
            (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]) /
            (_supp_pnts_y[y_interval_idx + 1] - _supp_pnts_y[y_interval_idx]);
        const double f_val_22 =
            _values_at_supp_pnts[x_interval_idx + 1 +
                                 x_supp_pnt_size * (y_interval_idx + 1)] *
            (pnt_x_to_interpolate - _supp_pnts_x[x_interval_idx]) *
            (pnt_y_to_interpolate - _supp_pnts_y[y_interval_idx]) /
            (_supp_pnts_x[x_interval_idx + 1] - _supp_pnts_x[x_interval_idx]) /
            (_supp_pnts_y[y_interval_idx + 1] - _supp_pnts_y[y_interval_idx]);

        return f_val_11 + f_val_21 + f_val_12 + f_val_22;
    }

private:
    std::vector<double> _supp_pnts_x = {
        0,    200,  400,  600,  800,  1000, 1200, 1400, 1600, 1800, 2000,
        2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000};
    std::vector<double> _supp_pnts_y = {
        0,    200,  400,  600,  800,  1000, 1200, 1400, 1600, 1800, 2000,
        2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000};
    std::vector<double> _porosity_at_supp_pnts = {
        0.076293302, 0.075791324, 0.075301303, 0.07481993,  0.074344971,
        0.073700221, 0.070649061, 0.067624029, 0.064498932, 0.060993032,
        0.057774372, 0.054663621, 0.05160242,  0.046826931, 0.045620505,
        0.042636186, 0.044636505, 0.04503369,  0.043701083, 0.043431782,
        0.043431781, 0.075531265, 0.07502836,  0.074537484, 0.074055334,
        0.073579688, 0.071264036, 0.068229574, 0.065208386, 0.06163246,
        0.05838122,  0.055259872, 0.05219524,  0.047852975, 0.046218374,
        0.043272937, 0.048120348, 0.04606729,  0.044782351, 0.044807159,
        0.044807159, 0.044807159, 0.074768968, 0.074265188, 0.073773495,
        0.073290641, 0.071890589, 0.068839676, 0.065812377, 0.062284361,
        0.058991398, 0.055857357, 0.05278887,  0.048883193, 0.046821535,
        0.044111917, 0.051617887, 0.047207129, 0.046186555, 0.046186555,
        0.046186555, 0.046186555, 0.046186555, 0.074006414, 0.073501827,
        0.073009352, 0.072542238, 0.069456436, 0.066420019, 0.062953537,
        0.059605924, 0.056456564, 0.053383697, 0.049916929, 0.04743586,
        0.045248495, 0.05214477,  0.048354136, 0.047569408, 0.047569409,
        0.047569409, 0.047569409, 0.047569409, 0.047569409, 0.073243621,
        0.072738246, 0.072245036, 0.07008379,  0.067032528, 0.063647002,
        0.06022616,  0.057058049, 0.052778919, 0.051008858, 0.048071735,
        0.047782076, 0.051055901, 0.049247964, 0.04895546,  0.048955461,
        0.048955461, 0.048955461, 0.048955461, 0.048955461, 0.048955461,
        0.072480581, 0.071974469, 0.070730773, 0.06765188,  0.06437551,
        0.060853959, 0.057662503, 0.052603071, 0.051614255, 0.048747333,
        0.051301556, 0.051455191, 0.050348775, 0.050344581, 0.050344581,
        0.050344582, 0.050344582, 0.050344582, 0.050344582, 0.050344582,
        0.050344582, 0.071717324, 0.071210524, 0.068281464, 0.065155424,
        0.061492082, 0.05827082,  0.052431241, 0.052226402, 0.049495774,
        0.054827183, 0.052624749, 0.051736694, 0.051736694, 0.051736695,
        0.051736695, 0.051736695, 0.051736695, 0.051736695, 0.051736695,
        0.051736695, 0.051736695, 0.070953851, 0.068928011, 0.065854499,
        0.062144582, 0.058884257, 0.052267996, 0.052852018, 0.050532838,
        0.058358548, 0.053803142, 0.053131745, 0.053131745, 0.053131745,
        0.053131745, 0.053131745, 0.053131745, 0.053131745, 0.053131745,
        0.053131745, 0.053131745, 0.053131745, 0.069608866, 0.066488949,
        0.062818034, 0.059504582, 0.052776835, 0.053501646, 0.051984567,
        0.0580906,   0.054928897, 0.054529681, 0.054529681, 0.054529681,
        0.054529681, 0.054529681, 0.054529681, 0.054529681, 0.054529681,
        0.054529681, 0.054529681, 0.054529681, 0.054529681, 0.067139475,
        0.063523259, 0.060134419, 0.053407303, 0.05419485,  0.054599675,
        0.057093085, 0.055964283, 0.055930455, 0.055930454, 0.055930454,
        0.055930454, 0.055930454, 0.055930454, 0.055930454, 0.055930454,
        0.055930454, 0.055930454, 0.055930454, 0.055930454, 0.055930454,
        0.064279005, 0.060777746, 0.054161328, 0.054967469, 0.058153676,
        0.058172756, 0.057334009, 0.05733401,  0.057334011, 0.057334011,
        0.057334011, 0.057334011, 0.057334011, 0.057334011, 0.057334011,
        0.057334011, 0.057334011, 0.057334011, 0.057334011, 0.057334011,
        0.057334011, 0.061440856, 0.05526379,  0.055882628, 0.061713876,
        0.059453967, 0.058740289, 0.05874029,  0.05874029,  0.05874029,
        0.05874029,  0.05874029,  0.05874029,  0.05874029,  0.05874029,
        0.05874029,  0.05874029,  0.05874029,  0.05874029,  0.05874029,
        0.05874029,  0.05874029,  0.056395643, 0.057045612, 0.065103651,
        0.060742421, 0.060149227, 0.060149227, 0.060149227, 0.060149227,
        0.060149227, 0.060149227, 0.060149227, 0.060149227, 0.060149227,
        0.060149227, 0.060149227, 0.060149227, 0.060149227, 0.060149227,
        0.060149227, 0.060149227, 0.060149227, 0.058988315, 0.064524085,
        0.06167538,  0.061560752, 0.061560752, 0.061560752, 0.061560752,
        0.061560752, 0.061560752, 0.061560752, 0.061560752, 0.061560752,
        0.061560752, 0.061560752, 0.061560752, 0.061560752, 0.061560752,
        0.061560752, 0.061560752, 0.061560752, 0.061560752, 0.064096878,
        0.06297479,  0.062974792, 0.062974792, 0.062974792, 0.062974792,
        0.062974792, 0.062974792, 0.062974792, 0.062974792, 0.062974792,
        0.062974792, 0.062974792, 0.062974792, 0.062974792, 0.062974792,
        0.062974792, 0.062974792, 0.062974792, 0.062974792, 0.062974792,
        0.064391265, 0.064391265, 0.064391265, 0.064391265, 0.064391265,
        0.064391265, 0.064391265, 0.064391265, 0.064391265, 0.064391265,
        0.064391265, 0.064391265, 0.064391265, 0.064391265, 0.064391265,
        0.064391265, 0.064391265, 0.064391265, 0.064391265, 0.064391265,
        0.064391265, 0.065810094, 0.065810094, 0.065810094, 0.065810094,
        0.065810094, 0.065810094, 0.065810094, 0.065810094, 0.065810094,
        0.065810094, 0.065810094, 0.065810094, 0.065810094, 0.065810094,
        0.065810094, 0.065810094, 0.065810094, 0.065810094, 0.065810094,
        0.065810094, 0.065810094, 0.067231195, 0.067231195, 0.067231195,
        0.067231195, 0.067231195, 0.067231195, 0.067231195, 0.067231195,
        0.067231195, 0.067231195, 0.067231195, 0.067231195, 0.067231195,
        0.067231195, 0.067231195, 0.067231195, 0.067231195, 0.067231195,
        0.067231195, 0.067231195, 0.067231195, 0.070435154, 0.070435154,
        0.070435154, 0.070435154, 0.070435154, 0.070435154, 0.070435154,
        0.070435154, 0.070435154, 0.070435154, 0.070435154, 0.070435154,
        0.070435154, 0.070435154, 0.070435154, 0.070435154, 0.070435154,
        0.070435154, 0.070435154, 0.070435154, 0.070435154, 0.074986518,
        0.074986524, 0.07498652,  0.074986521, 0.074986521, 0.074986518,
        0.074986521, 0.074986524, 0.074986521, 0.074986519, 0.074986519,
        0.074986524, 0.074986524, 0.07498652,  0.074986517, 0.074986531,
        0.074986524, 0.07498652,  0.07498653,  0.074986531, 0.074986524,
        0.074986522, 0.074986517, 0.07498653,  0.074986516, 0.07498652,
        0.07498652,  0.07498652,  0.07498653,  0.074986522, 0.074986526,
        0.074986522, 0.074986516, 0.074986516, 0.074986531, 0.07498652,
        0.074986531, 0.07498653,  0.074986516, 0.074986526, 0.07498652,
        0.074986516};
    std::vector<double> _pH_at_supp_pnt = {
        12.88155,  12.83001,  12.783924, 12.742414, 12.704776, 12.651733,
        12.437222, 12.289936, 12.163045, 12.047983, 11.925204, 11.791595,
        11.640398, 11.528473, 11.395379, 11.120054, 10.787758, 10.66685,
        10.347723, 9.9006945, 9.9006948, 12.879261, 12.828731, 12.783593,
        12.742995, 12.706253, 12.522938, 12.362857, 12.2312,   12.117826,
        11.997617, 11.866754, 11.717611, 11.623578, 11.469911, 11.139248,
        10.834449, 10.672162, 9.8428832, 9.9396281, 9.9396281, 9.9396281,
        12.876898, 12.827323, 12.783072, 12.743321, 12.620172, 12.436816,
        12.299533, 12.184135, 12.065512, 11.936331, 11.78773,  11.70894,
        11.528072, 11.091922, 10.862513, 10.698922, 9.9646645, 9.9646664,
        9.9646664, 9.9646664, 9.9646664, 12.874477, 12.825806, 12.782392,
        12.745391, 12.5147,   12.367628, 12.248182, 12.130136, 12.001506,
        11.851768, 11.786467, 11.565523, 10.970762, 10.822586, 10.719057,
        9.9824167, 9.9824168, 9.9824168, 9.9824168, 9.9824168, 9.9824168,
        12.872016, 12.8242,   12.781576, 12.600851, 12.437229, 12.311101,
        12.192599, 12.063267, 11.917321, 11.859163, 11.574474, 10.898752,
        10.83763,  10.540133, 9.9959195, 9.9959195, 9.9959195, 9.9959195,
        9.9959195, 9.9959195, 9.9959195, 12.869511, 12.822518, 12.703926,
        12.510424, 12.373926, 12.253953, 12.122519, 11.977844, 11.910323,
        11.544549, 10.913163, 10.747266, 10.028739, 10.006734, 10.006734,
        10.006734, 10.006734, 10.006734, 10.006734, 10.006734, 10.006734,
        12.866983, 12.820777, 12.59013,  12.4376,   12.315228, 12.180126,
        12.035871, 11.953279, 11.462657, 10.924729, 10.759198, 10.015733,
        10.015734, 10.015734, 10.015734, 10.015734, 10.015734, 10.015734,
        10.015734, 10.015734, 10.015734, 12.864432, 12.681352, 12.506659,
        12.377473, 12.236931, 12.091947, 11.985491, 11.322194, 10.934471,
        10.769898, 10.023441, 10.023441, 10.023441, 10.023441, 10.023441,
        10.023441, 10.023441, 10.023441, 10.023441, 10.023441, 10.023441,
        12.795851, 12.583257, 12.441797, 12.29378,  12.138641, 12.002224,
        11.109904, 10.883598, 10.706262, 10.030189, 10.030189, 10.030189,
        10.030189, 10.030189, 10.030189, 10.030189, 10.030189, 10.030189,
        10.030189, 10.030189, 10.030189, 12.669034, 12.50941,  12.351559,
        12.179058, 11.995388, 10.959919, 10.887308, 10.217923, 10.036197,
        10.036197, 10.036197, 10.036197, 10.036197, 10.036197, 10.036197,
        10.036197, 10.036197, 10.036197, 10.036197, 10.036197, 10.036197,
        12.581638, 12.411241, 12.217873, 11.951875, 10.967188, 10.795357,
        10.041616, 10.041616, 10.041616, 10.041616, 10.041616, 10.041616,
        10.041616, 10.041616, 10.041616, 10.041616, 10.041616, 10.041616,
        10.041616, 10.041616, 10.041616, 12.47397,  12.270825, 11.850557,
        10.97378,  10.800934, 10.046556, 10.046556, 10.046557, 10.046557,
        10.046557, 10.046557, 10.046557, 10.046557, 10.046557, 10.046557,
        10.046557, 10.046557, 10.046557, 10.046557, 10.046557, 10.046557,
        12.323641, 11.653496, 10.957089, 10.80607,  10.051096, 10.051096,
        10.051096, 10.051096, 10.051096, 10.051096, 10.051096, 10.051096,
        10.051096, 10.051096, 10.051096, 10.051096, 10.051096, 10.051096,
        10.051096, 10.051096, 10.051096, 11.376403, 10.914699, 10.459941,
        10.055296, 10.055296, 10.055296, 10.055296, 10.055296, 10.055296,
        10.055296, 10.055296, 10.055296, 10.055296, 10.055296, 10.055296,
        10.055296, 10.055296, 10.055296, 10.055296, 10.055296, 10.055296,
        10.814769, 10.059198, 10.059201, 10.059201, 10.059201, 10.059201,
        10.059201, 10.059201, 10.059201, 10.059201, 10.059201, 10.059201,
        10.059201, 10.059201, 10.059201, 10.059201, 10.059201, 10.059201,
        10.059201, 10.059201, 10.059201, 10.062844, 10.062846, 10.062846,
        10.062846, 10.062847, 10.062847, 10.062847, 10.062847, 10.062847,
        10.062847, 10.062847, 10.062847, 10.062847, 10.062847, 10.062847,
        10.062847, 10.062847, 10.062847, 10.062847, 10.062847, 10.062847,
        10.066268, 10.066268, 10.066268, 10.066268, 10.066268, 10.066268,
        10.066268, 10.066268, 10.066268, 10.066268, 10.066268, 10.066268,
        10.066268, 10.066268, 10.066268, 10.066268, 10.066268, 10.066268,
        10.066268, 10.066268, 10.066268, 10.069483, 10.069483, 10.069483,
        10.069483, 10.069483, 10.069483, 10.069483, 10.069483, 10.069483,
        10.069483, 10.069483, 10.069483, 10.069483, 10.069483, 10.069483,
        10.069483, 10.069483, 10.069483, 10.069483, 10.069483, 10.069483,
        8.3499436, 8.3499436, 8.3499436, 8.3499436, 8.3499437, 8.3499437,
        8.3499437, 8.3499437, 8.3499437, 8.3499437, 8.3499437, 8.3499437,
        8.3499437, 8.3499437, 8.3499437, 8.3499437, 8.3499437, 8.3499437,
        8.3499437, 8.3499437, 8.3499437, 6.1892969, 6.1892962, 6.1892966,
        6.1892966, 6.1892966, 6.1892968, 6.1892966, 6.1892962, 6.1892966,
        6.1892968, 6.1892968, 6.1892962, 6.1892962, 6.1892967, 6.1892969,
        6.1892955, 6.1892963, 6.1892967, 6.1892956, 6.1892955, 6.1892962,
        6.1892961, 6.1892966, 6.1892955, 6.1892967, 6.1892964, 6.1892964,
        6.1892964, 6.1892955, 6.1892962, 6.1892959, 6.1892962, 6.1892967,
        6.1892967, 6.1892954, 6.1892964, 6.1892955, 6.1892955, 6.1892967,
        6.1892959, 6.1892963, 6.1892966};
    std::vector<double> _quartz_rate_suppt_pnt = {
        -3.91256E-07, -3.66299E-07, -3.45079E-07, -3.26791E-07, -3.10841E-07,
        -2.90458E-07, -2.2536E-07,  -1.88907E-07, -1.62099E-07, -1.4099E-07,
        -1.21528E-07, -1.0343E-07,  -8.62229E-08, -7.52134E-08, -6.39546E-08,
        -4.59015E-08, -3.03816E-08, -2.56009E-08, -1.48896E-08, 1.10736E-14,
        1.85061E-15,  -3.90226E-07, -3.6576E-07,  -3.44948E-07, -3.27009E-07,
        -3.1137E-07,  -2.50431E-07, -2.06864E-07, -1.76544E-07, -1.53861E-07,
        -1.33021E-07, -1.13574E-07, -9.49014E-08, -8.45151E-08, -7.01537E-08,
        -4.71147E-08, -3.22914E-08, -2.57082E-08, 4.83842E-09,  3.29348E-15,
        7.94098E-16,  4.59299E-16,  -3.89166E-07, -3.65168E-07, -3.44741E-07,
        -3.27132E-07, -2.81992E-07, -2.26788E-07, -1.92309E-07, -1.67218E-07,
        -1.4484E-07,  -1.23908E-07, -1.03598E-07, -9.39037E-08, -7.55031E-08,
        -4.47103E-08, -3.3593E-08,  -2.67059E-08, 8.05755E-14,  1.10712E-15,
        3.68018E-16,  2.11449E-16,  1.46549E-16,  -3.88083E-07, -3.64531E-07,
        -3.44471E-07, -3.27913E-07, -2.49745E-07, -2.09414E-07, -1.81253E-07,
        -1.57109E-07, -1.34493E-07, -1.12296E-07, -1.03394E-07, -7.93255E-08,
        -3.8835E-08,  -3.20034E-08, -2.75302E-08, 5.41014E-15,  3.1266E-16,
        1.55185E-16,  1.00648E-16,  7.44277E-17,  5.90647E-17,  -3.86985E-07,
        -3.63857E-07, -3.44147E-07, -2.77644E-07, -2.28426E-07, -1.96202E-07,
        -1.69986E-07, -1.45403E-07, -1.21947E-07, -1.13209E-07, -8.06202E-08,
        -3.57839E-08, -3.27951E-08, -2.06463E-08, 1.1865E-15,   2.98891E-16,
        1.85333E-16,  1.31651E-16,  1.02176E-16,  8.35151E-17,  7.06377E-17,
        -3.85871E-07, -3.63154E-07, -3.14717E-07, -2.50185E-07, -2.12352E-07,
        -1.83676E-07, -1.56738E-07, -1.31651E-07, -1.20891E-07, -7.83004E-08,
        -3.66395E-08, -2.90538E-08, -1.08162E-09, 1.02888E-15,  4.45833E-16,
        2.87917E-16,  2.12771E-16,  1.68751E-16,  1.39855E-16,  1.19399E-16,
        1.04151E-16,  -3.8475E-07,  -3.62426E-07, -2.76067E-07, -2.30045E-07,
        -1.98442E-07, -1.6863E-07,  -1.41712E-07, -1.27868E-07, -7.15481E-08,
        -3.73914E-08, -2.96632E-08, 1.59607E-14,  1.278E-15,    6.50512E-16,
        4.40507E-16,  3.33007E-16,  2.67663E-16,  2.23773E-16,  1.92238E-16,
        1.685E-16,    1.49975E-16,  -3.83621E-07, -3.08677E-07, -2.50753E-07,
        -2.14625E-07, -1.81248E-07, -1.52193E-07, -1.33569E-07, -6.09903E-08,
        -3.80769E-08, -3.02411E-08, 3.33128E-15,  8.18059E-16,  4.92102E-16,
        3.44264E-16,  2.66151E-16,  2.16174E-16,  1.81999E-16,  1.57145E-16,
        1.38269E-16,  1.23439E-16,  1.11472E-16,  -3.54497E-07, -2.75697E-07,
        -2.32675E-07, -1.94813E-07, -1.61667E-07, -1.37037E-07, -4.76588E-08,
        -3.58045E-08, -2.76228E-08, 3.00758E-15,  1.09258E-15,  8.91208E-16,
        6.42763E-16,  5.02311E-16,  4.124E-16,    3.4981E-16,   3.03693E-16,
        2.6833E-16,   2.40331E-16,  2.17634E-16,  1.98851E-16,  -3.06328E-07,
        -2.53189E-07, -2.09608E-07, -1.70476E-07, -1.36797E-07, -4.00262E-08,
        -3.61631E-08, -8.52227E-09, 2.26275E-15,  1.82482E-15,  1.18177E-15,
        8.74245E-16,  6.93706E-16,  5.74965E-16,  4.90946E-16,  4.28338E-16,
        3.79889E-16,  3.41292E-16,  3.09815E-16,  2.83652E-16,  2.61567E-16,
        -2.76967E-07, -2.26013E-07, -1.79425E-07, -1.30858E-07, -4.06369E-08,
        -3.20071E-08, 3.05073E-14,  2.67392E-15,  1.30164E-15,  8.81092E-16,
        6.65937E-16,  5.35209E-16,  4.4739E-16,   3.84317E-16,  3.36846E-16,
        2.99806E-16,  2.70105E-16,  2.45765E-16,  2.25436E-16,  2.08218E-16,
        1.93449E-16,  -2.44549E-07, -1.9197E-07,  -1.17023E-07, -4.1223E-08,
        -3.24313E-08, 1.11783E-14,  2.89404E-15,  1.66326E-15,  1.16696E-15,
        8.98776E-16,  7.30826E-16,  6.15756E-16,  5.31986E-16,  4.68276E-16,
        4.18197E-16,  3.778E-16,    3.44527E-16,  3.16625E-16,  2.92906E-16,
        2.72496E-16,  2.54751E-16,  -2.05351E-07, -9.3544E-08,  -4.05753E-08,
        -3.2843E-08,  9.21273E-15,  3.43115E-15,  2.1147E-15,   1.5283E-15,
        1.19653E-15,  9.83121E-16,  8.34304E-16,  7.24624E-16,  6.40427E-16,
        5.73764E-16,  5.19665E-16,  4.74887E-16,  4.37219E-16,  4.05083E-16,
        3.77356E-16,  3.53175E-16,  3.31906E-16,  -6.79493E-08, -3.8636E-08,
        -1.82707E-08, 8.87625E-15,  4.10311E-15,  2.66846E-15,  1.97719E-15,
        1.57037E-15,  1.3024E-15,   1.11256E-15,  9.71022E-16,  8.6143E-16,
        7.7406E-16,   7.02789E-16,  6.4354E-16,   5.93497E-16,  5.50679E-16,
        5.1362E-16,   4.81239E-16,  4.52698E-16,  4.27353E-16,  -3.38386E-08,
        1.13281E-13,  9.33737E-15,  4.91662E-15,  3.33699E-15,  2.52561E-15,
        2.03163E-15,  1.69928E-15,  1.46039E-15,  1.28039E-15,  1.1399E-15,
        1.02719E-15,  9.34766E-16,  8.576E-16,    7.92197E-16,  7.36069E-16,
        6.87364E-16,  6.44715E-16,  6.07037E-16,  5.7352E-16,   5.43519E-16,
        9.89441E-14,  2.55624E-14,  1.4727E-14,   1.03432E-14,  7.97082E-15,
        6.48378E-15,  5.46435E-15,  4.72194E-15,  4.15712E-15,  3.71301E-15,
        3.35461E-15,  3.05932E-15,  2.81181E-15,  2.60135E-15,  2.4202E-15,
        2.26264E-15,  2.12435E-15,  2.00197E-15,  1.89294E-15,  1.79516E-15,
        1.70699E-15,  9.79437E-15,  6.02158E-15,  4.34731E-15,  3.40152E-15,
        2.79374E-15,  2.37022E-15,  2.0582E-15,   1.81879E-15,  1.62926E-15,
        1.4755E-15,   1.34827E-15,  1.24123E-15,  1.14994E-15,  1.07116E-15,
        1.00248E-15,  9.42079E-16,  8.88544E-16,  8.40764E-16,  7.97854E-16,
        7.59121E-16,  7.23969E-16,  3.40664E-15,  2.52905E-15,  2.01102E-15,
        1.66912E-15,  1.42659E-15,  1.24559E-15,  1.10536E-15,  9.93507E-16,
        9.0221E-16,   8.26285E-16,  7.62141E-16,  7.07235E-16,  6.59716E-16,
        6.18179E-16,  5.81562E-16,  5.49041E-16,  5.19964E-16,  4.9381E-16,
        4.70164E-16,  4.48678E-16,  4.29071E-16,  8.44179E-16,  6.83602E-16,
        5.74349E-16,  4.95205E-16,  4.3523E-16,   3.88214E-16,  3.50367E-16,
        3.19242E-16,  2.93196E-16,  2.7108E-16,   2.52065E-16,  2.35545E-16,
        2.21056E-16,  2.08247E-16,  1.9684E-16,   1.86619E-16,  1.77406E-16,
        1.6906E-16,   1.61464E-16,  1.54522E-16,  1.48151E-16,  9.30664E-17,
        4.07837E-17,  5.99993E-17,  5.09864E-17,  4.52325E-17,  4.1873E-17,
        3.69025E-17,  2.82936E-17,  3.11605E-17,  3.18221E-17,  2.96808E-17,
        2.11481E-17,  1.98914E-17,  2.48341E-17,  2.54637E-17,  2.37618E-18,
        1.78009E-17,  2.02896E-17,  2.83545E-18,  1.64939E-18,  9.75291E-18,
        6.89001E-17,  7.30387E-17,  7.15597E-18,  5.73125E-17,  4.39844E-17,
        3.95239E-17,  3.5884E-17,   4.85791E-18,  2.73689E-17,  2.09815E-17,
        2.36829E-17,  2.84297E-17,  2.67448E-17,  2.49231E-18,  2.0663E-17,
        2.16216E-18,  2.41142E-18,  2.06315E-17,  1.27078E-17,  1.6333E-17,
        1.81428E-17};
    std::vector<double> _fluid_volume_suppt_pnt = {
        7.62947E-05, 7.26806E-05, 6.90664E-05, 6.54522E-05, 6.18379E-05,
        5.82245E-05, 5.46175E-05, 5.10066E-05, 4.73559E-05, 4.35668E-05,
        3.98794E-05, 3.62298E-05, 3.25974E-05, 2.79434E-05, 2.63955E-05,
        2.27998E-05, 2.49507E-05, 2.60633E-05, 2.48393E-05, 2.42177E-05,
        2.42177E-05, 7.99077E-05, 7.62935E-05, 7.26793E-05, 6.90651E-05,
        6.54508E-05, 6.18437E-05, 5.82341E-05, 5.46229E-05, 5.08099E-05,
        4.71118E-05, 4.34593E-05, 3.98264E-05, 3.56687E-05, 3.36288E-05,
        3.00658E-05, 3.5491E-05,  3.39125E-05, 3.22782E-05, 3.23473E-05,
        3.23473E-05, 3.23473E-05, 8.35206E-05, 7.99064E-05, 7.62922E-05,
        7.2678E-05,  6.90687E-05, 6.54613E-05, 6.18509E-05, 5.80568E-05,
        5.4345E-05,  5.06889E-05, 4.70555E-05, 4.33981E-05, 4.08665E-05,
        3.75508E-05, 4.60421E-05, 4.19191E-05, 4.04783E-05, 4.04783E-05,
        4.04783E-05, 4.04783E-05, 4.04783E-05, 8.71335E-05, 8.35194E-05,
        7.99051E-05, 7.62907E-05, 7.26878E-05, 6.90789E-05, 6.53091E-05,
        6.15792E-05, 5.79189E-05, 5.42852E-05, 5.11304E-05, 4.81132E-05,
        4.53484E-05, 5.33373E-05, 4.99287E-05, 4.86107E-05, 4.86107E-05,
        4.86107E-05, 4.86107E-05, 4.86107E-05, 4.86107E-05, 9.07464E-05,
        8.71323E-05, 8.3518E-05,  7.99134E-05, 7.63066E-05, 7.25692E-05,
        6.88149E-05, 6.51496E-05, 6.04598E-05, 5.89279E-05, 5.53778E-05,
        5.47454E-05, 5.9067E-05,  5.74493E-05, 5.67442E-05, 5.67442E-05,
        5.67442E-05, 5.67442E-05, 5.67442E-05, 5.67442E-05, 5.67442E-05,
        9.43594E-05, 9.07452E-05, 8.71371E-05, 8.35338E-05, 7.98405E-05,
        7.60526E-05, 7.23812E-05, 6.70102E-05, 6.61649E-05, 6.26751E-05,
        6.53059E-05, 6.61662E-05, 6.48925E-05, 6.48788E-05, 6.48788E-05,
        6.48788E-05, 6.48788E-05, 6.48788E-05, 6.48788E-05, 6.48788E-05,
        6.48788E-05, 9.79723E-05, 9.4358E-05,  9.07602E-05, 8.71288E-05,
        8.32931E-05, 7.9614E-05,  7.35627E-05, 7.34066E-05, 7.00325E-05,
        7.58684E-05, 7.4185E-05,  7.30145E-05, 7.30145E-05, 7.30145E-05,
        7.30145E-05, 7.30145E-05, 7.30145E-05, 7.30145E-05, 7.30145E-05,
        7.30145E-05, 7.30145E-05, 0.000101585, 9.79849E-05, 9.43853E-05,
        9.05376E-05, 8.68488E-05, 8.01211E-05, 8.06591E-05, 7.76867E-05,
        8.64326E-05, 8.22074E-05, 8.11512E-05, 8.11512E-05, 8.11512E-05,
        8.11512E-05, 8.11512E-05, 8.11512E-05, 8.11512E-05, 8.11512E-05,
        8.11512E-05, 8.11512E-05, 8.11512E-05, 0.000105206, 0.000101614,
        9.77883E-05, 9.40862E-05, 8.7267E-05,  8.79314E-05, 8.57557E-05,
        9.29028E-05, 9.00908E-05, 8.9289E-05,  8.9289E-05,  8.9289E-05,
        8.9289E-05,  8.9289E-05,  8.9289E-05,  8.9289E-05,  8.9289E-05,
        8.9289E-05,  8.9289E-05,  8.9289E-05,  8.9289E-05,  0.000108841,
        0.000105048, 0.000101327, 9.45172E-05, 9.52402E-05, 9.51549E-05,
        9.86927E-05, 9.75293E-05, 9.74278E-05, 9.74278E-05, 9.74278E-05,
        9.74278E-05, 9.74278E-05, 9.74278E-05, 9.74278E-05, 9.74278E-05,
        9.74278E-05, 9.74278E-05, 9.74278E-05, 9.74278E-05, 9.74278E-05,
        0.000112324, 0.000108575, 0.000101901, 0.000102616, 0.000105724,
        0.000106504, 0.000105568, 0.000105568, 0.000105568, 0.000105568,
        0.000105568, 0.000105568, 0.000105568, 0.000105568, 0.000105568,
        0.000105568, 0.000105568, 0.000105568, 0.000105568, 0.000105568,
        0.000105568, 0.000115831, 0.000109673, 0.000110111, 0.000116294,
        0.00011459,  0.000113708, 0.000113708, 0.000113708, 0.000113708,
        0.000113708, 0.000113708, 0.000113708, 0.000113708, 0.000113708,
        0.000113708, 0.000113708, 0.000113708, 0.000113708, 0.000113708,
        0.000113708, 0.000113708, 0.000117464, 0.00011781,  0.000126657,
        0.00012268,  0.00012185,  0.00012185,  0.00012185,  0.00012185,
        0.00012185,  0.00012185,  0.00012185,  0.00012185,  0.00012185,
        0.00012185,  0.00012185,  0.00012185,  0.00012185,  0.00012185,
        0.00012185,  0.00012185,  0.00012185,  0.000126347, 0.000132762,
        0.00013021,  0.000129993, 0.000129993, 0.000129993, 0.000129993,
        0.000129993, 0.000129993, 0.000129993, 0.000129993, 0.000129993,
        0.000129993, 0.000129993, 0.000129993, 0.000129993, 0.000129993,
        0.000129993, 0.000129993, 0.000129993, 0.000129993, 0.00013906,
        0.000138136, 0.000138136, 0.000138136, 0.000138136, 0.000138136,
        0.000138136, 0.000138136, 0.000138136, 0.000138136, 0.000138136,
        0.000138136, 0.000138136, 0.000138136, 0.000138136, 0.000138136,
        0.000138136, 0.000138136, 0.000138136, 0.000138136, 0.000138136,
        0.000146281, 0.000146281, 0.000146281, 0.000146281, 0.000146281,
        0.000146281, 0.000146281, 0.000146281, 0.000146281, 0.000146281,
        0.000146281, 0.000146281, 0.000146281, 0.000146281, 0.000146281,
        0.000146281, 0.000146281, 0.000146281, 0.000146281, 0.000146281,
        0.000146281, 0.000154426, 0.000154426, 0.000154426, 0.000154426,
        0.000154426, 0.000154426, 0.000154426, 0.000154426, 0.000154426,
        0.000154426, 0.000154426, 0.000154426, 0.000154426, 0.000154426,
        0.000154426, 0.000154426, 0.000154426, 0.000154426, 0.000154426,
        0.000154426, 0.000154426, 0.000162572, 0.000162572, 0.000162572,
        0.000162572, 0.000162572, 0.000162572, 0.000162572, 0.000162572,
        0.000162572, 0.000162572, 0.000162572, 0.000162572, 0.000162572,
        0.000162572, 0.000162572, 0.000162572, 0.000162572, 0.000162572,
        0.000162572, 0.000162572, 0.000162572, 0.000171739, 0.000171739,
        0.000171739, 0.000171739, 0.000171739, 0.000171739, 0.000171739,
        0.000171739, 0.000171739, 0.000171739, 0.000171739, 0.000171739,
        0.000171739, 0.000171739, 0.000171739, 0.000171739, 0.000171739,
        0.000171739, 0.000171739, 0.000171739, 0.000171739, 0.000177038,
        0.000177038, 0.000177038, 0.000177038, 0.000177038, 0.000177038,
        0.000177038, 0.000177038, 0.000177038, 0.000177038, 0.000177038,
        0.000177038, 0.000177038, 0.000177038, 0.000177038, 0.000177038,
        0.000177038, 0.000177038, 0.000177038, 0.000177038, 0.000177038,
        0.000177038, 0.000177038, 0.000177038, 0.000177038, 0.000177038,
        0.000177038, 0.000177038, 0.000177038, 0.000177038, 0.000177038,
        0.000177038, 0.000177038, 0.000177038, 0.000177038, 0.000177038,
        0.000177038, 0.000177038, 0.000177038, 0.000177038, 0.000177038,
        0.000177038};
    std::vector<double> _flag_carbon_suppt_pnt = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1};
};

}  // end of namespace
}  // end of namespace

#include "TwoPhaseComponentialFlowLocalAssembler-impl.h"
