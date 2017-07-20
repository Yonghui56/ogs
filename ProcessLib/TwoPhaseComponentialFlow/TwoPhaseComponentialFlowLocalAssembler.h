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
          rho_mol_sio2_backfill(0.0),
          rho_mol_sio2_prev_backfill(0.0),
          porosity_backfill(0.095228012),
          porosity_prev_backfill(0.095228012),
          rho_mol_co2_cumul_total_backfill(0.0),
          rho_mol_co2_cumul_total_prev_backfill(0.0),
          fluid_volume_backfill(95.22765
          ),
          fluid_volume_prev_backfill(95.22765),
          porosity_waste(0.2),
          porosity_prev_waste(0.2),
          rho_mol_co2_cumul_total_waste(0.0),
          rho_mol_co2_cumul_total_prev_waste(0.0),
          fluid_volume_waste(9.0687704),
          fluid_volume_prev_waste(9.0687704)

    {
    }
    TwoPhaseComponentialFlowMaterialProperties& mat_property;
    double rho_mol_sio2_backfill;
    double rho_mol_sio2_prev_backfill;
    double porosity_backfill;
    double porosity_prev_backfill;
    double rho_mol_co2_cumul_total_backfill;
    double rho_mol_co2_cumul_total_prev_backfill;
    double fluid_volume_backfill;
    double fluid_volume_prev_backfill;
    double porosity_waste;
    double porosity_prev_waste;
    double rho_mol_co2_cumul_total_waste;
    double rho_mol_co2_cumul_total_prev_waste;
    double fluid_volume_waste;
    double fluid_volume_prev_waste;
    double integration_weight;
    NodalMatrixType massOperator;
    NodalMatrixType diffusionOperator;

    void pushBackState()
    {
        rho_mol_sio2_prev_backfill = rho_mol_sio2_backfill;
        porosity_prev_backfill = porosity_backfill;
        rho_mol_co2_cumul_total_prev_backfill =
            rho_mol_co2_cumul_total_backfill;
        fluid_volume_prev_backfill = fluid_volume_backfill;
        porosity_prev_waste = porosity_waste;
        rho_mol_co2_cumul_total_prev_waste = rho_mol_co2_cumul_total_waste;
        fluid_volume_prev_waste = fluid_volume_waste;
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

    virtual std::vector<double> const& getIntPtMolFracNonwetVapor(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtCO2Concentration(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtRhoMolCo2CumulTotalPrev(
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
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _co2_concentration(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _rho_mol_co2_cumulated_prev(
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

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
        double const /*t*/,
        double const /*delta_t*/) override
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

    std::vector<double> const& getIntPtCO2Concentration(
        std::vector<double>& /*cache*/) const override
    {
        assert(_co2_concentration.size() > 0);
        return _co2_concentration;
    }
    /*
    * used to store previous time step value of rho_mol_co2_cumul_total_prev
    */
    std::vector<double> const& getIntPtRhoMolCo2CumulTotalPrev(
        std::vector<double>& /*cache*/) const override
    {
        assert(_rho_mol_co2_cumulated_prev.size() > 0);
        return _rho_mol_co2_cumulated_prev;
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
    std::vector<double> _co2_concentration;
    std::vector<double> _rho_mol_co2_cumulated_prev;
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
    const double& R = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double Q_steel = 7.8591 * 4 / 3;  // generate H2
    const double para_slow = 401.55;
    const double para_fast = 191.47286;

private:
    const double get_P_sat(double T)
    {
        double P_sat(0.0);
        double T_0 = 373.15;
        double P_0 = 101325.0;
        double h_wg = 2258000.0;
        return P_0 * exp(((1 / T_0) - (1 / T)) * Water * h_wg / R);
    }
    const double get_x_nonwet_air_gp(double PG, double X1, double X2, double X3,
                                     double P_sat, double kelvin_term)
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
                                                  double const p_sat,
                                                  double const kelvin_term)
    {
        return (1 / (kelvin_term / p_sat - 1 / Hen_L_air)) * (-1 / pg / pg);
    }

    const double get_derivative_x_nonwet_h2o_d_x1(double const pg,
                                                  double const /*x1*/,
                                                  double const /*x2*/,
                                                  double const /*x3*/,
                                                  double const p_sat,
                                                  double const kelvin_term)
    {
        return (1 / (kelvin_term / p_sat - 1 / Hen_L_air)) *
               (1 / Hen_L_air - 1 / Hen_L_h);
    }
    const double get_derivative_x_nonwet_h2o_d_x2(double const pg,
                                                  double const /*x1*/,
                                                  double const /*x2*/,
                                                  double const /*x3*/,
                                                  double const p_sat,
                                                  double const kelvin_term)
    {
        return (1 / (kelvin_term / p_sat - 1 / Hen_L_air)) *
               (1 / Hen_L_air - 1 / Hen_L_c);
    }
    const double get_derivative_x_nonwet_h2o_d_x3(double const pg,
                                                  double const /*x1*/,
                                                  double const /*x2*/,
                                                  double const /*x3*/,
                                                  double const p_sat,
                                                  double const kelvin_term)
    {
        return (1 / (kelvin_term / p_sat - 1 / Hen_L_air)) *
               (1 / Hen_L_air - 1 / Hen_L_co2);
    }

    const double get_derivative_x_nonwet_h2o_d_kelvin(double const pg,
                                                      double const x1,
                                                      double const x2,
                                                      double const x3,
                                                      double const p_sat,
                                                      double const kelvin_term)
    {
        double const A = 1 - x1 - x2 - x3;
        double const B =
            1 - (x1 * pg / Hen_L_h + x2 * pg / Hen_L_c + x3 * pg / Hen_L_co2);
        return -(B / pg - A / Hen_L_air) /
               (kelvin_term / p_sat - 1 / Hen_L_air) /
               (kelvin_term / p_sat - 1 / Hen_L_air) / p_sat;
    }

private:
    const double bi_interpolation(
        double pnt_x_to_interpolate, /*x-coordinate*/
        double pnt_y_to_interpolate, /*y-coordinate*/
        std::vector<double>
            _values_at_supp_pnts /*value for the matrix*/)
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
        0,    500,  1000, 1500, 2000, 2500, 3000, 3500,
        4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500};  // cumulated SiO2
    std::vector<double> _supp_pnts_y = {
        0,    500,  1000, 1500, 2000, 2500, 3000, 3500,
        4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500};  // CUMULATED co2
    std::vector<double> _porosity_at_supp_pnts_backfill = {
        0.095228012, 0.085187413, 0.076921316, 0.069146561, 0.054203101,
        0.054411434, 0.04968754,  0.070527096, 0.063517643, 0.061614807,
        0.055984546, 0.056036837, 0.056036849, 0.05603684,  0.056036839,
        0.056036849, 0.087055207, 0.078610167, 0.070827897, 0.056986405,
        0.056444216, 0.054371328, 0.07463621,  0.06720369,  0.065309411,
        0.059586561, 0.059586567, 0.059586559, 0.059586551, 0.059586567,
        0.059586558, 0.059586567, 0.080356345, 0.072570845, 0.059824535,
        0.058785554, 0.061545443, 0.073113988, 0.07089605,  0.066657709,
        0.063136279, 0.063136267, 0.063136285, 0.063136285, 0.063136285,
        0.063136285, 0.063136286, 0.063136273, 0.074398846, 0.062735462,
        0.061610492, 0.070823191, 0.076473714, 0.074594181, 0.068455597,
        0.066685994, 0.066686003, 0.066686003, 0.066686003, 0.066686002,
        0.066686004, 0.066685983, 0.066686004, 0.066686004, 0.065744135,
        0.065234123, 0.080084279, 0.080170044, 0.078297401, 0.071161581,
        0.070235714, 0.07023572,  0.07023572,  0.070235721, 0.070235721,
        0.070235721, 0.070235721, 0.070235721, 0.07023572,  0.070235721,
        0.070069991, 0.089331148, 0.083871435, 0.08200523,  0.073645281,
        0.073785432, 0.073785437, 0.073785437, 0.073785438, 0.073785438,
        0.073785438, 0.073785438, 0.073785438, 0.073785439, 0.073785439,
        0.073785439, 0.093619744, 0.087577372, 0.085598064, 0.07733515,
        0.077335149, 0.077335152, 0.077335153, 0.077335154, 0.077335155,
        0.077335155, 0.077335155, 0.077335155, 0.077335155, 0.077335155,
        0.077335155, 0.077335155, 0.09128756,  0.084320422, 0.080884867,
        0.080884863, 0.080884868, 0.080884869, 0.08088487,  0.080884871,
        0.080884871, 0.080884871, 0.080884871, 0.080884871, 0.080884871,
        0.080884872, 0.080884871, 0.080884872, 0.087022805, 0.084434588,
        0.084434577, 0.084434582, 0.084434585, 0.084434584, 0.084434585,
        0.084434585, 0.084434584, 0.084434586, 0.084434586, 0.084434586,
        0.084434586, 0.084434586, 0.084434586, 0.084434587, 0.087984275,
        0.08798429,  0.087984295, 0.087984297, 0.087984297, 0.087984299,
        0.087984298, 0.087984299, 0.0879843,   0.087984299, 0.087984299,
        0.0879843,   0.0879843,   0.0879843,   0.087984301, 0.0879843,
        0.091534002, 0.091534006, 0.091534008, 0.09153401,  0.091534011,
        0.091534022, 0.091534011, 0.091534012, 0.091534013, 0.091534012,
        0.091534012, 0.091534012, 0.091534013, 0.091534012, 0.091534013,
        0.091534012, 0.095083716, 0.095083717, 0.09508372,  0.09508372,
        0.095083722, 0.095083721, 0.095083723, 0.095083722, 0.095083724,
        0.095083722, 0.095083724, 0.095083724, 0.095083723, 0.095083725,
        0.095083725, 0.095083725, 0.098633427, 0.098633429, 0.098633428,
        0.098633431, 0.098633429, 0.09863343,  0.098633432, 0.098633433,
        0.098633433, 0.098633433, 0.098633431, 0.098633434, 0.098633432,
        0.098633432, 0.098633434, 0.098633432, 0.11064294,  0.11064294,
        0.11064294,  0.11064294,  0.11064294,  0.11064294,  0.11064294,
        0.11064294,  0.11064294,  0.11064294,  0.11064294,  0.11064294,
        0.11064294,  0.11064294,  0.11064294,  0.11064294,  0.12647343,
        0.12647342,  0.12647343,  0.12647343,  0.12647342,  0.12647343,
        0.12647342,  0.12647342,  0.12647343,  0.12647342,  0.12647342,
        0.12647342,  0.12647343,  0.12647342,  0.12647342,  0.12647342,
        0.13009791,  0.1300979,   0.13009791,  0.13009791,  0.1300979,
        0.13009791,  0.13009791,  0.13009791,  0.13009791,  0.13009791,
        0.13009791,  0.1300979,   0.13009791,  0.13009789,  0.1300979,
        0.13009791};

    std::vector<double> _pH_at_supp_pnt_backfill = {
        12.97891,  12.885851, 12.777895, 12.644894, 12.494799, 12.260243,
        11.466019, 11.140468, 10.961412, 10.962806, 10.401687, 10.107074,
        10.107072, 10.107073, 10.107073, 10.107072, 12.928335, 12.812059,
        12.665889, 12.513661, 12.207812, 11.27934,  11.05757,  10.96196,
        10.963307, 10.107074, 10.107072, 10.107073, 10.107074, 10.107072,
        10.107073, 10.107072, 12.849316, 12.687447, 12.533165, 12.137228,
        11.153893, 11.039851, 10.962468, 10.889863, 10.107072, 10.107073,
        10.107072, 10.107072, 10.107072, 10.107072, 10.107072, 10.107073,
        12.709468, 12.553313, 12.043468, 11.15155,  10.961694, 10.962941,
        10.724599, 10.107072, 10.107072, 10.107072, 10.107072, 10.107072,
        10.107072, 10.107073, 10.107073, 10.107072, 12.574097, 11.923318,
        11.149432, 10.962176, 10.963382, 10.481562, 10.107072, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107072, 10.107072, 10.107072,
        10.107073, 10.107072, 11.788069, 11.147505, 10.962626, 10.963792,
        10.016801, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        11.061197, 10.963046, 10.938225, 10.107072, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.963441, 10.87238,
        10.107072, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.55986,  10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 10.107073, 10.107073,
        10.107073, 10.107073, 10.107073, 10.107073, 8.3863666, 8.3863666,
        8.3863666, 8.3863666, 8.3863666, 8.3863666, 8.3863666, 8.3863666,
        8.3863666, 8.3863666, 8.3863666, 8.3863666, 8.3863666, 8.3863666,
        8.3863666, 8.3863666, 7.3567287, 7.3567285, 7.3567287, 7.3567288,
        7.3567286, 7.3567289, 7.3567286, 7.3567286, 7.3567289, 7.3567286,
        7.3567285, 7.3567285, 7.3567287, 7.3567285, 7.3567286, 7.3567285,
        6.3585847, 6.3585848, 6.3585847, 6.3585847, 6.358585,  6.3585847,
        6.3585847, 6.3585845, 6.3585847, 6.3585846, 6.3585847, 6.3585848,
        6.3585847, 6.3585855, 6.3585851, 6.3585847};
    std::vector<double> _quartz_rate_suppt_pnt_backfill = {
        -2.87E-07, -2.51E-07, -2.16E-07, -1.81E-07, -1.48E-07, -1.09E-07,
        -4.19E-08, -2.76E-08, -2.10E-08, -2.04E-08, -6.20E-09, 2.44E-15,
        2.98E-16,  -2.57E-15, -4.17E-16, 7.49E-16,  -2.71E-07, -2.31E-07,
        -1.90E-07, -1.55E-07, -1.06E-07, -3.46E-08, -2.55E-08, -2.17E-08,
        -2.11E-08, -4.32E-14, 5.53E-16,  -3.67E-15, 1.95E-15,  4.62E-16,
        -8.96E-16, -7.72E-17, -2.47E-07, -2.00E-07, -1.63E-07, -1.00E-07,
        -3.06E-08, -2.56E-08, -2.23E-08, -1.95E-08, 1.63E-14,  5.39E-15,
        3.64E-17,  1.26E-16,  -7.98E-17, -6.15E-17, -9.95E-17, -1.47E-17,
        -2.10E-07, -1.71E-07, -9.22E-08, -3.13E-08, -2.36E-08, -2.30E-08,
        -1.52E-08, 1.38E-14,  1.09E-16,  -9.80E-17, 7.35E-18,  3.67E-16,
        -1.15E-16, 9.86E-16,  -1.11E-15, -7.73E-17, -1.80E-07, -8.23E-08,
        -3.21E-08, -2.43E-08, -2.37E-08, -9.17E-09, 5.34E-15,  3.63E-16,
        7.04E-17,  1.96E-17,  4.27E-18,  -5.23E-17, -7.64E-17, -4.49E-17,
        5.32E-17,  1.49E-17,  -7.21E-08, -3.29E-08, -2.50E-08, -2.44E-08,
        2.43E-09,  3.14E-15,  4.79E-16,  1.17E-16,  7.29E-17,  8.95E-17,
        -5.19E-18, 1.48E-17,  9.62E-18,  -1.59E-17, -1.64E-17, -1.64E-17,
        -3.01E-08, -2.57E-08, -2.41E-08, 1.00E-14,  2.49E-15,  7.44E-16,
        3.21E-16,  1.68E-16,  1.10E-16,  7.18E-17,  5.97E-17,  4.09E-17,
        2.70E-17,  2.49E-17,  2.05E-17,  1.78E-17,  -2.63E-08, -2.24E-08,
        4.54E-15,  3.07E-15,  8.89E-16,  5.33E-16,  3.38E-16,  2.12E-16,
        1.83E-16,  1.46E-16,  1.08E-16,  1.09E-16,  8.85E-17,  7.04E-17,
        6.93E-17,  5.67E-17,  -1.27E-08, 1.42E-15,  3.23E-15,  1.37E-15,
        6.93E-16,  5.55E-16,  4.16E-16,  3.60E-16,  2.97E-16,  2.31E-16,
        2.01E-16,  1.77E-16,  1.58E-16,  1.43E-16,  1.41E-16,  1.02E-16,
        1.32E-14,  3.59E-15,  1.78E-15,  1.14E-15,  8.20E-16,  6.36E-16,
        5.64E-16,  4.35E-16,  3.75E-16,  3.59E-16,  3.19E-16,  2.88E-16,
        2.40E-16,  2.40E-16,  2.03E-16,  2.05E-16,  4.05E-15,  2.25E-15,
        1.52E-15,  1.14E-15,  9.04E-16,  2.20E-16,  6.93E-16,  5.54E-16,
        4.90E-16,  4.78E-16,  4.33E-16,  3.95E-16,  3.34E-16,  3.37E-16,
        2.88E-16,  2.93E-16,  2.79E-15,  2.13E-15,  1.51E-15,  1.32E-15,
        1.03E-15,  9.55E-16,  7.73E-16,  7.46E-16,  6.19E-16,  6.11E-16,
        5.16E-16,  4.77E-16,  4.80E-16,  4.13E-16,  3.87E-16,  3.64E-16,
        2.49E-15,  1.95E-15,  1.73E-15,  1.36E-15,  1.27E-15,  1.12E-15,
        9.27E-16,  8.39E-16,  7.65E-16,  7.03E-16,  7.03E-16,  6.06E-16,
        6.12E-16,  5.74E-16,  5.01E-16,  5.12E-16,  1.04E-16,  8.32E-17,
        7.20E-17,  6.34E-17,  5.67E-17,  5.12E-17,  4.43E-17,  4.30E-17,
        3.98E-17,  4.32E-17,  3.46E-17,  3.25E-17,  3.02E-17,  2.89E-17,
        2.74E-17,  2.61E-17,  2.03E-17,  2.12E-17,  1.59E-17,  1.19E-17,
        1.47E-17,  8.66E-18,  1.29E-17,  1.24E-17,  7.47E-18,  1.09E-17,
        1.09E-17,  1.01E-17,  7.97E-18,  9.91E-18,  8.39E-18,  8.47E-18,
        -2.00E-18, -3.33E-19, -1.65E-18, -1.26E-18, 9.03E-19,  -1.42E-18,
        -1.41E-18, -2.09E-18, -1.22E-18, -1.67E-18, -1.16E-18, -2.63E-19,
        -7.07E-19, 2.96E-18,  9.93E-19,  -8.19E-19};
    std::vector<double> _fluid_volume_suppt_pnt_backfill = {
        95.22765,  85.331506, 76.085842, 67.019768, 51.660815, 51.416037,
        44.355274, 66.601954, 60.069912, 59.584932, 52.81459,  52.015426,
        52.015434, 52.015436, 52.015431, 52.015432, 103.52836, 94.241774,
        85.206075, 71.115295, 69.923249, 64.706386, 86.817996, 80.468606,
        79.986833, 72.369165, 72.36916,  72.369162, 72.369148, 72.369159,
        72.369159, 72.369161, 112.42432, 103.43373, 90.604349, 88.681898,
        88.139755, 101.91107, 100.86965, 98.103461, 92.722879, 92.722873,
        92.722887, 92.722887, 92.722887, 92.722887, 92.722888, 92.722881,
        121.71971, 110.1397,  107.83343, 114.50672, 121.74913, 121.27286,
        116.1576,  113.0766,  113.07661, 113.07661, 113.07661, 113.07661,
        113.07661, 113.0766,  113.07662, 113.07661, 129.73776, 127.59087,
        140.87441, 142.15166, 141.67796, 134.80626, 133.43033, 133.43034,
        133.43034, 133.43034, 133.43034, 133.43034, 133.43034, 133.43034,
        133.43034, 133.43034, 148.05889, 167.24264, 162.55609, 162.08477,
        153.58609, 153.78406, 153.78406, 153.78406, 153.78406, 153.78406,
        153.78406, 153.78406, 153.78406, 153.78407, 153.78407, 153.78407,
        187.86055, 182.9622,  182.31373, 174.13779, 174.13779, 174.13779,
        174.13779, 174.13779, 174.13779, 174.13779, 174.13779, 174.13779,
        174.13779, 174.13779, 174.13779, 174.13779, 203.3699,  198.19923,
        194.49151, 194.49151, 194.49151, 194.49151, 194.49151, 194.49151,
        194.49151, 194.49151, 194.49152, 194.49152, 194.49152, 194.49152,
        194.49152, 194.49152, 217.05873, 214.84524, 214.84523, 214.84524,
        214.84524, 214.84524, 214.84524, 214.84524, 214.84524, 214.84524,
        214.84524, 214.84524, 214.84524, 214.84524, 214.84524, 214.84524,
        235.19895, 235.19896, 235.19896, 235.19896, 235.19896, 235.19896,
        235.19896, 235.19896, 235.19896, 235.19896, 235.19896, 235.19896,
        235.19896, 235.19896, 235.19896, 235.19896, 255.55268, 255.55268,
        255.55268, 255.55268, 255.55268, 255.55269, 255.55268, 255.55269,
        255.55269, 255.55269, 255.55269, 255.55269, 255.55269, 255.55269,
        255.55269, 255.55269, 275.9064,  275.9064,  275.90641, 275.90641,
        275.90641, 275.90641, 275.90641, 275.90641, 275.90641, 275.90641,
        275.90641, 275.90641, 275.90641, 275.90641, 275.90641, 275.90641,
        296.26013, 296.26013, 296.26013, 296.26013, 296.26013, 296.26013,
        296.26013, 296.26013, 296.26013, 296.26013, 296.26013, 296.26013,
        296.26013, 296.26013, 296.26013, 296.26013, 322.70038, 322.70038,
        322.70038, 322.70038, 322.70038, 322.70038, 322.70038, 322.70038,
        322.70038, 322.70038, 322.70038, 322.70038, 322.70038, 322.70038,
        322.70038, 322.70038, 346.23258, 346.23257, 346.23258, 346.23258,
        346.23258, 346.23258, 346.23257, 346.23257, 346.23258, 346.23257,
        346.23257, 346.23257, 346.23258, 346.23257, 346.23257, 346.23257,
        349.66438, 349.66437, 349.66437, 349.66437, 349.66437, 349.66438,
        349.66438, 349.66438, 349.66438, 349.66438, 349.66438, 349.66437,
        349.66437, 349.66437, 349.66437, 349.66438};
    std::vector<double> _saturation_index_suppt_pnts_backfill = {
        0.000101319, 0.000214774, 0.000490402, 0.001150398,
        0.002369789, 0.00416683,  0.018833242, 0.032584237,
        0.066964862, 0.066969149, 0.44112056,  1.0000003,
        1,           0.99999967,  0.99999995,  1.0000001,
        0.000188382, 0.000459403, 0.001158823, 0.002364584,
        0.004783349, 0.026147744, 0.04513495,  0.066966541,
        0.066970695, 0.99999472,  1.0000001,   0.99999955,
        1.0000002,   1.0000001,   0.99999989,  0.99999999,
        0.000425271, 0.001168962, 0.002358595, 0.005646063,
        0.032363714, 0.048986827, 0.066968094, 0.089372858,
        1.0000019,   1.0000006,   1,           1,
        0.99999999,  0.99999999,  0.99999999,  1,
        0.001181343, 0.002351677, 0.006899904, 0.032403405,
        0.06696573,  0.066969566, 0.16379294,  1.0000016,
        1,           0.99999999,  1,           1,
        0.99999999,  1.0000001,   0.99999987,  0.99999999,
        0.002343624, 0.008746205, 0.03243885,  0.066967201,
        0.066970926, 0.3530189,   1.0000006,   1,
        1,           1,           1,           0.99999999,
        0.99999999,  1,           1,           1,
        0.011193448, 0.032470718, 0.066968594, 0.066972199,
        1.2840608,   1.0000003,   1.0000001,   1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        0.045207964, 0.066969888, 0.074193434, 1.000001,
        1.0000003,   1.0000001,   1,           1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        0.066971109, 0.095278844, 1.0000004,   1.0000003,
        1.0000001,   1.0000001,   1,           1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        0.28364658,  1.0000001,   1.0000003,   1.0000001,
        1.0000001,   1.0000001,   1,           1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        1.0000012,   1.0000003,   1.0000002,   1.0000001,
        1.0000001,   1.0000001,   1.0000001,   1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        1.0000004,   1.0000002,   1.0000001,   1.0000001,
        1.0000001,   1,           1.0000001,   1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        1.0000002,   1.0000002,   1.0000001,   1.0000001,
        1.0000001,   1.0000001,   1.0000001,   1.0000001,
        1.0000001,   1.0000001,   1,           1,
        1,           1,           1,           1,
        1.0000002,   1.0000002,   1.0000001,   1.0000001,
        1.0000001,   1.0000001,   1.0000001,   1.0000001,
        1.0000001,   1.0000001,   1.0000001,   1.0000001,
        1.0000001,   1,           1,           1,
        1.0000001,   1,           1,           1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        1,           1,           1,           1,
        0.99999999,  1,           0.99999999,  0.99999999,
        1,           0.99999999,  0.99999999,  0.99999999,
        0.99999999,  0.99999999,  0.99999999,  1,
        1,           1,           1,           1};
    // look up table for waste matrix
    const double piecewiselinear_interpolation(
        double pnt_to_interpolate, std::vector<double> _values_at_supp_pnts)
    {
        // search interval that has the point inside
        if (pnt_to_interpolate <= _supp_pnts.front())
        {
            return _values_at_supp_pnts[0];
        }

        if (_supp_pnts.back() <= pnt_to_interpolate)
        {
            return _values_at_supp_pnts[_supp_pnts.size() - 1];
        }

        auto const& it(std::lower_bound(_supp_pnts.begin(), _supp_pnts.end(),
                                        pnt_to_interpolate));
        std::size_t const interval_idx =
            std::distance(_supp_pnts.begin(), it) - 1;

        // support points.
        double const x = _supp_pnts[interval_idx];
        double const x_r = _supp_pnts[interval_idx + 1];

        // values.
        double const f = _values_at_supp_pnts[interval_idx];
        double const f_r = _values_at_supp_pnts[interval_idx + 1];

        // compute linear interpolation polynom: y = m * (x - support[i]) +
        // value[i]
        const double m = (f_r - f) / (x_r - x);

        return m * (pnt_to_interpolate - x) + f;
    }
    std::vector<double> _supp_pnts = {50,  100, 150, 200, 250,
                                      300, 350, 400, 450, 500};
    std::vector<double> _pH_at_supp_pnt_waste = {
        10.065252, 10.07599,  10.080844, 10.085368, 10.089579,
        10.093479, 6.4492598, 6.2057855, 6.2057854, 6.2057854};
    std::vector<double> _porosity_change_at_supp_pnt_waste = {
        0,          0.00070226, 0.00105532, 0.00140958, 0.001765,
        0.00212163, 0.00335214, 0.00335031, 0.00335031, 0.00335031};
    std::vector<double> _quartz_rate_suppt_pnt_waste = {
        1.63E-18, 4.16E-18,  1.16E-17, 1.04E-17,  2.04E-18,
        7.30E-19, -4.20E-21, 2.86E-21, -1.16E-21, -1.21E-21};
    std::vector<double> _fluid_volume_suppt_pnt_waste = {
        9.0687704, 13.1384,  15.173906, 17.209846, 19.246203,
        21.283195, 23.61377, 23.629577, 23.629577, 23.629577};
};

}  // end of namespace
}  // end of namespace

#include "TwoPhaseComponentialFlowLocalAssembler-impl.h"