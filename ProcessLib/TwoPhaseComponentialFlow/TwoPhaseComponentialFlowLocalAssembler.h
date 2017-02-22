/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_TWOPHASECOMPONENTIALFLOWLOCALASSEMBLER_H
#define OGS_TWOPHASECOMPONENTIALFLOWLOCALASSEMBLER_H

#include <vector>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "MaterialLib/TwoPhaseModels/TwoPhaseFlowWithPPMaterialProperties.h"
#include "TwoPhaseComponentialFlowProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseComponentialFlow
{
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

    // Note: currently only isothermal case is considered, so the temperature is
    // assumed to be const
    // the variation of temperature will be taken into account in future
    double _temperature = 293.15;
    const double T_0 = 303.15;
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
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

    const double mu_L = 3.171e-11;
    const double mu_G = 2.8539e-13;  // viscosity
                                     // const double D_G = 31.536;
                                     // const double D_L= 9.47E-2;

    const double Q_steel = 7.8591;           // generate H2
    double Q_organic_fast_co2_ini = 0.035;   //
    double Q_organic_fast_ch4_ini = 0.035;   // 6.70155
    double Q_organic_slow_co2_ini = 0.0019;  //
    double Q_organic_slow_ch4_ini = 0.0019;  // 0.76294573

    const double para_slow = 401.55;
    const double para_fast = 191.47286;

    std::vector<double> _saturation;
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
    const double get_P_sat(double T)
    {
        // Here unit of T is Celsius;
        double P_sat(0.0);
        double T_0 = 373.15;
        double P_0 = 101325.0;
        double h_wg = 2258000.0;
        P_sat = P_0 * exp(((1 / T_0) - (1 / T)) * M_L * h_wg / R);
        return P_sat;
    }
    const double get_X_G_air_gp(double PG, double X1, double X2, double X3,
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
    const double bi_interpolation(double P_, double T_)
    {
        std::vector<double> p_supp_pnts = {2, 3, 4, 5, 6};
        std::vector<double> T_supp_pnts = {20, 30, 40, 50, 60};
        auto const& it_P(
            std::lower_bound(p_supp_pnts.begin(), p_supp_pnts.end(), P_));
        std::size_t const P_interval_idx =
            std::distance(p_supp_pnts.begin(), it_P) - 1;
        auto const& it_T(
            std::lower_bound(T_supp_pnts.begin(), T_supp_pnts.end(), T_));
        std::size_t const T_interval_idx =
            std::distance(T_supp_pnts.begin(), it_T) - 1;

        const double k[5][5] = {{1, 2, 3, 4, 5},
                                {6, 7, 8, 9, 10},
                                {11, 12, 13, 14, 15},
                                {16, 17, 18, 19, 20},
                                {21, 22, 23, 24, 25}};
        const double f_r1 =
            k[P_interval_idx][T_interval_idx] *
                (p_supp_pnts[P_interval_idx + 1] - P_) /
                (p_supp_pnts[P_interval_idx + 1] -
                 p_supp_pnts[P_interval_idx]) +
            k[P_interval_idx + 1][T_interval_idx] *
                (P_ - p_supp_pnts[P_interval_idx]) /
                (p_supp_pnts[P_interval_idx + 1] - p_supp_pnts[P_interval_idx]);
        const double f_r2 =
            k[P_interval_idx][T_interval_idx + 1] *
                (p_supp_pnts[P_interval_idx + 1] - P_) /
                (p_supp_pnts[P_interval_idx + 1] -
                 p_supp_pnts[P_interval_idx]) +
            k[P_interval_idx + 1][T_interval_idx + 1] *
                (P_ - p_supp_pnts[P_interval_idx]) /
                (p_supp_pnts[P_interval_idx + 1] - p_supp_pnts[P_interval_idx]);
        const double f_p =
            f_r1 * (T_supp_pnts[T_interval_idx + 1] - T_) /
                (T_supp_pnts[T_interval_idx + 1] -
                 T_supp_pnts[T_interval_idx]) +
            f_r2 * (T_ - T_supp_pnts[T_interval_idx]) /
                (T_supp_pnts[T_interval_idx + 1] - T_supp_pnts[T_interval_idx]);
        const double f_val_11 =
            k[P_interval_idx][T_interval_idx] *
            (p_supp_pnts[P_interval_idx + 1] - P_) *
            (T_supp_pnts[T_interval_idx + 1] - T_) /
            (p_supp_pnts[P_interval_idx + 1] - p_supp_pnts[P_interval_idx]) /
            (T_supp_pnts[T_interval_idx + 1] - T_supp_pnts[T_interval_idx]);
        const double f_val_21 =
            k[P_interval_idx + 1][T_interval_idx] *
            (P_ - p_supp_pnts[P_interval_idx]) *
            (T_supp_pnts[T_interval_idx + 1] - T_) /
            (p_supp_pnts[P_interval_idx + 1] - p_supp_pnts[P_interval_idx]) /
            (T_supp_pnts[T_interval_idx + 1] - T_supp_pnts[T_interval_idx]);
        const double f_val_12 =
            k[P_interval_idx][T_interval_idx + 1] *
            (p_supp_pnts[P_interval_idx + 1] - P_) *
            (T_ - T_supp_pnts[T_interval_idx]) /
            (p_supp_pnts[P_interval_idx + 1] - p_supp_pnts[P_interval_idx]) /
            (T_supp_pnts[T_interval_idx + 1] - T_supp_pnts[T_interval_idx]);
        const double f_val_22 =
            k[P_interval_idx + 1][T_interval_idx + 1] *
            (P_ - p_supp_pnts[P_interval_idx]) *
            (T_ - T_supp_pnts[T_interval_idx]) /
            (p_supp_pnts[P_interval_idx + 1] - p_supp_pnts[P_interval_idx]) /
            (T_supp_pnts[T_interval_idx + 1] - T_supp_pnts[T_interval_idx]);

        return f_val_11 + f_val_21 + f_val_12 + f_val_22;
    }
};

}  // end of namespace
}  // end of namespace

#include "TwoPhaseComponentialFlowLocalAssembler-impl.h"

#endif /* TWOPHASECOMPONENTIALFLOWLOCALASSEMBLER_H */
