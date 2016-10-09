/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_TWOPHASECOMPONENTIAL_FEM_H_
#define PROCESS_LIB_TWOPHASECOMPONENTIAL_FEM_H_

#include <array>
#include <iostream>
#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MeshLib/CoordinateSystem.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "TwoPhaseComponentialMaterialModel.h"
#include "TwoPhaseComponentialProcessData.h"
namespace ProcessLib
{
namespace TwoPhaseComponential
{
const unsigned NUM_NODAL_DOF = 5;
template <typename ShapeMatrixType>
struct SecondaryData
{
    // std::vector<ShapeMatrixType> N;
    std::vector<double> _saturation;
    std::vector<double> _liquid_pressure;  // =
    // std::vector<double>(ShapeFunction::NPOINTS);
};

class TwoPhaseComponentialLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSat(
        std::vector<double>& /*cache*/) const = 0;
    virtual std::vector<double> const& getIntPtPL(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public TwoPhaseComponentialLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    /// The thermal_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       TwoPhaseComponentialProcessData const& process_data)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method))
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
        const MeshLib::CoordinateSystem coordsystem(element);
        // const MeshLib::ElementCoordinatesMappingLocal
        // ele_local_coord(element, coordsystem);// wp.getCoords());
        dim = coordsystem.getDimension();
        _secondary_data._saturation.resize(
            _integration_method.getNumberOfPoints());
        _secondary_data._liquid_pressure.resize(
            _integration_method.getNumberOfPoints());
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = local_x.size();
        // auto const n_nodes = _element.getNodes();
        auto const n_nodes = ShapeFunction::NPOINTS;
        const int n_dof = local_x.size();
        const int N_size = n_dof / n_nodes;
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        Eigen::MatrixXd mass_mat_coeff =
            Eigen::MatrixXd::Zero(NUM_NODAL_DOF, NUM_NODAL_DOF);
        Eigen::MatrixXd K_mat_coeff = Eigen::MatrixXd::Zero(N_size, N_size);
        Eigen::VectorXd H_vec_coeff = Eigen::VectorXd::Zero(N_size);
        Eigen::VectorXd F_vec_coeff = Eigen::VectorXd::Zero(N_size);
        Eigen::MatrixXd localMass_tmp = Eigen::MatrixXd::Zero(
            ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

        Eigen::MatrixXd localDispersion_tmp = Eigen::MatrixXd::Zero(
            ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
        Eigen::VectorXd localGravity_tmp =
            Eigen::VectorXd::Zero(ShapeFunction::NPOINTS);
        Eigen::VectorXd localSource_tmp =
            Eigen::VectorXd::Zero(ShapeFunction::NPOINTS);

        auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);  // mass matrix
        auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_K_data, local_matrix_size,
            local_matrix_size);  // laplace matrix
        auto local_b = MathLib::createZeroedVector<NodalVectorType>(
            local_b_data, local_matrix_size);
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition pos;
        size_t ele_id = _element.getID();
        pos.setElementID(_element.getID());
        /*primary variable value on each gauss point*/
        double PG_int_pt = 0.0;
        double X1_int_pt = 0.0;
        double X2_int_pt = 0.0;
        double X3_int_pt = 0.0;
        double PC_int_pt = 0.0;
        /*paramters*/
        double P_sat_gp = 0.0;  // vapor pressure
        double X_L_h_gp = 0.0;
        double X_L_c_gp = 0.0;
        double X_L_co2_gp = 0.0;
        double X_G_air_gp = 0.0;
        double X_G_h2o_gp = 0.0;
        double dXgairdPG = 0.0;
        double dXgh2odPG = 0.0;
        double dxgairdX1 = 0.0;
        double dxgairdX2 = 0.0;
        double dxgairdX3 = 0.0;
        double dxgh2odX1 = 0.0;
        double dxgh2odX2 = 0.0;
        double dxgh2odX3 = 0.0;
        double rho_mol_G_gp = 0.0;
        double rho_mol_L_gp = 0.0;
        double rho_mass_G_gp = 0.0;
        double rho_mass_L_gp = 0.0;
        double Kr_L_gp = 0.0;
        double Kr_G_gp = 0.0;
        double lambda_L = 0.0;
        double lambda_G = 0.0;
        /*secondary variables*/
        double S_G_gp = 0.0;
        double PL_gp = 0.0;
        double dSgdPC = 0.0;
        double Q_organic_slow_co2 = 0.0;
        double Q_organic_fast_co2 = 0.0;
        _secondary_data._saturation.resize(n_integration_points);
        _secondary_data._liquid_pressure.resize(n_integration_points);

        MathLib::PiecewiseLinearInterpolation const& interQ_slow =
            *_process_data.curves.at("curveA");
        MathLib::PiecewiseLinearInterpolation const& interQ_fast =
            *_process_data.curves.at("curveB");
        /*Loop for each gauss point*/
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            if (ele_id == 5758)
                double test = 1e-15;
            F_vec_coeff.setZero();
            pos.setIntegrationPoint(ip);
            auto const& sm = _shape_matrices[ip];
            auto const& wp = _integration_method.getWeightedPoint(ip);
            auto const K = _process_data.intrinsic_permeability(
                t, pos)[0];  // permeability
            auto const mat_id = _process_data.mat_id(t, pos)[0];
            auto const poro = _process_data.porosity(t, pos)[0];  // porosity
            auto const D_G = _process_data.viscosity_gas(t, pos)[0];
            auto const D_L = _process_data.viscosity_liquid(t, pos)[0];
            auto const rho_l_std =
                _process_data.water_density(t, pos)[0];  // water density

            NumLib::shapeFunctionInterpolate(local_x, sm.N, PG_int_pt,
                                             X1_int_pt, X2_int_pt, X3_int_pt,
                                             PC_int_pt);
            X_L_h_gp = PG_int_pt * X1_int_pt / Hen_L_h;  // Henry law
            X_L_c_gp = PG_int_pt * X2_int_pt / Hen_L_c;
            X_L_co2_gp = PG_int_pt * X3_int_pt / Hen_L_co2;

            P_sat_gp = get_P_sat(T_0);
            double K_G_w = PG_int_pt / P_sat_gp;
            double K_G_air = PG_int_pt / Hen_L_air;
            double L = 1 - (X_L_h_gp + X_L_c_gp + X_L_co2_gp);
            double G = 1 - X1_int_pt - X2_int_pt - X3_int_pt;
            X_G_air_gp = get_X_G_air_gp(PG_int_pt, X1_int_pt, X2_int_pt,
                                        X3_int_pt, P_sat_gp);
            X_G_h2o_gp = get_X_G_h2o_gp(PG_int_pt, X1_int_pt, X2_int_pt,
                                        X3_int_pt, P_sat_gp);
            // double X_L_h2o_gp = X_G_h2o_gp *K_G_w;
            // double dXgairdPG_test = (_EOS->get_X_G_air_gp(PG_int_pt  + eps,
            // X1_int_pt, X2_int_pt, X3_int_pt , P_sat_gp ) - X_G_air_gp ) /
            // eps;
            double dLdPG = -X1_int_pt / Hen_L_h - X2_int_pt / Hen_L_c -
                           X3_int_pt / Hen_L_co2;
            dXgairdPG = ((G / P_sat_gp - dLdPG) * (K_G_w - K_G_air) -
                         (K_G_w * G - L) * (1 / P_sat_gp - 1 / Hen_L_air)) /
                        (K_G_w - K_G_air) / (K_G_w - K_G_air);
            dXgh2odPG = -((G / Hen_L_air - dLdPG) * (K_G_w - K_G_air) +
                          (L - K_G_air * G) * (1 / P_sat_gp - 1 / Hen_L_air)) /
                        (K_G_w - K_G_air) / (K_G_w - K_G_air);
            dxgairdX1 = ((-K_G_w + PG_int_pt / Hen_L_h) / (K_G_w - K_G_air));
            dxgairdX2 = ((-K_G_w + PG_int_pt / Hen_L_c) / (K_G_w - K_G_air));
            dxgairdX3 = ((-K_G_w + PG_int_pt / Hen_L_co2) / (K_G_w - K_G_air));
            dxgh2odX1 = ((K_G_air - PG_int_pt / Hen_L_h) / (K_G_w - K_G_air));
            dxgh2odX2 = ((K_G_air - PG_int_pt / Hen_L_c) / (K_G_w - K_G_air));
            dxgh2odX3 = ((K_G_air - PG_int_pt / Hen_L_co2) / (K_G_w - K_G_air));

            S_G_gp = getSat_byPC(PC_int_pt, 4);
            _secondary_data._saturation[ip] =
                S_G_gp;  // store the secondary variable
            PL_gp = PG_int_pt - PC_int_pt;
            _secondary_data._liquid_pressure[ip] =
                PL_gp;  // store the secondary variable
            dSgdPC = get_diriv_PC(PC_int_pt, 4);  // dSgdPC

            rho_mol_G_gp = PG_int_pt / R / T_0;
            rho_mol_L_gp = rho_l_std / M_L;
            rho_mass_G_gp =
                rho_mol_G_gp *
                (X1_int_pt * M_H + X2_int_pt * M_C + X3_int_pt * M_CO2 +
                 X_G_air_gp * M_AIR + X_G_h2o_gp * M_L);
            rho_mass_L_gp = rho_l_std;

            Q_organic_slow_co2_ini =
                interQ_slow.getValue(t);  // read from curves
            Q_organic_fast_co2_ini =
                interQ_fast.getValue(t);  // read from curves

            /*calculate each entry of mass matrix*/
            mass_mat_coeff(0, 0) = poro * (S_G_gp * X1_int_pt / R / T_0 +
                                           (1 - S_G_gp) * rho_mol_L_gp *
                                               X1_int_pt / Hen_L_h);  // dPG
            mass_mat_coeff(0, 1) = poro * (rho_mol_G_gp * S_G_gp +
                                           (1 - S_G_gp) * rho_mol_L_gp *
                                               PG_int_pt / Hen_L_h);  // dX1 h2
            mass_mat_coeff(0, 2) = 0.0;                               // Dx2 ch4
            mass_mat_coeff(0, 3) = 0.0;                               // dx3 CO2
            mass_mat_coeff(0, 4) =
                poro * (rho_mol_G_gp * X1_int_pt * dSgdPC -
                        rho_mol_L_gp * X_L_h_gp * dSgdPC);  // dPC
            // CH4
            mass_mat_coeff(1, 0) = poro * (S_G_gp * X2_int_pt / R / T_0 +
                                           (1 - S_G_gp) * rho_mol_L_gp *
                                               X2_int_pt / Hen_L_c);  // dPG
            mass_mat_coeff(1, 1) = 0.0;                               // dX1 h2
            mass_mat_coeff(1, 2) = poro * (rho_mol_G_gp * S_G_gp +
                                           (1 - S_G_gp) * rho_mol_L_gp *
                                               PG_int_pt / Hen_L_c);  // Dx2 ch4
            mass_mat_coeff(1, 3) = 0.0;                               // dx3 CO2
            mass_mat_coeff(1, 4) =
                poro * (rho_mol_G_gp * X2_int_pt * dSgdPC -
                        rho_mol_L_gp * X_L_c_gp * dSgdPC);  // dPC
            // co2
            mass_mat_coeff(2, 0) = poro * (S_G_gp * X3_int_pt / R / T_0 +
                                           (1 - S_G_gp) * rho_mol_L_gp *
                                               X3_int_pt / Hen_L_co2);  // dPG
            mass_mat_coeff(2, 1) = 0.0;  // dX1 h2
            mass_mat_coeff(2, 2) = 0.0;  // Dx2 ch4
            mass_mat_coeff(2, 3) =
                poro * (rho_mol_G_gp * S_G_gp +
                        (1 - S_G_gp) * rho_mol_L_gp * PG_int_pt /
                            Hen_L_co2);  // dx3 CO2
            mass_mat_coeff(2, 4) =
                poro * (rho_mol_G_gp * X3_int_pt * dSgdPC -
                        rho_mol_L_gp * X_L_co2_gp * dSgdPC);  // dPC
            // air
            mass_mat_coeff(3, 0) =
                poro *
                (S_G_gp * X_G_air_gp / R / T_0 +
                 S_G_gp * rho_mol_G_gp * dXgairdPG +
                 (1 - S_G_gp) * rho_mol_L_gp * dXgairdPG * K_G_air +
                 (1 - S_G_gp) * rho_mol_L_gp * X_G_air_gp / Hen_L_air);  // dPG
            mass_mat_coeff(3, 1) = poro * (S_G_gp * rho_mol_G_gp * dxgairdX1 +
                                           (1 - S_G_gp) * rho_mol_L_gp *
                                               K_G_air * dxgairdX1);  // dX1 h2

            mass_mat_coeff(3, 2) = poro * (S_G_gp * rho_mol_G_gp * dxgairdX2 +
                                           (1 - S_G_gp) * rho_mol_L_gp *
                                               K_G_air * dxgairdX2);  // dX2 CH4

            mass_mat_coeff(3, 3) = poro * (S_G_gp * rho_mol_G_gp * dxgairdX3 +
                                           (1 - S_G_gp) * rho_mol_L_gp *
                                               K_G_air * dxgairdX3);  // dX3 co2

            mass_mat_coeff(3, 4) =
                poro * (rho_mol_G_gp * X_G_air_gp * dSgdPC -
                        rho_mol_L_gp * K_G_air * X_G_air_gp * dSgdPC);  // dPC
            // h20
            mass_mat_coeff(4, 0) =
                poro *
                (S_G_gp * X_G_h2o_gp / R / T_0 +
                 S_G_gp * rho_mol_G_gp * dXgh2odPG +
                 (1 - S_G_gp) * rho_mol_L_gp * dXgh2odPG * K_G_w +
                 (1 - S_G_gp) * rho_mol_L_gp * X_G_h2o_gp / P_sat_gp);  // dPG
            mass_mat_coeff(4, 1) = poro * (S_G_gp * rho_mol_G_gp * dxgh2odX1 +
                                           (1 - S_G_gp) * rho_mol_L_gp * K_G_w *
                                               dxgh2odX1);  // dX1 h2

            mass_mat_coeff(4, 2) = poro * (S_G_gp * rho_mol_G_gp * dxgh2odX2 +
                                           (1 - S_G_gp) * rho_mol_L_gp * K_G_w *
                                               dxgh2odX2);  // dX2 CH4

            mass_mat_coeff(4, 3) = poro * (S_G_gp * rho_mol_G_gp * dxgh2odX3 +
                                           (1 - S_G_gp) * rho_mol_L_gp * K_G_w *
                                               dxgh2odX3);  // dX3 co2

            mass_mat_coeff(4, 4) =
                poro * (rho_mol_G_gp * X_G_h2o_gp * dSgdPC -
                        rho_mol_L_gp * K_G_w * X_G_h2o_gp * dSgdPC);  // dPC
            //-------------debugging------------------------
            // std::cout << "mass_mat_coeff=" << std::endl;
            // std::cout << mass_mat_coeff << std::endl;
            //--------------end debugging-------------------
            for (int ii = 0; ii < NUM_NODAL_DOF; ii++)
            {
                for (int jj = 0; jj < NUM_NODAL_DOF; jj++)
                {
                    localMass_tmp.setZero();
                    localMass_tmp.noalias() = sm.N.transpose() *
                                              mass_mat_coeff(ii, jj) * sm.N *
                                              sm.detJ * wp.getWeight();
                    local_M.block(n_nodes * ii, n_nodes * jj, n_nodes, n_nodes)
                        .noalias() += localMass_tmp;
                }
            }
            // std::cout << local_M << std::endl;
            Kr_L_gp = getKr_L_bySg(S_G_gp, 4);  // change the Relative
                                                // permeability calculation into
                                                // porous media file
            Kr_G_gp = getKr_g_bySg(S_G_gp, 4);
            lambda_L = K * Kr_L_gp / mu_L;
            lambda_G = K * Kr_G_gp / mu_G;

            // Calc each entry of the Laplace Matrix
            // h2
            K_mat_coeff(0, 0) =
                lambda_G * rho_mol_G_gp * X1_int_pt +
                lambda_L * rho_mol_L_gp * X_L_h_gp +
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * X1_int_pt / Hen_L_h;
            K_mat_coeff(0, 1) =
                poro * D_G * S_G_gp * rho_mol_G_gp +
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * PG_int_pt / Hen_L_h;
            K_mat_coeff(0, 2) = 0.0;
            K_mat_coeff(0, 3) = 0.0;
            K_mat_coeff(0, 4) = -lambda_L * rho_mol_L_gp * X_L_h_gp;
            // ch4
            K_mat_coeff(1, 0) =
                lambda_G * rho_mol_G_gp * X2_int_pt +
                lambda_L * rho_mol_L_gp * X_L_c_gp +
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * X2_int_pt / Hen_L_c;
            K_mat_coeff(1, 1) = 0.0;
            K_mat_coeff(1, 2) =
                poro * D_G * S_G_gp * rho_mol_G_gp +
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * PG_int_pt / Hen_L_c;
            K_mat_coeff(1, 3) = 0.0;
            K_mat_coeff(1, 4) = -lambda_L * rho_mol_L_gp * X_L_c_gp;
            // co2
            K_mat_coeff(2, 0) = lambda_G * rho_mol_G_gp * X3_int_pt +
                                lambda_L * rho_mol_L_gp * X_L_co2_gp +
                                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp *
                                    X3_int_pt / Hen_L_co2;
            K_mat_coeff(2, 1) = 0.0;
            K_mat_coeff(2, 2) = 0.0;
            K_mat_coeff(2, 3) = poro * D_G * S_G_gp * rho_mol_G_gp +
                                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp *
                                    PG_int_pt / Hen_L_co2;
            K_mat_coeff(2, 4) = -lambda_L * rho_mol_L_gp * X_L_co2_gp;
            // air
            K_mat_coeff(3, 0) =
                lambda_G * rho_mol_G_gp * X_G_air_gp +
                lambda_L * rho_mol_L_gp * K_G_air * X_G_air_gp +
                poro * D_G * S_G_gp * rho_mol_G_gp * dXgairdPG +
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * dXgairdPG * K_G_air;
            K_mat_coeff(3, 1) =
                poro * D_G * S_G_gp * rho_mol_G_gp *
                    ((-K_G_w + PG_int_pt / Hen_L_h) / (K_G_w - K_G_air)) +
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * dxgairdX1 * K_G_air;
            K_mat_coeff(3, 2) =
                poro * D_G * S_G_gp * rho_mol_G_gp *
                    ((-K_G_w + PG_int_pt / Hen_L_c) / (K_G_w - K_G_air)) +
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * dxgairdX2 * K_G_air;
            K_mat_coeff(3, 3) =
                poro * D_G * S_G_gp * rho_mol_G_gp *
                    ((-K_G_w + PG_int_pt / Hen_L_co2) / (K_G_w - K_G_air)) +
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * dxgairdX3 * K_G_air;
            K_mat_coeff(3, 4) = -lambda_L * rho_mol_L_gp * K_G_air * X_G_air_gp;
            // h2o
            K_mat_coeff(4, 0) =
                lambda_G * rho_mol_G_gp * X_G_h2o_gp +
                lambda_L * rho_mol_L_gp * K_G_w * X_G_h2o_gp -
                poro * D_G * S_G_gp * rho_mol_G_gp * dXgairdPG -
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp *
                    (X1_int_pt / Hen_L_h + X2_int_pt / Hen_L_c +
                     X3_int_pt / Hen_L_co2) -
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * dXgairdPG * K_G_air;
            K_mat_coeff(4, 1) =
                -poro * D_G * S_G_gp * rho_mol_G_gp -
                poro * D_G * S_G_gp * rho_mol_G_gp * dxgairdX1 -
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * PG_int_pt / Hen_L_h -
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * dxgairdX1 * K_G_air;
            K_mat_coeff(4, 2) =
                -poro * D_G * S_G_gp * rho_mol_G_gp -
                poro * D_G * S_G_gp * rho_mol_G_gp * dxgairdX2 -
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * PG_int_pt / Hen_L_c -
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * dxgairdX2 * K_G_air;
            K_mat_coeff(4, 3) =
                -poro * D_G * S_G_gp * rho_mol_G_gp -
                poro * D_G * S_G_gp * rho_mol_G_gp * dxgairdX3 -
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * PG_int_pt /
                    Hen_L_co2 -
                poro * D_L * (1 - S_G_gp) * rho_mol_L_gp * dxgairdX3 * K_G_air;
            K_mat_coeff(4, 4) = -lambda_L * rho_mol_L_gp * K_G_w * X_G_h2o_gp;

            //-------------debugging------------------------
            // std::cout << "K_mat_coeff=" << std::endl;
            // std::cout << K_mat_coeff << std::endl;
            //--------------end debugging-------------------

            for (int ii = 0; ii < NUM_NODAL_DOF; ii++)
            {
                for (int jj = 0; jj < NUM_NODAL_DOF; jj++)
                {
                    localDispersion_tmp.setZero();
                    localDispersion_tmp.noalias() =
                        sm.dNdx.transpose() * K_mat_coeff(ii, jj) * sm.dNdx *
                        sm.detJ * wp.getWeight();
                    local_K.block(n_nodes * ii, n_nodes * jj, n_nodes, n_nodes)
                        .noalias() += localDispersion_tmp;
                }
            }
            // std::cout << local_K << std::endl;
            /*Here is the gravity term*/
            H_vec_coeff(0) =
                -lambda_G * rho_mol_G_gp * X1_int_pt * rho_mass_G_gp -
                lambda_L * rho_mol_L_gp * X1_int_pt * PG_int_pt *
                    rho_mass_L_gp / Hen_L_h;
            H_vec_coeff(1) =
                -lambda_G * rho_mol_G_gp * X2_int_pt * rho_mass_G_gp -
                lambda_L * rho_mol_L_gp * X2_int_pt * PG_int_pt *
                    rho_mass_L_gp / Hen_L_c;
            H_vec_coeff(2) =
                -lambda_G * rho_mol_G_gp * X3_int_pt * rho_mass_G_gp -
                lambda_L * rho_mol_L_gp * X3_int_pt * PG_int_pt *
                    rho_mass_L_gp / Hen_L_co2;
            H_vec_coeff(3) =
                -lambda_G * rho_mol_G_gp * X_G_air_gp * rho_mass_G_gp -
                lambda_L * rho_mol_L_gp * X_G_air_gp * PG_int_pt *
                    rho_mass_L_gp / Hen_L_air;
            H_vec_coeff(4) =
                -lambda_G * rho_mol_G_gp * X_G_h2o_gp * rho_mass_G_gp -
                lambda_L * rho_mol_L_gp * X_G_h2o_gp * PG_int_pt *
                    rho_mass_L_gp / P_sat_gp;

            //-------------debugging------------------------
            // std::cout << "H_vec_coeff=" << std::endl;
            // std::cout << H_vec_coeff << std::endl;
            //--------------end debugging-------------------

            if (_process_data.has_gravity)
            {
                for (int idx = 0; idx < NUM_NODAL_DOF; idx++)
                {
                    typename ShapeMatricesType::GlobalDimVectorType vec_g;
                    vec_g.resize(dim);
                    vec_g(dim - 1) = 9.81;
                    // since no primary vairable involved
                    // directly assemble to the Right-Hand-Side
                    // F += dNp^T * K * gz
                    localGravity_tmp.setZero();
                    localGravity_tmp.noalias() = sm.dNdx.transpose() *
                                                 H_vec_coeff(idx) * vec_g *
                                                 sm.detJ * wp.getWeight();
                    local_b.block(n_nodes * idx, 0, n_nodes, 1).noalias() +=
                        localGravity_tmp;
                }
            }  // end of if hasGravityEffect
            // std::cout << t << "  " << local_b << "\n";

            // load the source term
            if (mat_id == 1)
            {
                F_vec_coeff(0) += Q_steel;
                F_vec_coeff(4) -= Q_steel;
                // curve_slow->eval(time.getTime(), Q_organic_slow_co2);
                Q_organic_slow_co2 = Q_organic_slow_ch4_ini * para_slow;
                // curve_fast->eval(time.getTime(), Q_organic_fast_co2);
                Q_organic_fast_co2 = Q_organic_fast_co2_ini * para_fast;
                F_vec_coeff(1) += Q_organic_slow_co2 * 31 /
                                  19;  // CH4 organic degradation slow
                F_vec_coeff(2) +=
                    Q_organic_slow_co2;  // co2 organic degradation slow
                F_vec_coeff(4) -= 2 * Q_organic_slow_co2;  // h2o consume slow
                F_vec_coeff(1) +=
                    Q_organic_fast_co2;  // co2 organic degradation fast
                F_vec_coeff(2) +=
                    Q_organic_fast_co2;  // ch4 organic degradation fast
                F_vec_coeff(4) -= Q_organic_fast_co2 / 3;  // h2o consume fast
            }
            else if (mat_id == 0)
            {
                // curve_slow->eval(time.getTime(), Q_organic_slow_co2);
                Q_organic_slow_co2 = Q_organic_slow_co2_ini * para_slow;
                // curve_fast->eval(time.getTime(), Q_organic_fast_co2);
                Q_organic_fast_co2 = Q_organic_fast_co2_ini * para_fast;
                if (X3_int_pt > 0.01)
                {
                    F_vec_coeff(2) -=
                        Q_organic_slow_co2;  // carbonation consume co2
                    F_vec_coeff(4) += Q_organic_slow_co2;
                }

                F_vec_coeff(4) -= 2.57635;  // concrete degradation
            }
            // std::cout << F_vec_coeff << std::endl;
            for (int idx = 0; idx < 5; idx++)
            {
                // tmp(0, 0) = F_vec_coeff(idx);
                localSource_tmp.setZero();
                // fe->integrateDWxvec_g(j, tmp, localGravity_tmp, vec_g);
                localSource_tmp = sm.N.transpose() * F_vec_coeff(idx) *
                                  sm.detJ * wp.getWeight();
                local_b.block(n_nodes * idx, 0, n_nodes, 1).noalias() +=
                    localSource_tmp;
            }
            // std::cout << t << "  " << local_b << "\n";
        }
        if (true)
        {
            for (int idx_ml = 0; idx_ml < local_M.cols(); idx_ml++)
            {
                double mass_lump_val;
                mass_lump_val = local_M.col(idx_ml).sum();
                local_M.col(idx_ml).setZero();
                local_M(idx_ml, idx_ml) = mass_lump_val;
            }
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSat(
        std::vector<double>& /*cache*/) const override
    {
        return _secondary_data._saturation;
    }

    std::vector<double> const& getIntPtPL(
        std::vector<double>& /*cache*/) const override
    {
        return _secondary_data._liquid_pressure;
    }

private:
    MeshLib::Element const& _element;
    unsigned dim;
    TwoPhaseComponentialProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    // std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
    // _shape_matrices;
    std::vector<ShapeMatrices> _shape_matrices;
    SecondaryData<typename ShapeMatrices::ShapeType>
        _secondary_data;  // secondary variables

    NodalMatrixType _localA;
    NodalMatrixType _localM;
    NodalVectorType _localRhs;

    // std::vector<std::vector<double>> _saturation
    //= std::vector<std::vector<double>>(1,
    // std::vector<double>(ShapeFunction::NPOINTS));

    /*std::vector<std::vector<double>> _liquid_pressure
        = std::vector<std::vector<double>>(1,
            std::vector<double>(ShapeFunction::NPOINTS));*/

    const double T_0 = 303.15;
    const double Hen_L_h = 7.26e+9;     // Henry constant in [Pa]
    const double Hen_L_c = 4.13e+9;     // Henry constant in [Pa]
    const double Hen_L_air = 9.077e+9;  // Henry constant in [Pa]
    const double Hen_L_co2 = 0.163e+9;  // Henry constant in [Pa]
    const double rho_l_std = 1000.0;

    const double M_H = 0.002;
    const double M_L = 0.02;
    const double M_C = 0.016;
    const double M_AIR = 0.029;
    const double M_CO2 = 0.044;

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
};

}  // namespace TwoPhaseComponential
}  // namespace ProcessLib

#endif  // PROCESS_LIB_TwoPhaseComponential_FEM_H_
