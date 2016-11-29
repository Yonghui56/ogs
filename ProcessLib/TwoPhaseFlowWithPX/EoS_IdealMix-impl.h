/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_TWOPHASEFLOWWITHPX_EOS_IDEALMIX_IMPL_H_
#define OGS_TWOPHASEFLOWWITHPX_EOS_IDEALMIX_IMPL_H_

#include "NewtonRaphson.h"
#include "EoS_IdealMix.h"
namespace ProcessLib
{
namespace TwoPhaseFlowWithPX
{
bool EoS_IdealMix::computeConstitutiveRelation(
    double const t,
    ProcessLib::SpatialPosition const& x,
	double const PG,
	double const X,
	double& Sw,
	double& X_m,
	double& dsw_dpg,
	double& dsw_dX,
	double& dxm_dpg,
	double& dxm_dX)
{
    using LocalJacobianMatrix =
        Eigen::Matrix<double, 2, 2,
                      Eigen::RowMajor>;
	LocalJacobianMatrix J_loc;
    // Linear solver for the newton loop is required after the loop with the
    // same matrix. This saves one decomposition.
    Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(2);

    // Different solvers are available for the solution of the local system.
    // TODO Make the following choice of linear solvers available from the
    // input file configuration:
    //      K_loc.partialPivLu().solve(-res_loc);
    //      K_loc.fullPivLu().solve(-res_loc);
    //      K_loc.householderQr().solve(-res_loc);
    //      K_loc.colPivHouseholderQr().solve(res_loc);
    //      K_loc.fullPivHouseholderQr().solve(-res_loc);
    //      K_loc.llt().solve(-res_loc);
    //      K_loc.ldlt().solve(-res_loc);

    {  // Local Newton solver
        using LocalResidualVector =
            Eigen::Matrix<double, 2, 1>;
		using LocalUnknownVector = Eigen::Matrix<double, 2, 1>;
        LocalJacobianMatrix J_loc;
		//LocalUnknownVector sec_var_unknown;
		//sec_var_unknown(0) = Sw;
		//sec_var_unknown(1) = X_m;
        auto const update_residual = [&](LocalResidualVector& residual) {
            calculateResidual(PG, X, Sw, X_m, residual);
        };

        auto const update_jacobian = [&](LocalJacobianMatrix& jacobian) {
            calculateJacobian(
                t, x, PG,X,jacobian, Sw, X_m);  // for solution dependent Jacobians
        };

        auto const update_solution = [&](LocalResidualVector const& increment) {
            // increment solution vectors
			Sw += increment[0];
			X_m += increment[1];
        };

        // TODO Make the following choice of maximum iterations and convergence
        // criteria available from the input file configuration:
        const int maximum_iterations(20);
        const double tolerance(1.e-15);

        auto newton_solver =
            NewtonRaphson<decltype(linear_solver), LocalJacobianMatrix,
                          decltype(update_jacobian), LocalResidualVector,
                          decltype(update_residual), decltype(update_solution)>(
                linear_solver, update_jacobian, update_residual,
                update_solution, maximum_iterations, tolerance);

        auto const success_iterations = newton_solver.solve(J_loc);

        if (!success_iterations)
            return false;

        // If the Newton loop didn't run, the linear solver will not be
        // initialized.
        // This happens usually for the first iteration of the first timestep.
        if (*success_iterations == 0)
            linear_solver.compute(J_loc);
    }
	dsw_dpg = Calculate_dSwdP(PG, Sw, X_m);
	dsw_dX = Calculate_dSwdX(PG, X, Sw,X_m);
	dxm_dpg=Calculate_dX_mdP(PG, Sw, X_m);
	dxm_dX = Calculate_dX_mdX(PG, Sw, X_m);
    return true;
}

void EoS_IdealMix::calculateResidual(double const PG, double const X, double Sw, double X_m,
	ResidualVector& res)
{
	// getting unknowns
	const double X_M = 1.0;

	// calculating residual
	res(0) = Calc_equili_molfrac_X_m(PG, Sw, X_m);
	res(1) = Calc_Saturation(PG, X, Sw, X_m, X_M, 303.15);
}

void EoS_IdealMix::calculateJacobian(double const t,
	ProcessLib::SpatialPosition const& x,
	double const PG, double const X,
	JacobianMatrix& Jac, double Sw, double X_m)
{
	// getting unknowns
	double X_M = 1.0;
	double const x_equili_m = PG * Hen / (PG * Hen+rho_mol_h2o);
	const double N_G = PG / R / 303.15;
	const double N_L = rho_mol_h2o / (1 - X_m);
	const double RT= R * 303.15;
	// evaluate J
	Jac.setZero();
	if ((1 - Sw) < (x_equili_m - X_m))
	{
		Jac(0, 0) = -1;
		Jac(0, 1) = 0.0;
	}
	else
	{
		Jac(0, 0) = 0.0;
		Jac(0, 1) = - 1;
	}
	/*double test1 = ((-N_G * (X - X_M) + N_L * (X - X_m)) *
		((1 - Sw) * N_G + Sw * N_L) -
		((1 - Sw) * N_G * (X - X_M) + Sw * N_L * (X - X_m)) *
		(N_G - N_L)) /
		pow(((1 - Sw) * N_G + Sw * N_L), 2);*/
	//double test2 = (rho_mol_h2o* PG * RT*pow(-1 + X_m, 2)) / pow((rho_mol_h2o* RT * Sw + PG*(-1 + Sw) * (-1 + X_m)), 2);
	Jac(1, 0) = (rho_mol_h2o* PG * RT*pow(-1 + X_m, 2)) / pow((rho_mol_h2o* RT * Sw + PG*(-1 + Sw) * (-1 + X_m)), 2);
	//Jac(1, 1) = -Sw*N_L + Sw*(X - X_m)*rho_mol_h2o / std::pow(1 - X_m, 2);
	Jac(1, 1) = -(pow(rho_mol_h2o, 2) * pow(RT, 2) * pow(Sw, 2)) / pow((rho_mol_h2o*RT*Sw + PG*(-1 + Sw) * (-1 + X_m)), 2);
	//std::cout << Jac << std::endl;
}

}  // namespace Solids
}  // namespace MaterialLib

#endif  // 
