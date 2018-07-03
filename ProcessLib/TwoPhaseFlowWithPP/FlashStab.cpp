/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FlashStab.h"
#include <math.h>
#include <malloc.h>

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
    static int stab_itr = 0;
    static double stab_solve_time = 0.;
    static double stab_pred_time = 0.;

    double * flash_calculation_estimate_K(double *K, double P, double T)
    {
        int i, ncomp = 4;
        double P, T;
        
        if (K == NULL) {
            K = (double *)malloc(ncomp * sizeof(*K));
        }


        for (i = 0; i < ncomp; i++) {
            K[i] = comp_PC[i] / P
                * exp(5.37 * (1.0 + comp_AC[i]) * (1.0 - comp_TC[i] / T));
        }

        return K;
    }
    /* ### Composition Initial Guess List
    # The initial guess is from the paper "General Strategy for Stability Testing and Phase-Split Calculation in Two and Three Phases" by Zhidong Li and Abbas Firoozabadi in SPE Journal, 2012.
    #
    # In total, there are $N_c + 4$ sets of initial guess:
    # 1. Wilson's equation: $X_i = K_i x_i$
    # 2. Inverse Wilson's equation: $X_i = x_i / K_i$
    # 3. $X_i = K_i^{\frac{1}{3}} x_i$
    # 4. $X_i = x_i / K_i^{\frac{1}{3}}$
    # 5. Pure component:
    # $$
    # X_i =
    # \left\{
    # \begin{array}[c]
    #  0.9 \quad i = j \\
    #  0.1 / (N_c - 1) \quad otherwise
    # \end{array}
    # \right.
    # $$ for $j = 1, \cdots, N_c$
    */
    double * flash_calculation_stability_analysis_initial_estimate(double P, double T)
    {
        int ncomp = 4, n_guess, i, j;
        double *K, Xi, *est;

        est = (double *)malloc((ncomp + 4) * ncomp * sizeof(*est));
        n_guess = 0;

        K = (double *)malloc(ncomp * sizeof(*K));
        flash_calculation_estimate_K(K, P,T);

        /* Wilson correlation */
        for (i = 0; i < ncomp; i++) {
            Xi = K[i] * z_molarfraction[i];
            *(est + n_guess * ncomp + i) = Xi;
        }
        n_guess += 1;

        /* Inverse Wilson correlation */
        for (i = 0; i < ncomp; i++) {
            Xi = z_molarfraction[i] / K[i];
            *(est + n_guess * ncomp + i) = Xi;
        }
        n_guess += 1;

        /* Wilson correlation power(1./3.) */
        for (i = 0; i < ncomp; i++) {
            Xi = pow(K[i], 1.0 / 3.0) * z_molarfraction[i];
            *(est + n_guess * ncomp + i) = Xi;
        }
        n_guess += 1;

        /* Inverse Wilson correlation power(1./3.) */
        for (i = 0; i < ncomp; i++) {
            Xi = z_molarfraction[i] / pow(K[i], 1.0 / 3.0);
            *(est + n_guess * ncomp + i) = Xi;
        }
        n_guess += 1;

        /* A pure phase */
        for (i = 0; i < ncomp; i++) {
            for (j = 0; j < ncomp; j++) {
                if (i == j) {
                    Xi = 0.9;
                }
                else {
                    Xi = 0.1 / (ncomp - 1);
                }
                *(est + n_guess * ncomp + j) = Xi;
            }
            n_guess += 1;
        }


        /* A hypothetical idea gas: TODO */
        free(K);

        return est;
    }

    /* ### Calculate compositions using X
    # $$
    # x_i = \frac{X_i}{\sum_j X_j}
    y_trial_i=Y_i/{\sum_j Y_j}
    # $$
    */

    void flash_calculation_calculate_trial_phase_composition(double *X_t, double *x)
    {
        double Xs = 0.0;
        int i;
        int ncomp = 4;
        for (i = 0; i < ncomp; i++) {
            Xs += X_t[i];
        }

        for (i = 0; i < ncomp; i++) {
            x[i] = X_t[i] / Xs;
        }
    }

    /* ### Update X using SS method
    # $$
    # X_{t,i} = \exp(\log{x_i} + \log{\phi_i} - \log{\phi_{t,i}})
    # $$
    */

    void flash_calculation_SS_method_update_X(PHASE *phase, PHASE *phase_t, double *X_t)
    {
        int ncomp = 4, i;

        for (i = 0; i < ncomp; i++) {
            X_t[i] = exp(log(z_molarfraction[i]) + log(phase->phi[i]) - log(phase_t->phi[i]));
        }
    }

    /* ### Calculate derivatives of x with respecte to X
    # From $x_i = \frac{X_i}{\sum{X_j}}$, we have
    # $$
    # \frac{\partial x_i}{\partial X_j} = \frac{\sigma_{i,j}}{\sum{X_j}}
    #     - \frac{X_i}{(\sum{X_j})^2}
    # $$
    # where
    # $$
    # \sigma_{i,j} =
    # \left\{
    # \begin{array}{cc}
    # 1 & i = j\\
    # 0 & \text{otherwise}
    # \end{array}
    # \right.
    # $$
    */

    void flash_calculation_calculate_trial_phase_composition_derivative(double *X_t, double *dx_t)
    {
        double sum_X = 0.0, sigma = 0.0;
        int i, j, ncomp=4;

        for (i = 0; i < ncomp; i++) {
            sum_X += X_t[i];
        }

        for (i = 0; i < ncomp; i++) {
            for (j = 0; j < ncomp; j++) {
                if (i == j) {
                    sigma = 1.0;
                }
                else {
                    sigma = 0.0;
                }

                dx_t[i * ncomp + j] = sigma / sum_X - X_t[i] / (sum_X * sum_X);
            }
        }
    }

    /* ### Calculate the values of equilibrium equations
    # $$
    #     D_i = \log{X_i} - \log{z_i} + \log{\phi_i(x)} - \log{\phi_i(z)}
    # $$
    */

    void flash_calculation_calculate_stability_equilibrium_equation(PHASE *phase, PHASE *phase_t, double *X_t, double *D)
    {
        int ncomp = 4, i;
        for (i = 0; i < ncomp; i++) {
            D[i] = log(X_t[i]) - log(phase->mf[i]) + log(phase_t->phi[i]) - log(phase->phi[i]);
        }
    }

    /* # ### Calculate the derivatives of equilibrium equations with respect to X
    # From the equations
    # $$
    #     D_i = \log{X_i} - \log{z_i} + \log{\phi_i(x)} - \log{\phi_i(z)},
    # $$
    # we have the following derivatives
    # $$
    # \frac{\partial D_i}{\partial X_j} = \frac{\sigma_{i,j}}{X_i}
    #                             + \frac{1}{\phi_i(x)}
    #                             \sum_k{\frac{\partial \phi_i(x)}{\partial x_k} \frac{\partial x_k}{\partial X_j}}
    # $$
    */

    void flash_calculation_calculate_stability_equilibrium_equation_derivative(PHASE *phase_t, double *dx_t,
        double *X_t, double *dD)
    {
        int ncomp = 4, i, j, k;
        double sigma = 0.0, temp;

        for (i = 0; i < ncomp; i++) {
            for (j = 0; j < ncomp; j++) {
                if (i == j) {
                    sigma = 1.0;
                }
                else {
                    sigma = 0.0;
                }

                dD[i * ncomp + j] = 1.0 / X_t[i] * sigma;

                temp = 0.0;
                for (k = 0; k < ncomp; k++) {
                    temp += phase_t->dphi_dx[i * ncomp + k] * dx_t[k * ncomp + j];
                }

                dD[i * ncomp + j] += 1.0 / phase_t->phi[i] * temp;
            }
        }
    }

    double flash_calculation_calculate_stability_residual(PHASE *phase, PHASE *phase_t,
        double *X_t, double *res)
    {
        int ncomp = 4, i;
        double tol = 0.0, tmp1, tmp2;

        for (i = 0; i < ncomp; i++) {
            res[i] = log(X_t[i]) + log(phase_t->phi[i]) - log(phase->mf[i]) - log(phase->phi[i]);

            tmp1 = res[i] * res[i];
            tmp2 = log(phase->mf[i]) + log(phase->phi[i]);
            tmp2 = tmp2 * tmp2;

            if ((tmp1 / tmp2) > tol) {
                tol = tmp1 / tmp2;
            }
        }

        return tol;
    }

    /* ### Update X using QNSS method
    # Solve $\delta X$ through the following equation
    # $$
    # J \delta X = - D,
    # $$
    # where $D_i$ is the value of equilibrium equation,
    # $J$ is Jacobian matrix and $J_{i,j} = \frac{\partial D_i}{\partial X_j}$,
    # and update X by
    # $$
    # X = X + \delta X
    # $$
    */

    void flash_calculation_QNSS_method_update_X(double *dD, double *D, double *X_t, int ncomp)
    {
        int i;
        double *x;

        x = (double *)malloc(ncomp * sizeof(*x));

        for (i = 0; i < ncomp; i++) {
            D[i] = -D[i];
        }

        flash_calculation_solve_dense_linear_system(dD, D, x, ncomp);

        for (i = 0; i < ncomp; i++) {
            X_t[i] += x[i];
        }

        free(x);
    }

    /* ### Check Stability using X
    # 1. If $\sum_i{(\log{\frac{X_i}{z_i}})^2} < \epsilon$, the solution is trivial, try next initial guess.
    # 2. If $\sum_i X_i < 1.0$, the phase is stable; otherwise, it is unstable.
    */

    int flash_calculation_check_stability(double *X_t, double *z, int ncomp)
    {
        int i;
        double sum_Y = 0.0, sum_K = 0.0;

        for (i = 0; i < ncomp; i++) {
            sum_Y += X_t[i];

            sum_K += log(X_t[i] / z[i]) * log(X_t[i] / z[i]);
        }

        if (sum_K < 1e-2)
            return -1;

        if (sum_Y < 1.0 + 1e-8) {
            return 1;
        }
        else {
            return 0;
        }

        return -1;
    }

    /* ### QNSS method to solve stability analysis
    # The tangent-plane criterion of stability analysis of a phase with
    #         compositions z results in solving the following set of equations
    #         (Nghiem and Li, 1984):
    # $$
    #             D_i = \log{X_i} - \log{z_i} + \log{\phi_i(x)} - \log{\phi_i(z)}
    # $$
    # where,
    # $$
    #             x_i = X_i / \sum{X_j}
    # $$
    # for the primary variables X.
    #         The phase will be stable if (1) $\sum{X_i} < 1$
    #                                     (2) $\sum{X_i} = 1$ and $X_i \neq z_i$
    #         otherwise the phase is unstable.
    #         Equations D are solved using the SS method first and then QNSS
    #         method.
    #         To solve the above equations, several sets of initial guesses
    #         should be used to avoid trivial solutions.
    */

    int flash_calculation_stability_analysis_QNSS(PHASE *phase, double *K, double tol)
    {
        int ncomp = 4, i, j, n_guess, itr = 0;
        double *D, *dD, *res, *est;
        double *x_t, *dx_t, *X_t, error;
        int system_status = -1;
        PHASE *phase_t;
        D = (double *)malloc(ncomp * sizeof(*D));
        dD = (double *)malloc(ncomp * ncomp * sizeof(*dD));
        res = (double *)malloc(ncomp * sizeof(*res));

        n_guess = ncomp + 4;
        
        est = flash_calculation_stability_analysis_initial_estimate(eos_pressure,eos_temp);

        x_t = (double *)malloc(ncomp * sizeof(*x_t));
        X_t = (double *)malloc(ncomp * sizeof(*X_t));
        dx_t = (double *)malloc(ncomp * ncomp * sizeof(*dx_t));

        phase_t = flash_calculation_phase_new(x_t);//x

        flash_calculation_compute_phase_parameter(phase);
        flash_calculation_calculate_compressibility_factor(phase);
        flash_calculation_calculate_fugacity(phase);

        for (i = 0; i < n_guess; i++) {
            for (j = 0; j < ncomp; j++) {
                X_t[j] = est[i * ncomp + j];
            }

            itr = 0;

            while (1) {
                /* Calculate trial phase composition */
                flash_calculation_calculate_trial_phase_composition(X_t, phase_t->mf);

                /* Calculation trial phase compostion derivative */
                flash_calculation_calculate_trial_phase_composition_derivative(X_t, dx_t);

                /* Calculate the compressibility factor and
                fugacity coefficient for the trial phase
                with compositions x_t */
                flash_calculation_compute_phase_parameter(phase_t);
                flash_calculation_calculate_compressibility_factor(phase_t);
                flash_calculation_calculate_fugacity(phase_t);

                /* Calculate the residual */
                error = flash_calculation_calculate_stability_residual(phase, phase_t, X_t, res);

                /* Check if stop */
                if (error < tol)
                    break;

                /* Update X_t */
                if (error > 1e-5) {
                    flash_calculation_SS_method_update_X(phase, phase_t, X_t);
                }
                else {
                    /* Calculate the equilibrim equation and its derivative */
                    flash_calculation_calculate_stability_equilibrium_equation(phase, phase_t, X_t, D);
                    flash_calculation_calculate_stability_equilibrium_equation_derivative(phase_t, dx_t, X_t, dD);

                    /* Update X_t by X_t += - dD^(-1) * D */
                    flash_calculation_QNSS_method_update_X(dD, D, X_t, ncomp);
                }

                /* Maximum itrations */
                itr += 1;
                if (itr > 1000) {
                    //printf("##### WARNING: Stability_analysis_QNSS reach the maximum iterations!\n");
                    break;
                }
            }

            /* Check stability based on Sum(X_t) */
            system_status = flash_calculation_check_stability(X_t, phase->mf, ncomp);

            /* If the solution is trivial, try the next initial guess;
            otherwise break
            if system_status is 'Unstable' or system_status is 'Stable':
            break */
            if (system_status != -1) {
                break;
            }
        }

        stab_itr += itr;

        /* if K is not None, we output K = X_t / z */
        if (K != NULL) {
            for (i = 0; i < ncomp; i++) {
                K[i] = X_t[i] / phase->mf[i];
            }
        }

        if (system_status == -1) {
            system_status = 1;
        }


        free(x_t);
        free(X_t);
        free(dx_t);
        free(D);
        free(dD);
        free(res);
        free(est);

        //flash_calculation_phase_free(&phase_t);

        return system_status;
    }

    /*
    # ### Calculate phase parameters
    # #### Calculate parameters $a_i$ and $b_i$ of each component
    # For PR EOS:
    # $$
    # b_i = \frac{0.07780 R T_{c,i}}{P_{c,i}}, \\
    # a_i = \frac{0.45724 R^2 T^2_{c,i}}{P_{c,i}} (1 + \lambda_i (1 - (\frac{T}{T_{c,i}})^{0.5})) \\
    # \lambda_i =
    # \left\{
    # \begin{array}{ccc}
    # 0.37464 + 1.5432 \omega_i - 0.26992 \omega_i^2 & \text{when} &\omega_i < 0.49\\
    # 0.3796 + 1.485 \omega_i - 0.1644 \omega_i^2 + 0.01666 \omega_i^3  & \text{when} & \omega_i >= 0.49
    # \end{array}
    # \right.
    # $$
    #
    # For SRK EOS:
    # $$
    # b_i = \frac{0.08664 R T_{c,i}}{P_{c,i}}, \\
    # a_i = \frac{0.42747 R^2 T^2_{c,i}}{P_{c,i}} (1 + \lambda_i (1 - (\frac{T}{T_{c,i}})^{0.5})) \\
    # \lambda_i = 0.48 + 1.574 \omega_i - 0.176 \omega_i^2
    # $$
    # #### Calculate phase parameters and their derivatives
    # $$
    # \left\{
    # \begin{array}{l}
    # a = \sum_{i = 1}^{N_c} x_i S_i \\
    # S_i = \sqrt{a_i} \sum_{i = 1}^{N_c} x_j (1 - k_{i,j}) \sqrt{a_j} \\
    # b = \sum_{i = 1}^{N_c} x_i b_i
    # \end{array}
    # \right.
    # $$
    # and
    # $$
    # \left\{
    # \begin{array}{l}
    # \frac{\partial a}{\partial x_i} = 2 S_i, \quad S_i = \sqrt{a_i} \sum_{i = 1}^{N_c} x_j (1 - k_{i,j}) \sqrt{a_j} \\
    # \frac{\partial b}{\partial x_i} = b_i
    # \end{array}
    # \right.
    # $$
    # where $k_{i,j}$ is the binary interaction parameter.
    # $$
    # \left\{
    # \begin{array}{l}
    # A = \frac{a P}{R^2 T^2}\\
    # B = \frac{b P}{R T}
    # \end{array}
    # \right.
    # $$
    */

    void flash_calculation_compute_phase_parameter(PHASE *phase)
    {
        double R, T, P;
        int ncomp=4, i, j;
        double a, b, da_dT, da_dT2;

        ncomp = 4;
        T = 310.95;
        P = 76;
        R = 8.314;
        int const type = 1;
        if (type == 0) {
            for (i = 0; i < ncomp; i++) {
                double AC, TC, PC, lambda_i, tmp, alpha_i;
                double dtmp, dalpha_i, ddtmp, ddalpha_i;

                AC = comp_AC[i];
                TC = comp_TC[i];
                PC = comp_PC[i];

                lambda_i = 0.48 + 1.574 * AC - 0.176 * AC * AC;

                tmp = 1.0 + lambda_i * (1.0 - sqrt(T / TC));
                alpha_i = tmp * tmp;

                phase->ai[i] = 0.42747 * alpha_i * R * R * TC * TC / PC;
                phase->bi[i] = 0.08664 * R * TC / PC;

                dtmp = lambda_i * (-1.0) * 0.5 / sqrt(T / TC) / TC;
                dalpha_i = 2.0 * tmp * dtmp;
                phase->dai_dT[i] = 0.42747 * dalpha_i * R * R * TC * TC / PC;

                ddtmp = lambda_i * (-1.0) * 0.5 / TC * (-0.5) * pow(T / TC, -1.5) / TC;
                ddalpha_i = 2.0 * (dtmp * dtmp + tmp * ddtmp);
                phase->dai_dT2[i] = 0.42747 * ddalpha_i * R * R * TC * TC / PC;
            }
        }
        else if (type == 1) {// here is the peng robinson
            for (i = 0; i < ncomp; i++) {
                double AC, TC, PC, lambda_i, tmp, alpha_i;
                double dtmp, dalpha_i, ddtmp, ddalpha_i;
                phase->mf[i] = z_molarfraction[i];
                AC = comp_AC[i];//w^i
                TC = comp_TC[i];
                PC = comp_PC[i];

                if (AC < 0.49) {
                    lambda_i = 0.37464 + 1.54226 * AC - 0.26992 * AC * AC;
                }
                else {
                    lambda_i = 0.3796 + 1.485 * AC - 0.1644 * AC * AC + 0.01666 * AC * AC * AC;
                }

                tmp = 1.0 + lambda_i * (1.0 - sqrt(T / TC));//alpha function
                alpha_i = 1.0 + lambda_i * (1.0 - sqrt(T / TC));
                alpha_i = alpha_i * alpha_i;

                phase->ai[i] = 0.45724 * alpha_i * R * R * TC * TC / PC;
                phase->bi[i] = 0.077796 * R * TC / PC;
                //calculate the derivatives 
                dtmp = lambda_i * (-1.0) * 0.5 / sqrt(T / TC) / TC;
                dalpha_i = 2.0 * tmp * dtmp;
                phase->dai_dT[i] = 0.45724 * dalpha_i * R * R * TC * TC / PC;
                //calculate the second order derivatives
                ddtmp = lambda_i * (-1.0) * 0.5 / TC * (-0.5) * pow(T / TC, -1.5) / TC;
                ddalpha_i = 2.0 * (dtmp * dtmp + tmp * ddtmp);
                phase->dai_dT2[i] = 0.45724 * ddalpha_i * R * R * TC * TC / PC;
            }
        }


        a = 0.0;
        b = 0.0;
        da_dT = 0.0;
        da_dT2 = 0.0;
        for (i = 0; i < ncomp; i++) {
            double da, dda_dT, dda_dT2;
            b += phase->mf[i] * phase->bi[i];
            phase->db[i] = phase->bi[i];

            da = 0.0;
            dda_dT = 0.0;
            dda_dT2 = 0.0;

            for (j = 0; j < ncomp; j++) {
                a += phase->mf[i] * phase->mf[j] * (1.0 - comp_binary[i][j]) * sqrt(phase->ai[i] * phase->ai[j]);
                da += phase->mf[j] * (1.0 - comp_binary[i][j]) * sqrt(phase->ai[j]);

                da_dT += phase->mf[i] * phase->mf[j] * (1.0 - comp_binary[i][j]) * 0.5 / sqrt(phase->ai[i] * phase->ai[j])
                    * (phase->ai[j] * phase->dai_dT[i] + phase->ai[i] * phase->dai_dT[j]);
                da_dT2 += phase->mf[i] * phase->mf[j] * (1.0 - comp_binary[i][j]) * 0.5 / sqrt(phase->ai[i] * phase->ai[j])
                    * (phase->ai[j] * phase->dai_dT2[i] + phase->dai_dT[j] * phase->dai_dT[i] + phase->ai[i] * phase->dai_dT2[j] + phase->dai_dT[i] * phase->dai_dT[j])
                    + phase->mf[i] * phase->mf[j] * (1.0 - comp_binary[i][j]) * 0.5 * (-0.5) * pow(phase->ai[i] * phase->ai[j], -1.5)
                    * pow(phase->ai[j] * phase->dai_dT[i] + phase->ai[i] * phase->dai_dT[j], 2.0);


                dda_dT += phase->mf[j] * (1.0 - comp_binary[i][j]) * 0.5 / sqrt(phase->ai[j]) * phase->dai_dT[j];
                dda_dT2 += phase->mf[j] * (1.0 - comp_binary[i][j]) * 0.5 * (-0.5) * pow(phase->ai[j], -1.5) * pow(phase->dai_dT[j], 2.0)
                    + phase->mf[j] * (1.0 - comp_binary[i][j]) * 0.5 / sqrt(phase->ai[j]) * phase->dai_dT2[j];
            }

            phase->da[i] = 2.0 * sqrt(phase->ai[i]) * da;

            phase->dda_dT[i] = 2.0 * sqrt(phase->ai[i]) * dda_dT
                + 2.0 * 0.5 / sqrt(phase->ai[i]) * phase->dai_dT[i] * da;
            phase->dda_dT2[i] = 2.0 * sqrt(phase->ai[i]) * dda_dT2
                + 2.0 * 0.5 / sqrt(phase->ai[i]) * phase->dai_dT[i] * dda_dT
                + 2.0 * 0.5 / sqrt(phase->ai[i]) * phase->dai_dT[i] * da_dT
                + 2.0 * 0.5 / sqrt(phase->ai[i]) * phase->dai_dT2[i] * da
                + 2.0 * 0.5 * (-0.5) * pow(phase->ai[i], -1.5) * pow(phase->dai_dT[i], 2.0) * da;
        }

        phase->a = a;
        phase->b = b;
        phase->da_dT = da_dT;
        phase->da_dT2 = da_dT2;

        /* A, B */

        phase->A = a * P / (R * R * T * T);
        phase->dAp = phase->A / P;

        phase->dA_dT = da_dT * phase->A / a + a * P / (R * R) * (-2.0) / (T * T * T);
        phase->dA_dT2 = da_dT2 * phase->A / a
            + da_dT * phase->dA_dT / a
            + da_dT * phase->A * (-1.0) / (a * a) * da_dT
            + da_dT * P / (R * R) * (-2.0) / (T * T * T)
            + a * P / (R * R) * (-2.0) * (-3.0) / (T * T * T * T);

        for (i = 0; i < ncomp; i++) {
            phase->dAx[i] = phase->A / phase->a * phase->da[i];
        }

        phase->B = b * P / (R * T);
        phase->dBp = phase->B / P;

        phase->dB_dT = b * P / R * (-1.0) / (T * T);
        phase->dB_dT2 = b * P / R * (-1.0) * (-2.0) / (T * T * T);

        for (i = 0; i < ncomp; i++) {
            phase->dBx[i] = phase->B / phase->b * phase->db[i];
        }
    }


/*
# ### Calculate compressibility factore Z and its derivatives
# Solve the equation
# $$Z^3 + c_2 Z^2 + c_1 Z + c_0 = 0,$$
# where
# $$
# \left\{
# \begin{array}{l}
# c_2 = (u - 1) B - 1 \\
# c_1 = A + (w - u) B^2 - u B \\
# c_0 = - A B - w B^2 - w B^3
# \end{array}
# \right.
# $$
# to get the compressibility factore Z.
# The derivatives are
# $$
# \left\{
# \begin{array}{l}
# \frac{\partial Z}{\partial P} = - \frac{\frac{\partial c_2}{\partial P} Z^2
#                                     + \frac{\partial c_1}{\partial P} Z + \frac{\partial c_0}{\partial P}}
#                                 {3 Z^2 + 2 c_2 Z + c_1} \\
# \frac{\partial Z}{\partial x_i} = - \frac{\frac{\partial c_2}{\partial x_i} Z^2
#                                     + \frac{\partial c_1}{\partial x_i} Z + \frac{\partial c_0}{\partial x_i}}
#                                 {3 Z^2 + 2 c_2 Z + c_1}
# \end{array}
# \right.
# $$
# #### Root selection
# ##### Three real roots
# If three roots are obtained, the middle one is ignored. From the rest two roots, the one that results in the lowest Gibb's free energy will be selected. Let $Z_A$ and $Z_B$ be the two real roots resulting in free energy $G_A$ and $G_B$ respectively. Since free energy
# $$
# G = \sum_i x_i \log{f_i} \\
# G_A - G_B = \log{\frac{Z_B - B}{Z_A - B}} + \frac{1}{\sigma_2 - \sigma_1} \frac{A}{B} \log{\frac{Z_B + \sigma_2 B}{Z_A + \sigma_2 B} \frac{Z_A + \sigma_1 B}{Z_B + \sigma_1 B}} - (Z_B - Z_A)
# $$
# If $G_A - G_B > 0$, $Z_B$ will be selected and vice versa. For single-phase fluids, if the above scheme selects the largest $Z$ root, the fluid is said to be vapour. Similarly, if the smallest $Z$ root is chosen, the fluid is said to be liquid.
# ##### One real root
# If only one real root is obtained, we use the method introduced in the paper "**Comparison of Phase Identification Methods Used in Oil Industry Flow Simulations**" to identify phase.
# The thermal expansion coefficent $\alpha$ is defined by:
# $$
# \alpha = (\frac{\partial \log{V}}{\partial T})_P = \frac{1}{V} (\frac{\partial V}{\partial T})_P
# $$
# If $\frac{\partial \alpha}{\partial T} > 0$, the fluid is liquid; otherwise, it is vapour.
*/

void flash_calculation_calculate_compressibility_factor(PHASE *phase)
{
    int ncomp, i, nroot;
    double u, w, A, B, dAp, dBp;
    double c0, c1, c2, dc0, dc1, dc2,
        *dc0_dx, *dc1_dx, *dc2_dx;
    double dA_dT, dB_dT, dc2_dT, dc1_dT, dc0_dT;
    double dA_dT2, dB_dT2, dc2_dT2, dc1_dT2, dc0_dT2;
    double Q, J, D, Z1, Z2, Z3;
    double Z_h, Z_l;
    double top, down;
    ncomp =4;

    u = 2;
    w = -1;
    A = phase->A;
    B = phase->B;
    dAp = phase->dAp;
    dBp = phase->dBp;

    /* Solve the Cubit EOS and get one or three real roots */
    c2 = (u - 1.0) * B - 1.0;
    c1 = A + (w - u) * B * B - u * B;
    c0 = -A * B - w * B * B - w * B * B * B;

    dc2 = (u - 1.0) * dBp;
    dc1 = dAp + (w - u) * 2.0 * B * dBp - u * dBp;
    dc0 = -dAp * B - A * dBp - w * 2.0 * B * dBp - w * 3.0 * B * B * dBp;

    dc2_dx = (double *)malloc(ncomp * sizeof(*dc2_dx));
    dc1_dx = (double *)malloc(ncomp * sizeof(*dc1_dx));
    dc0_dx = (double *)malloc(ncomp * sizeof(*dc0_dx));

    for (i = 0; i < ncomp; i++) {
        dc2_dx[i] = (u - 1.0) * phase->dBx[i];
        dc1_dx[i] = phase->dAx[i] + (w - u) * 2.0 * B * phase->dBx[i] - u * phase->dBx[i];
        dc0_dx[i] = -phase->dAx[i] * B - A * phase->dBx[i]
            - w * 2.0 * B * phase->dBx[i] - w * 3.0 * B * B * phase->dBx[i];
    }


    dA_dT = phase->dA_dT;
    dB_dT = phase->dB_dT;
    dc2_dT = (u - 1.0) * dB_dT;
    dc1_dT = dA_dT + (w - u) * 2.0 * B * dB_dT - u * dB_dT;
    dc0_dT = -dA_dT * B - A * dB_dT - w * 2.0 * B * dB_dT
        - w * 3.0 * B * B * dB_dT;

    dA_dT2 = phase->dA_dT2;
    dB_dT2 = phase->dB_dT2;
    dc2_dT2 = (u - 1.0) * dB_dT2;
    dc1_dT2 = dA_dT2 + (w - u) * 2.0 * (dB_dT * dB_dT
        + B * dB_dT2) - u * dB_dT2;
    dc0_dT2 = -dA_dT2 * B - dA_dT * dB_dT
        - dA_dT * dB_dT - A * dB_dT2
        - w * 2.0 * (dB_dT * dB_dT + B * dB_dT2)
        - w * 3.0 * (2.0 * B * dB_dT * dB_dT
            + B * B * dB_dT2);

    Q = (3.0 * c1 - c2 * c2) / 9.0;
    J = (9.0 * c2 * c1 - 27.0 * c0 - 2.0 * c2 * c2 * c2) / 54.0;
    D = Q * Q * Q + J * J;

    Z1 = 0.0;
    Z2 = 0.0;
    Z3 = 0.0;

    if (D > 0.0) {
        double tmp1, tmp2;

        nroot = 1;
        tmp1 = J + sqrt(D);
        tmp2 = J - sqrt(D);

        if (tmp1 > 0) {
            Z1 += pow(tmp1, 1.0 / 3.0);
        }
        else {
            Z1 += -pow(-tmp1, 1.0 / 3.0);
        }

        if (tmp2 > 0) {
            Z1 += pow(tmp2, 1.0 / 3.0);
        }
        else {
            Z1 += -pow(-tmp2, 1.0 / 3.0);
        }

        Z1 += -c2 / 3.0;
    }
    else if (D < 0.0) {
        double theta;

        nroot = 3;
        theta = acos(J / sqrt(-Q * Q * Q));
        Z1 = 2.0 * sqrt(-Q) * cos(theta / 3.0) - c2 / 3.0;
        Z2 = 2.0 * sqrt(-Q) * cos(theta / 3.0 + 2.0 * M_PI / 3.0) - c2 / 3.0;
        Z3 = 2.0 * sqrt(-Q) * cos(theta / 3.0 + 4.0 * M_PI / 3.0) - c2 / 3.0;
    }
    else {
        double tmp;

        nroot = 3;

        if (J > 0) {
            tmp = pow(J, 1.0 / 3.0);
        }
        else {
            tmp = -pow(-J, 1.0 / 3.0);
        }

        Z1 = 2.0 * tmp - c2 / 3.0;
        Z2 = Z3 = -tmp - c2 / 3.0;
    }

    phase->nroot = nroot;

    Z_h = -1e20;
    Z_l = 1e20;
    if (nroot == 3) {
        double sigma_1, sigma_2, dG;
        //return the minimum z-factor if multiple roots are found.
        //for the liquid phase
        if (Z_l > Z1) {
            Z_l = Z1;
        }
        if (Z_l > Z2) {
            Z_l = Z2;
        }
        if (Z_l > Z3) {
            Z_l = Z3;
        }
        //return the max z-factor
        //for the vapor/gas phase
        if (Z_h < Z1) {
            Z_h = Z1;
        }
        if (Z_h < Z2) {
            Z_h = Z2;
        }
        if (Z_h < Z3) {
            Z_h = Z3;
        }

        sigma_1 = sqrt(2)-1;//need conform
        sigma_2 = eos->para_sigma2;

        dG = (Z_h - Z_l) + log((Z_l - B) / (Z_h - B))
            - A / (B * (sigma_2 - sigma_1))
            * log((Z_l + sigma_1 * B) / (Z_l + sigma_2 * B)
                * (Z_h + sigma_2 * B) / (Z_h + sigma_1 * B));

        if (dG > 0.0) {
            phase->Z = Z_l;
            phase->phase_no = 0;
        }
        else {
            phase->Z = Z_h;
            phase->phase_no = 1;
        }
    }
    else {
        double Tc, Dc;

        Tc = 0;
        Dc = 0;
        for (i = 0; i < ncomp; i++) {
            Tc += z_molarfraction[i] * comp_TC[i] * comp[i].VC;
            Dc += z_molarfraction[i] * comp[i].VC;
        }
        Tc = Tc / Dc;
        if (eos_temp < Tc) {
            phase->phase_no = 0;
        }
        else {
            phase->phase_no = 1;
        }
        phase->Z = Z1;
    }

    phase->dZ = -(dc2 * phase->Z * phase->Z + dc1 * phase->Z + dc0)
        / (3.0 * phase->Z * phase->Z + 2.0 * c2 * phase->Z + c1);

    for (i = 0; i < ncomp; i++) {
        phase->dZ_dx[i] = -(dc2_dx[i] * phase->Z * phase->Z
            + dc1_dx[i] * phase->Z + dc0_dx[i])
            / (3.0 * phase->Z * phase->Z + 2.0 * c2 * phase->Z + c1);
    }

    phase->dZ_dT = -(dc2_dT * phase->Z * phase->Z + dc1_dT * phase->Z + dc0_dT)
        / (3.0 * phase->Z * phase->Z + 2.0 * c2 * phase->Z + c1);


    top = -(dc2_dT * phase->Z * phase->Z + dc1_dT * phase->Z + dc0_dT);
    down = 3.0 * phase->Z * phase->Z + 2.0 * c2 * phase->Z + c1;

    phase->dZ_dT2 = -(dc2_dT2 * phase->Z * phase->Z + dc2_dT * 2.0 * phase->Z * phase->dZ_dT
        + dc1_dT2 * phase->Z + dc1_dT * phase->dZ_dT + dc0_dT2) / down
        + top * (-1.0) / (down * down) * (3.0 * 2.0 * phase->Z * phase->dZ_dT
            + 2.0 * c2 * phase->dZ_dT + 2.0 * dc2_dT * phase->Z
            + dc1_dT);

    if (nroot > 0) {
        double P, Z, T, R, V, dV_dT, dV_dT2, dalpha_dT;

        P = eos_pressure;
        T = eos_temp;
        Z = phase->Z;
        R = phase->R;

        V = Z * R * T / P;

        dV_dT = phase->dZ_dT * R * T / P + Z * R / P;
        dV_dT2 = phase->dZ_dT2 * R * T / P + 2.0 * phase->dZ_dT * R / P;

        dalpha_dT = dV_dT2 / V - 1.0 / (V * V) * dV_dT * dV_dT;

        if (dalpha_dT > 0.0) {
            /* liquid */
            phase->phase_no = 0;
        }
        else {
            /* vapor */
            phase->phase_no = 1;
        }
    }

    free(dc2_dx);
    free(dc1_dx);
    free(dc0_dx);
}

/* ### Calculate component fugacity
# $$
# f_i = P x_i \phi_i \\
# \log{\phi_i} = \frac{b_i}{b}(Z - 1) - \log{(Z-B)}
#             - \frac{A}{B v}(\frac{b_i}{b} - \frac{1}{a} \frac{\partial a}{\partial x_i})
#                 \log{(\frac{2 Z + B(u+v)}{2 Z + B (u-v)})} \\
# v = \sqrt{u^2 - 4 w}
# $$
*/

void flash_calculation_calculate_fugacity(PHASE *phase)
{
    int ncomp = phase->ncomp, i, j;
    double A, B, Z, a, b, u, w, v, p;

    A = phase->A;
    B = phase->B;
    Z = phase->Z;
    a = phase->a;
    b = phase->b;
    u = 2;
    w = -1;
    v = sqrt(u * u - 4.0 * w);
    p = eos_pressure;

    for (i = 0; i < ncomp; i++) {
        phase->phi[i] = phase->bi[i] / b * (Z - 1.0) - log(Z - B);
        phase->phi[i] += A / (B * v) * (phase->bi[i] / b - phase->da[i] / a)
            * log((2.0 * Z + B * (u + v)) / (2.0 * Z + B * (u - v)));
        phase->phi[i] = exp(phase->phi[i]);
        phase->fug[i] = p * phase->mf[i] * phase->phi[i];

        phase->dphi[i] = phase->bi[i] / b * phase->dZ;
        phase->dphi[i] += -(phase->dZ - phase->dBp) / (Z - B);
        phase->dphi[i] += (phase->bi[i] / b - phase->da[i] / a) / v
            * ((phase->dAp / B - A / (B * B) * phase->dBp)
                * log((2.0 * Z + B * (u + v)) / (2.0 * Z + B * (u - v)))
                + A / B * ((2.0 * phase->dZ + phase->dBp * (u + v)) / (2.0 * Z + B * (u + v))
                    - (2.0 * phase->dZ + phase->dBp * (u - v)) / (2.0 * Z + B * (u - v))));
        phase->dphi[i] *= phase->phi[i];
        phase->dfug[i] = phase->fug[i] / p + phase->mf[i] * p * phase->dphi[i];

        /*print("phi[%d]: %lf, pres: %f, mf: %f, phi: %f"%(i, phase.phi[i], phase.cubic_eos->pres, \
        phase.mf[i], phase.phi[i])) */

        for (j = 0; j < ncomp; j++) {
            phase->dphi_dx[i * ncomp + j] = phase->bi[i] / b * phase->dZ_dx[j]
                - phase->bi[i] / (b * b) * (Z - 1.0) * phase->db[j]
                - (phase->dZ_dx[j] - phase->dBx[j]) / (Z - B)
                + (phase->bi[i] / b - phase->da[i] / a) / v
                * ((phase->dAx[j] / B - phase->dBx[j] * A / (B * B))
                    * log((2.0 * Z + B * (u + v)) / (2.0 * Z + B * (u - v)))
                    + A / B * ((2.0 * phase->dZ_dx[j] + phase->dBx[j] * (u + v))
                        / (2.0 * Z + B * (u + v))
                        - (2.0 * phase->dZ_dx[j] + phase->dBx[j] * (u - v))
                        / (2.0 * Z + B * (u - v))))
                + A / (B * v) * log((2.0 * Z + B * (u + v))
                    / (2.0 * Z + B * (u - v)))
                * (-phase->bi[i] / (b * b) * phase->db[j]
                    + phase->da[i] * phase->da[j] / (a * a)
                    - 2.0 / a * sqrt(phase->ai[i] * phase->ai[j]) * (1.0 - comp_binary[i][j]));
            phase->dphi_dx[i * ncomp + j] *= phase->phi[i];
        }

        phase->dphi_dT[i] = phase->bi[i] / b * phase->dZ_dT;
        phase->dphi_dT[i] += -(phase->dZ_dT - phase->dB_dT) / (Z - B);
        phase->dphi_dT[i] += (phase->bi[i] / b - phase->da[i] / a) / v
            * ((phase->dA_dT / B - A / (B * B) * phase->dB_dT)
                * log((2.0 * Z + B * (u + v)) / (2.0 * Z + B * (u - v)))
                + A / B * ((2.0 * phase->dZ_dT + phase->dB_dT * (u + v)) / (2.0 * Z + B * (u + v))
                    - (2.0 * phase->dZ_dT + phase->dB_dT * (u - v)) / (2.0 * Z + B * (u - v))));

        phase->dphi_dT[i] += A / B / v * log((2.0 * Z + B * (u + v)) / (2.0 * Z + B * (u - v)))
            * (phase->da[i] / (a * a) * phase->da_dT - phase->dda_dT[i] / a);
        phase->dphi_dT[i] *= phase->phi[i];
    }
}


int flash_calculation_solve_dense_linear_system_LU(int n, double *a0, int pvt[])
{
    int i, j, k;
    double d, r;

#define a(i,j)	(a0[(i) * n + (j)])

    for (i = 0; i < n; i++) {
        d = fabs(a(i, i));
        k = i;
        for (j = i + 1; j < n; j++) {
            if ((r = fabs(a(j, i))) > d) {
                d = r;
                k = j;
            }
        }

        if (d == 0.0)
            return 0;

        pvt[i] = k;
        if (k != i) { /* exchange row i and row k */
            for (j = i; j < n; j++) {
                d = a(i, j);
                a(i, j) = a(k, j);
                a(k, j) = d;
            }
        }

        if ((d = a(i, i)) != 1.0) {
            a(i, i) = d = 1.0 / d;
            for (j = i + 1; j < n; j++)
                a(i, j) *= d;
        }

        for (j = i + 1; j < n; j++) {
            if ((d = a(j, i)) == 0.0) continue;

            for (k = i + 1; k < n; k++) a(j, k) -= d * a(i, k);
        }
    }

#undef a

    return 1;
}

/* solves L*U*X = B. */
void flash_calculation_solve_dense_linear_system_SV(int n, double *a0,
    int pvt[], int m, double *b0)
{
    int i, j, k;
    double d;

#define a(i,j)	(a0[(i) * n + (j)])
#define b(i,j)	(b0[(i) * m + (j)])

    for (i = 0; i < n; i++) {
        k = pvt[i];
        if (k != i) {
            /* exchange row i with row k */
            for (j = 0; j < m; j++) {
                d = b(i, j);
                b(i, j) = b(k, j);
                b(k, j) = d;
            }
        }

        if ((d = a(i, i)) != 1.0) {
            for (j = 0; j < m; j++) b(i, j) *= d;
        }

        for (j = i + 1; j < n; j++) {
            if ((d = a(j, i)) == 0.0) continue;
            for (k = 0; k < m; k++) b(j, k) -= d * b(i, k);
        }
    }

    for (i = n - 2; i >= 0; i--)
        for (j = i + 1; j < n; j++) {
            d = a(i, j);
            for (k = 0; k < m; k++) b(i, k) -= d * b(j, k);
        }
#undef a
#undef b
    return;
}

int flash_calculation_solve_dense_linear_system(double *M, double *b,
    double *x, int n)
{
    int *pvt, i, LU;

    pvt = (int *)malloc(n * sizeof(*pvt));

    LU = flash_calculation_solve_dense_linear_system_LU(n, M, pvt);

    if (LU) {
        flash_calculation_solve_dense_linear_system_SV(n, M,
            pvt, 1, b);

        for (i = 0; i < n; i++) {
            x[i] = b[i];
        }

        free(pvt);
    }
    else {
        free(pvt);

        return LU;
    }

    return 1;
}
PHASE * flash_calculation_phase_new(double *mf)
{
    PHASE *phase;
    int ncomp;

    phase = (PHASE *)malloc(sizeof*(phase));

    ncomp = 4;
    phase->ncomp = ncomp;
    phase->mf = mf;
    phase->R = 82.05736;
    phase->density = 0.0;

    phase->A = 0.0;
    phase->dAp = 0.0;
    phase->dA_dT = 0.0;
    phase->dA_dT2 = 0.0;
    phase->dAx = (double *)malloc(ncomp * sizeof(*(phase->dAx)));

    phase->B = 0.0;
    phase->dBp = 0.0;
    phase->dB_dT = 0.0;
    phase->dB_dT2 = 0.0;
    phase->dBx = (double *)malloc(ncomp * sizeof(*(phase->dAx)));

    phase->nroot = 0;
    phase->Z = -1;
    phase->dZ = 0;
    phase->dZ_dT = 0.0;
    phase->dZ_dT2 = 0.0;
    phase->dZ_dx = (double *)malloc(ncomp * sizeof(*(phase->dZ_dx)));
    phase->phase_no = -1;

    phase->fug = (double *)malloc(ncomp * sizeof(*(phase->fug)));
    phase->phi = (double *)malloc(ncomp * sizeof(*(phase->phi)));
    phase->dphi = (double *)malloc(ncomp * sizeof(*(phase->dphi)));
    phase->dphi_dT = (double *)malloc(ncomp * sizeof(*(phase->dphi_dT)));
    phase->dphi_dx = (double *)malloc(ncomp * ncomp * sizeof(*(phase->dphi_dx)));
    phase->dfug = (double *)malloc(ncomp * sizeof(*(phase->dfug)));

    phase->ai = (double *)malloc(ncomp * sizeof(*(phase->ai)));
    phase->dai_dT = (double *)malloc(ncomp * sizeof(*(phase->dai_dT)));
    phase->dai_dT2 = (double *)malloc(ncomp * sizeof(*(phase->dai_dT2)));
    phase->a = 0.0;
    phase->da_dT = 0.0;
    phase->da = (double *)malloc(ncomp * sizeof(*(phase->da)));
    phase->dda_dT = (double *)malloc(ncomp * sizeof(*(phase->dda_dT)));
    phase->dda_dT2 = (double *)malloc(ncomp * sizeof(*(phase->dda_dT2)));

    phase->bi = (double *)malloc(ncomp * sizeof(*(phase->dZ_dx)));
    phase->b = 0.0;
    phase->db = (double *)malloc(ncomp * sizeof(*(phase->dZ_dx)));

    return phase;
}





    double const comp_PC[4] = { 46,73.8,89.4,220.5 };
    double const comp_TC[4] = { 190.6,204.2,373.2,647.3 };
    double const comp_AC[4] = { 0.008,0.225,0.1,0.344 };
    double const comp_binary[4][4] = {
        {0, 0, 0, 0},
    {0.1005,0,0,0},
    {0.0755,0.0999,0,0},
    {0.4928,0,0.04,0}
    };
    double const eos_temp = 310.95;
    double const eos_pressure = 76;
    double const z_molarfraction[4] = { 0.1488,0.2991,0.0494,0.5027 };
    //define a struct "phase"
    

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
