/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowProcess.cpp
 *
 * Created on August 19, 2016, 1:38 PM
 */

#include "FlashSplit.h"

#include <cassert>

#include "MeshLib/PropertyVector.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
    static int verb = 0;
    static int split_failure = 0;
    static int split_itr = 0;
    static double split_solve_time = 0.;
    static double split_pred_time = 0.;
    /* ## 3. Two-phase Flash Calculation
    # The following code is designed for two-phase flash calculations.
    */

    /* ### Calculate compositions for liquid and vapour phases
    # $$
    # x_{l,i} = \frac{z_i}{1 + (K_i - 1) F_v} \\
    # x_{v,i} = \frac{K_i z_i}{1 + (K_i - 1) F_v}
    # $$
    # where $z_i$ is the feed composition, $F_v$ is mole fraction of vapour phase.
    */

    void flash_calculation_calculate_composition(double *K, double *z, double F_v,
        double *x_l, double *x_v, int ncomp)
    {
        int i;

        for (i = 0; i < ncomp; i++) {
            x_l[i] = z[i] / (1.0 + (K[i] - 1.0) * F_v);
            x_v[i] = K[i] * z[i] / (1.0 + (K[i] - 1.0) * F_v);
        }
    }

    /* ### Calculate derivatives of compositions for liquid and vapour phases with respect to $K_i$
    # $$
    # \frac{\partial x_{l,i}}{\partial K_i} = - \frac{z_i F_v}{(1 + (K_i - 1) F_v)^2}  \\
    # \frac{\partial x_{v,i}}{\partial K_i} = \frac{z_i}{1 + (K_i - 1) F_v} - \frac{z_i K_i F_v}{(1 + (K_i - 1) F_v)^2}
    # $$
    */

    void flash_calculation_calculate_composition_derivative(double *K, double *z, double F_v,
        double *dx_l, double *dx_v, int ncomp)
    {
        int i;
        double temp;

        for (i = 0; i < ncomp; i++) {
            temp = 1.0 + (K[i] - 1.0) * F_v;
            dx_l[i] = -z[i] / (temp * temp) * F_v;
            dx_v[i] = z[i] / temp - z[i] * K[i] / (temp * temp) * F_v;
        }
    }

    /* ### Calculate equilibrium equation
    # $$
    # G_i = \log{x_{v,i}} + \log{\phi_{v,i}(x_v)}
    #             - \log{x_{l,i}} - \log{\phi_{l,i}(x_l)}
    # $$
    */

    double flash_calculation_calculate_equilibrium_equation(PHASE *phase_L, PHASE *phase_V, double *G)
    {
        int ncomp = phase_L->ncomp, i;
        double error = 0.0;

        for (i = 0; i < ncomp; i++) {
            G[i] = log(phase_V->mf[i]) + log(phase_V->phi[i]);
            G[i] += -log(phase_L->mf[i]) - log(phase_L->phi[i]);

            error += G[i] * G[i];
        }

        return error;
    }

    /* ### Calculate derivatives of equilibrium equations with respect to K
    # $$
    # \frac{\partial G_i}{\partial K_j} = \frac{1}{x_{v,i}} \frac{\partial x_{v,i}}{\partial K_i} \sigma_{i,j}
    #                 + \frac{1}{\phi_{v,i}} \frac{\partial \phi_{v,i}}{\partial x_{v,j}} \frac{\partial x_{v,j}}{\partial K_j}
    #                - \frac{1}{x_{l,i}} \frac{\partial x_{l,i}}{\partial K_i} \sigma_{i,j}
    #                - \frac{1}{\phi_{l,i}} \frac{\partial \phi_{l,i}}{\partial x_{l,j}} \frac{\partial x_{l,j}}{\partial K_j}
    # $$
    */

    void flash_calculation_calculate_equilibrium_equation_derivative(PHASE *phase_L, PHASE *phase_V,
        double *dx_l, double *dx_v, double *dG)
    {
        int ncomp = phase_L->ncomp, i, j;

        for (i = 0; i < ncomp; i++) {
            for (j = 0; j < ncomp; j++) {
                if (i == j) {
                    dG[i * ncomp + j] = 1.0 / phase_V->mf[i] * dx_v[i] - 1.0 / phase_L->mf[i] * dx_l[i];
                }
                else {
                    dG[i * ncomp + j] = 0.0;
                }

                dG[i * ncomp + j] += 1.0 / phase_V->phi[i] * phase_V->dphi_dx[i * ncomp + j] * dx_v[j]
                    - 1.0 / phase_L->phi[i] * phase_L->dphi_dx[i * ncomp + j] * dx_l[j];
            }
        }
    }

    /* ### Update K using QNSS method
    # Solve
    # $$
    # dG \delta K = - G,
    # $$
    # and update $K$ by
    # $$
    # K = K + \delta K
    # $$
    */

    void flash_calculation_QNSS_method_update_K(double *dG, double *G, double *K, int ncomp)
    {
        int i;
        double *x;

        x = (double *)malloc(ncomp * sizeof(*x));

        for (i = 0; i < ncomp; i++) {
            G[i] = -G[i];
        }

        flash_calculation_solve_dense_linear_system(dG, G, x, ncomp);

        for (i = 0; i < ncomp; i++) {
            K[i] += x[i];
        }

        free(x);
    }

    /* ### Calculate the value of the Rachford-Rice equation
    # $$
    # V = \sum_i{\frac{z_i (K_i - 1)}{1 + (K_i - 1) F_v}}
    # $$
    */

    double flash_calculation_calculate_RachfordRice_equation_value(double *K, double *z, double n_V, int ncomp)
    {
        int i;
        double value = 0.0;

        for (i = 0; i < ncomp; i++) {
            value += z[i] * (K[i] - 1.0) / (1.0 + (K[i] - 1.0) * n_V);
        }

        return value;
    }

    /* ### Calculate the derivative of Rachford-Rice equation
    # $$
    # \frac{\partial V}{\partial F_v} = - \sum_i{\frac{z_i (K_i - 1)^2}{(1 + (K_i - 1)F_v)^2}}
    # $$
    */

    double flash_calculation_calculate_RachfordRice_equation_derivative(double *K, double *z, double n_V, int ncomp)
    {
        int i;
        double value = 0.0;

        for (i = 0; i < ncomp; i++) {
            value += -z[i] * (K[i] - 1.0) * (K[i] - 1.0)
                / (1.0 + (K[i] - 1.0) * n_V)
                / (1.0 + (K[i] - 1.0) * n_V);
        }

        return value;
    }

    /* ### Solve the Rachford-Rice equation
    # Solve $F_v$ using an iterative method:
    # 1. Solve $\delta F_v$
    # $$
    # \frac{\partial V}{\partial F_v} \delta F_v = - V
    # $$
    # 2. Update $F_v$
    # $$
    # F_v = F_v + \delta F_v
    # $$
    # 3. Re-calculate $V$ and $\frac{\partial V}{\partial F_v}$
    */

    double flash_calculation_solve_RachfordRice_equation(double *K, double *z, double n_V0, int ncomp)
    {
        int itr = 0;
        double n_V = n_V0, F, J, d;

        while (1) {
            F = flash_calculation_calculate_RachfordRice_equation_value(K, z, n_V, ncomp);

            if (fabs(F) < 1e-10)
                break;

            J = flash_calculation_calculate_RachfordRice_equation_derivative(K, z, n_V, ncomp);

            d = -F / J;
            n_V += d;

            if (n_V < 0.0) {
                n_V -= d;
                n_V *= 0.5;
            }

            if (n_V > 1.0) {
                n_V -= d;
                n_V = (1.0 + n_V) * 0.5;
            }

            itr += 1;
            if (itr > 100)
                break;
        }

        return n_V;
    }

    /* ### Update K using SS method
    # $$
    # K_i = K_i \frac{\phi_{l,i}}{\phi_{v,i}}
    # $$
    */
    void flash_calculation_SS_method_update_K(double *fug_L, double *fug_V, double *K, int ncomp)
    {
        int i;

        for (i = 0; i < ncomp; i++) {
            K[i] = K[i] * fug_L[i] / fug_V[i];
        }
    }


    /*
    At range [0, 1.0] find a value which makes Rachford-Rich is zero
    using bisection method
    */
    double flash_calculation_two_phase_flash_calculation_calculate_initial_F(double *K, double *z, int ncomp)
    {
        double tol = 1e-3, F_min, F_max,
            V_min, V_max, F_mid, V_mid;
        int itr;

        F_min = 0.0;
        F_max = 1.0;
        V_min = flash_calculation_calculate_RachfordRice_equation_value(K, z, F_min, ncomp);
        V_max = flash_calculation_calculate_RachfordRice_equation_value(K, z, F_max, ncomp);
        F_mid = 0.0;
        itr = 0;

        while (1) {
            /* calculation the value at bisection position */
            F_mid = (F_min + F_max) * 0.5;
            V_mid = flash_calculation_calculate_RachfordRice_equation_value(K, z, F_mid, ncomp);

            /* if Value is zero, break */
            if (fabs(V_mid) < 1e-10)
                break;

            /* # If V_mid has different signal from V_min,
            #  set F_mid as F_max and V_mid as V_max */
            if (V_mid * V_min < 0.0) {
                V_max = V_mid;
                F_max = F_mid;
            }

            /* # If V_mid has different signal from V_max,
            #  set F_mid as F_min and V_mid as V_min */
            if (V_mid * V_max < 0.0) {
                V_min = V_mid;
                F_min = F_mid;
            }

            if (fabs(F_min - F_max) < tol)
                break;

            itr += 1;
            if (itr > 100)
                break;
        }

        return F_mid;
    }

    /* ### QNSS method for Two-phase flash calculation
    # Two-phase flahs calculation requires the solution of the following
    # equilibrium and material balance equaions:
    # Equilibrium equation:
    # $$
    #             G_i = \log{x_{v,i}} + \log{\phi_{v,i}(x_v)}
    #             - \log{x_{l,i}} - \log{\phi_{l,i}(x_l)}
    # $$
    # Material balance equation:
    # $$
    # V = \sum_i{\frac{z_i (K_i - 1)}{1 + (K_i - 1) F_v}}
    # $$
    # We take $K_i$ and $F_v$ as the primary variables.
    #
    # K-values are solved vis the equilibrium equations, then solve the material balance equation to get the $F_v$. Repeat the process until converge.
    */
    /*
    Two-phase flahs calculation requires the solution of the following
    equilibrium and material balance equaions:
    Equilibrium equation:
    G[i] = ln(x_v[i]) + ln(phi_v[i]) - ln(x_l[i]) - ln(phi_l[i])
    for i = 1,...,Nc
    Material balance equation:
    G[Nc + 1] = Sum_k((z[k] * (K[k] - 1.0)) / (1.0 + F_v * (K[k] - 1.0)))
    We take K[i] and F_v as the primary variables.

    */
    double flash_calculation_two_phase_flash_Calculation_QNSS(double *z,
        double *K, double Fv, double tol)
    {
        int i, ncomp = 4, itr;
        double F_v, sum_K, error, *x_l, *x_v;
        double *G, *dG, *dx_l, *dx_v, *K0;
        PHASE *phase_L, *phase_V;

        /* Initial estimate K */
        K0 = (double *)malloc(ncomp * sizeof(*K0));
        if (Fv <= 0.0) {
            flash_calculation_estimate_K(K0,eos_pressure,eos_temp);
            F_v = 0.5;

#if 0
            printf("KI: \n");
            for (i = 0; i < ncomp; i++) {
                printf("%e ", K0[i]);
            }
            printf("\n");
#endif
        }
        else {
            sum_K = 0.0;
            for (i = 0; i < ncomp; i++) {
                sum_K += log(K[i]) * log(K[i]);
            }

            if (sum_K < 1e-5) {
                flash_calculation_estimate_K(K0, eos_pressure, eos_temp);
                F_v = 0.5;
            }
            else {
                for (i = 0; i < ncomp; i++)
                    K0[i] = K[i];
                F_v = Fv;
            }
        }

#if 0
        {
            double *KK;

            KK = malloc(ncomp * sizeof(*KK));
            flash_calculation_estimate_K(eos, KK);

            printf("Estimate K:\n");
            for (i = 0; i < ncomp; i++) {
                printf("%e ", KK[i]);
            }
            printf("\n");

            free(KK);
        }
#endif

        /* Initial compositions x_v and x_l */
        x_l = (double *)malloc(ncomp * sizeof(*x_l));
        x_v = (double *)malloc(ncomp * sizeof(*x_v));
        flash_calculation_calculate_composition(K0, z, F_v, x_l, x_v, ncomp);

        phase_L = flash_calculation_phase_new(x_l);
        phase_V = flash_calculation_phase_new(x_v);

        G = (double *)malloc(ncomp * sizeof(*G));
        dG = (double *)malloc(ncomp * ncomp * sizeof(*dG));
        dx_l = (double *)malloc(ncomp * sizeof(*dx_l));
        dx_v = (double *)malloc(ncomp * sizeof(*dx_v));

        itr = 0;

        while (1) {
            /* Calculate liquid phase fugacity */
            flash_calculation_compute_phase_parameter(phase_L);
            flash_calculation_calculate_compressibility_factor(phase_L);
            flash_calculation_calculate_fugacity(phase_L);

            /* Calculate vapour phase fugacity */
            flash_calculation_compute_phase_parameter(phase_V);
            flash_calculation_calculate_compressibility_factor(phase_V);
            flash_calculation_calculate_fugacity(phase_V);

            /* Calculate G error */
            error = flash_calculation_calculate_equilibrium_equation(phase_L, phase_V, G);

            /* Check convergence */
            if (error < tol)
                break;

            if (error > 1e-5) {
                flash_calculation_SS_method_update_K(phase_L->fug, phase_V->fug, K0, ncomp);
            }
            else {
                /* Calculate the derivatives */
                flash_calculation_calculate_composition_derivative(K0, z, F_v, dx_l, dx_v, ncomp);
                flash_calculation_calculate_equilibrium_equation_derivative(phase_L, phase_V, dx_l, dx_v, dG);

                /* Update K */
                flash_calculation_QNSS_method_update_K(dG, G, K0, ncomp);
            }

            /* ## Solve Rachford-Rice equation and get F_v */
            F_v = flash_calculation_solve_RachfordRice_equation(K0, z, F_v, ncomp);

            /* ## Calculate compositions */
            flash_calculation_calculate_composition(K0, z, F_v, x_l, x_v, ncomp);

            itr += 1;
            if (itr > 1000) {
                //printf("##### WARNING: two-phase flash calculation reach maximum iterations!\n");
                break;
            }
        }

        if (verb) {
            printf("Pres: %e, Temp: %e, Itr: %d\n", eos_pressure,eos_temp, itr);
        }

        if (fabs(F_v) < 1e-5 || fabs(F_v - 1.0) < 1e-5) {
            split_failure++;
        }
        split_itr += itr;

        flash_calculation_calculate_phase_density(phase_V);
        flash_calculation_calculate_phase_density(phase_L);

        if (phase_V->density > phase_L->density) {
            F_v = 1.0 - F_v;
            for (i = 0; i < ncomp; i++) {
                K[i] = 1.0 / K0[i];
            }
        }
        else {
            for (i = 0; i < ncomp; i++) {
                K[i] = K0[i];
            }
        }

        free(x_l);
        free(x_v);
        free(G);
        free(dG);
        free(dx_l);
        free(dx_v);
        free(K0);

        //flash_calculation_phase_free(&phase_L);
        //flash_calculation_phase_free(&phase_V);
        /* printf("##### Two-phase flash calculation iterations: %d" %itr); */

        return F_v;
    }
    // some functions to return the calculation time
    double flash_calculation_split_time_cost()
    {
        return split_solve_time;
    }

    int flash_calculation_split_iteration_number()
    {
        return split_itr;
    }

    int flash_calculation_split_failure_number()
    {
        return split_failure;
    }

    double flash_calculation_split_pred_time_cost()
    {
        return split_pred_time;
    }



}  // end of namespace
}  // end of namespace
