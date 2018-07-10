/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once
#include <math.h>

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
    typedef struct PHASE_ {
        int ncomp;
        double *mf;

        double R;

        double density;

        double A;
        double dAp;
        double dA_dT;
        double dA_dT2;
        double *dAx;

        double B;
        double dBp;
        double dB_dT;
        double dB_dT2;
        double *dBx;

        int nroot;
        double Z;
        double dZ;
        double dZ_dT;
        double dZ_dT2;
        double *dZ_dx;

        int phase_no;

        double *fug;
        double *phi;
        double *dphi;
        double *dphi_dT;
        double *dphi_dx;
        double *dfug;

        double *ai;
        double *dai_dT;
        double *dai_dT2;
        double a;
        double da_dT;
        double da_dT2;
        double *da;
        double *dda_dT;
        double *dda_dT2;

        double *bi;
        double b;
        double *db;
    } PHASE;
    double * flash_calculation_estimate_K(double *K, double P, double T);
    double * flash_calculation_stability_analysis_initial_estimate(double P, double T);
    void flash_calculation_calculate_trial_phase_composition(double *X_t, double *x);
    void flash_calculation_SS_method_update_X(PHASE *phase, PHASE *phase_t, double *X_t);
    void flash_calculation_calculate_trial_phase_composition_derivative(double *X_t, double *dx_t);
    void flash_calculation_calculate_stability_equilibrium_equation(PHASE *phase, PHASE *phase_t, double *X_t, double *D);
    void flash_calculation_calculate_stability_equilibrium_equation_derivative(PHASE *phase_t, double *dx_t,
        double *X_t, double *dD);
    int flash_calculation_stability_analysis_QNSS(PHASE *phase, double *K, double tol);
    void flash_calculation_compute_phase_parameter(PHASE *phase);

    void flash_calculation_calculate_compressibility_factor(PHASE *phase);
    void flash_calculation_calculate_fugacity(PHASE *phase);
    void flash_calculation_calculate_compressibility_factor(PHASE *phase);

    int flash_calculation_solve_dense_linear_system_LU(int n, double *a0, int pvt[]);
    void flash_calculation_solve_dense_linear_system_SV(int n, double *a0,
        int pvt[], int m, double *b0);
    int flash_calculation_solve_dense_linear_system(double *M, double *b,
        double *x, int n);
    PHASE * flash_calculation_phase_new(double *mf);
    void flash_calculation_calculate_phase_density(PHASE *phase);
    double const M_PI = 3.1415926;
    double const comp_PC[4] = { 46,73.8,89.4,220.5 };
    double const comp_TC[4] = { 190.6,204.2,373.2,647.3 };
    double const comp_AC[4] = { 0.008,0.225,0.1,0.344 };
    double const comp_binary[4][4] = {
        { 0, 0, 0, 0 },
    { 0.1005,0,0,0 },
    { 0.0755,0.0999,0,0 },
    { 0.4928,0,0.04,0 }
    };
    double const comp_molarweight[4] = {0.016,0.044,0.034,0.018};
    double const eos_temp = 310.95;
    double const eos_pressure = 76;
    double const z_molarfraction[4] = { 0.1488,0.2991,0.0494,0.5027 };
    //define a struct "phase"
    
}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
