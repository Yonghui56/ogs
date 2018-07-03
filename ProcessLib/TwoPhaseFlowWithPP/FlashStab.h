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
    double const M_PI = 3.1415926;
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
}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
