/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   FlashSplit.h
 *
 * Created on August 19, 2016, 1:38 PM
 */

#pragma once

#include "ProcessLib\TwoPhaseFlowWithPP\FlashStab.h"
#include <memory>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Process.h"

namespace MeshLib
{
class Element;
class Mesh;
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
    void flash_calculation_calculate_composition(double *K, double *z, double F_v,
        double *x_l, double *x_v, int ncomp);
    void flash_calculation_calculate_composition_derivative(double *K, double *z, double F_v,
        double *dx_l, double *dx_v, int ncomp);
    double flash_calculation_calculate_equilibrium_equation(PHASE *phase_L, PHASE *phase_V, double *G);
    void flash_calculation_calculate_equilibrium_equation_derivative(PHASE *phase_L, PHASE *phase_V,
        double *dx_l, double *dx_v, double *dG);
    void flash_calculation_QNSS_method_update_K(double *dG, double *G, double *K, int ncomp);
    double flash_calculation_calculate_RachfordRice_equation_value(double *K, double *z, double n_V, int ncomp);
    double flash_calculation_calculate_RachfordRice_equation_derivative(double *K, double *z, double n_V, int ncomp);
    double flash_calculation_solve_RachfordRice_equation(double *K, double *z, double n_V0, int ncomp);
    void flash_calculation_SS_method_update_K(double *fug_L, double *fug_V, double *K, int ncomp);
    double flash_calculation_two_phase_flash_calculation_calculate_initial_F(double *K, double *z, int ncomp);
    double flash_calculation_two_phase_flash_Calculation_QNSS(double *z,
        double *K, double Fv, double tol);


    double flash_calculation_split_time_cost(void);
    int flash_calculation_split_iteration_number(void);
    int flash_calculation_split_failure_number(void);
    double flash_calculation_split_pred_time_cost(void);




}  // end of namespace
}  // end of namespace
