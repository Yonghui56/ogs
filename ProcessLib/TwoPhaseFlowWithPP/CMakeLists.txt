AddTest(
    NAME 2D_TwoPhase_PP_Lia_quad
    PATH Parabolic/TwoPhaseFlowPP/Liakopoulos
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_Lia_quad_short.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-1 RELTOL 1e-2
    DIFF_DATA
    Lia_20.vtu twophaseflow_pcs_0_ts_218_t_20.000000.vtu capillary_pressure capillary_pressure
    Lia_20.vtu twophaseflow_pcs_0_ts_218_t_20.000000.vtu gas_pressure gas_pressure
    Lia_20.vtu twophaseflow_pcs_0_ts_218_t_20.000000.vtu saturation saturation
)
AddTest(
    NAME LARGE_2D_TwoPhase_PP_Lia_quad
    PATH Parabolic/TwoPhaseFlowPP/Liakopoulos
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_Lia_quad_large.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-1 RELTOL 1e-2
    DIFF_DATA
    h2_Liako_1198.vtu twophaseflow_pcs_0_ts_1198_t_1000.000000.vtu capillary_pressure capillary_pressure
    h2_Liako_1198.vtu twophaseflow_pcs_0_ts_1198_t_1000.000000.vtu gas_pressure gas_pressure
    h2_Liako_1198.vtu twophaseflow_pcs_0_ts_1198_t_1000.000000.vtu saturation saturation
)
