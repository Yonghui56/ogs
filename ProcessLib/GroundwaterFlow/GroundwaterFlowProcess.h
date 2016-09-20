/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
#define PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_

#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"
#include "GroundwaterFlowFEM.h"
#include "GroundwaterFlowProcessData.h"


namespace ProcessLib
{
namespace GroundwaterFlow
{
class GroundwaterFlowProcess final : public Process
{
    using Base = Process;

public:
    GroundwaterFlowProcess(
        MeshLib::Mesh& mesh,
        Base::NonlinearSolver& nonlinear_solver,
        std::unique_ptr<Base::TimeDiscretization>&& time_discretization,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        GroundwaterFlowProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        ProcessOutput&& process_output);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return true; }
    //! @}

private:
    using GlobalAssembler = NumLib::VectorMatrixAssembler<
        GroundwaterFlowLocalAssemblerInterface,
        NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    using ExtrapolatorInterface =
        NumLib::Extrapolator<IntegrationPointValue,
                             GroundwaterFlowLocalAssemblerInterface>;
    using ExtrapolatorImplementation =
        NumLib::LocalLinearLeastSquaresExtrapolator<
            IntegrationPointValue, GroundwaterFlowLocalAssemblerInterface>;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    GroundwaterFlowProcessData _process_data;

    std::unique_ptr<GlobalAssembler> _global_assembler;
    std::vector<std::unique_ptr<GroundwaterFlowLocalAssemblerInterface>>
        _local_assemblers;

    std::unique_ptr<ExtrapolatorInterface> _extrapolator;
};

}   // namespace GroundwaterFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
