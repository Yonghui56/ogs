/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TwoPhaseComponentialFlowProcess.h"

#include <cassert>
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

#include "TwoPhaseComponentialFlowLocalAssembler.h"

namespace ProcessLib
{
namespace TwoPhaseComponentialFlow
{
    TwoPhaseComponentialFlowProcess::TwoPhaseComponentialFlowProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    TwoPhaseComponentialFlowProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    BaseLib::ConfigTree const& /*config*/,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
    /*curves*/)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
    DBUG("Create Two Phase Componential Flow Process.");
}

void TwoPhaseComponentialFlowProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::ProcessVariable const& pv = getProcessVariables()[0];
    ProcessLib::createLocalAssemblers<TwoPhaseComponentialFlowLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    _secondary_variables.addSecondaryVariable(
        "saturation", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtSaturation));

    _secondary_variables.addSecondaryVariable(
        "pressure_wetting", 1,
        makeExtrapolator(getExtrapolator(), _local_assemblers,
                         &TwoPhaseComponentialFlowLocalAssemblerInterface::
                             getIntPtWettingPressure));
    _secondary_variables.addSecondaryVariable(
        "mol_frac_nonwet_vapor", 1,
        makeExtrapolator(getExtrapolator(), _local_assemblers,
            &TwoPhaseComponentialFlowLocalAssemblerInterface::
            getIntPtMolFracNonwetVapor));
    _secondary_variables.addSecondaryVariable(
        "co2_concentration", 1,
        makeExtrapolator(getExtrapolator(), _local_assemblers,
            &TwoPhaseComponentialFlowLocalAssemblerInterface::
            getIntPtCO2Concentration));
    _secondary_variables.addSecondaryVariable(
        "porosity", 1,
        makeExtrapolator(getExtrapolator(), _local_assemblers,
            &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtPorosityValue));
    _secondary_variables.addSecondaryVariable(
        "pH_value", 1,
        makeExtrapolator(getExtrapolator(), _local_assemblers,
            &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtpHValue));
    _secondary_variables.addSecondaryVariable(
        "co2_cumulated_prev", 1,
        makeExtrapolator(getExtrapolator(), _local_assemblers,
            &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtRhoMolCo2CumulTotalPrev));

}

void TwoPhaseComponentialFlowProcess::assembleConcreteProcess(const double t,
                                                          GlobalVector const& x,
                                                          GlobalMatrix& M,
                                                          GlobalMatrix& K,
                                                          GlobalVector& b,
                                                          StaggeredCouplingTerm
                                                          const& coupling_term)
{
    DBUG("Assemble TwoPhaseFlowWithPrhoProcess.");
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b, coupling_term);
}

void TwoPhaseComponentialFlowProcess::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac,
    StaggeredCouplingTerm const& coupling_term)
{
    DBUG("AssembleWithJacobian TwoPhaseFlowWithPrhoProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac, coupling_term);
}

}  // end of namespace
}  // end of namespace