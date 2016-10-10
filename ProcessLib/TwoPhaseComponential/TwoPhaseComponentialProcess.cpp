/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TwoPhaseComponentialProcess.h"

#include <cassert>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace TwoPhaseComponential
{
TwoPhaseComponentialProcess::TwoPhaseComponentialProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
	TwoPhaseComponentialProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              std::move(process_variables), std::move(secondary_variables),
              std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
}

void TwoPhaseComponentialProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    //ProcessLib::createLocalAssemblers<LocalAssemblerData>(
       // mesh.getDimension(), mesh.getElements(), dof_table, integration_order,
       // _local_assemblers, _process_data);

	ProcessLib::createLocalAssemblers<LocalAssemblerData>(
		mesh.getDimension(), mesh.getElements(), dof_table, _local_assemblers,
		mesh.isAxiallySymmetric(), integration_order, _process_data);

	_secondary_variables.addSecondaryVariable(
		"saturation", 1,
		makeExtrapolator(getExtrapolator(), _local_assemblers,
			&TwoPhaseComponentialLocalAssemblerInterface::getIntPtSat));

	_secondary_variables.addSecondaryVariable(
		"liquid_pressure", 1,
		makeExtrapolator(getExtrapolator(), _local_assemblers,
			&TwoPhaseComponentialLocalAssemblerInterface::getIntPtPL));


    
}

void TwoPhaseComponentialProcess::assembleConcreteProcess(const double t,
                                                    GlobalVector const& x,
                                                    GlobalMatrix& M,
                                                    GlobalMatrix& K,
                                                    GlobalVector& b)
{
    DBUG("Assemble TwoPhaseComponentialProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b);
}

void TwoPhaseComponentialProcess::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian TwoPhaseComponentialProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac);
}

}  // namespace TwoPhaseComponential
}  // namespace ProcessLib
