/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_GENERICNATURALBOUNDARYCONDITIONLOCALASSEMBLER_H
#define PROCESSLIB_GENERICNATURALBOUNDARYCONDITIONLOCALASSEMBLER_H

#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
class GenericNaturalBoundaryConditionLocalAssemblerInterface
{
public:
    virtual ~GenericNaturalBoundaryConditionLocalAssemblerInterface() = default;

    virtual void assemble(
        std::size_t const id,
        NumLib::LocalToGlobalIndexMap const& dof_table_boundary, double const t,
        const GlobalVector& x, GlobalMatrix& K, GlobalVector& b) = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class GenericNaturalBoundaryConditionLocalAssembler
    : public GenericNaturalBoundaryConditionLocalAssemblerInterface
{
protected:
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;

public:
    GenericNaturalBoundaryConditionLocalAssembler(
        MeshLib::Element const& e, unsigned const integration_order)
        : _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              e, integration_order)),
          _integration_order(integration_order)
    {
    }

protected:
    std::vector<typename ShapeMatricesType::ShapeMatrices> _shape_matrices;
    unsigned const _integration_order;
};

}  // ProcessLib

#endif  // PROCESSLIB_GENERICNATURALBOUNDARYCONDITIONLOCALASSEMBLER_H
