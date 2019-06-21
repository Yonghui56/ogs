/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   FinesTransportLocalAssemblerInterface.h
 *  Created on October 11, 2017, 1:35 PM
 */

#pragma once

#include <Eigen/Dense>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
struct CoupledSolutionsForStaggeredScheme;

namespace FinesTransport
{
template <typename NodalRowVectorType,
          typename GlobalDimNodalMatrixType>
struct IntegrationPointData final
{
    IntegrationPointData(NodalRowVectorType N_,
                         GlobalDimNodalMatrixType dNdx_,
                         double const& integration_weight_
                        )
        : N(std::move(N_)),
          dNdx(std::move(dNdx_)),
          integration_weight(integration_weight_),
          porosity_curr(0.3),
          porosity_prev(0.3),
          permeability_curr(4.4e-10),
          permeability_prev(4.4e-10)
    {
    }

    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;
    double porosity_curr, porosity_prev;
    double permeability_curr, permeability_prev;
    double pore_throat_concentration;
    double pore_body_concentration;
    double permeability_ratio;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    void push_back()
    {
        porosity_prev = porosity_curr;
        permeability_prev = permeability_curr;
    }
};

class FinesTransportLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    FinesTransportLocalAssemblerInterface() = default;
    void setStaggeredCoupledSolutions(
        std::size_t const /*mesh_item_id*/,
        CoupledSolutionsForStaggeredScheme* const coupling_term)
    {
        _coupled_solutions = coupling_term;
    }

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtPoreThroatConcentration(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtPoreBodyConcentration(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtPermeabilityRatio(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

    Eigen::Vector3d getFlux(MathLib::Point3d const& pnt_local_coords,
                            double const t,
                            std::vector<double> const& local_x) const override =
        0;

protected:
    // TODO: remove _coupled_solutions or move integration point data from
    // local assembler class to a new class to make local assembler unique
    // for each process.
    /** Pointer to CoupledSolutionsForStaggeredScheme that is set in a
     * member of Process class,
     * setCoupledTermForTheStaggeredSchemeToLocalAssemblers. It is used for
     * calculate the secondary variables like velocity for coupled
     * processes.
     */
    CoupledSolutionsForStaggeredScheme* _coupled_solutions{nullptr};
};

}  // namespace FinesTransport
}  // namespace ProcessLib
