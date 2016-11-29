/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_TWOPHASEFLOWWITHPX_EOSBASE_H_
#define OGS_TWOPHASEFLOWWITHPX_EOSBASE_H_

#include <memory>

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
	class SpatialPosition;
namespace TwoPhaseFlowWithPX
{
/// Interface for mechanical solid material models. Provides updates of the
/// stress for a given current state and also a tangent at that position. If the
/// implemented material model stores an internal state, the nested
/// MaterialStateVariables class should be used; it's only responsibility is to
/// provide state's push back possibility.
class EoSBase
{
public:
    /// constitutive relation compute function.
    /// Returns false in case of errors in the computation if Newton iterations
    /// did not converge, for example.
    virtual bool computeConstitutiveRelation(
        double const t,
        SpatialPosition const& x,
        double const PG,
		double const X,
		double& Sw,
		double& X_m,
		double& dsw_dpg,
		double& dsw_dX,
		double& dxm_dpg,
		double& dxm_dX) = 0;

    virtual ~EoSBase() = default;
};

}  // namespace TwoPhaseFlowWithPX
}  // namespace MaterialLib

#endif  // OGS_TWOPHASEFLOWWITHPX_EOSBASE_H_
