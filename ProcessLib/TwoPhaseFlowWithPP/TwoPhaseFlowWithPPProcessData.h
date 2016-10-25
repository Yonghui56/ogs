/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_TWOPHASEFLOWWITHPP_TWOPHASEFLOWWITHPPPROCESSDATA_H
#define PROCESSLIB_TWOPHASEFLOWWITHPP_TWOPHASEFLOWWITHPPPROCESSDATA_H

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace TwoPhaseFlowWithPP
{
struct TwoPhaseFlowWithPPProcessData
{
	TwoPhaseFlowWithPPProcessData(
        Parameter<double> const& specific_body_force_,
        bool const has_gravity_,
        bool const has_mass_lumping_)
        : specific_body_force(specific_body_force_),
          has_gravity(has_gravity_),
          has_mass_lumping(has_mass_lumping_)
    {
    }

	TwoPhaseFlowWithPPProcessData(TwoPhaseFlowWithPPProcessData&& other)
        : specific_body_force(other.specific_body_force),
          has_gravity(other.has_gravity),
          has_mass_lumping(other.has_mass_lumping)
    {
    }

    //! Copies are forbidden.
	TwoPhaseFlowWithPPProcessData(TwoPhaseFlowWithPPProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseFlowWithPPProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseFlowWithPPProcessData&&) = delete;

    Parameter<double> const& specific_body_force;
    bool const has_gravity;
    bool const has_mass_lumping;
};

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib

#endif  // PROCESSLIB_TwoPhaseFlowWithPP_TwoPhaseFlowWithPPPROCESSDATA_H
