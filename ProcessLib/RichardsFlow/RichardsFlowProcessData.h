/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_RICHARDSFLOW_RICHARDSFLOWPROCESSDATA_H
#define PROCESSLIB_RICHARDSFLOW_RICHARDSFLOWPROCESSDATA_H

namespace MeshLib
{
    class Element;
}


namespace ProcessLib
{

template <typename T>
struct Parameter;

namespace RichardsFlow
{

struct RichardsFlowProcessData
{
    RichardsFlowProcessData(
        Parameter<double> const& intrinsic_permeability_,
		Parameter<double> const& porosity_,
		Parameter<double> const& viscosity_,
		bool const has_gravity_,
		std::map<std::string,
		std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
		curves_)
        : intrinsic_permeability(intrinsic_permeability_)
		, porosity(porosity_)
		, viscosity(viscosity_)
		, has_gravity(has_gravity_)
		, curves(curves_)
    {}

    RichardsFlowProcessData(RichardsFlowProcessData&& other)
        : intrinsic_permeability(other.intrinsic_permeability)
		, porosity(other.porosity)
		, viscosity(other.viscosity)
		, has_gravity(other.has_gravity)
		, curves(other.curves)
    {}

    //! Copies are forbidden.
    RichardsFlowProcessData(RichardsFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(RichardsFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(RichardsFlowProcessData&&) = delete;

    Parameter<double> const& intrinsic_permeability;
	Parameter<double> const& porosity;
	Parameter<double> const& viscosity;
	std::map<std::string,
		std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
		curves;
	bool const has_gravity;
};

} // namespace RichardsFlow
} // namespace ProcessLib

#endif // PROCESSLIB_RichardsFlow_RichardsFlowPROCESSDATA_H
