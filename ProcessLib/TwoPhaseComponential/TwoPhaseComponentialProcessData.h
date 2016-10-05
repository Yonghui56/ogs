/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_TWOPHASECOMPONENTIAL_TWOPHASECOMPONENTIALPROCESSDATA_H
#define PROCESSLIB_TWOPHASECOMPONENTIAL_TWOPHASECOMPONENTIALPROCESSDATA_H

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace TwoPhaseComponential
{

struct TwoPhaseComponentialProcessData
{
    TwoPhaseComponentialProcessData(Parameter<double> const& porosity_,
                              Parameter<double> const& intrinsic_permeability_,
                              Parameter<double> const& water_density_,
		                      Parameter<double> const& mat_id_,
		                      Parameter<double> const& viscosity_gas_, 
		                      Parameter<double> const& viscosity_liquid_,
		                      bool const has_gravity_,
		                      std::map<std::string,
		                               std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
		                                curves_)
        : porosity(porosity_),
          intrinsic_permeability(intrinsic_permeability_),
          water_density(water_density_)
		, mat_id(mat_id_)
		, viscosity_gas(viscosity_gas_)
		, viscosity_liquid(viscosity_liquid_)
		, has_gravity(has_gravity_)
		, curves(curves_)
    {
    }

    TwoPhaseComponentialProcessData(TwoPhaseComponentialProcessData&& other)
        : porosity(other.porosity)
		, intrinsic_permeability(other.intrinsic_permeability)
		, water_density(other.water_density)
		, mat_id(other.mat_id)
		, viscosity_gas(other.viscosity_gas)
		, viscosity_liquid(other.viscosity_liquid)
		, has_gravity(other.has_gravity)
		, curves(other.curves)
    {
    }

    //! Copies are forbidden.
    TwoPhaseComponentialProcessData(TwoPhaseComponentialProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseComponentialProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseComponentialProcessData&&) = delete;

    Parameter<double> const& porosity;
    Parameter<double> const& intrinsic_permeability;
    Parameter<double> const& water_density;
	Parameter<double> const& mat_id;
	Parameter<double> const& viscosity_gas;
	Parameter<double> const& viscosity_liquid;
	std::map<std::string,
		std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
		curves;
	bool const has_gravity;
};

}  // namespace TwoPhaseComponential
}  // namespace ProcessLib

#endif  // PROCESSLIB_TwoPhaseComponential_TwoPhaseComponentialPROCESSDATA_H
