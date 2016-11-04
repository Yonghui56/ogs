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
        bool const has_mass_lumping_,
		Parameter<double> const& henry_const_,
		Parameter<double> const& diffusion_coeff_componentb_,
		Parameter<double> const& diffusion_coeff_componenta_,
        std::unique_ptr<MaterialLib::TwoPhaseFlowWithPP::
                            TwoPhaseFlowWithPPMaterialProperties>&& material_)
        : _specific_body_force(specific_body_force_),
          _has_gravity(has_gravity_),
          _has_mass_lumping(has_mass_lumping_),
		  _henry_const(henry_const_),
		  _diffusion_coeff_componentb(diffusion_coeff_componentb_),
		  _diffusion_coeff_componenta(diffusion_coeff_componenta_),
          _material(std::move(material_))
    {
    }

    TwoPhaseFlowWithPPProcessData(TwoPhaseFlowWithPPProcessData&& other)
        : _specific_body_force(other._specific_body_force),
          _has_gravity(other._has_gravity),
          _has_mass_lumping(other._has_mass_lumping),
		_henry_const(other._henry_const),
		_diffusion_coeff_componentb(other._diffusion_coeff_componentb),
		_diffusion_coeff_componenta(other._diffusion_coeff_componenta),
          _material(std::move(other._material))
    {
    }

    //! Copies are forbidden.
    TwoPhaseFlowWithPPProcessData(TwoPhaseFlowWithPPProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseFlowWithPPProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseFlowWithPPProcessData&&) = delete;

    Parameter<double> const& _specific_body_force;
    bool const _has_gravity;
    bool const _has_mass_lumping;
	Parameter<double> const& _henry_const;
	Parameter<double> const& _diffusion_coeff_componentb;
	Parameter<double> const& _diffusion_coeff_componenta;

    std::unique_ptr<
        MaterialLib::TwoPhaseFlowWithPP::TwoPhaseFlowWithPPMaterialProperties>
        _material;

    
};

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib

#endif  // PROCESSLIB_TwoPhaseFlowWithPP_TwoPhaseFlowWithPPPROCESSDATA_H
