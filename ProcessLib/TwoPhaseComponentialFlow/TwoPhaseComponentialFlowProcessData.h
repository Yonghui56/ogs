/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TwoPhaseComponentialFlowMaterialProperties.h"
namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace TwoPhaseComponentialFlow
{
struct TwoPhaseComponentialFlowProcessData
{
    TwoPhaseComponentialFlowProcessData(
        Eigen::VectorXd const specific_body_force_,
        bool const has_gravity_,
        bool const has_mass_lumping_,
        Parameter<double> const& diffusion_coeff_component_b_,
        Parameter<double> const& diffusion_coeff_component_a_,
        Parameter<double> const& temperature_,
        std::unique_ptr<TwoPhaseComponentialFlowMaterialProperties>&& material_,
        MathLib::PiecewiseLinearInterpolation const& interpolated_Q_slow_,
        MathLib::PiecewiseLinearInterpolation const& interpolated_Q_fast_,
        MathLib::PiecewiseLinearInterpolation const& interpolated_kinetic_rate_
        )
        : _specific_body_force(specific_body_force_),
          _has_gravity(has_gravity_),
          _has_mass_lumping(has_mass_lumping_),
          _diffusion_coeff_component_b(diffusion_coeff_component_b_),
          _diffusion_coeff_component_a(diffusion_coeff_component_a_),
          _temperature(temperature_),
          _material(std::move(material_)),
          _interpolated_Q_slow(interpolated_Q_slow_),
          _interpolated_Q_fast(interpolated_Q_fast_),
          _interpolated_kinetic_rate(interpolated_kinetic_rate_)


    {
    }

    TwoPhaseComponentialFlowProcessData(TwoPhaseComponentialFlowProcessData&& other)
        : _specific_body_force(other._specific_body_force),
          _has_gravity(other._has_gravity),
          _has_mass_lumping(other._has_mass_lumping),
          _diffusion_coeff_component_b(other._diffusion_coeff_component_b),
          _diffusion_coeff_component_a(other._diffusion_coeff_component_a),
          _temperature(other._temperature),
          _material(std::move(other._material)),
        _interpolated_Q_slow(other._interpolated_Q_slow),
        _interpolated_Q_fast(other._interpolated_Q_fast),
        _interpolated_kinetic_rate(other._interpolated_kinetic_rate),
        _dt{ other._dt }
    {
    }

    //! Copies are forbidden.
    TwoPhaseComponentialFlowProcessData(TwoPhaseComponentialFlowProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseComponentialFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseComponentialFlowProcessData&&) = delete;
    Eigen::VectorXd const _specific_body_force;

    bool const _has_gravity;
    bool const _has_mass_lumping;
    Parameter<double> const& _diffusion_coeff_component_b;
    Parameter<double> const& _diffusion_coeff_component_a;
    Parameter<double> const& _temperature;
    std::unique_ptr<TwoPhaseComponentialFlowMaterialProperties> _material;
    MathLib::PiecewiseLinearInterpolation const& _interpolated_Q_slow;
    MathLib::PiecewiseLinearInterpolation const& _interpolated_Q_fast;
    MathLib::PiecewiseLinearInterpolation const& _interpolated_kinetic_rate;
    double _dt = 0;
};

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib