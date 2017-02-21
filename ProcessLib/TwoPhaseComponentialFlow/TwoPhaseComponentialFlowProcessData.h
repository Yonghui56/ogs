/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_TWOPHASECOMPONENTIALFLOW_TWOPHASECOMPONENTIALFLOWPROCESSDATA_H
#define PROCESSLIB_TWOPHASECOMPONENTIALFLOW_TWOPHASECOMPONENTIALFLOWPROCESSDATA_H
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
        std::unique_ptr<TwoPhaseComponentialFlowMaterialProperties>&& material_,
        Parameter<double> const& diffusion_coeff_componentb_,
        Parameter<double> const& diffusion_coeff_componenta_,
        MathLib::PiecewiseLinearInterpolation const& interpolated_Q_slow_,
        MathLib::PiecewiseLinearInterpolation const& interpolated_Q_fast_)
        : _specific_body_force(specific_body_force_),
          _has_gravity(has_gravity_),
          _has_mass_lumping(has_mass_lumping_),
          _material(std::move(material_)),
          _diffusion_coeff_componentb(diffusion_coeff_componentb_),
          _diffusion_coeff_componenta(diffusion_coeff_componenta_),
          _interpolated_Q_slow(interpolated_Q_slow_),
          _interpolated_Q_fast(interpolated_Q_fast_)
    {
    }

    TwoPhaseComponentialFlowProcessData(
        TwoPhaseComponentialFlowProcessData&& other)
        : _specific_body_force(other._specific_body_force),
          _has_gravity(other._has_gravity),
          _has_mass_lumping(other._has_mass_lumping),
          _material(std::move(other._material)),
          _diffusion_coeff_componentb(other._diffusion_coeff_componentb),
          _diffusion_coeff_componenta(other._diffusion_coeff_componenta),
          _interpolated_Q_slow(other._interpolated_Q_slow),
          _interpolated_Q_fast(other._interpolated_Q_fast),
          _dt{other._dt}
    {
    }

    //! Copies are forbidden.
    TwoPhaseComponentialFlowProcessData(
        TwoPhaseComponentialFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseComponentialFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseComponentialFlowProcessData&&) = delete;

    //! Specific body forces applied to solid and fluid.
    //! It is usually used to apply gravitational forces.
    //! A vector of displacement dimension's length.
    Eigen::VectorXd const _specific_body_force;

    bool const _has_gravity;

    //! Enables lumping of the mass matrix.
    bool const _has_mass_lumping;
    std::unique_ptr<TwoPhaseComponentialFlowMaterialProperties> _material;
    Parameter<double> const& _diffusion_coeff_componentb;
    Parameter<double> const& _diffusion_coeff_componenta;
    MathLib::PiecewiseLinearInterpolation const& _interpolated_Q_slow;
    MathLib::PiecewiseLinearInterpolation const& _interpolated_Q_fast;
    double _dt = 0;
};

}  // namespace TwoPhaseComponentialFlow
}  // namespace ProcessLib

#endif  // PROCESSLIB_TwoPhaseComponentialFlow_TwoPhaseComponentialFlowPROCESSDATA_H
