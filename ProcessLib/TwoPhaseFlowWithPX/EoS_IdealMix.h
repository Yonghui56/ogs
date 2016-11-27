/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_TWOPHASEWITHPX_EOS_IDEALMIX_H_
#define OGS_TWOPHASEWITHPX_EOS_IDEALMIX_H_

#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"
#include "EoSBase.h"
#include "MaterialLib/PhysicalConstant.h"

namespace ProcessLib
{
    class SpatialPosition;
namespace TwoPhaseFlowWithPX
{
class EoS_IdealMix final : public EoSBase
{
public:
	//
	// Variables specific to the material model.
	//
	struct EoS_IdealMix_Properties
	{
		using P = ProcessLib::Parameter<double>;
		EoS_IdealMix_Properties(P const& henry_const_)
			: _henry_const(henry_const_)
		{
		}

		// basic eos parameters
		P const& _henry_const;
	};

public:

    static int const JacobianResidualSize = 2;
    using ResidualVector = Eigen::Matrix<double, JacobianResidualSize, 1>;
    using JacobianMatrix = Eigen::Matrix<double,
                                         JacobianResidualSize,
                                         JacobianResidualSize,
                                         Eigen::RowMajor>;
	using UnknownVector= Eigen::Matrix<double, JacobianResidualSize, 1>;

public:
	explicit EoS_IdealMix(EoS_IdealMix_Properties& eos_properties)
		: _eos_prop(eos_properties)
	{
	}

    bool computeConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x_position,
        double const PG,
        double const X,
        double& Sw,
        double& X_m
        ) override;

private:
    /// Calculates the 18x1 residual vector.
    void calculateResidual(double const PG, double const X, UnknownVector vec_unknown,
                                  ResidualVector& res);

    /// Calculates the 18x18 Jacobian.
    void calculateJacobian(double const t,
                                  ProcessLib::SpatialPosition const& x,
                                  double const PG, double const X,
                                  JacobianMatrix& Jac,UnknownVector& vec_unknown);

    /**
    * Complementary condition 1
    * for calculating molar fraction of light component in the liquid phase
    */
    const double Calc_equili_molfrac_X_m(double const PG, double const Sw, double const X_m)
    {	
        const double N_L = rho_mol_h2o / (1 - X_m);
        double const x_equili_m = PG * Hen / N_L;
        return std::min(1 - Sw, x_equili_m - X_m);
    }
	/**
	* Complementary condition 3
	* for calculating the saturation
	*/
	const double Calc_Saturation(double PG, double X, double Sw, double X_m,
		double X_M, double T)
	{
		const double N_G = PG / 8.314 / T;
		const double N_L = rho_mol_h2o / (1 - X_m);
		return ((1 - Sw) * N_G * (X - X_M) + Sw * N_L * (X - X_m)) / ((1 - Sw) * N_G + Sw * N_L);
	}

private:
	EoS_IdealMix_Properties _eos_prop;
    double const Hen = 7.65e-6;  // mol/Pa./m3
                                 /**
                                 * Molar mass of water
                                 */
    double const molar_mass_h2o =
        MaterialLib::PhysicalConstant::MolarMass::Water;
    /**
    * mass density of water
    */
    double const rho_mass_h20 = 1000;
    /**
    * molar density of water
    */
    double const rho_mol_h2o = rho_mass_h20 / molar_mass_h2o;
};



}  // namespace Solids
}  // namespace MaterialLib

#include "EoS_IdealMix-impl.h"

#endif  // OGS_TWOPHASEWITHPX_EOS_IDEALMIX_H_
