/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_TWOPHASEWITHPRHO_EOS_IDEALMIX_H_
#define OGS_TWOPHASEWITHPRHO_EOS_IDEALMIX_H_

#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"
#include "EoSBase.h"
#include "MaterialLib/PhysicalConstant.h"

namespace ProcessLib
{
    class SpatialPosition;
namespace TwoPhaseFlowWithPrho
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
        double& X_m,
		double& dsw_dpg,
		double& dsw_dX,
		double& dxm_dpg,
		double& dxm_dX
        ) override;

private:
    /// Calculates the 18x1 residual vector.
    void calculateResidual(double const PG, double const X, double Sw, double X_m,
                                  ResidualVector& res);

    /// Calculates the 18x18 Jacobian.
    void calculateJacobian(double const t,
                                  ProcessLib::SpatialPosition const& x,
                                  double const PG, double const X,
                                  JacobianMatrix& Jac, double Sw, double X_m);

    /**
    * Complementary condition 1
    * for calculating molar fraction of light component in the liquid phase
    */
    const double Calc_equili_rho_m(double const PG, double const Sw, double const rho_h2_wet)
    {	
        double const rho_equili_h2_wet = PG * Hen * molar_mass_h2;
        return std::min(1 - Sw, rho_equili_h2_wet - rho_h2_wet);
    }
	/**
	* Complementary condition 2
	* for calculating the saturation
	*/
	const double Calc_Saturation(double PG, double X, double Sw, double rho_h2_wet,
		double rho_h2_nonwet, double T)
	{
		return X - (Sw*rho_h2_wet + (1-Sw)*rho_h2_nonwet);
	}
	/**
	* Calculate the derivatives using the analytical way
	*/
	const double Calculate_dSwdP(double PG, double S, double rho_h2_wet,
		double T = 303.15)
	{
		double const rho_equili_h2_wet = PG * Hen * molar_mass_h2;
		const double rho_h2_nonwet = PG*molar_mass_h2 / R / 303.15;
		if ((1 - S) < (rho_equili_h2_wet - rho_h2_wet))
		{
			return 0.0;
		}
		else {
			double const drhoh2wet_dpg = Hen*molar_mass_h2;
			double const drhoh2nonwet_dpg = molar_mass_h2 / R / T;
			return -(S*drhoh2wet_dpg + (1 - S)*drhoh2nonwet_dpg) / (rho_h2_wet - rho_h2_nonwet);
		}
	}
	/*
	* Calculate the derivative using the analytical way
	*/
	const double Calculate_dSwdX(double PG, double X, double S,double rho_h2_wet,
		double T = 303.15)
	{
		double const rho_equili_h2_wet = PG * Hen * molar_mass_h2;
		const double rho_h2_nonwet = PG*molar_mass_h2 / R / 303.15;
		if ((1 - S) < (rho_equili_h2_wet - rho_h2_wet))
		{
			return 0.0;
		}
		else{
		
		return 1/ (rho_h2_wet - rho_h2_nonwet);
		}
	}
	/*
	* Calculate the derivative using the analytical way
	*/
	const double Calculate_dX_mdX(double PG, double Sw, double rho_h2_wet)
	{
		double const rho_equili_h2_wet = PG * Hen * molar_mass_h2;
		if ((1 - Sw) < (rho_equili_h2_wet - rho_h2_wet))
			return 1.0;
		return 0.0;
	}
	/*
	* Calculate the derivative using the analytical way
	*/
	const double Calculate_dX_mdP(double PG, double Sw, double rho_h2_wet)
	{
		double const rho_equili_h2_wet = PG * Hen * molar_mass_h2;
		if ((1 - Sw) < (rho_equili_h2_wet - rho_h2_wet))
			return 0.0;
		return Hen * molar_mass_h2;
	}
private:
	EoS_IdealMix_Properties _eos_prop;
    double const Hen = 7.65e-6;  // mol/Pa./m3
                                 /**
                                 * Molar mass of water
                                 */
    double const molar_mass_h2o =
        MaterialLib::PhysicalConstant::MolarMass::Water;
	double const molar_mass_h2 = MaterialLib::PhysicalConstant::MolarMass::H2;
    /**
    * mass density of water
    */
    double const rho_mass_h20 = 1000;
    /**
    * molar density of water
    */
    double const rho_mol_h2o = rho_mass_h20 / molar_mass_h2o;
	double const R = MaterialLib::PhysicalConstant::IdealGasConstant;
};



}  // namespace 
}  // namespace ProcessLib

#include "EoS_IdealMix-impl.h"

#endif  // OGS_TWOPHASEWITHPRHO_EOS_IDEALMIX_H_
