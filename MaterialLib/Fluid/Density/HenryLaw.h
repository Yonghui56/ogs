/*!
   \file  IdealGasLow.h
   \brief Declaration of class IdealGasLow for fluid density by the ideal gas
          law depending on one variable linearly.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
 */
#ifndef HENRY_LAW_H_
#define HENRY_LAW_H_

#include <cassert>
#include <string>

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialLib
{
namespace Fluid
{
/// Fluid density by ideal gas law
class HenryLaw final : public FluidProperty
{
public:
    ///   \param molar_mass Molar mass of the gas phase.
    explicit HenryLaw(const double molar_mass, const std::string component)
        : FluidProperty(), _molar_mass(molar_mass), _component(component)
    {
    }

    /// Get density model name.
    std::string getName() const override { return "Henry's Law"; }
    /// Get density value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.
	/// rho_L^{air} = H(T) * P_G^{air} *M_{air}
	/// where H(T) is in mol.Pa^{-1}.m^{-3}
	/// molar fraction x_L^{air} = H(T) * P_G^{air} 
	/// where H(T) is in Pa^{-1}
    double getValue(const ArrayType& var_vals) const override
    {

		return _molar_mass * getHenryconst(var_vals[static_cast<int>(PropertyVariableType::T)]) *
			var_vals[static_cast<int>(PropertyVariableType::pg)] *
			var_vals[static_cast<int>(PropertyVariableType::molarg)]*
			1000/0.018;
    }

    /// Get the partial differential of the density with respect to temperature
    /// or gas pressure.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariableType.
    /// \param var       Variable type.
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::pg)];
		const double molar_gas = var_vals[static_cast<int>(PropertyVariableType::molarg)];
        switch (var)
        {
            case PropertyVariableType::T:
                return dHenryLaw_dT(T, p, molar_gas);
            case PropertyVariableType::pg:
                return dHenryLaw_dp(T, p, molar_gas);
			case PropertyVariableType::molarg:
				return dHenryLaw_dmol(T, p, molar_gas);
            default:
                return 0.;
        }
    }

private:
    /// Molar mass of gas phase.
    const double _molar_mass;
	const std::string _component;

    /// Get the partial differential of density with the respect to temperature
    /// \param T  Temperature in K.
    /// \param pg Gas phase pressure in Pa.
    double dHenryLaw_dT(const double T, const double pg, const double molarg) const
    {
		return _molar_mass*molarg*pg*getdHenryconstdT(T);
    }

    /// Get the partial differential of density with the respect to pressure
    /// \param T  Temperature in K.
    /// \param pg Gas phase pressure in Pa.
    double dHenryLaw_dp(const double T, const double /* pg */, const double molarg) const
    {
        return _molar_mass * getHenryconst(T)*molarg;
    }
	double dHenryLaw_dmol(const double T, const double pg , const double molarg) const
	{
		return _molar_mass*getHenryconst(T)*pg;
	}
	double getHenryconst(const double T) const
	{
		double henry_const;
		if (_component == "air")
		{
			henry_const= (0.8942 + 1.47*exp(-0.04394*(T - 273.15)))*1.E-10;
			return henry_const;
		}
		else
		{
			OGS_FATAL("This component has not been implemented yet");
		}
	}
	double getdHenryconstdT(const double T) const
	{
		double dhenryconstdT;
		double dhenryconstdT_ana;
		const double eps = 1e-6;
		if (_component == "air")
		{
			dhenryconstdT = (getHenryconst(T + eps*T) - getHenryconst(T- eps*T)) / eps/T/2;
			dhenryconstdT_ana = 1.47*1.E-10*(-0.04394)*exp(-0.04394*(T - 273.15));
			return dhenryconstdT;
		}
		else
		{
			OGS_FATAL("The derivative for this component has not been implemented yet");
		}
	}
};
}  // end namespace
}  // end namespace
#endif
