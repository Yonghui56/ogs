/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   VanGenuchtenCapillaryPressureSaturation.h
 *
 *  Created on October 28, 2016, 6:05 PM
 */

#ifndef OGS_VAN_GENUCHTEN_CAPILLARY_PRESSURE_SATURATION_H
#define OGS_VAN_GENUCHTEN_CAPILLARY_PRESSURE_SATURATION_H

#include <limits>

#include "CapillaryPressureSaturation.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
 *   \brief van Genuchten water retention model
 *
 *   \f[p_c=p_b (S_e^{-1/m}-1)^{1-m}\f]
 *   with
 *   \f[S_e=\frac{S-S_r}{S_{\mbox{max}}-S_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &p_b&            \mbox{ entry pressure,}\\
 *       &S_r&            \mbox{ residual saturation,}\\
 *       &S_{\mbox{max}}& \mbox{ maximum saturation,}\\
 *       &m(<=1) &        \mbox{ exponent.}\\
 *    \f}
 *
 *    Note:
 *     \f[m=1/(1-n)\f].
 *
 *    If \f$\alpha\f$ instead of \f$p_b\f$ is available, \f$p_b\f$ can be
 * calculated
 * as
 *    \f[p_b=\rho g/\alpha\f]
 */
class VanGenuchtenCapillaryPressureSaturation final
    : public CapillaryPressureSaturation
{
public:
    /**
     * @param pb     Entry pressure, \f$ p_b \f$
     * @param Sr     Residual saturation, \f$ S_r \f$
     * @param Smax   Maximum saturation, \f$ S_{\mbox{max}} \f$
     * @param m      Exponent, \f$ m \f$
     * @param Pc_max Maximum capillary pressure, \f$ P_c^{\mbox{max}}\f$
     */
    VanGenuchtenCapillaryPressureSaturation(const double pb, const double Sr,
                                            const double Smax, const double m,
                                            const double Pc_max)
        : CapillaryPressureSaturation(Sr, Smax, Pc_max), _pb(pb), _m(m)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "van Genuchten water retention model.";
    }

    /// Get capillary pressure.
    double getCapillaryPressure(const double saturation) const override;
	double getRegularizedCapillaryPressure(const double saturation) const override;
    /// Get saturation.
    double getSaturation(const double capillary_pressure) const override;

    /// Get the derivative of the capillary pressure with respect to saturation
    double getdPcdS(const double saturation) const override;
	double getRegularizedPcdS(const double saturation) const override;
private:
    const double _pb;  ///< Entry pressure.
    const double _m;   ///< Exponent m, m in [0,1]. n=1/(1-m).
	const double xi = 1e-5;
private:
	/**
	* regularized van Genuchten capillary pressure-saturation Model
	*/
	double getPc_bar_vG_Sg(double Sg) const
	{
		double const S_gr = 0.0;
		double const S_bar = getS_bar(Sg);
		double const Pc_bar_vG = getPc_vG_Sg(S_bar) - getPc_vG_Sg(S_gr + (1 - S_gr - _saturation_r)*xi / 2);
		return Pc_bar_vG;
	}
	/**
	* regularized van Genuchten capillary pressure-saturation Model
	*/
	double getS_bar(double Sg) const
	{

		double const S_gr = 0.0;
		return  S_gr + (1 - xi)*(Sg - S_gr) + 0.5*xi*(1 - S_gr - _saturation_r);
	}
	/**
	*  van Genuchten capillary pressure-saturation Model
	*/
	double getPc_vG_Sg(double Sg) const
	{
		double const S_gr = 0.0;
		//effective saturation
		double const S_le = (1 - Sg - _saturation_r) / (1 - S_gr - _saturation_r);
		//Pc_vG = P_r*(S_le. ^ (-1 / m) - 1). ^ (1 / n);
		return _pb * std::pow(std::pow(S_le, (-1.0 / _m)) - 1.0, 1.0 - _m);
	}
	/**
	* derivative dPCdS based on regularized van Genuchten capillary pressure-saturation Model
	*/
	double get_dPCdS_vG_bar(double Sg) const
	{
		double dPCdS(0.0);
		double S_bar = 0.0;
		S_bar = getS_bar(Sg);
		dPCdS = get_dPCdS_vG(S_bar)*(1 - xi);
		return dPCdS;
	}
	/**
	* derivative dPCdS based on standard van Genuchten capillary pressure-saturation Model
	*/
	virtual double get_dPCdS_vG(double Sg) const
	{
		double dPcdSg = 0.0;
		double const S_gr = 0.0;
		double const _nn = 1 / (1 - _m);
		double S_le = (1 - Sg - _saturation_r) / (1 - S_gr - _saturation_r);
		return _pb*(1 / (_m*_nn))*(1 / (1 - _saturation_r - S_gr))*pow(pow(S_le, (-1 / _m)) - 1, (1 / _nn) - 1)*pow(S_le, (-1 / _m)) / S_le;

	}
};

}  // end namespace
}  // end namespace
#endif /* OGS_VAN_GENUCHTEN_CAPILLARY_PRESSURE_SATURATION_H */
