#ifndef PROCESS_LIB_TWOPHASECOMPONENTIALMATERIALMODEL_H_
#define PROCESS_LIB_TWOPHASECOMPONENTIALMATERIALMODEL_H_

#include <math.h>
#include <cmath>

namespace ProcessLib
{

namespace TwoPhaseComponential {
	const double M_H = 0.002;
	const double M_L = 0.02;
	const double M_AIR = 0.029;
	const double M_C = 0.016;
	const double M_CO2 = 0.044;
	const double R = 8.314;
	const double Hen_L_h = 7.26e+9;
	const double Hen_L_c = 4.13e+9;
	const double Hen_L_co2 = 0.163e+9;
	const double Hen_L_air = 9.077e+9;
	/*van Genuchten parameters*/
	const double S_gr = 0.0;
	const double S_lr = 0.3;
	const double n_van_Genuchten = 2.0;
	const double Pb_van_Genuchten = 1e+5;
	/*Brooks corey parameters*/
	const double P_entry_Brook = 1e+5;
	const double Lambda_Brook = 2.0;
static double CalcIceVolFrac(double T_in_dC, double freezing_sigmoid_coeff, double porosity)
{
   double phi_i = 0.0;
   
   phi_i = porosity* (1.0 - 1.0 / (1.0 + std::exp(-1.0 * freezing_sigmoid_coeff * T_in_dC)));

   return phi_i;
}

static double get_P_sat(double T)
{
	// Here unit of T is Celsius;
	double P_sat(0.0);
	double T_0 = 373.15;
	double P_0 = 101325.0;
	double h_wg = 2258000.0;
	P_sat = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R);

	return P_sat;
}
static double get_X_G_air_gp(double PG, double X1, double X2, double X3, double P_sat)
{
	double K_G_w = PG / P_sat;
	double K_G_air = PG / Hen_L_air;
	double L = 1 - (X1*PG / Hen_L_h + X2*PG / Hen_L_c + X3*PG / Hen_L_co2);
	double G = 1 - X1 - X2 - X3;
	double X_G_air = (K_G_w*G - L) / (K_G_w - K_G_air);
	return X_G_air;
}
static double get_X_G_h2o_gp(double PG, double X1, double X2, double X3, double P_sat)
{
	double K_G_w = PG / P_sat;
	double K_G_air = PG / Hen_L_air;
	double L = 1 - (X1*PG / Hen_L_h + X2*PG / Hen_L_c + X3*PG / Hen_L_co2);
	double G = 1 - X1 - X2 - X3;
	double X_G_h2o = (L - K_G_air*G) / (K_G_w - K_G_air);
	return X_G_h2o;
}

/**
* calculate the effective saturation of gas phase by Gas saturation
*/
static double getEffectSat_gbySg(double Sg)
{
	double EffectSat_g = 0.0;
	// this->res_saturation_g->eval(Sg, res_S_g);

	EffectSat_g = (Sg - S_gr) / (1 - S_gr - S_lr);

	return EffectSat_g;
}
static double getEffectSat_lbySg(double Sg)
{
	double EffectSat_l = 0.0;
	// this->res_saturation_g->eval(Sg, res_S_g);

	if (Sg < S_gr)
		EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
	else if (Sg>(1 - S_lr))
		EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
	else
		EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);

	return EffectSat_l;
}
static double getKr_L_bySg(double Sg/*gas phase saturation*/, size_t capp_sat_model)//double m = 1.0 - 1.0 / 1.49, double res_S_l = 0.40, double res_S_g = 0.0 
{
	double Kr_L = 0.0;
	double EffectSat_l = 0.0;
	double m = 0.0;
	double Kr_L2 = 0.0;
	EffectSat_l = getEffectSat_lbySg(Sg);//for MoMaS benchmark 2

	switch (capp_sat_model)
	{
	case 4:
		m = 1.0 - 1.0 / n_van_Genuchten;
		Kr_L = sqrt(EffectSat_l)*pow(1 - pow(1 - pow(EffectSat_l, 1 / m), m), 2);
		if (Sg > 1)
			Kr_L = 0;
		else if (Sg < 0)
			Kr_L = 1;
		break;
	case 6:
		// krw = pow(se,3.0+2.0/m)
		// Se  = (sl - slr) / (slm - slr)
		//Kr_L = pow(EffectSat_l, 3.0 + 2.0 / Lambda_Brook);
		Kr_L = pow(EffectSat_l, 3);
		if (Sg > 1)
			Kr_L = 0;
		else if (Sg < 0)
			Kr_L = 1;
		break;
	case 8:
		Kr_L = pow(EffectSat_l, 3);
		if (Sg > 1 - S_lr)
			Kr_L = 0;
		else if (Sg < 0)
			Kr_L = 1;
		break;
	default:
		ERR("Error in getKr_L_bySg: No valid relative permeability vs saturation model! ");
		break;
	}
	//if (Kr_L < 1e-9)
	//Kr_L = 1e-9;

	return Kr_L;
}

/**
* get the relative permeability of gas phase based on Sg
*/
static double getKr_g_bySg(double Sg/*gas phase saturation*/, size_t capp_sat_model)
{
	double Kr_G = 0.0;

	double EffectSat_g = 0.0;
	double EffectSat_l = 0.0;
	EffectSat_g = getEffectSat_gbySg(Sg);
	EffectSat_l = getEffectSat_lbySg(Sg);
	double m = 0.0;
	switch (capp_sat_model)
	{
	case 4:// index 4 stands for the van Genuchten Model
		m = 1.0 - 1.0 / n_van_Genuchten;
		Kr_G = sqrt(EffectSat_g)*pow(1 - pow(1 - EffectSat_g, 1 / m), 2 * m);
		if (Sg > 1)
			Kr_G = 1;
		else if (Sg < 0)
			Kr_G = 0;
		break;
	case 6: //index 6 stands for the Brooks corey  model
			// krg = pow(1.0-se,2)*(1.0-pow(se,1.0+2.0/m))
			// Se  = (1-Sg - S_lr) / (1-S_lr - S_gr)
			//Kr_G = pow(1.0 - EffectSat_l, 2)*(1.0 - pow(EffectSat_l, 1.0 + 2.0 / Lambda_Brook));
		Kr_G = pow(1 - EffectSat_l, 3);
		if (Sg > 1)
			Kr_G = 1.0;
		else if (Sg < 0)
			Kr_G = 0.0;
		break;
	case 8:
		Kr_G = pow(1 - EffectSat_l, 3);
		if (Sg > 1 - S_lr)
			Kr_G = 1.0;
		else if (Sg < 0)
			Kr_G = 0.0;
		break;
	default:
		ERR("Error in getKr_g_bySg: No valid relative permeability vs saturation model! ");
		break;
	}
	return Kr_G;
}

static double getSat_byPC(double P_c, size_t capp_sat_model)//for MoMaS benchmark2 n = 1.49
{

	double Sg = 0.0;
	double m = 0.0;
	switch (capp_sat_model)
	{
	case 4:// index 4 stands for the van Genuchten Model

		m = 1.0 - 1.0 / n_van_Genuchten;
		if (P_c > 0) {
			Sg = 1 - (((1 - S_gr - S_lr) / pow((pow((P_c / Pb_van_Genuchten), n_van_Genuchten) + 1), m)) + S_lr);
		}
		else
			Sg = 0.0;//here need to pay attention to!!!!
		break;
	case 6: // index 6 stands for the Brooks Corey Model
		if (P_c<P_entry_Brook)
		{
			Sg = (P_c - P_entry_Brook)*Lambda_Brook / P_entry_Brook;
			//Sg = 1 - (S_lr + (1 - S_lr - S_gr)*pow(P_entry_Brook / P_c, Lambda_Brook));
		}
		else if (P_c > 4 * P_entry_Brook)
		{
			Sg = 1 - 1 / pow(4, Lambda_Brook) + Lambda_Brook*(P_c - 4 * P_entry_Brook) / P_entry_Brook / pow(4, 1 + Lambda_Brook);
		}

		else
		{
			Sg = 1 - (S_lr + (1 - S_lr - S_gr)*pow(P_entry_Brook / P_c, Lambda_Brook));// 1 - pow(P_entry_Brook / P_c, Lambda_Brook);
		}
		break;
	default:
		ERR("Error in getSat_byPC: No valid capillary pressure vs saturation model! ");
		break;
	}
	//if (Sg < 1e-5)
	//Sg = 1e-5;
	return Sg;

}

static double get_diriv_PC(double P_c, size_t capp_sat_model)//for MoMaS benchmark2 n = 1.49
{
	double dPC(0.0);
	double m = 0.0;
	switch (capp_sat_model)
	{
	case 4:
		m = 1.0 - 1.0 / n_van_Genuchten;
		dPC = m*n_van_Genuchten*(1 - S_gr - S_lr)*(1 / Pb_van_Genuchten)*(pow((P_c / Pb_van_Genuchten), (n_van_Genuchten - 1)))*pow(((pow((P_c / Pb_van_Genuchten), n_van_Genuchten)) + 1.), (-m - 1));
		if (P_c <= 0)
		{
			dPC = 0.0;//just for test
		}
		break;
	case 6:
		//dPC = (1 - S_lr - S_gr)*Lambda_Brook*pow(P_entry_Brook, Lambda_Brook)*pow(P_c, -Lambda_Brook - 1);
		if (P_c <= P_entry_Brook)
		{
			dPC = Lambda_Brook / P_entry_Brook;
			//Sg = 1 - (S_lr + (1 - S_lr - S_gr)*pow(P_entry_Brook / P_c, -Lambda_Brook));
		}
		else if (P_c > 4 * P_entry_Brook)
		{
			dPC = Lambda_Brook / P_entry_Brook / pow(4, 1. + Lambda_Brook);
		}
		else
		{
			dPC = (1 - S_lr - S_gr)* Lambda_Brook*pow(P_entry_Brook, Lambda_Brook)*pow(P_c, -Lambda_Brook - 1.);

		}
		break;
	default:
		ERR("Error in get_diriv_PC: No valid inverse capillary pressure vs saturation model! ");
		break;

	}
	//if (dPC < 1e-5)
	//dPC = 1e-5;
	return dPC;

}

}  // twophasecomponential

}  // ProcessLib


#endif  // PROCESS_LIB_FREEZINGMATERIALMODEL_H_
