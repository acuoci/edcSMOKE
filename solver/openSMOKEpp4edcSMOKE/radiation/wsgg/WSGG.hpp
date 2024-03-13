/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

namespace OpenSMOKE
{
	
	// Truelove model
	void WSGG_Truelove(const double T, const double PCO2_Pa, const double PH2O_Pa, double* kappa, double* a)
	{
		const double PH2O = PH2O_Pa / 101325.;
		const double PCO2 = PCO2_Pa / 101325.;
		const double P = PH2O + PCO2;

		a[0] = 0.5880 - 0.2401 / 1000.*T - 0.1650 + 0.2834 / 1000.*T;
		a[1] = 0.4120 - 0.1665 / 1000.*T - 0.1270 + 0.2178 / 1000.*T;
		a[2] = 0.2375 - 0.0941 / 1000.*T - 0.0105 + 0.0265 / 1000.*T;
		a[3] = 0.0585 - 0.0243 / 1000.*T + 0.0065 - 0.0027 / 1000.*T;

		kappa[0] = 0.;
		kappa[1] = 0.89*P;
		kappa[2] = 15.5*P;
		kappa[3] = 239.*P;

		// Test
		//kappa[0] = 0; a[0] = 0.;
		//kappa[1] = 1; a[1] = 0.4;
		//kappa[2] = 1; a[2] = 0.3;
		//kappa[3] = 1; a[3] = 0.3;

	}

	double SmithEmissivity(const double b1, const double b2, const double b3, const double b4, const double T)
	{
		return ( (b1 / 10.) + (b2/1.e4)*T + (b3/1.e7)*T*T + (b4/1.e11)*T*T*T );
	}

	void WSGG_Smith(const double T, const double PCO2_Pa, const double PH2O_Pa, double* kappa, double* a)
	{
		const double Pw = PH2O_Pa / 101325.;
		const double Pc = PCO2_Pa / 101325.;
		const double P = Pw + Pc;

		if (Pw <= 0.5*Pc)	//use Pc → 0 atm table 
		{
			kappa[1] = 0.3966;	a[1] = SmithEmissivity(0.4334, 2.620, -1.560, 2.565, T);
			kappa[2] = 15.64;	a[2] = SmithEmissivity(-0.4814, 2.822, -1.794, 3.274, T);
			kappa[3] = 394.3;	a[3] = SmithEmissivity(0.5492, 0.1087, -0.3500, 0.9123, T);
		}
		else if (Pw <= 1.5*Pc) //use(Pw / Pc) = 1 table 
		{
			kappa[1] = 0.4303;	a[1] = SmithEmissivity(5.150, -2.303, 0.9779, -1.494, T);
			kappa[2] = 7.055;	a[2] = SmithEmissivity(0.7749, 3.399, -2.297, 3.770, T);
			kappa[3] = 178.1;	a[3] = SmithEmissivity(1.907, -1.824, 0.5608, -0.5122, T);
		}
		else if (Pw <= 2.5*Pc) // use(Pw / Pc) = 2 table 
		{
			kappa[1] = 0.4201;	a[1] = SmithEmissivity(6.508, -5.551, 3.029, -5.353, T);
			kappa[2] = 6.516;	a[2] = SmithEmissivity(-0.2504, 6.112, -3.882, 6.528, T);
			kappa[3] = 131.9;	a[3] = SmithEmissivity(2.718, -3.118, 1.221, -1.612, T);
		}
		else if (Pw <= 0.5) // use Pw → 0 atm table 
		{
			kappa[1] = 0.4098;	a[1] = SmithEmissivity(5.977, -5.119, 3.042, -5.564, T);
			kappa[2] = 6.325;	a[2] = SmithEmissivity(0.5677, 3.333, -1.967, 2.718, T);
			kappa[3] = 120.5;	a[3] = SmithEmissivity(1.800, -2.334, 1.008, -1.454, T);
		}
		else // use Pw = 1 atm table 
		{
			kappa[1] = 0.4496;	a[1] = SmithEmissivity(6.324, -8.358, 6.135, -13.03, T);
			kappa[2] = 7.113;	a[2] = SmithEmissivity(-0.2016, 7.145, -5.212, 9.868, T);
			kappa[3] = 119.7;	a[3] = SmithEmissivity(3.500, -5.040, 2.425, -3.888, T);
		}
	
		// Spectral window
		a[0] = 1. - a[1] - a[2] - a[3];
		kappa[0] = 0.;

		// Correction
		for (unsigned int j = 0; j < 4; j++)
			kappa[j] *= P;
	}

	void WSGG_Smith_Soot_Singleband(const double T, const double PCO2_Pa, const double PH2O_Pa, const double soot_fv, double* kappa, double* a)
	{
		const double Pw = PH2O_Pa / 101325.;
		const double Pc = PCO2_Pa / 101325.;
		const double P = Pw + Pc;

		double *kappa_gas = new double[4];
		double *a_gas = new double[4];

		double *kappa_soot = new double[2];
		double *a_soot = new double[2];

		double *kappa_mix = new double[8];
		double *a_mix = new double[8];

		kappa[0] = 0;
		a[0] = 1;

		// Absorbing gas coefficients
		if (Pw <= 0.5*Pc)	//use Pc → 0 atm table 
		{
			kappa_gas[1] = 0.3966;	a_gas[1] = SmithEmissivity(0.4334, 2.620, -1.560, 2.565, T);
			kappa_gas[2] = 15.64;	a_gas[2] = SmithEmissivity(-0.4814, 2.822, -1.794, 3.274, T);
			kappa_gas[3] = 394.3;	a_gas[3] = SmithEmissivity(0.5492, 0.1087, -0.3500, 0.9123, T);
		}
		else if (Pw <= 1.5*Pc) //use(Pw / Pc) = 1 table 
		{
			kappa_gas[1] = 0.4303;	a_gas[1] = SmithEmissivity(5.150, -2.303, 0.9779, -1.494, T);
			kappa_gas[2] = 7.055;	a_gas[2] = SmithEmissivity(0.7749, 3.399, -2.297, 3.770, T);
			kappa_gas[3] = 178.1;	a_gas[3] = SmithEmissivity(1.907, -1.824, 0.5608, -0.5122, T);
		}
		else if (Pw <= 2.5*Pc) // use(Pw / Pc) = 2 table 
		{
			kappa_gas[1] = 0.4201;	a_gas[1] = SmithEmissivity(6.508, -5.551, 3.029, -5.353, T);
			kappa_gas[2] = 6.516;	a_gas[2] = SmithEmissivity(-0.2504, 6.112, -3.882, 6.528, T);
			kappa_gas[3] = 131.9;	a_gas[3] = SmithEmissivity(2.718, -3.118, 1.221, -1.612, T);
		}
		else if (Pw <= 0.5) // use Pw → 0 atm table 
		{
			kappa_gas[1] = 0.4098;	a_gas[1] = SmithEmissivity(5.977, -5.119, 3.042, -5.564, T);
			kappa_gas[2] = 6.325;	a_gas[2] = SmithEmissivity(0.5677, 3.333, -1.967, 2.718, T);
			kappa_gas[3] = 120.5;	a_gas[3] = SmithEmissivity(1.800, -2.334, 1.008, -1.454, T);
		}
		else // use Pw = 1 atm table 
		{
			kappa_gas[1] = 0.4496;	a_gas[1] = SmithEmissivity(6.324, -8.358, 6.135, -13.03, T);
			kappa_gas[2] = 7.113;	a_gas[2] = SmithEmissivity(-0.2016, 7.145, -5.212, 9.868, T);
			kappa_gas[3] = 119.7;	a_gas[3] = SmithEmissivity(3.500, -5.040, 2.425, -3.888, T);
		}

		a_gas[0] = 1. - a_gas[1] - a_gas[2] - a_gas[3];
		kappa_gas[0] = 0.;

		//Soot coefficients
		kappa_soot[0] = 1.00802e6 * soot_fv;	a_soot[0] = SmithEmissivity(1.420 * 10, -0.77942e-3 * 1e4, -0.38408e-7 * 1e7, 0.24166e-10 * 1e11, T);
		kappa_soot[1] = 3.23520e6 * soot_fv;	a_soot[1] = SmithEmissivity(-0.420 * 10, 0.77942e-3 * 1e4, 0.38408e-7 * 1e7, -0.24166e-10 * 1e11, T);

		unsigned int n_band = 0;
		for(unsigned int i = 0; i < 2; i++)
			for(unsigned int j = 0; j < 4; j++)
			{
				kappa_mix[n_band] = kappa_soot[i] + kappa_gas[j] * P;
				a_mix[n_band] = a_soot[i] * a_gas[j];
				n_band++;
			}

		//Mixture absorption coefficient
		for(unsigned int p = 0; p < n_band; p++)
			kappa[0] += kappa_mix[p] * a_mix[p];
                
                if(kappa[0] < 0 && kappa[0] > -1.e-20)
                  kappa[0] = 0;
                
                delete [] kappa_gas;
                delete [] a_gas;
                
                delete [] kappa_soot;
                delete [] a_soot;
                
                delete [] kappa_mix;
                delete [] a_mix;
	}

	void WSGG_Smith_Soot(const double T, const double PCO2_Pa, const double PH2O_Pa, const double soot_fv, double* kappa, double* a)
	{
		const double Pw = PH2O_Pa / 101325.;
		const double Pc = PCO2_Pa / 101325.;
		const double P = Pw + Pc;

		double *kappa_gas = new double[4];
		double *a_gas = new double[4];

		double *kappa_soot = new double[2];
		double *a_soot = new double[2];

		// Absorbing gas coefficients
		if (Pw <= 0.5*Pc)	//use Pc → 0 atm table 
		{
			kappa_gas[1] = 0.3966;	a_gas[1] = SmithEmissivity(0.4334, 2.620, -1.560, 2.565, T);
			kappa_gas[2] = 15.64;	a_gas[2] = SmithEmissivity(-0.4814, 2.822, -1.794, 3.274, T);
			kappa_gas[3] = 394.3;	a_gas[3] = SmithEmissivity(0.5492, 0.1087, -0.3500, 0.9123, T);
		}
		else if (Pw <= 1.5*Pc) //use(Pw / Pc) = 1 table 
		{
			kappa_gas[1] = 0.4303;	a_gas[1] = SmithEmissivity(5.150, -2.303, 0.9779, -1.494, T);
			kappa_gas[2] = 7.055;	a_gas[2] = SmithEmissivity(0.7749, 3.399, -2.297, 3.770, T);
			kappa_gas[3] = 178.1;	a_gas[3] = SmithEmissivity(1.907, -1.824, 0.5608, -0.5122, T);
		}
		else if (Pw <= 2.5*Pc) // use(Pw / Pc) = 2 table 
		{
			kappa_gas[1] = 0.4201;	a_gas[1] = SmithEmissivity(6.508, -5.551, 3.029, -5.353, T);
			kappa_gas[2] = 6.516;	a_gas[2] = SmithEmissivity(-0.2504, 6.112, -3.882, 6.528, T);
			kappa_gas[3] = 131.9;	a_gas[3] = SmithEmissivity(2.718, -3.118, 1.221, -1.612, T);
		}
		else if (Pw <= 0.5) // use Pw → 0 atm table 
		{
			kappa_gas[1] = 0.4098;	a_gas[1] = SmithEmissivity(5.977, -5.119, 3.042, -5.564, T);
			kappa_gas[2] = 6.325;	a_gas[2] = SmithEmissivity(0.5677, 3.333, -1.967, 2.718, T);
			kappa_gas[3] = 120.5;	a_gas[3] = SmithEmissivity(1.800, -2.334, 1.008, -1.454, T);
		}
		else // use Pw = 1 atm table 
		{
			kappa_gas[1] = 0.4496;	a_gas[1] = SmithEmissivity(6.324, -8.358, 6.135, -13.03, T);
			kappa_gas[2] = 7.113;	a_gas[2] = SmithEmissivity(-0.2016, 7.145, -5.212, 9.868, T);
			kappa_gas[3] = 119.7;	a_gas[3] = SmithEmissivity(3.500, -5.040, 2.425, -3.888, T);
		}

		//Spectral window
		a_gas[0] = 1. - a_gas[1] - a_gas[2] - a_gas[3];
		kappa_gas[0] = 0.;

		//Soot coefficients
		kappa_soot[0] = 1.00802e6 * soot_fv;	a_soot[0] = SmithEmissivity(1.420 * 10, -0.77942e-3 * 1e4, -0.38408e-7 * 1e7, 0.24166e-10 * 1e11, T);
		kappa_soot[1] = 3.23520e6 * soot_fv;	a_soot[1] = SmithEmissivity(-0.420 * 10, 0.77942e-3 * 1e4, 0.38408e-7 * 1e7, -0.24166e-10 * 1e11, T);

		unsigned int n_band = 0;
		for(unsigned int i = 0; i < 2; i++)
			for(unsigned int j = 0; j < 4; j++)
			{
				kappa[n_band] = kappa_soot[i] + kappa_gas[j] * P;
				a[n_band] = a_soot[i] * a_gas[j];
				n_band++;
			}
                
                delete [] kappa_gas;
                delete [] a_gas;
                
                delete [] kappa_soot;
                delete [] a_soot;
                
	}

	double YinEmissivity(const double b1, const double b2, const double b3, const double b4, const double T)
	{
		const double t = T / 1200.;
		return b1 + t*(b2 + t*(b3 + t*b4));
	}

	// Yin model
	void WSGG_Yin(const double T, const double PCO2_Pa, const double PH2O_Pa, double* kappa, double* a)
	{
		const double Pw = PH2O_Pa / 101325.;
		const double Pc = PCO2_Pa / 101325.;
		const double P = Pw + Pc;

		if (Pw <= 0.01*Pc)		// use Pc → 0 atm table
		{
			kappa[1] = 0.163233;	a[1] = YinEmissivity(0.204623, -0.378060, 0.666639, -0.203453, T);
			kappa[2] = 13.096584;	a[2] = YinEmissivity(-0.020227, 0.256006, -0.195201, 0.040493, T);
			kappa[3] = 175.474735;	a[3] = YinEmissivity(0.044221, 0.003850, -0.020175, 0.004919, T);
			kappa[4] = 1310.847307;	a[4] = YinEmissivity(0.039311, -0.054832, 0.025370, -0.003891, T);
		}
		else if (Pw <= 0.5*Pc)	// use(Pw / Pc) = 0.005 table
		{
			kappa[1] = 0.352505;	a[1] = YinEmissivity(0.315106, 0.023475, -0.057930, 0.008408, T);
			kappa[2] = 8.210621;	a[2] = YinEmissivity(0.092474, 0.109146, -0.121000, 0.027145, T);
			kappa[3] = 137.410012;	a[3] = YinEmissivity(0.031702, 0.037396, -0.040731, 0.008742, T);
			kappa[4] = 1269.710976;	a[4] = YinEmissivity(0.046138, -0.061392, 0.027164, -0.003996, T);
		}
		else if (Pw <= 1.5*Pc)	// use(Pw / Pc) = 1 table
		{
			kappa[1] = 0.261021;	a[1] = YinEmissivity(0.500119, -0.447068, 0.286878, -0.059165, T);
			kappa[2] = 3.147817;	a[2] = YinEmissivity(0.071592, 0.508252, -0.384253, 0.073477, T);
			kappa[3] = 54.265868;	a[3] = YinEmissivity(0.155320, -0.104294, 0.014096, 0.001643, T);
			kappa[4] = 482.900353;	a[4] = YinEmissivity(0.072615, -0.100601, 0.046681, -0.007224, T);
		}
		else if (Pw <= 2.5*Pc)	// use(Pw / Pc) = 2 table
		{
			kappa[1] = 0.179160;	a[1] = YinEmissivity(0.542458, -0.658411, 0.466444, -0.100186, T);
			kappa[2] = 2.388971;	a[2] = YinEmissivity(0.101734, 0.518429, -0.386151, 0.073453, T);
			kappa[3] = 28.415805;	a[3] = YinEmissivity(0.146066, -0.008745, -0.058325, 0.015984, T);
			kappa[4] = 253.059089;	a[4] = YinEmissivity(0.129511, -0.187993, 0.090709, -0.014493, T);
		}
		else if (Pw <= 0.01)	// use Pw → 0 atm table
		{
			kappa[1] = 0.085523;	a[1] = YinEmissivity(0.966357, -0.790165, -0.050144, 0.115202, T);
			kappa[2] = 0.475777;	a[2] = YinEmissivity(0.662059, -2.262877, 2.309473, -0.572895, T);
			kappa[3] = 8.549733;	a[3] = YinEmissivity(0.060870, 0.436788, -0.395493, 0.085146, T);
			kappa[4] = 201.906503;	a[4] = YinEmissivity(0.103568, -0.153135, 0.074910, -0.012091, T);
		}
		else if (Pw <= 0.2)		// use Pw = 0.05 atm table
		{
			kappa[1] = 0.232724;	a[1] = YinEmissivity(0.340618, -0.105469, 0.068051, -0.017828, T);
			kappa[2] = 2.134299;	a[2] = YinEmissivity(0.175818, -0.063466, 0.086631, -0.026581, T);
			kappa[3] = 9.266065;	a[3] = YinEmissivity(0.044325, 0.288376, -0.258205, 0.054333, T);
			kappa[4] = 134.988332;	a[4] = YinEmissivity(0.126628, -0.186480, 0.090755, -0.014569, T);
		}
		else					// use Pw = 1 atm table
		{
			kappa[1] = 0.065411;	a[1] = YinEmissivity(-0.077336, 0.661776, -0.362515, 0.053534, T);
			kappa[2] = 0.696552;	a[2] = YinEmissivity(0.506777, -0.758948, 0.516146, -0.102909, T);
			kappa[3] = 4.862610;	a[3] = YinEmissivity(-0.079989, 0.851078, -0.604264, 0.113500, T);
			kappa[4] = 60.255980;	a[4] = YinEmissivity(0.373898, -0.540887, 0.258923, -0.040957, T);
		}

		// Spectral window
		a[0] = 1. - a[1] - a[2] - a[3] - a[4];
		kappa[0] = 0.;

		// Correction
		for (unsigned int j = 0; j < 5; j++)
			kappa[j] *= P;
	}
        
        void WSGG_Cassol(const double T, const double PCO2_Pa, const double PH2O_Pa, const double soot_fv, double* kappa, double* a)
        {
		const double Pw = PH2O_Pa / 101325.;
		const double Pc = PCO2_Pa / 101325.;
		const double P = Pw + Pc;
                
                double *kappa_soot = new double[4];
		double *a_soot = new double[4];
                double alfa = 6.3; //Cassol et al. (2014)
                
                double *kappa_CO2 = new double[5];
		double *a_CO2 = new double[5];
                
                double *kappa_H2O = new double[5];
		double *a_H2O = new double[5];
                
                //Coefficients for H2O
                {
                  kappa_H2O[1] = 0.171 * Pw;    a_H2O[1] = CassolEmissivity( 0.06617,  55.48, -48.41,  22.27, -40.17, T);
                  kappa_H2O[2] = 1.551 * Pw;    a_H2O[2] = CassolEmissivity( 0.11045,  0.576,  24.00, -17.01,  30.96, T);
                  kappa_H2O[3] = 5.562 * Pw;    a_H2O[3] = CassolEmissivity(-0.04915,  70.63, -70.12,  26.07, -34.94, T);
                  kappa_H2O[4] = 49.159 * Pw;   a_H2O[4] = CassolEmissivity( 0.23675, -18.91, -0.907,  4.082, -8.778, T);
                }
                
                //Coefficients for CO2
                {
                  kappa_CO2[1] = 0.138 * Pc;    a_CO2[1] = CassolEmissivity( 0.09990,  64.41, -86.94,  41.27, -67.74, T);
                  kappa_CO2[2] = 1.895 * Pc;    a_CO2[2] = CassolEmissivity( 0.00942,  10.36, -2.277, -2.134,  6.497, T);
                  kappa_CO2[3] = 13.301 * Pc;   a_CO2[3] = CassolEmissivity( 0.14511, -30.73,  37.65, -18.41,  30.16, T);
                  kappa_CO2[4] = 340.811 * Pc;  a_CO2[4] = CassolEmissivity(-0.02915,  25.23, -26.10,  9.965, -13.26, T);
                }
                
                //Coefficients for soot
                {
                  kappa_soot[0] = 2875.86 * soot_fv * alfa;     a_soot[0] = CassolEmissivity( 0.00129, -0.545,  0.123*10, -0.847,  1.6807, T);
                  kappa_soot[1] = 39234.9 * soot_fv * alfa;     a_soot[1] = CassolEmissivity( 1.26110, -319.2,  27.72*10, -100.5,  132.80, T);
                  kappa_soot[2] = 160748.0 * soot_fv * alfa;    a_soot[2] = CassolEmissivity(-0.25757,  362.1, -40.12*10,  154.9, -207.80, T);
                  kappa_soot[3] = 495898.0 * soot_fv * alfa;    a_soot[3] = CassolEmissivity( 0.07980, -72.08,  15.87*10, -70.89,  97.690, T);
                }
                
                a_H2O[0] = 1. - a_H2O[1] - a_H2O[2] - a_H2O[3] - a_H2O[4];
                kappa_H2O[0] = 0.;
                
                a_CO2[0] = 1. - a_CO2[1] - a_CO2[2] - a_CO2[3] - a_CO2[4];
                kappa_CO2[0] = 0.;
                
                unsigned int n_band = 0;
                for(unsigned int i = 0; i < 4; i++)
                  for(unsigned int j = 0; j < 5; j++)
                    for(unsigned int k = 0; k < 5; k++)
                      {
                        kappa[n_band] = kappa_soot[i] + kappa_H2O[j] + kappa_CO2[k];
                        a[n_band] = a_soot[i] * a_H2O[j] * a_CO2[k];
                        n_band++;
                      }
                
                
                
                delete [] kappa_soot;
                delete [] a_soot;
                
                delete [] kappa_CO2;
                delete [] a_CO2;
                
                delete [] kappa_H2O;
                delete [] a_H2O;
        }
        
        double CassolEmissivity(const double b1, const double b2, const double b3, const double b4, const double b5, const double T)
        {
          return (b1 + b2*1e-5*T + b3*1e-8*T*T + b4*1e-11*T*T*T + b5*1e-15*T*T*T*T);
        }
}
