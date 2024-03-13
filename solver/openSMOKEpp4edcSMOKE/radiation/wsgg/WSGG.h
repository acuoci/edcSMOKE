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

#ifndef OpenSMOKE_WSGG_H
#define	OpenSMOKE_WSGG_H

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp> 
#include <boost/math/special_functions/pow.hpp> 
#include "utilities/grids/multi-purpose/Grid1DMultiPurpose.h"

namespace OpenSMOKE
{

	enum ExtinctionCoefficientModel {	EXTINCTION_COEFFICIENT_GRAY, EXTINCTION_COEFFICIENT_WSGG_SMITH, EXTINCTION_COEFFICIENT_WSGG_TRUELOVE, EXTINCTION_COEFFICIENT_WSGG_YIN,
										EXTINCTION_COEFFICIENT_WSGG_SMITH_SOOT, EXTINCTION_COEFFICIENT_WSGG_SMITH_SOOT_SINGLEBAND, EXTINCTION_COEFFICIENT_WSGG_CASSOL};

	// Truelove model (Modest, Radiative Heat Transfer, 2nd Edition, Table 19.3)
	void WSGG_Truelove(const double T, const double PCO2_Pa, const double PH2O_Pa, double* kappa, double* a);

	// Smith model (Smith et al., 1982)
	double SmithEmissivity(const double b1, const double b2, const double b3, const double b4, const double T);
	void WSGG_Smith(const double T, const double PCO2_Pa, const double PH2O_Pa, double* kappa, double* a);
	
	// Yin model (Yin et al. Energy Fuels 2010, 24, p. 6275–6282, DOI:10.1021/ef101211p)
	double YinEmissivity(const double b1, const double b2, const double b3, const double b4, const double T);
	void WSGG_Yin(const double T, const double PCO2_Pa, const double PH2O_Pa, double* kappa, double* a);

	// Smith model accounting for soot (Smith et al., 1987)
	void WSGG_Smith_Soot(const double T, const double PCO2_Pa, const double PH2O_Pa, const double soot_fv, double* kappa, double* a);

	// Smith model accounting for soot, with a single weighted coefficient (Smith et al., 1987)
	void WSGG_Smith_Soot_Singleband(const double T, const double PCO2_Pa, const double PH2O_Pa, const double soot_fv, double* kappa, double* a);
	
	// Cassol model including soot (Cassol et al., 2014)
    void WSGG_Cassol(const double T, const double PCO2_Pa, const double PH2O_Pa, const double soot_fv, double* kappa, double* a);
    double CassolEmissivity(const double b1, const double b2, const double b3, const double b4, const double b5, const double T);
}

#include "WSGG.hpp"

#endif	/* OpenSMOKE__WSGG_H */
