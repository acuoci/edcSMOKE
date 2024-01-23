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
|   License                                                               |
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

#include "math/OpenSMOKEUtilities.h"

namespace OpenSMOKE
{
	void TransportPropertiesMap::Thermodiffusion(const double* moleFractions, const double* cpspecies, const unsigned int Mmain, const std::vector<unsigned int>& kreorder)
	{
		solve_thermodiffusion_onePlusM(moleFractions, cpspecies, Mmain, kreorder);
	}

	double TransportPropertiesMap::ThermalConductivity(const double* moleFractions)
	{
		lambda();
		return lambdaMix(moleFractions);
	}

	double TransportPropertiesMap::ThermalConductivity(const double* moleFractions, const unsigned int Mmain, const std::vector<unsigned int>& kreorder)
	{
		return lambdaMult(moleFractions, Mmain, kreorder);
	}

	double TransportPropertiesMap::DynamicViscosity(const double* moleFractions)
	{
		eta();
		return etaMix(moleFractions);
	}

	double TransportPropertiesMap::DynamicViscosity(const double* moleFractions, const unsigned int Mmain, const std::vector<unsigned int>& kreorder)
	{
		eta(Mmain, kreorder);
		return etaMix(moleFractions, Mmain, kreorder);
	}

	void TransportPropertiesMap::MassDiffusionCoefficients(double* gammamix, const double* moleFractions, const bool bundling)
	{
		if (bundling == false)
		{
			gamma();
			gammaMix(gammamix, moleFractions);
		}
		else
		{
			bundling_gamma();
			bundling_gammaMix(gammamix, moleFractions);
		}
	}

	void TransportPropertiesMap::MassDiffusionCoefficients(double* gammamult, const double* moleFractions, const unsigned int Mmain, const std::vector<unsigned int>& kreorder)
	{
		gamma();
		gammaMult(gammamult, moleFractions, Mmain, kreorder);
	}

	void TransportPropertiesMap::BinaryDiffusionCoefficients(const bool bundling)
	{
		if (bundling == false)
			gamma();
		else
			bundling_gamma();
	}

	void TransportPropertiesMap::ThermalDiffusionRatios(double* tetamix, const double* moleFractions)
	{
		teta();
		tetaMix(tetamix, moleFractions);
	}

	void TransportPropertiesMap::ThermalDiffusionRatios(double* tetamix, const double* moleFractions, const unsigned int Mmain, const std::vector<unsigned int>& kreorder)
	{
		tetaMult(tetamix, moleFractions, Mmain, kreorder);
	}

}
