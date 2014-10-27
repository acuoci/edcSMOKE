/*----------------------------------------------------------------------*\
|                  _       _____ __  __  ____  _  ________                |
|                 | |     / ____|  \/  |/ __ \| |/ /  ____|               |
|          ___  __| | ___| (___ | \  / | |  | | ' /| |__                  |
|         / _ \/ _` |/ __|\___ \| |\/| | |  | |  < |  __|                 |
|        |  __/ (_| | (__ ____) | |  | | |__| | . \| |____                |
|         \___|\__,_|\___|_____/|_|  |_|\____/|_|\_\______|               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, M.R. Malik, A. Parente                             |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|   Contacts: Mohammad Rafi Malik                                         |
|   Aero-Thermo-Mechanical Department                                     |
|   Université Libre de Bruxelles                                         |
|   Avenue F. D. Roosevelt 50, 1050 Bruxelles (Belgium)                   |
|                                                                         |
|   Contacts: Alessandro Parente                                          |
|   email: alessandro.parente@ulb.ac.be                                   |
|   Aero-Thermo-Mechanical Department                                     |
|   Université Libre de Bruxelles                                         |
|   Avenue F. D. Roosevelt 50, 1050 Bruxelles (Belgium)                   |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of edcSMOKE solver.                                 |
|                                                                         |
|	License                                                           |
|                                                                         |
|   Copyright(C) 2014 A. Cuoci, M.R. Malik, A. Parente                    |
|   edcSMOKE is free software: you can redistribute it and/or modify      |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   edcSMOKE is distributed in the hope that it will be useful,           |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with edcSMOKE. If not, see <http://www.gnu.org/licenses/>.      |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef ODE_PFR_H
#define ODE_PFR_H

namespace OpenSMOKE
{
	template<typename type>
	class ThermodynamicsMap_CHEMKIN;

	template<typename type>
	class KineticsMap_CHEMKIN;
}

class ODE_PFR
{
public:

	ODE_PFR(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& 		thermodynamicsMapXML, 
			OpenSMOKE::KineticsMap_CHEMKIN<double>& 		kineticsMapXML);

	void Set(const double P0_Pa)
	{
		P_Pa_ = P0_Pa;
	}

	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };
	unsigned int NumberOfEquations() const { return number_of_equations_; }

	int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);
	int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

private:

	unsigned int number_of_equations_;
	unsigned int number_of_gas_species_;

	double P_Pa_;
	double MWStar_;
	double cTotStar_;
	double rhoStar_;
	double cpStar_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& 	thermodynamicsMapXML_;
	OpenSMOKE::KineticsMap_CHEMKIN<double>& 	kineticsMapXML_;

	OpenSMOKE::OpenSMOKEVectorDouble omegaStar_;
	OpenSMOKE::OpenSMOKEVectorDouble xStar_;
	OpenSMOKE::OpenSMOKEVectorDouble cStar_;
	OpenSMOKE::OpenSMOKEVectorDouble RStar_;

	bool checkMassFractions_;
};

#include "ODE_PFR.C"

#endif // ODE_PFR_H

