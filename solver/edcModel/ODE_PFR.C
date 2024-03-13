/*-----------------------------------------------------------------------*\
|                  _       _____ __  __  ____  _  ________                |
|                 | |     / ____|  \/  |/ __ \| |/ /  ____|               |
|          ___  __| | ___| (___ | \  / | |  | | ' /| |__                  |
|         / _ \/ _` |/ __|\___ \| |\/| | |  | |  < |  __|                 |
|        |  __/ (_| | (__ ____) | |  | | |__| | . \| |____                |
|         \___|\__,_|\___|_____/|_|  |_|\____/|_|\_\______|               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, M.R. Malik, Z. Li, A. Parente                      |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|   Contacts: Mohammad Rafi Malik, Zhiyi Li, Alessandro Parente           |
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
|   Copyright(C) 2017-2014 A. Cuoci, A. Parente                           |
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

// OpenSMOKE
#include "OpenSMOKE_Definitions.h"
#include <string>
#include <iostream>
#include <numeric>
#include <Eigen/Dense>

// Base classes
#include "kernel/thermo/ThermoPolicy_CHEMKIN.h"
#include "kernel/kinetics/ReactionPolicy_CHEMKIN.h"
#include "math/PhysicalConstants.h"
#include "math/OpenSMOKEUtilities.h"

// DRG
#include "DRG.H"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"


ODE_PFR::ODE_PFR(
	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXML, 
	OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapXML) :
	thermodynamicsMapXML_(thermodynamicsMapXML),
	kineticsMapXML_(kineticsMapXML)
{
	number_of_gas_species_ = thermodynamicsMapXML_.NumberOfSpecies();
	number_of_reactions_ = kineticsMapXML_.NumberOfReactions();
	number_of_equations_ = number_of_gas_species_ + 1 + 2;	// species and temperature + 2 dummy variables

	ChangeDimensions(number_of_gas_species_, &omegaStar_, true);
	ChangeDimensions(number_of_gas_species_, &xStar_, true);
	ChangeDimensions(number_of_gas_species_, &cStar_, true);
	ChangeDimensions(number_of_gas_species_, &RStar_, true);
	ChangeDimensions(number_of_reactions_, 	 &rStar_, true);

	checkMassFractions_ = false;	
	drgAnalysis_ = false;
	energyEquation_ = true;
	debug_ = false;
}

int ODE_PFR::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	if (drgAnalysis_ == false)
	{
		// Recover mass fractions
		if (checkMassFractions_ == true)
		{	for(unsigned int i=1;i<=number_of_gas_species_;++i)
				omegaStar_[i] = min(1., max(y[i], 0.)); 
		}
		else
		{
			for(unsigned int i=1;i<=number_of_gas_species_;++i)
				omegaStar_[i] = y[i];
		}
		// Recover temperature
		const double TStar_ = y[number_of_gas_species_+1];

		// Recover dummy variables
		// There are 2 additional dummy variables (not needed to recover them)
	
		// Calculates the pressure and the concentrations of species
		thermodynamicsMapXML_.MoleFractions_From_MassFractions(xStar_.GetHandle(), MWStar_, omegaStar_.GetHandle());
		cTotStar_ = P_Pa_/(PhysicalConstants::R_J_kmol * TStar_);
		rhoStar_ = cTotStar_*MWStar_;
		Product(cTotStar_, xStar_, &cStar_);

		// Calculates thermodynamic properties
		thermodynamicsMapXML_.SetTemperature(TStar_);
		thermodynamicsMapXML_.SetPressure(P_Pa_);
		cpStar_ = thermodynamicsMapXML_.cpMolar_Mixture_From_MoleFractions(xStar_.GetHandle());
		cpStar_/=MWStar_;
	
		// Calculates kinetics
		kineticsMapXML_.SetTemperature(TStar_);
		kineticsMapXML_.SetPressure(P_Pa_);
		kineticsMapXML_.ReactionEnthalpiesAndEntropies();
		kineticsMapXML_.KineticConstants();
		kineticsMapXML_.ReactionRates(cStar_.GetHandle());
		kineticsMapXML_.FormationRates(RStar_.GetHandle());

		// Recovering residuals
		for (unsigned int i=1;i<=number_of_gas_species_;++i)	
			dy[i] = thermodynamicsMapXML_.MW(i-1)*RStar_[i]/rhoStar_;

		if (energyEquation_ == true)
		{	
			const double Q = 0.; // radiation contribution
			const double QRStar_ = kineticsMapXML_.HeatRelease(RStar_.GetHandle());
			dy[number_of_gas_species_+1] = (QRStar_ - Q)/(rhoStar_*cpStar_);
		}
		else		
		{
			dy[number_of_gas_species_+1] = 0.;
		}

		// Dummy equations
		dy[number_of_gas_species_+2] = 0.;
		dy[number_of_gas_species_+3] = 0.;

		return 0;
	}
	else
	{
		// Recover mass fractions 
		if (checkMassFractions_ == true)
		{	
			for (unsigned int i=0;i<drg_->number_important_species();++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				omegaStar_[j] = max(y[i+1], 0.);
			}	
		}
		else
		{
			for (unsigned int i=0;i<drg_->number_important_species();++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				omegaStar_[j] = y[i+1];
			}
		}

		// Recover temperature
		unsigned int index_TStar = drg_->number_important_species()+1;
		const double TStar_ = y[index_TStar];
	
		// Calculates the pressure and the concentrations of species
		thermodynamicsMapXML_.MoleFractions_From_MassFractions(xStar_.GetHandle(), MWStar_, omegaStar_.GetHandle());
		cTotStar_ = P_Pa_/(PhysicalConstants::R_J_kmol * TStar_);
		rhoStar_ = cTotStar_*MWStar_;
		Product(cTotStar_, xStar_, &cStar_);

		// Calculates thermodynamic properties
		thermodynamicsMapXML_.SetTemperature(TStar_);
		thermodynamicsMapXML_.SetPressure(P_Pa_);
		cpStar_ = thermodynamicsMapXML_.cpMolar_Mixture_From_MoleFractions(xStar_.GetHandle());
		cpStar_/=MWStar_;
	
		// Calculates kinetics
		kineticsMapXML_.SetTemperature(TStar_);
		kineticsMapXML_.SetPressure(P_Pa_);
		kineticsMapXML_.ReactionEnthalpiesAndEntropies();
		kineticsMapXML_.KineticConstants();
		kineticsMapXML_.ReactionRates(cStar_.GetHandle());

		// Remove useless reactions
		for (unsigned int i=0;i<drg_->indices_unimportant_reactions().size();++i)
			rStar_[drg_->indices_unimportant_reactions()[i]+1] = 0.;

		// Formation rates
		kineticsMapXML_.stoichiometry().FormationRatesFromReactionRates(RStar_.GetHandle(), rStar_.GetHandle());
		

		// Recovering residuals
		for (unsigned int i=0;i<drg_->number_important_species();++i)	
		{
			const unsigned int j = drg_->indices_important_species()[i]+1;
			dy[i+1] = thermodynamicsMapXML_.MW(j-1)*RStar_[j]/rhoStar_;
		}

		if (energyEquation_ == true)
		{	
			const double Q = 0.; // radiation contribution
			const double QRStar_ = kineticsMapXML_.HeatRelease(RStar_.GetHandle());
			dy[index_TStar] = (QRStar_ - Q)/(rhoStar_*cpStar_);
		}
		else
		{
			dy[index_TStar] = 0.;
		}

		return 0;
	}
}

int ODE_PFR::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
	//std::cout << t << std::endl;
	return 0;
}

