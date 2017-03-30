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



NLS_PSR::NLS_PSR(
	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXML, 
	OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapXML) :
	thermodynamicsMapXML_(thermodynamicsMapXML),
	kineticsMapXML_(kineticsMapXML)
{
	number_of_gas_species_ = thermodynamicsMapXML_.NumberOfSpecies();
	number_of_reactions_ = kineticsMapXML_.NumberOfReactions();
	number_of_equations_ = number_of_gas_species_-1;

	ChangeDimensions(number_of_gas_species_, &omegaSurr_, true);
	ChangeDimensions(number_of_gas_species_, &omegaStar_, true);
	ChangeDimensions(number_of_gas_species_, &xStar_, true);
	ChangeDimensions(number_of_gas_species_, &cStar_, true);
	ChangeDimensions(number_of_gas_species_, &RStar_, true);
	ChangeDimensions(number_of_gas_species_, &RStar_, true);
	ChangeDimensions(number_of_reactions_, 	 &rStar_, true);

	checkMassFractions_ = false;	
	drgAnalysis_ = false;
	index_max_species_ = 0;
}

int NLS_PSR::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	if (drgAnalysis_ == false)
	{
		// Recover mass fractions
		{
			double sum = 0.;
			unsigned int count = 1;
			if (checkMassFractions_ == true)
			{	
				for(unsigned int i=1;i<index_max_species_;++i)
				{
					omegaStar_[i] = max(y[count++], 0.);
					sum += omegaStar_[i];
				}
				for(unsigned int i=index_max_species_+1;i<=number_of_gas_species_;++i)
				{
					omegaStar_[i] = max(y[count++], 0.);
					sum += omegaStar_[i];
				}
			}
			else
			{
				for(unsigned int i=1;i<index_max_species_;++i)
				{
					omegaStar_[i] = y[count++];
					sum += omegaStar_[i];
				}
				for(unsigned int i=index_max_species_+1;i<=number_of_gas_species_;++i)
				{
					omegaStar_[i] = y[count++];
					sum += omegaStar_[i];
				}
			}

			omegaStar_[index_max_species_] = max(1.-sum,0.);
		}

		// Composition of fine structure
		thermodynamicsMapXML_.MoleFractions_From_MassFractions(xStar_.GetHandle(), MWStar_, omegaStar_.GetHandle());

		// Fine structure temperature from composition
		const double hStar = hMean_;
		TStar_ = thermodynamicsMapXML_.GetTemperatureFromEnthalpyAndMoleFractions(hStar*MWStar_, P_Pa_, xStar_.GetHandle(), TStar_);
		TStar_ = ( TStar_<250.  ) ?  250. : TStar_;
		TStar_ = ( TStar_>4000. ) ? 4000. : TStar_;

		// Calculates the pressure and the concentrations of species
		cTotStar_ = P_Pa_/(PhysicalConstants::R_J_kmol * TStar_);
		rhoStar_ = cTotStar_*MWStar_;
		Product(cTotStar_, xStar_, &cStar_);

		// Calculates kinetics
		thermodynamicsMapXML_.SetTemperature(TStar_);
		thermodynamicsMapXML_.SetPressure(P_Pa_);
		kineticsMapXML_.SetTemperature(TStar_);
		kineticsMapXML_.SetPressure(P_Pa_);
		kineticsMapXML_.ReactionEnthalpiesAndEntropies();
		kineticsMapXML_.KineticConstants();
		kineticsMapXML_.FormationRates(RStar_.GetHandle());

		// Surrounding environment composition
		for(unsigned int i=1;i<=number_of_gas_species_;i++)
			omegaSurr_[i] = (omegaMean_[i] - omegaStar_[i]*gammaStar_)/(1.-gammaStar_);

		// Mass balance equations
		unsigned int count = 1;
		for(unsigned int i=1;i<index_max_species_;++i)
			dy[count++] = thermodynamicsMapXML_.MW(i-1)*RStar_[i]/rhoStar_ + mDotStar_*(omegaSurr_[i]-omegaStar_[i]);
		for(unsigned int i=index_max_species_+1;i<=number_of_gas_species_;++i)
			dy[count++] = thermodynamicsMapXML_.MW(i-1)*RStar_[i]/rhoStar_ + mDotStar_*(omegaSurr_[i]-omegaStar_[i]);

		return 0;
	}
	else
	{
		//Recover mass fractions 
		{
			double sum = 0.;
			for (unsigned int i=0;i<drg_->number_unimportant_species();++i)	
			{
				const unsigned int j = drg_->indices_unimportant_species()[i]+1;
				sum += omegaStar_[j];
			}	

			if (checkMassFractions_ == true)
			{	
				unsigned int count = 1;
				for (unsigned int i=0;i<index_max_species_;++i)	
				{
					const unsigned int j = drg_->indices_important_species()[i]+1;
					omegaStar_[j] = max(y[count++], 0.);
					sum += omegaStar_[j];
				}
				for (unsigned int i=index_max_species_+1;i<drg_->number_important_species();++i)	
				{
					const unsigned int j = drg_->indices_important_species()[i]+1;
					omegaStar_[j] = max(y[count++], 0.);
					sum += omegaStar_[j];
				}
			}
			else
			{
				unsigned int count = 1;
				for (unsigned int i=0;i<index_max_species_;++i)	
				{
					const unsigned int j = drg_->indices_important_species()[i]+1;
					omegaStar_[j] = y[count++];
					sum += omegaStar_[j];
				}
				for (unsigned int i=index_max_species_+1;i<drg_->number_important_species();++i)	
				{
					const unsigned int j = drg_->indices_important_species()[i]+1;
					omegaStar_[j] = y[count++];
					sum += omegaStar_[j];
				}
			}
			omegaStar_[drg_->indices_important_species()[index_max_species_]+1] = max(1.-sum,0.);
		}

		// Composition of fine structure
		thermodynamicsMapXML_.MoleFractions_From_MassFractions(xStar_.GetHandle(), MWStar_, omegaStar_.GetHandle());

		// Fine structure temperature from composition
		const double hStar = hMean_;
		TStar_ = thermodynamicsMapXML_.GetTemperatureFromEnthalpyAndMoleFractions(hStar*MWStar_, P_Pa_, xStar_.GetHandle(), TStar_);
		TStar_ = ( TStar_<250.  ) ?  250. : TStar_;
		TStar_ = ( TStar_>4000. ) ? 4000. : TStar_;	

		// Calculates the pressure and the concentrations of species
		cTotStar_ = P_Pa_/(PhysicalConstants::R_J_kmol * TStar_);
		rhoStar_ = cTotStar_*MWStar_;
		Product(cTotStar_, xStar_, &cStar_);

		// Calculates kinetics
		thermodynamicsMapXML_.SetTemperature(TStar_);
		thermodynamicsMapXML_.SetPressure(P_Pa_);
		kineticsMapXML_.SetTemperature(TStar_);
		kineticsMapXML_.SetPressure(P_Pa_);
		kineticsMapXML_.ReactionEnthalpiesAndEntropies();
		kineticsMapXML_.KineticConstants();
		kineticsMapXML_.ReactionRates(cStar_.GetHandle());
		kineticsMapXML_.GiveMeReactionRates(rStar_.GetHandle());

		// Remove useless reactions
		for (unsigned int i=0;i<drg_->indices_unimportant_reactions().size();++i)
			rStar_[drg_->indices_unimportant_reactions()[i]+1] = 0.;

		// Formation rates
		kineticsMapXML_.stoichiometry().FormationRatesFromReactionRates(RStar_.GetHandle(), rStar_.GetHandle());

		// Composition of surrouding environment
		for(unsigned int i=0;i<drg_->number_important_species();++i)
		{
			const unsigned int j = drg_->indices_important_species()[i]+1;
			omegaSurr_[j] = (omegaMean_[j] - omegaStar_[j]*gammaStar_)/(1.-gammaStar_);
		}

		// Recovering residuals
		{
			unsigned int count = 1;
			for (unsigned int i=0;i<index_max_species_;++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				dy[count++] = thermodynamicsMapXML_.MW(j-1)*RStar_[j]/rhoStar_ + mDotStar_*(omegaSurr_[j]-omegaStar_[j]);
			}
			for (unsigned int i=index_max_species_+1;i<drg_->number_important_species();++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				dy[count++] = thermodynamicsMapXML_.MW(j-1)*RStar_[j]/rhoStar_ + mDotStar_*(omegaSurr_[j]-omegaStar_[j]);
			}
		}

	}
}

void NLS_PSR::ReconstructData(const Eigen::VectorXd& y, double& TStar, double& omegaReconstructed)
{
	if (drgAnalysis_ == false)
	{
		unsigned int count = 0;
		double sum = 0.;
		for(unsigned int i=1;i<index_max_species_;++i)
		{
			omegaStar_[i] = y(count++);
			sum += omegaStar_[i];
		}
		for(unsigned int i=index_max_species_+1;i<=number_of_gas_species_;++i)
		{
			omegaStar_[i] = y(count++);
			sum += omegaStar_[i];
		}

		omegaStar_[index_max_species_] = max(1.-sum,0.);
		omegaReconstructed = omegaStar_[index_max_species_];
	}
	else
	{
		double sum = 0.;
		for (unsigned int i=0;i<drg_->number_unimportant_species();++i)	
		{
			const unsigned int j = drg_->indices_unimportant_species()[i]+1;
			sum += omegaStar_[j];
		}	

		if (checkMassFractions_ == true)
		{	
			unsigned int count = 0;
			for (unsigned int i=0;i<index_max_species_;++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				omegaStar_[j] = max(y(count++), 0.);
				sum += omegaStar_[j];
			}
			for (unsigned int i=index_max_species_+1;i<drg_->number_important_species();++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				omegaStar_[j] = max(y(count++), 0.);
				sum += omegaStar_[j];
			}
		}
		else
		{
			unsigned int count = 0;
			for (unsigned int i=0;i<index_max_species_;++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				omegaStar_[j] = y(count++);
				sum += omegaStar_[j];
			}
			for (unsigned int i=index_max_species_+1;i<drg_->number_important_species();++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				omegaStar_[j] = y(count++);
				sum += omegaStar_[j];
			}
		}

		omegaStar_[drg_->indices_important_species()[index_max_species_]+1] = max(1.-sum,0.);
		omegaReconstructed = omegaStar_[drg_->indices_important_species()[index_max_species_]+1];
	}

	// Composition of fine structure
	thermodynamicsMapXML_.MoleFractions_From_MassFractions(xStar_.GetHandle(), MWStar_, omegaStar_.GetHandle());

	// Fine structure temperature from composition
	const double hStar = hMean_;
	TStar = thermodynamicsMapXML_.GetTemperatureFromEnthalpyAndMoleFractions(hStar*MWStar_, P_Pa_, xStar_.GetHandle(), TStar_);
}

int NLS_PSR::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
	return 0;
}

