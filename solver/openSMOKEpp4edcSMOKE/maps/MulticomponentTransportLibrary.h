/*-----------------------------------------------------------------------*\
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
|   Copyright(C) 2022  Alberto Cuoci                                      |
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

#ifndef OpenSMOKEpp_MulticomponentTransportLibrary_H
#define OpenSMOKEpp_MulticomponentTransportLibrary_H

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// MuTLib
#include "maps/mutlib/MuTLib_Transport.h"

// Grammar utilities
#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"


namespace OpenSMOKE
{
	//!  A class to manage the multicomponent transport properties based on the MuTLib
	/*!
	A class to manage the multicomponent transport properties based on the MuTLib
	*/

	class MulticomponentTransportLibrary
	{

	public:

		/**
		*@brief Initialize
		*@param file_transport path to the transport file in CHEMKIN format
		*@param file_thermo path to the thermodynamic file in CHEMKIN format
		*@param names list of names species
		*@param M vector of molecular weights of species (in kg/kmol)
		*/
		void Initialize(const std::string file_transport, const std::string file_thermo, const std::vector<std::string> names, const double* M);


		/**
		*@brief Initialize
		*@param transportMapXML transport map
		*@param thermodynamicsMapXML thermodynamic map
		*/
		void Initialize(const OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMapXML, const OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXML);


		/**
		*@brief Setup from a dictionary
		*@param dictionary name of the dictionary
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);


		/**
		*@brief Set the operating conditions
		*@param T temperature (in K)
		*@param P pressure (in Pa)
		*@param X mole fractions
		*/
		void SetupMixture(const double T, const double P, const double* X);


		/**
		*@brief Calculates the Fick's fluxes
		*@param nablaX gradient of mole fractions (in 1/m)
		*@param jfick Fick's fluxes (in kg/m2/s)
		*/
		void GetFickFluxes(const Eigen::VectorXd& nablaX, Eigen::VectorXd& jfick);


		/**
		*@brief Calculates the Soret's fluxes
		*@param nablaLogT gradient of temperature divided by temperature (in 1/m)
		*@param jsoret Soret's fluxes (in kg/m2/s)
		*/
		void GetSoretFluxes(const double nablaLogT, Eigen::VectorXd& jsoret);

		/**
		*@brief Returns the thermal conductivity
		*@return the thermal conductivity (in W/m/K)
		*/
		double GetThermalConductivity() { return transport_.GetThCond(); }


		/**
		*@brief Returns true if the multicomponent (accurate) calculation of thermal conductivity is enabled
		*/
		bool UseThermalConductivity() const { return use_thermal_conductivity_; }


		/**
		*@brief Sets the exact calculation of Fick's and Soret's fluxes
		*/
		void SetExact(const bool flag);


		/**
		*@brief Sets the number of Neumann's steps for evaluation of Fick's fluxes
		*/
		void SetFickNeumannSteps(const unsigned int n);


		/**
		*@brief Sets the number of Neumann's steps for evaluation of Soret's fluxes
		*/
		void SetSoretNeumannSteps(const unsigned int n);


		/**
		*@brief Sets the swap policy for selection of species on which to force the closure
		*@param policy the swap policy: -2=most-abundant species, -1=last species, 0...N=fixed species
		*/
		void SetSwapPolicy(const int policy);


		/**
		*@brief Enables/Disables the multicomponent (accurate) calculation of thermal conductivity
		*/
		void SetUserThermalConductivity(const bool flag);


	private:

		MuTLib::MuTLib_Species     species_;		// MuTLib: main container combustion state of the mixture
		MuTLib::MuTLib_ThermoData  thermodata_;		// MuTLib: main container to JANAF thermodynamic properties interpolation
		MuTLib::MuTLib_KineticData kineticdata_;	// MuTLib: main container of transport parameters
		MuTLib::MuTLib_Transport   transport_;		// MuTLib: main container to coefficients, fluxes and conductivity members

		int N_;										// number of species (of the mixture)
		bool exact_;								// exact Fick diffusion coefficients (default: true)
		unsigned int Fickr_;						// Fick diffusion Neumann steps (default: 1)
		int FickM_;									// M number of species in model, 1+M for Fick diffusion coefficients (default: 0)
		double FickThr_;							// molar fraction threshold in model 1+M for Fick diffusion coefficients (default: 0.1)
		unsigned int Soretr_;						// Soret diffusion Neumann steps (default: 1)
		int SwapPolicy_;							// Swap policy
		bool use_thermal_conductivity_;				// use thermal conductivity
	};


	class Grammar_MulticomponentTransportLibrary : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{

	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SwapPolicy",
				OpenSMOKE::SINGLE_STRING,
				"Type of swap policy: most-abundant (default) | last | species name",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FickNeumannSteps",
				OpenSMOKE::SINGLE_INT,
				"Number of Neumann steps for Fick coefficients (default: 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SoretNeumannSteps",
				OpenSMOKE::SINGLE_INT,
				"Number of Neumann steps for Soret coefficients (default: 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Exact",
				OpenSMOKE::SINGLE_BOOL,
				"If true, the calculations are done exactly (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@UseMulticomponentThermalConductivity",
				OpenSMOKE::SINGLE_BOOL,
				"If true, the thermal conductivity is calculated using the multicomponent transport approach, instead of polynomial fitting (default: true)",
				false));
		}
	};

}

#include "MulticomponentTransportLibrary.hpp"

#endif /* OpenSMOKEpp_MulticomponentTransportLibrary_H */
