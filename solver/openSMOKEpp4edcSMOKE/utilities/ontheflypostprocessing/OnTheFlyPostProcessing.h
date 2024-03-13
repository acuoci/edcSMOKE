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

#ifndef OpenSMOKE_OnTheFlyPostProcessing_H
#define	OpenSMOKE_OnTheFlyPostProcessing_H

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Reactor utilities
#include "idealreactors/utilities/Utilities"
#include "utilities/Utilities"

namespace OpenSMOKE
{
	//!  A class for performing on the fly post-processing on kinetics
	/*!
	A class for performing on the fly post-processing on kinetics
	*/

	class OnTheFlyPostProcessing
	{
	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap map containing the thermodynamic data
		*@param kineticsMap map containing the kinetic mechanism
		*@param path_output path pointing to the folder where the output file will be written
		*/
		OnTheFlyPostProcessing(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								const boost::filesystem::path path_output);

		/**
		*@brief Setup from Dictionary
		*@param dictionary name of dictionary
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Returns true is the CEMA is requested
		*/
		bool is_active() { return is_active_; }

		/**
		*@brief Prepares the output files
		*/
		void PrepareOutputFiles();

		/**
		*@brief Prepares the output files with additional variables
		*/
		void PrepareOutputFiles(const std::vector<std::string>& additional);

		/**
		*@brief Writes the output file
		*@param t time [s]
		*@param x coordinate [m]
		*@param y coordinate [m]
		*@param z coordinate [m]
		*@param T temperature [K]
		*@param P_Pa pressure [Pa]
		*@param omega mass fractions
		*/
		void WriteOnFile(const double t, const double x, const double y, const double z, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omega);

		/**
		*@brief Writes the output file
		*@param t time [s]
		*@param x coordinate [m]
		*@param y coordinate [m]
		*@param z coordinate [m]
		*@param T temperature [K]
		*@param P_Pa pressure [Pa]
		*@param omega mass fractions
		*@param additional variables
		*/
		void WriteOnFile(const double t, const double x, const double y, const double z, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omega, const std::vector<double>& additional);
		
		/**
		*@brief Closes the output files
		*/
		void CloseOutputFiles();

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&			kineticsMap_;			//!< kinetic map
		boost::filesystem::path					path_output_;			//!< path pointing to the folder where the output file will be written							

		bool is_active_;					//!< true if the CEMA is active
		bool print_formation_rates_;		//!< true if the formation rates are required
		bool print_reaction_rates_;			//!< true if the reaction rates are required
		bool formation_rates_moles_;		//!< true if the formation rates are required in [kmol/m3s/] instead of [kg/m3/s]
		bool class_grouping_;				//!< true if the class grouping is required

		unsigned int NS_;					//!< number of species
		unsigned int NR_;					//!< number of reactions

		OpenSMOKE::OutputFileColumns fFormationRates_;		//!< output file for formation rates [kmol/m3/s or kg/m3/s]
		OpenSMOKE::OutputFileColumns fReactionRates_;		//!< output file for reaction rates [kmol/m3/s]
		OpenSMOKE::OutputFileColumns fClassGrouping_;		//!< output file for class grouping

		std::vector<unsigned int> indices_of_reaction_rates_;				//!< indices of species for which the profiles are written on the ASCII file
		std::vector<unsigned int> indices_of_formation_rates_species_;		//!< indices of species for which the profiles are written on the ASCII file
		
		OpenSMOKE::OpenSMOKEVectorDouble R_;	//!< formation rates of species [kmol/m3/s]
		OpenSMOKE::OpenSMOKEVectorDouble x_;	//!< mole fractions
		OpenSMOKE::OpenSMOKEVectorDouble c_;	//!< concentrations [kmol/m3]

		std::vector<std::string> class_grouping_names_families_;
		std::vector<unsigned int> class_grouping_indices_species_;
		std::vector<unsigned int> class_grouping_indices_families_;

	private:

		/**
		*@brief Memory allocation
		*/
		void MemoryAllocation();
	};
}

#include "OnTheFlyPostProcessing.hpp"

#endif	/* OpenSMOKE_OnTheFlyPostProcessing_H */

