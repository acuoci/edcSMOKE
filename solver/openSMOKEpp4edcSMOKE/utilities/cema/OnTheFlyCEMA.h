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

#ifndef OpenSMOKE_OnTheFlyCEMA_H
#define	OpenSMOKE_OnTheFlyCEMA_H

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Reactor utilities
#include "idealreactors/utilities/Utilities"

namespace OpenSMOKE
{
	//!  A class for performing Chemical Explosive Mode Analysis (CEMA)
	/*!
	A class for performing Chemical Explosive Mode Analysis (CEMA)
	*/

	class OnTheFlyCEMA
	{
	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap map containing the thermodynamic data
		*@param kineticsMap map containing the kinetic mechanism
		*@param path_output path pointing to the folder where the output file will be written
		*/
		OnTheFlyCEMA(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
						OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
						const boost::filesystem::path path_output);

		/**
		*@brief Setup CEMA (Chemical Explosive Mode Analysis) from Dictionary
		*@param lapack_mode Lapack libraries used
		*@param additional_conservative_modes user defined additional conservative modes
		*@param output_species list of output species (ALL, NONE, or list of names)
		*@param output_reactions list of output reactions (ALL, NONE, or list of numbers, from 1)
		*/
		void Setup(const bool lapack_mode, const int additional_conservative_modes, const std::vector<std::string>& output_species, const std::vector<std::string> output_reactions);

		/**
		*@brief Setup CEMA (Chemical Explosive Mode Analysis) from Dictionary
		*@param dictionary name of dictionary
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Returns true is the CEMA is requested
		*/
		bool is_active() { return is_active_; }

		/**
		*@brief CEMA (Chemical Explosive Mode Analysis)
		*@param T temperature [K]
		*@param P_Pa pressure [Pa]
		*@param c concentrations [kmol/m3]
		*@param J dense Jacobian matrix
		*/
		void Calculate(const double t, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEMatrixDouble& J);

		/**
		*@brief Prepares the output files
		*/
		void PrepareOutputFiles();

		/**
		*@brief Writes the output file
		*@param t time [s]
		*@param T temperature [K]
		*@param P_PA pressure [Pa]
		*/
		void WriteOnFile(const double t, const double T, const double P_Pa);

		/**
		*@brief Closes the output files
		*/
		void CloseOutputFiles();

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&			kineticsMap_;			//!< kinetic map

		boost::filesystem::path path_output_;	//!< path pointing to the folder where the output file will be written							

		bool is_active_;					//!< true if the CEMA is active
		bool is_print_on_file_active_;		//!< true if printing on file is active
		bool print_ei_;						//!< true if the Explosion Indices (EI) have to be written on file
		bool print_pi_;						//!< true if the Participation Indices (PI) have to be written on file
		bool lapack_mode_;					//!< true if the Lapack libraries will be used to carry out CEMA (instead of Eigen C++)

		unsigned int NE_;					//!< number of elements
		unsigned int NS_;					//!< number of species
		unsigned int NR_;					//!< number of reactions
		unsigned int N_;					//!< number of non-conservative modes
		unsigned int Nstar_;				//!< total number of equations (i.e. eigenvalues)
		int N_add_conservative_modes_;		//!< user-defined additional conservative modes

		Eigen::MatrixXd Jomega_;			//!< Jacobian matrix
		Eigen::MatrixXd JomegaTranspose_;	//!< transpose Jacobian matrix

		Eigen::VectorXcd ae_;				//!< right eigenvector associated to CEM
		Eigen::VectorXcd be_;				//!< left eigenvector associated to CEM

		double cem_lambda_;					//!< Chemical Explosive Mode (CEM)
		unsigned int cem_index_;			//!< index of CEM (from 0)
		Eigen::VectorXd EI_;				//!< vector of Explosion Indices (PI)
		Eigen::VectorXd PI_;				//!< vector of Participation Indices (PI)

		std::ofstream fCEMA_;				//!< output file

		std::vector<unsigned int> indices_of_output_reactions_;		//!< indices of species for which the profiles are written on the ASCII file
		std::vector<unsigned int> indices_of_output_species_;		//!< indices of species for which the profiles are written on the ASCII file
		std::vector<unsigned int> widths_of_output_species_;        //!< width   of species for which the profiles are written on the ASCII file
		
	private:

		/**
		*@brief Memory allocation
		*/
		void MemoryAllocation();

		/**
		*@brief Search for conservative modes by looking at the eigenvalues of the Jacobian matrix
		*@param lambda_real eigenvalues (real part)
		*@param lambda_mod eigenvalues (module)
		*@param lambda_real_cleaned cleaned eigenevalues
		*@param lambda_original_associations references to the original indices
		*/
		void SearchConservativeModes(const std::vector<double>& lambda_real, const std::vector<double>& lambda_mod, std::vector<double>& lambda_real_cleaned, std::vector<size_t>& lambda_original_associations);
		
		/**
		*@brief Sorts the eigenvalues in decreasing order
		*@param lambda_real_cleaned cleaned eigenevalues
		*@param lambda_original_associations references to the original indices
		*@param lambda_original_associations_sorted references to the original indices after sorting
		*/
		void SortVectors(const std::vector<double>& lambda_real_cleaned, const std::vector<size_t>& lambda_original_associations, std::vector<size_t>& lambda_original_associations_sorted);
		
		/**
		*@brief Calculates the chemical explosive mode (CEM)
		*@param lambda_indices_sorted references to the original indices after sorting
		*@param lambda_real eigenvalues (real part)
		*/
		void CalculateChemicalExplosiveMode(const std::vector<size_t>& lambda_indices_sorted, const std::vector<double>& lambda_real);
		
		/**
		*@brief Calculates the Explosion Indices (EI)
		*/
		void CalculateExplosionIndices();

		/**
		*@brief Calculates the Participation Indices (PI)
		*@param T temperature [K]
		*@param P_Pa pressure [Pa]
		*@param c concentrations [kmol/m3]
		*/
		void CalculateParticipationIndices(const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c);

		/**
		*@brief CEMA (Chemical Explosive Mode Analysis) carried out through Lapack libraries
		*@param T temperature [K]
		*@param P_Pa pressure [Pa]
		*@param c concentrations [kmol/m3]
		*@param J dense Jacobian matrix
		*/
		void CalculateLapack(const double t, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEMatrixDouble& J);
		
		/**
		*@brief CEMA (Chemical Explosive Mode Analysis) carried out through Eigen C++ libraries
		*@param T temperature [K]
		*@param P_PA pressure [Pa]
		*@param c concentrations [kmol/m3]
		*@param J dense Jacobian matrix
		*/
		void CalculateEigen(const double t, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEMatrixDouble& J);
	};
}

#include "OnTheFlyCEMA.hpp"

#endif	/* OpenSMOKE_OnTheFlyCEMA_H */

