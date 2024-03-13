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

#ifndef OpenSMOKE_OnTheFlyROPA_Liquid_H
#define	OpenSMOKE_OnTheFlyROPA_Liquid_H

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Reactor utilities
#include "idealreactors/utilities/Utilities"

namespace OpenSMOKE
{
	class OnTheFlyROPA_Liquid
	{
	public:

		OnTheFlyROPA_Liquid(	OpenSMOKE::ThermodynamicsMap_Liquid_CHEMKIN& thermodynamicsMap,
						OpenSMOKE::KineticsMap_Liquid_CHEMKIN& kineticsMap);

		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, boost::filesystem::path path_kinetics_output);

		bool is_active() { return is_active_; }

		/**
		*@brief Perform the ROPA analysis for liquid-phase:
		*@param fOut			the file stream were ROPAL is pirnted
		*@param current_step	the current iteration step
		*@param t				the physical time [s]
		*@param T				the liquid-phase temperature [K] 
		*@param P				the liquid-phase pressure [Pa] 
		*@param c				the vector of concentrations for gas-phase species in liquid-phase (to be better implemented) [kmol/m3] 
		*@param cL				the vector of concentrations of liquid-phase species [kmol/m3] 
		*@param massL			the vector of masses of liquid-phase species [kg/kg0]
		*@param omega0			the vector of initial mass fractions of liquid-phase species [-] 
		*@param mL_ev			the vector of molar evaporation rates per unit volume [kmol/m3/s], liquid-index based 
		*/
		void Analyze(std::ofstream& fOut, const int current_step, const double t, const double T, const double P,
			const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEVectorDouble& cL,
			const OpenSMOKE::OpenSMOKEVectorDouble& massL, const OpenSMOKE::OpenSMOKEVectorDouble& omega0,
			const OpenSMOKE::OpenSMOKEVectorDouble& mL_ev);

		/**
		*@brief Prints information on the ROPAL file for the total liquid mass variation:
		*@param fOut			the file stream were ROPAL is pirnted
		*@param mL_ev			the vector of mass evaporation rates per unit volume [kmol/m3/s], liquid-index based
		*/
		void Print_Evaporation(std::ofstream& fOut, const OpenSMOKE::OpenSMOKEVectorDouble& mL_ev);

		/**
		*@brief Prints information on the ROPAL file for the total liquid mass variation:
		*@param fOut					the file stream were ROPAL is pirnted
		*@param mL_ev					the vector of mass evaporation rates per unit volume [kg/m3/s], liquid-index based
		*@param Omega_Gas_from_Liquid	the vector of mass formation rates per unit volume of gas-species formed from liquid-phase reactions [kg/m3/s], gas-index based
		*/
		void Print_Evaporation(std::ofstream& fOut, const OpenSMOKE::OpenSMOKEVectorDouble& mL_ev, 
			const OpenSMOKE::OpenSMOKEVectorDouble& Omega_Gas_from_Liquid);

		void WriteHead(std::ofstream& fOut, const std::string reactor_type);

		void SetThreshold(const double threshold) { threshold_contributions_ = threshold; }
		void SetCompactMode(const bool compact_mode) { compact_mode_ = compact_mode; }
		void SetMergeForwardAndBackwardReactions(const bool flag) { merge_reaction_rates_ = flag; }
		void SetSpecies(const std::vector<std::string> list_of_species) { list_species_ = list_of_species; }
		void SetActive(const bool flag) { is_active_ = flag; }
		void SetReferenceSpecies(const std::string name) { reference_species_ = name; }
		void SetNumberOfSteps(const unsigned int number_steps) { number_steps_ = number_steps; next_step_ = number_steps; }
		void Setup(const boost::filesystem::path& path_kinetics_output);

		void SetIndexLiquidSpecies(const Eigen::VectorXi& index_fuel_species, const Eigen::VectorXi& index_Lmixture_Lkin);
		
		/**
		* @brief Returns the name of the reference species (without "(L)")
		*/
		std::string	GetRefSpecies() { return reference_species_; };

		/**
		* @brief Returns flag if at the current iteration the ROPAL analysis was performed
		*/
		bool IsAnalyzing() { return am_i_analyzing; };

	private:

		enum writing_policy_type { FIXED_TIME_STEPS_, LIST_CONVERSIONS, LIST_TIMES };
		bool CheckForROPA(const int current_step, const double time, OpenSMOKE::OpenSMOKEVectorDouble& conversions);

	private:
		
		OpenSMOKE::ThermodynamicsMap_Liquid_CHEMKIN& thermodynamicsMap_;
		OpenSMOKE::KineticsMap_Liquid_CHEMKIN& kineticsMap_;

		Eigen::VectorXi index_of_fuel_species_;
		Eigen::VectorXi index_Lmixture_Lkin_;

		std::vector<std::string> reaction_names_;
		std::vector<std::string> list_species_;
		std::vector<double> list_conversions_;
		std::vector<double> list_times_;
		int number_steps_;
		double threshold_contributions_;
		double threshold_formation_rate_;

		writing_policy_type writing_policy_;
		bool merge_reaction_rates_;
		std::string reference_species_;

		int next_step_;
		int next_index_conversion_;
		int next_index_time_;
		int index_ref_species_;
		bool compact_mode_;

		bool is_active_;
		bool total_liquid_conversion;
		bool am_i_analyzing;
		
		OpenSMOKE::OpenSMOKEVectorDouble Rg_;
		OpenSMOKE::OpenSMOKEVectorDouble Rl_;
		OpenSMOKE::OpenSMOKEVectorDouble P_;
		OpenSMOKE::OpenSMOKEVectorDouble D_;
		OpenSMOKE::OpenSMOKEVectorDouble rf_;
		OpenSMOKE::OpenSMOKEVectorDouble rb_;
	};
}

#include "OnTheFlyROPA_Liquid.hpp"

#endif	/* OpenSMOKE_OnTheFlyROPA_Liquid_H */

