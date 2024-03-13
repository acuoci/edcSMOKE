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

#ifndef OpenSMOKE_OnTheFlyROPA_H
#define	OpenSMOKE_OnTheFlyROPA_H

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Reactor utilities
#include "idealreactors/utilities/Utilities"

namespace OpenSMOKE
{
	class OnTheFlyROPA
	{
	public:

		OnTheFlyROPA(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
						OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, boost::filesystem::path path_kinetics_output);

		bool is_active() { return is_active_; }

		void Analyze(std::ofstream& fOut, const int current_step, const double t, const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEVectorDouble& omega, const OpenSMOKE::OpenSMOKEVectorDouble& omega0);

		void WriteHead(std::ofstream& fOut, const std::string reactor_type);
		void WriteHeadXML(const std::string reactor_type);

		void SetThreshold(const double threshold) { threshold_contributions_ = threshold; }
		void SetCompactMode(const bool compact_mode) { compact_mode_ = compact_mode; }
		void SetMergeForwardAndBackwardReactions(const bool flag) { merge_reaction_rates_ = flag; }
		void SetSpecies(const std::vector<std::string> list_of_species) { list_species_ = list_of_species; }
		void SetActive(const bool flag) { is_active_ = flag; }
		void SetReferenceSpecies(const std::string name) { reference_species_ = name; }
		void SetNumberOfSteps(const unsigned int number_steps) { number_steps_ = number_steps; next_step_ = number_steps; }
		void Setup(const boost::filesystem::path& path_kinetics_output);
		
		/**
		* @Allows to set CheckForROPA true through the set flag (allows to link Liquid and Gas ROPA)
		*/
		void SetAnalyze(const bool flag) { do_the_analysis_ = flag; };
		
		void Close();

	private:

		enum writing_policy_type { FIXED_TIME_STEPS_, LIST_CONVERSIONS, LIST_TIMES };
		bool CheckForROPA(const int current_step, const double time, OpenSMOKE::OpenSMOKEVectorDouble& conversions);

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&			kineticsMap_;			//!< kinetic map

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
		bool compact_mode_;
		int current_step_;
		double told_;

		bool is_active_;
		OpenSMOKE::OpenSMOKEVectorDouble R_;
		OpenSMOKE::OpenSMOKEVectorDouble P_;
		OpenSMOKE::OpenSMOKEVectorDouble D_;
		OpenSMOKE::OpenSMOKEVectorDouble rf_;
		OpenSMOKE::OpenSMOKEVectorDouble rb_;
		
		bool do_the_analysis_;
		
		bool is_write_xml_;
		std::ofstream xml_string_;
	};
}

#include "OnTheFlyROPA.hpp"

#endif	/* OpenSMOKE_OnTheFlyROPA_H */

