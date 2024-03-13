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
|   Copyright(C) 2019  Alberto Cuoci                                      |
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

#ifndef OpenSMOKE_KineticsModifier_H
#define	OpenSMOKE_KineticsModifier_H

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

namespace OpenSMOKE
{
	class KineticsModifier
	{
	public:

		KineticsModifier();

		/**
		*@brief Setup from dictionary
		*@param dictionary dictionary name
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Setup from indices (frequency factors only)
		*@param indices reaction indices (0-based)
		*@param coefficients correction coefficients
		*/
		void SetupFromIndices(const std::vector<unsigned int>& indices, const std::vector<double>& coefficients);
		
		/**
		*@brief Return true if the kinetic modifier is turned on
		*/
		bool is_active() { return is_active_; }

		/**
		*@brief Returns the indices (0-based) of species to be kept frozen during the simulation
		*/
		const std::vector<unsigned int>& list_frozen_species_index() const { return list_frozen_species_index_; }

		/**
		*@brief Modifies the kinetic parameters directly on the kinetic map
		*@param thermodynamicsMap thermodynamic map
		*@param kineticsMap kinetic map
		*/
		void Setup(ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

	private:

		bool is_active_;
		bool is_verbose_;

		std::vector<unsigned int>	list_A_index_;
		std::vector<double>			list_A_value_;

		std::vector<unsigned int>	list_n_index_;
		std::vector<double>			list_n_value_;

		std::vector<unsigned int>	list_EoverR_index_;
		std::vector<double>			list_EoverR_value_;

		std::vector<unsigned int>	list_ThirdBody_index_;
		std::vector<std::string>	list_ThirdBody_name_;
		std::vector<double>			list_ThirdBody_value_;

		std::vector<unsigned int>	list_A_multiplier_index_;
		std::vector<double>			list_A_multiplier_value_;

		std::vector<std::string>	list_frozen_species_name_;
		std::vector<unsigned int>	list_frozen_species_index_;
	};
}

#include "KineticsModifier.hpp"

#endif	/* OpenSMOKE_KineticsModifier_H */

