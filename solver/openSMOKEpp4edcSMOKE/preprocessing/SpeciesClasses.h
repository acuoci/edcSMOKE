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
|   Copyright(C) 20202  Alberto Cuoci                                     |
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

#ifndef OpenSMOKE_SpeciesClasses_H
#define	OpenSMOKE_SpeciesClasses_H

namespace OpenSMOKE
{
	//!  Class for managing the classes of species in the Polimi kinetic mechanism
	/*!
		This class provides the tools for managing the classes of species in the Polimi kinetic mechanism
	*/

	class SpeciesClasses 
	{
	public:

		/**
		* Default constructor
		*/
		SpeciesClasses();

		/**
		* Default constructor
		*/
		SpeciesClasses(const std::vector<std::string>& species_class_names, const std::vector<int>& species_class_lines_abs);

		/**
		* Default constructor
		*/
		SpeciesClasses(const boost::filesystem::path& file_name);

		/**
		* Checks
		*/
		bool Checking();

		/**
		* Prints a summary on a stream
		*/
		void Summary(std::ostream& out);

		/**
		* Returns true if species classes exist
		*/
		bool is_active() const { return is_active_;  }

		/**
		*@brief Check if species j at line i is part of one of the species classes
		*@param j species index (from 0)
		*@param i species line (absolute value, from 1)
		*@return true if the species is part of a species class
		*/
		bool filter_species(const unsigned int j, const unsigned int i);

		/**
		*@brief Check if first species of each species class is among the species in the kinetic model 
		*@param names_species species in the kinetic model
		*@param firstspecies_names_ first species of each species class
		*@return true if the species is part of a species class
		*/
		bool filter_species_nobili(const std::vector<std::string>& names_species,
			const std::vector<std::string>& firstspecies_names_);

		/**
		*@brief Writes the species classes in the XML file
		*@param fOutput stream corresponding to the XML file
		*/
		void WriteXMLFile(std::stringstream& fOutput) const;

		/**
		*@brief Reads the species classes from the XML file
		*@param file_name file name corresponding to the XML file
		*/
		void ReadXMLFile(const boost::filesystem::path& file_name);

		const std::vector<std::string>& species_class_names() const  { return species_class_names_;  }

		// Da vedere forse sono gli species indices
		const std::vector<unsigned int>& species_indices(const unsigned int k) const { return species_indices_[k]; }
		
		unsigned int index_from_species_name(const std::string name);

    
	private:

	private:

		bool is_active_;
		unsigned int n_classes_;
		std::vector<std::string> species_class_names_;
		std::vector<int> species_class_lines_abs_;
		std::vector<std::vector<unsigned int>> species_indices_; // ?index based

	};

}

#include "SpeciesClasses.hpp"

#endif	/* OpenSMOKE_SpeciesClasses_H */

