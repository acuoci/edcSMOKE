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

#ifndef OpenSMOKE_PolimiSootClasses_H
#define	OpenSMOKE_PolimiSootClasses_H

namespace OpenSMOKE
{
	//!  Class for managing the classes of soot reactions in the Polimi kinetic mechanism
	/*!
		 This class provides the tools for managing the classes of soot reactions in the Polimi kinetic mechanism
	*/

	class PolimiSootClasses 
	{
	public:

		/**
		* Default constructor
		*/
		PolimiSootClasses();

		/**
		* Default constructor
		*/
		PolimiSootClasses(const std::vector<std::string>& soot_class_names, const std::vector<int>& soot_class_lines_abs);

		/**
		* Default constructor
		*/
		PolimiSootClasses(const boost::filesystem::path& file_name);

		/**
		* Checks
		*/
		bool Checking();

		/**
		* Prints a summary on a stream
		*/
		void Summary(std::ostream& out);

		/**
		* Returns true if soot classes exist
		*/
		bool is_active() const { return is_active_;  }

		/**
		*@brief Check if reaction j at line i is part of one of the soot classes
		*@param j reaction index (from 0)
		*@param i reaction line (absolute value, from 1)
		*@return true if the reaction is part of a soot class
		*/
		bool filter_reaction(const unsigned int j, const unsigned int i);

		/**
		*@brief Writes the soot classes in the XML file
		*@param fOutput stream corresponding to the XML file
		*/
		void WriteXMLFile(std::stringstream& fOutput) const;

		/**
		*@brief Reads the soot classes from the XML file
		*@param file_name file name corresponding to the XML file
		*/
		void ReadXMLFile(const boost::filesystem::path& file_name);

		const std::vector<std::string>& soot_class_names() const  { return soot_class_names_;  }
		const std::vector<int>& reaction_indices(const unsigned int k) const { return reaction_indices_[k]; }
		unsigned int index_from_class_name(const std::string name);

    
	private:

	private:

		bool is_active_;
		unsigned int n_classes_;
		std::vector<std::string> soot_class_names_;
		std::vector<int> soot_class_lines_abs_;
		std::vector<std::vector<int>> reaction_indices_; // ?index based

	};

}

#include "PolimiSootClasses.hpp"

#endif	/* OpenSMOKE_PolimiSootClasses_H */

