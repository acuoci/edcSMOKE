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

#ifndef OpenSMOKE_InputFileCHEMKIN_H
#define	OpenSMOKE_InputFileCHEMKIN_H

#include <boost/filesystem.hpp>

namespace OpenSMOKE
{
	//!  A class for storing input files in CHEMKIN format
	/*!
		 This class basically provides the possibility to store in memory a CHEMKIN input file
		 by distinguishing among the blank lines, the comment lines and the so-called 
		 good lines, i.e. the lines which are useful to be interpreted
	*/

	class InputFileCHEMKIN {
	public:

		/**
		*@brief Constructor: default
		*/
		InputFileCHEMKIN();

		/**
		*@brief Constructor: from a file
		*/
		InputFileCHEMKIN(const std::string file_name);

		/**
		*@brief Default copy constructor
		*/
		InputFileCHEMKIN(const InputFileCHEMKIN& orig);

		/**
		*@brief Default destructor
		*/
		virtual ~InputFileCHEMKIN();

	public:
    
		/**
		*@brief Returns the number of lines in the file
		*/
		int number_of_lines() const { return number_of_lines_; }

		/**
		*@brief Returns the number of blank lines in the file
		*/
		int number_of_blank_lines() const { return number_of_blank_lines_; }

		/**
		*@brief Returns the number of useful lines in the file (i.e. lines which
				have to be interpreted). Of course blank lines and comment lines are
				automatically excluded from this list
		*/
		int number_of_good_lines() const { return number_of_good_lines_; }

		/**
		*@brief Returns the number of string comment (only) lines
		*/
		int number_of_strong_comment_lines() const { return number_of_strong_comment_lines_; }
    
		/**
		*@brief Returns the indices of good lines (i.e. lines to be interpreted)
		*/
		const std::vector<int>& indices_of_good_lines() const { return indices_of_good_lines_; }

		/**
		*@brief Returns the indices of blanck lines
		*/
		const std::vector<int>& indices_of_blank_lines() const { return indices_of_blank_lines_; }
    
		/**
		*@brief Returns the good lines (i.e. lines to be interpreted)
		*/
		const std::vector<std::string>& good_lines() const { return good_lines_; }

		/**
		*@brief Returns the indices of lines containing strong comments only 
		*/
		const std::vector<int>& indices_of_strong_comment_lines() const { return indices_of_strong_comment_lines_; }

		/**
		*@brief Returns the lines containing strong comments only 
		*/
		const std::vector<std::string>& strong_comment_lines() const { return strong_comment_lines_; }

		/**
		*@brief Returns the blank lines (i.e. lines to be interpreted)
		*/
		const std::vector<std::string>& blank_lines() const { return blank_lines_; }

		/**
		*@brief Returns the strong comments 
		*/
		const std::vector<std::string>& strong_comments() const { return strong_comments_; }

		/**
		*@brief Returns the names of soot classes
		*/
		const std::vector<std::string>& soot_class_names() const { return soot_class_names_; }

		/**
		*@brief Returns the indices of lines (absolute values) from which the different soot classes start
		*/
		const std::vector<int>& soot_class_lines_abs() const { return soot_class_lines_abs_; }

		/**
		*@brief Returns the names of species classes
		*/
		const std::vector<std::string>& species_class_names() const { return species_class_names_; }

		/**
		*@brief Returns the names of the first species in each species classes
		*/
		const std::vector<std::string>& firstspecies_names() const { return firstspecies_names_; }

		/**
		*@brief Returns the indices of lines (absolute values) from which the different soot classes start
		*/
		const std::vector<int>& species_class_lines_abs() const { return species_class_lines_abs_; }

		/**
		*@brief Returns the names of reaction classes
		*/
		const std::vector<std::string>& reaction_mainclass_names() const { return reaction_mainclass_names_; }

		/**
		*@brief Returns the indices of lines (absolute values) from which the different reaction classes start
		*/
		const std::vector<int>& reaction_mainclass_lines_abs() const { return reaction_mainclass_lines_abs_; }

		/**
		*@brief Returns the names of reaction classes
		*/
		const std::vector<std::string>& reaction_subclass_names() const { return reaction_subclass_names_; }

		/**
		*@brief Returns the indices of lines (absolute values) from which the different reaction classes start
		*/
		const std::vector<int>& reaction_subclass_lines_abs() const { return reaction_subclass_lines_abs_; }

		/**
		*@brief Returns the names of endreaction classes
		*/
		const std::vector<std::string>& endreaction_mainclass_names() const { return endreaction_mainclass_names_; }

		/**
		*@brief Returns the indices of lines (absolute values) from which the different reaction classes end
		*/
		const std::vector<int>& endreaction_mainclass_lines_abs() const { return endreaction_mainclass_lines_abs_; }

		/**
		*@brief Returns the names of endreaction classes
		*/
		const std::vector<std::string>& endreaction_subclass_names() const { return endreaction_subclass_names_; }

		/**
		*@brief Returns the indices of lines (absolute values) from which the different reaction classes end
		*/
		const std::vector<int>& endreaction_subclass_lines_abs() const { return endreaction_subclass_lines_abs_; }
    
		/**
		*@brief Returns the name of the file
		*/
		boost::filesystem::path file_name() const { return file_name_->filename(); }

		/**
		*@brief Returns the path of the file (without the name of the file)
		*/
		boost::filesystem::path folder_path() const { return file_name_->parent_path(); }

		/**
		*@brief Returns the size of the file in bytes
		*/
		boost::uintmax_t size() const { return boost::filesystem::file_size(*file_name_); }

		/**
		*@brief Writes information about the status of the file on output stream
		*/
		void Status(std::ostream &fOut) const;

		void ConvertGoodLinesIntoBlankLines(const std::vector<unsigned int> indices);

		/**
		*@brief Returns the good lines (i.e. lines to be interpreted)
		*/
		void ReplaceGoodLine(const unsigned int i, const std::string& line) { good_lines_[i] = line; }
    
	private:
    
		int number_of_lines_;							//!< total number of lines
		int number_of_good_lines_;						//!< total number of good lines
		int number_of_blank_lines_;						//!< total number of blank lines
		int number_of_strong_comment_lines_;					//!< total number of strong comment (only) lines
		std::vector<int> indices_of_good_lines_;		//!< indices of good lines
		std::vector<int> indices_of_blank_lines_;		//!< indices of blank lines
		std::vector<std::string> good_lines_;			//!< good lines
		std::vector<std::string> blank_lines_;			//!< blank lines
		std::vector<std::string> strong_comments_;		//!< strong comments
		std::vector<std::string> strong_comment_lines_;		//!< strong comment lines
		std::vector<int> indices_of_strong_comment_lines_;	//!< indices of strong comment lines
		
		std::vector<std::string> soot_class_names_;		//!< soot class names (if any)
		std::vector<int> soot_class_lines_abs_;			//!< soot class line indices (absolute, if any)

		std::vector<std::string> species_class_names_;		//!< species class names (if any)
		std::vector<int> species_class_lines_abs_;			//!< species class line indices (absolute, if any)
		std::vector<std::string> firstspecies_names_;		//!< name of the first species in each species class


		std::vector<std::string> reaction_mainclass_names_;		//!< reaction class names (if any)
		std::vector<int> reaction_mainclass_lines_abs_;			//!< reaction class line indices (absolute, if any)
		std::vector<std::string> reaction_subclass_names_;		//!< reaction class names (if any)
		std::vector<int> reaction_subclass_lines_abs_;			//!< reaction class line indices (absolute, if any)

		std::vector<std::string> endreaction_mainclass_names_;		//!< end reaction class names (if any)
		std::vector<int> endreaction_mainclass_lines_abs_;			//!< end reaction class line indices (absolute, if any)
		std::vector<std::string> endreaction_subclass_names_;		//!< end reaction class names (if any)
		std::vector<int> endreaction_subclass_lines_abs_;			//!< end reaction class line indices (absolute, if any)

		boost::filesystem::path* file_name_;			//!< full name (path + file name)
    
	};
}

#include "InputFileCHEMKIN.hpp"

#endif	/* OpenSMOKE_InputFileCHEMKIN_H */

