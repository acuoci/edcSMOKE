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
|   Copyright(C) 2023  Alberto Cuoci                                      |
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

#ifndef OpenSMOKE_OutputFileColumns_H
#define OpenSMOKE_OutputFileColumns_H

namespace OpenSMOKE
{
	class OutputFileColumns
	{
	public:

		/**
		*@brief Default constructor
		*/
		OutputFileColumns();

		/**
		*@brief Default constructor
		*@param file_name file name
		*/
		OutputFileColumns(const boost::filesystem::path& file_name);
	
		/**
		*@brief Default constructor
		*@param file_name file name
		*/
		OutputFileColumns(const std::string file_name);
	
		/**
		*@brief Open the file
		*@param file_name file name
		*/
		void Open(const boost::filesystem::path& file_name);
	
		/**
		*@brief Open the file
		*@param file_name file name
		*/
		void Open(const std::string file_name);

		/**
		*@brief Close the file
		*/
		void Close();
		
		/**
		*@brief Add a column
		*@param label label of the column
		*@param precision precision (i.e., number of significant digits)
		*/
		void AddColumn(const std::string label, const unsigned int precision);

		/**
		*@brief Complete the column initialization
		*/
		void Complete();

		/**
		*@brief New row of data
		*/
		void NewRow();

		/**
		*@brief Write a new numerical value
		*@param number the value
		*/
		template<typename T>
		void operator<<(const T& number);

	private:

		/**
		*@brief Set the default values
		*/
		void Default();

	private:

		std::ofstream fOut_;				// !< the output stream

		unsigned int n_;				// !< number of columns
		unsigned int iterator_;				// !< current iterator (from 0 to n_-1)

		unsigned int tab_;				// !< tab size
		unsigned int minimum_;				// !< minimum width

		std::vector<std::string> tag_;			// !< column tags
		std::vector<unsigned int> width_;		// !< column widths
		std::vector<unsigned int> precision_;		// !< column precisions	

		bool status_completed_;				// true if the declarations of columns has been completed
		boost::filesystem::path file_name_;		// file name
	};

}

#include "OutputFileColumns.hpp"

#endif	// OpenSMOKE_OutputFileColumns_H
