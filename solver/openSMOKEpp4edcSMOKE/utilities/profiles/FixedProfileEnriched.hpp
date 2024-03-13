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

#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <cctype>
#include <algorithm>
#include <boost/lexical_cast.hpp>

namespace OpenSMOKE
{
	FixedProfileEnriched::FixedProfileEnriched(const std::string file_name)
	{
		// Check file existence
		if (!boost::filesystem::exists(file_name))
			OpenSMOKE::FatalErrorMessage("The " + file_name + " file does not exist");

		// Check delimiter
		char delimiter = ';';
		{
			std::ifstream iFile(file_name.c_str(), std::ios::in);
			std::string line;
			std::getline(iFile, line);
			iFile.close();

			if (std::count(line.begin(), line.end(), ',') != 0)
			{
				delimiter = ',';
			}
			else
			{
				if (std::count(line.begin(), line.end(), ';') == 0)
					OpenSMOKE::FatalErrorMessage("The " + file_name + " does not appear as a CSV file. Accepted delimiters: , | ;");
			}
		}
		
		// Read file
		std::vector<std::vector<std::string> > parsedCSV;
		{
			std::ifstream iFile(file_name.c_str(), std::ios::in);

			std::string line;

			while (std::getline(iFile, line))
			{
				std::stringstream lineStream(line);
				std::string cell;
				std::vector<std::string> parsedRow;
				while (std::getline(lineStream, cell, delimiter))
				{
					parsedRow.push_back(cell);
				}

				parsedCSV.push_back(parsedRow);
			}

			iFile.close();
		}

		// Search for filter
		int line_filter = -1;
		for (unsigned int i = 0; i < parsedCSV.size(); i++)
			for (unsigned int j = 0; j < parsedCSV[i].size(); j++)
				if (parsedCSV[i][j] == "filter")
				{
					line_filter = i;
					break;
				}

		// Indices of lines
		unsigned int line_profile = 4;
		if (line_filter != -1)
			line_profile = 5;

		// Total number of lines
		if (parsedCSV.size() < 7)
			OpenSMOKE::FatalErrorMessage(file_name + " file: The minimum number of lines must be 7");

		// Check number of elements
		{
			for (unsigned int i = 0; i < line_profile; i++)
				if (parsedCSV[i].size() != 2)
					OpenSMOKE::FatalErrorMessage(file_name + " file: error in header lines (wrong number of columns)");

			if (parsedCSV[line_profile].size() != 1)
				OpenSMOKE::FatalErrorMessage(file_name + " file: error header lines (wrong number of columns)");

			for (unsigned int i = line_profile + 1; i < parsedCSV.size(); i++)
				if (parsedCSV[i].size() != 2)
					OpenSMOKE::FatalErrorMessage(file_name + " file: error in (x,y) profile lines (wrong number of columns)");

			if (parsedCSV[line_profile][0] != "profile")
				OpenSMOKE::FatalErrorMessage(file_name + " file: missing 'profile' keyword");
		}

		// Read x variable type and unit
		double x_conversion = 1.;
		{
			unsigned int line_profile_x = 2;

			std::string x_units;
			for (unsigned int j = 0; j < parsedCSV[line_profile_x].size(); j++)
			{
				if (j == 0)	x_type_ = parsedCSV[line_profile_x][j];
				if (j == 1)	x_units = parsedCSV[line_profile_x][j];
			}

			x_type_.erase(std::remove_if(x_type_.begin(), x_type_.end(), ::isspace), x_type_.end());
			x_units.erase(std::remove_if(x_units.begin(), x_units.end(), ::isspace), x_units.end());

			if (x_type_ == "time")
			{
				if (x_units == "s")			x_conversion = 1.;
				else if (x_units == "ms")	x_conversion = 1.e-3;
				else if (x_units == "min")	x_conversion = 60.;
				else if (x_units == "h")	x_conversion = 3600.;
				else OpenSMOKE::FatalErrorMessage(file_name + " file: Unknown time units. Allowed time units: s | ms | min | h");
			}
			else if (x_type_ == "length")
			{
				if (x_units == "m")			x_conversion = 1.;
				else if (x_units == "cm")	x_conversion = 1.e-2;
				else if (x_units == "mm")	x_conversion = 1.e-3;
				else OpenSMOKE::FatalErrorMessage(file_name + " file: Unknown length units. Allowed length units: m | cm | mm");
			}
			else OpenSMOKE::FatalErrorMessage(file_name + " file: Unknown x type. Allowed types: time | length");
		}

		// Read y variable type and unit
		double y_conversion = 1.;
		{
			unsigned int line_profile_y = 3;

			std::string y_units;
			for (unsigned int j = 0; j < parsedCSV[line_profile_y].size(); j++)
			{
				if (j == 0)	y_type_ = parsedCSV[line_profile_y][j];
				if (j == 1)	y_units = parsedCSV[line_profile_y][j];
			}

			y_type_.erase(std::remove_if(y_type_.begin(), y_type_.end(), ::isspace), y_type_.end());
			y_units.erase(std::remove_if(y_units.begin(), y_units.end(), ::isspace), y_units.end());

			if (y_type_ == "temperature")
			{
				if (y_units == "K")			y_conversion = 1.;
				else OpenSMOKE::FatalErrorMessage(file_name + " file: Unknown temperature units. Allowed temperature units: K");
			}
			else if (y_type_ == "pressure")
			{
				if (y_units == "Pa")		y_conversion = 1.;
				else if (y_units == "atm")	y_conversion = 101325.;
				else if (y_units == "bar")	y_conversion = 100000.;
				else OpenSMOKE::FatalErrorMessage(file_name + " file: Unknown pressure units. Allowed pressure units: Pa | atm | bar");
			}
			else if (y_type_ == "volume")
			{
				if (y_units == "m3")		y_conversion = 1.;
				else if (y_units == "l")	y_conversion = 1.e-3;
				else if (y_units == "dm3")	y_conversion = 1.e-3;
				else if (y_units == "cm3")	y_conversion = 1.e-6;
				else if (y_units == "mm3")	y_conversion = 1.e-9;
				else OpenSMOKE::FatalErrorMessage(file_name + " file: Unknown volume units. Allowed volume units: m3 | l | dm3 | cm3 | mm3");
			}
			else OpenSMOKE::FatalErrorMessage(file_name + " file: Unknown y type. Allowed types: temperature | pressure | volume");
		}

		// Read operating conditions: temperature and pressure
		{
			{
				unsigned int line_variable_1 = 0;
				
				std::string line = parsedCSV[line_variable_1][0];
					line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

				if (line != "temperature")
					OpenSMOKE::FatalErrorMessage(file_name + " file: Missing 'temperature' keyword in first line");

				std::vector<std::string> pieces;
				boost::split(pieces, parsedCSV[line_variable_1][1], boost::is_any_of(" "), boost::token_compress_on);       //Split data line
				if (pieces.size() != 2)
					OpenSMOKE::FatalErrorMessage(file_name + " file: Error in defining the temperature in the first line");

				{
					if (pieces[1] == "K")		temperature_ = boost::lexical_cast<double>(pieces[0]);
					else if (pieces[1] == "C")	temperature_ = boost::lexical_cast<double>(pieces[0]) + 273.15;
					else OpenSMOKE::FatalErrorMessage(file_name + " file: Unknown temperature units. Allowed temperature units: K | C");
				}
			}

			{
				unsigned int line_variable_2 = 1;

				std::string line = parsedCSV[line_variable_2][0];
				line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

				if (line != "pressure")
					OpenSMOKE::FatalErrorMessage(file_name + " file: Missing 'pressure' keyword in second line");

				std::vector<std::string> pieces;
				boost::split(pieces, parsedCSV[line_variable_2][1], boost::is_any_of(" "), boost::token_compress_on);       //Split data line
				if (pieces.size() != 2)
					OpenSMOKE::FatalErrorMessage(file_name + " file: Error in defining the pressure in the second line");

				{
					if (pieces[1] == "Pa")			pressure_ = boost::lexical_cast<double>(pieces[0]);
					else if (pieces[1] == "atm")	pressure_ = boost::lexical_cast<double>(pieces[0])*101325.;
					else if (pieces[1] == "bar")	pressure_ = boost::lexical_cast<double>(pieces[0])*100000.;
					else OpenSMOKE::FatalErrorMessage(file_name + " file: Unknown pressure units. Allowed pressure units: Pa | atm | bar");
				}
			}
		}

		// Recover profile
		std::vector<double> x;
		std::vector<double> y;

		// Recover profile
		for (unsigned int i = line_profile+1; i < parsedCSV.size(); i++)
			for (unsigned int j = 0; j < parsedCSV[i].size(); j++)
			{
				if (j == 0)	x.push_back(boost::lexical_cast<double>(parsedCSV[i][j])*x_conversion);
				if (j == 1)	y.push_back(boost::lexical_cast<double>(parsedCSV[i][j])*y_conversion);
			}
		

		// Read filter width
		if (line_filter != -1)
		{
			std::string line = parsedCSV[line_filter][0];
			line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

			if (line != "filter")
				OpenSMOKE::FatalErrorMessage(file_name + " file: Missing 'filter' keyword in first line");

			std::vector<std::string> pieces;
			boost::split(pieces, parsedCSV[line_filter][1], boost::is_any_of(" "), boost::token_compress_on);       //Split data line
			if (pieces.size() != 2)
				OpenSMOKE::FatalErrorMessage(file_name + " file: Error in defining the filter width in the first line");

			double conversion = 1.;

			if (x_type_ == "time")
			{
				if (pieces[1] == "s")			conversion = 1.;
				else if (pieces[1] == "ms")		conversion = 1.e-3;
				else if (pieces[1] == "min")	conversion = 60.;
				else if (pieces[1] == "h")		conversion = 3600.;
				else OpenSMOKE::FatalErrorMessage(file_name + " file: Unknown time units. Allowed time units: s | ms | min | h");
			}
			else if (x_type_ == "length")
			{
				if (pieces[1] == "m")		conversion = 1.;
				else if (pieces[1] == "cm")	conversion = 1.e-2;
				else if (pieces[1] == "mm")	conversion = 1.e-3;
				else OpenSMOKE::FatalErrorMessage(file_name + " file: Unknown length units. Allowed length units: m | cm | mm");
			}

			double filter_width = boost::lexical_cast<double>(pieces[0])*conversion;

			// Filtering
			{
				OpenSMOKE::OpenSMOKEVectorDouble xo(static_cast<int>(x.size()));
				OpenSMOKE::OpenSMOKEVectorDouble yo(static_cast<int>(x.size()));
				for (unsigned int i = 1; i <= x.size(); i++)
				{
					xo[i] = x[i - 1];
					yo[i] = y[i - 1];
				}

				OpenSMOKE::OpenSMOKEVectorDouble x_new(static_cast<int>(x.size()));
				OpenSMOKE::OpenSMOKEVectorDouble y_new(static_cast<int>(x.size()));
				OpenSMOKE::ApplyFilter(filter_width, 5, x_new, y_new, xo, yo);

				for (unsigned int i = 1; i <= x.size(); i++)
				{
					x[i - 1] = x_new[i];
					y[i-1] = y_new[i];
				}
			}
		}

		fixed_profile_ = new FixedProfile(static_cast<unsigned int>(x.size()), x.data(), y.data());
	}
}
