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

namespace OpenSMOKE
{
	OutputFileColumns::OutputFileColumns()
	{
		Default();
	}

	OutputFileColumns::OutputFileColumns(const boost::filesystem::path& file_name)
	{
		Default();
		Open(file_name);
		file_name_ = file_name;
	}

	OutputFileColumns::OutputFileColumns(const std::string file_name)
	{
		Default();
		Open(file_name);
		file_name_ = file_name;
	}

	void OutputFileColumns::Default()
	{
		n_ = 0;
		iterator_ = 0;
		minimum_ = 16;
		tab_ = 4;
		status_completed_ = false;
		file_name_ = "";
	}

	void OutputFileColumns::Open(const boost::filesystem::path& file_name)
	{
		file_name_ = file_name;
		fOut_.open(file_name.c_str(), std::ios::out);
		fOut_.setf(std::ios::scientific);
		fOut_.setf(std::ios::left);
	}

	void OutputFileColumns::Open(const std::string file_name)
	{
		file_name_ = file_name;
		fOut_.open(file_name.c_str(), std::ios::out);
		fOut_.setf(std::ios::scientific);
		fOut_.setf(std::ios::left);
	}

	void OutputFileColumns::Close()
	{
		if (fOut_.is_open())
			fOut_.close();
		Default();
	}

	void OutputFileColumns::NewRow()
	{
		if (iterator_ != n_ && iterator_ != 0)
		{
			std::cout << "Filename: " << file_name_ << std::endl;
			OpenSMOKE::FatalErrorMessage("The number of written columns is not equal to the number of declared columns");
		}

		fOut_ << std::endl;
		iterator_ = 0;
	}

	void OutputFileColumns::AddColumn(const std::string label, const unsigned int precision)
	{
		std::stringstream number;

		iterator_++;
		number << iterator_;
		std::string tag = label + "(" + number.str() + ")";
		tag_.push_back(tag);
		width_.push_back( std::max( static_cast<unsigned int>(tag.size()+tab_), std::max(precision+9, minimum_) ) );
		precision_.push_back(precision);	
	}

	void OutputFileColumns::Complete()
	{
		if (status_completed_ == true)
		{
			std::cout << "Filename: " << file_name_ << std::endl;
			OpenSMOKE::FatalErrorMessage("In OutputFileColumns object, the columns have been already added.");
		}

		if (iterator_ == 0)
		{
			std::cout << "Filename: " << file_name_ << std::endl;
			OpenSMOKE::FatalErrorMessage("In OutputFileColumns object, no columns have been added.");
		}

		status_completed_ = true;
		n_ = iterator_;
		iterator_ = 0;
		
		for (unsigned int i=0;i<n_;i++)
			fOut_ << std::setprecision(precision_[i]) << std::setw(width_[i]) << tag_[i];
		fOut_ << std::endl;
	}

	template<typename T>
	void OutputFileColumns::operator<<(const T& number)
	{
		fOut_ << std::setprecision(precision_[iterator_]) << std::setw(width_[iterator_]) << number;
		iterator_++;
	}

}
