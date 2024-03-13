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

#ifndef OpenSMOKE_FixedProfileEnriched_H
#define OpenSMOKE_FixedProfileEnriched_H

#include "FixedProfile.h"

namespace OpenSMOKE
{
	class FixedProfileEnriched
	{
	public:

		/**
		*@brief Default constructor
		*@param file_name file name
		*/
		FixedProfileEnriched(const std::string file_name);

		/**
		*@brief Returns the nominal temperature (in K) corresponding to the given profile
		*/
		double temperature() const { return temperature_; }

		/**
		*@brief Returns the nominal pressure (in Pa) corresponding to the given profile
		*/
		double pressure() const { return pressure_; }

		/**
		*@brief Returns the type of x variable in the given profile (time | length)
		*/
		std::string x_type() const { return x_type_; }

		/**
		*@brief Returns the type of y variable in the given profile (temperature | pressure | volume)
		*/
		std::string y_type() const { return y_type_; }

		/**
		*@brief Returns the (x,y) profile
		*/
		const FixedProfile& profile() const { return *fixed_profile_; }

	private:

		FixedProfile*	fixed_profile_;		//!< profile
		double temperature_;				//!< nominal temperature (in K)
		double pressure_;					//!< nominal pressure (in Pa)

		std::string x_type_;				//!< x-type variable (time | length)
		std::string y_type_;				//!< y-type variable (temperature | pressure | volume)
	};

}

#include "FixedProfileEnriched.hpp"

#endif	// OpenSMOKE_FixedProfileEnriched_H
