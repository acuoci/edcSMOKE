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

#ifndef OpenSMOKE_FixedProfile_H
#define OpenSMOKE_FixedProfile_H

namespace OpenSMOKE
{
	class FixedProfile
	{
	public:

		/**
		*@brief Default constructor
		*/
		FixedProfile();

		/**
		*@brief Default constructor
		*@param n number of grid points
		*@param x vector of abscissas
		*@param x vector of ordinates
		*/
		FixedProfile(const unsigned int n, const double* x, const double* y);

		/**
		*@brief Sets the profile
		*@param n number of grid points
		*@param x vector of abscissas
		*@param x vector of ordinates
		*/
		void operator() (const unsigned int n, const double* x, const double* y);

		/**
		*@brief Returns the abscissas
		*/
		const Eigen::VectorXd& x() const { return x_; }

		/**
		*@brief Retrurns the ordinates
		*/
		const Eigen::VectorXd& y() const { return y_; }

		/**
		*@brief Perform the interpolation based on the stored profile
		*@param x vector of points where to perform the interpolation
		*@param y interpolated ordinates
		*/
		void Interpolate(const Eigen::VectorXd& x, Eigen::VectorXd& y);

		/**
		*@brief Perform the interpolation based on the stored profile (without checks on the boundaries)
		*@param x vector of points where to perform the interpolation
		*@param y interpolated ordinates
		*/
		void InterpolateWithoutChecks(const Eigen::VectorXd& x, Eigen::VectorXd& y);

		/**
		*@brief Perform the interpolation based on the stored profile
		*@param x point where to perform the interpolation
		*@param y interpolated ordinate
		*/
		double Interpolate(const double x);

		/**
		*@brief Perform the interpolation based on the stored profile (without checks on the boundaries)
		*@param x point where to perform the interpolation
		*@param y interpolated ordinate
		*/
		double InterpolateWithoutChecks(const double x);

	private:

		Eigen::VectorXd x_;		//!< vector of abscissas
		Eigen::VectorXd y_;		//!< vector of ordinates
		double xstart_;			//!< first abscissa
		double xend_;			//!< last abscissa
		unsigned int n_;		//!< number of points
	};

}

#include "FixedProfile.hpp"

#endif	// OpenSMOKE_FixedProfile_H
