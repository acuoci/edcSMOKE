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
|	License                                                               |
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

#ifndef OpenSMOKE_GridAdapter_H
#define OpenSMOKE_GridAdapter_H

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"
#include <Eigen/Dense>

namespace OpenSMOKE
{
	enum Adapter_Grid1D_Status { REGRID_FAILURE, REGRID_SUCCESS, NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED, NEW_POINTS_ARE_NEEDED, MAXIMUM_NUMBER_POINTS };

	//!  A class to manage refinement in 1D grids
	/*!
	This class provides the tools to manage adaptive refinement on 1D grids
	*/

	class Adapter_Grid1D
	{
	public:

		/**
		*@brief Initialize the internal parameters using default values
		*/
		Adapter_Grid1D();
		
		/**
		*@brief Set the maximum number of points
		*@param max_points maximum number of points (default: 300)
		*/
		void SetMaxPoints(const unsigned int max_points) { max_points_ = max_points; }
		
		/**
		*@brief Set the maximum number of points which can be added in a time
		*@param max_points_to_be_added maximum number of points which can be added in a time (default: 10)
		*/
		void SetMaxPointsToBeAdded(const unsigned int max_points_to_be_added) { max_points_to_be_added_ = max_points_to_be_added; }
		
		/**
		*@brief Set the coefficient governing the addition of points to resolve high gradient regions
		*@param coefficient_gradient coefficient (default: 0.1)
		*/
		void SetCoefficientGradient(const double coefficient_gradient) { coeff_grad_ = coefficient_gradient; }

		/**
		*@brief Set the coefficient governing the addition of points to resolve high curvature regions
		*@param coefficient_curvature coefficient (default: 0.5)
		*/
		void SetCoefficientCurvature(const double coefficient_curvature) { coeff_curv_ = coefficient_curvature; }
		
		/**
		*@brief Set the minimum absolute value for variable to be active in the refinement process
		*@param tolerance minimum value (default: 1e-8)
		*/
		void SetThreshold(const double threshold) { threshold_ = threshold; }

		/**
		*@brief Set the number of points to be used in the regrid procedure
		*@param regrid_points number of points to be used in the regrid procedure (default: 30)
		*/
		void SetNumberOfRegridPoints(const unsigned int regrid_points) { regrid_points_ = regrid_points; }

		/**
		*@brief Set the fraction of active points to be considered for regridding process
		*@param fraction_active_points fraction of active points to be considered for regridding process (default: 0.60)
		*/
		void SetFractionActivePoints(const double fraction_active_points) { fraction_active_points_ = fraction_active_points; }

		/**
		*@brief Set the ratio gradient over curvature to be considered for regridding process
		*@param ratio_gradient_curvature ratio gradient over curvature to be considered for regridding process (default: 1.50)
		*/
		void SetRatioGradientCurvature(const double ratio_gradient_curvature) { ratio_gradient_curvature_ = ratio_gradient_curvature; }

		/**
		*@brief Refine the grid
		*@param x the original coordinates (i.e. the grid to be refined)
		*@param phi the variables to look at for the refinement
		*@param x_new the new coordinates (i.e. the new refined grid)
		*/
		Adapter_Grid1D_Status Refine(const Eigen::VectorXd& x, const std::vector<Eigen::VectorXd>& phi, Eigen::VectorXd& x_new);

		/**
		*@brief Regrid operation
		*@param x the original coordinates (i.e. the grid to be regridded)
		*@param phi the variable to look at for calculating the weights
		*@param x_new the new coordinates (i.e. the new refined grid)
		*/
		Adapter_Grid1D_Status Regrid(const Eigen::VectorXd& x, const Eigen::VectorXd& phi, Eigen::VectorXd& x_new);

		/**
		*@brief Returns the number of points on which the regridding operations will be performed
		*/
		unsigned int regrid_points() const { return regrid_points_;  }

	private:

		double coeff_grad_;						//!< coefficient governing grid refinement (gradient)
		double coeff_curv_;						//!< coefficient governing grid refinement (curvarure)
		unsigned int max_points_;				//!< maximum number of points 
		unsigned int max_points_to_be_added_;	//!< maximum number of points which can be added in a time
		double threshold_;						//!< minimum value for a variable to be active
		double fraction_active_points_;			//!< fraction of active points to be considered for regridding process
		double ratio_gradient_curvature_;		//!< ratio gradient over curvature to be considered for regridding process
		unsigned int regrid_points_;			//!< number of points to be used in the regrid procedure
	};
}

#include "Adapter_Grid1D.hpp"

#endif // OpenSMOKE_GridAdapter_H
