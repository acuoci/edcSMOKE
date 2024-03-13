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

#ifndef OpenSMOKE_Grid1D_H
#define OpenSMOKE_Grid1D_H

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"
#include <Eigen/Dense>
#include "Adapter_Grid1D.h"
#include "Grammar_Grid1D.h"

namespace OpenSMOKE
{
	enum derivative_type { DERIVATIVE_1ST_UPWIND, DERIVATIVE_1ST_CENTERED, DERIVATIVE_1ST_BACKWARD, DERIVATIVE_1ST_FORWARD };

	//!  A class to manage 1D meshes with non-uniform step size
	/*!
	This class provides the tools to manage 1D meshes with non-uniform step size and automatic grid refinement
	*/

	class Grid1D
	{
	public:

		/**
		*@brief Creates a 1D mesh with uniform step size
		*@param Np number of points
		*@param x0 initial coordinate
		*@param xF final coordinate
		*/
		Grid1D(const int Np, const double x0, const double xF);

		/**
		*@brief Creates a 1D mesh from the disctionary
		*@param dictionary name of dictionary
		*/
		Grid1D(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, Eigen::VectorXd& w);

		/**
		*@brief Creates a stretched 1D mesh
		*@param Np number of points
		*@param x0 initial coordinate
		*@param xF final coordinate
		*@param alpha stretching factor
		*/
		Grid1D(const int Np, const double x0, const double xF, const double alpha);

		/**
		*@brief Creates a 1D mesh directly from a vector of coordinates
		*@param coordinates vector of coordinates
		*@param x_fixed coordinate of fixed point temperature
		*/
		Grid1D(const Eigen::VectorXd& coordinates);

		/**
		*@brief Creates a 1D mesh directly from a vector of coordinates
		*@param coordinates vector of coordinates
		*@param x_fixed coordinate of fixed point temperature
		*/
		Grid1D(const std::vector<double>& coordinates);

		/**
		*@brief Set a fixed point (if any)
		*@param x_point coordinate of fixed point (if any)
		*/
		void SetFixedPoint(const double x_point);

		/**
		*@brief Reset fixed point (if any)
		*/
		void ResetFixedPoint();

		/**
		*@brief Refine the grid
		*@param phi the variables to look at for the refinement
		*@param phi_new the updated variables, linearly interpolated on the new grid
		*/
		Adapter_Grid1D_Status Refine(const std::vector<Eigen::VectorXd>& phi, std::vector<Eigen::VectorXd>& phi_new);

		/**
		*@brief Double the number of steps in the grid
		*@param phi the variables to look at
		*@param phi_new the updated variables, linearly interpolated on the new grid
		*/
		void Double(const std::vector<Eigen::VectorXd>& phi, std::vector<Eigen::VectorXd>& phi_new);

		/**
		*@brief Double the number of steps only in the specified interval
		*@param xA left side
		*@param xB right side
		*@param phi the variables to look at
		*@param phi_new the updated variables, linearly interpolated on the new grid
		*/
		void Refine(const double xA, const double xB, const std::vector<Eigen::VectorXd>& phi, std::vector<Eigen::VectorXd>& phi_new);

		/**
		*@brief Regrid operation
		*@param index_psi the index of variable used for evaluationg the weights during the regrid operation
		*@param phi the variables to look at for the refinement
		*@param phi_new the updated variables, linearly interpolated on the new grid
		*/
		Adapter_Grid1D_Status Regrid(const unsigned int index_psi, const std::vector<Eigen::VectorXd>& phi, std::vector<Eigen::VectorXd>& phi_new);

		/**
		*@brief Update the grid on the basis of a new set of coordinates
		*@param coordinates the new set of coordinates
		*/
		void Update(const Eigen::VectorXd& coordinates);

		/**
		*@brief Calculate the derivative (first order) of variable phi (vector formulation)
		*@param derivative_type the type of derivative (upwind, centered, backward, forward)
		*@param u the velocity vector
		*@param phi the vector containing the variable to be differentiated
		*@param dphi the result of derivation
		*/
		void Derivative(const enum derivative_type, const Eigen::VectorXd &u, const Eigen::VectorXd &phi, Eigen::VectorXd *dphi);

		/**
		*@brief Calculate the derivative (first order) of vector of variables phi (matrix formulation)
		*@param derivative_type the type of derivative (upwind, centered, backward, forward)
		*@param u the velocity vector
		*@param phi the vector containing the vectors of variables to be differentiated
		*@param dphi the result of derivation
		*/
		void Derivative(const enum derivative_type type, const Eigen::VectorXd &u, const std::vector<Eigen::VectorXd> &phi, std::vector<Eigen::VectorXd>* dphi);

		/**
		*@brief Calculate the derivative (second order) of variable phi (without coefficient)
		*@param phi the vector containing the variable to be differentiated
		*@param d2phi the result of derivation
		*/
		void SecondDerivative(const Eigen::VectorXd &phi, Eigen::VectorXd* d2phi);

		/**
		*@brief Calculate the derivative (second order) of variable phi (without coefficient)
		*@param phi the vector containing the variable to be differentiated
		*@param d2phi the result of derivation
		*/
		void SecondDerivative(const std::vector<Eigen::VectorXd> &phi, std::vector<Eigen::VectorXd>* d2phi);

		/**
		*@brief Calculate the derivative (second order) of variable phi (without coefficient)
		*@param coeff the (diffusion) coefficient 
		*@param phi the vector containing the variable to be differentiated
		*@param d2phi the result of derivation
		*/
		void SecondDerivative(const Eigen::VectorXd &coeff, const Eigen::VectorXd &phi, Eigen::VectorXd* d2phi);
		
		/**
		*@brief Returns the number of points of the grid
		*/
		int Np() const { return Np_; }

		/**
		*@brief Returns the number of internal points of the grid
		*/
		int Ni() const { return Ni_; }

		/**
		*@brief Returns the coordinate of first point of the grid
		*/
		double x0() const { return x0_; }

		/**
		*@brief Returns the coordinate of the last point of the grid
		*/
		double xF() const { return xF_; }

		/**
		*@brief Returns the total length of the grid
		*/
		double L() const { return L_; }

		/**
		*@brief Returns the coordinates of the points of the grid
		*/
		const Eigen::VectorXd& x() const { return x_; }

		/**
		*@brief Returns the interval size (west side)
		*/
		const Eigen::VectorXd& dxw() const { return dxw_; }

		/**
		*@brief Returns the interval size (east side)
		*/
		const Eigen::VectorXd& dxe() const { return dxe_; }

		/**
		*@brief Returns the interval size (centred, east+west sides)
		*/
		const Eigen::VectorXd& dxc() const { return dxc_; }

		/**
		*@brief Returns the following quantity: [x(j+1/2) - x(j-1/2)]
		*/
		const Eigen::VectorXd& dxc_over_2() const { return dxc_over_2_; }
		
		/**
		*@brief Returns the index of fixed point coordinate
		*/
		unsigned int fixed_point() const { return fixed_point_; }
		
		/**
		*@brief Returns the coordinate of fixed point
		*/
		double x_fixed_point() const { return x_fixed_point_; }

		/**
		*@brief Returns the object governing the grid refinement
		*/
		Adapter_Grid1D& grid_adapter() { return *grid_adapter_; }
		
	private:

		Adapter_Grid1D*	grid_adapter_;	//!< the object governing the grid refinement

		unsigned int Np_;	//!< total number of points
		unsigned int Ni_;	//!< total number of intervals
		double x0_;			//!< coordinate of first point
		double xF_;			//!< coordinate of last point
		double L_;			//!< length

		unsigned int fixed_point_;	//!< index of fixed point (if any)
		double x_fixed_point_;		//!< coordinate of fixed point (if any)

		Eigen::VectorXd x_, dxe_, dxw_;
		Eigen::VectorXd x2_, ux_, udxe_, udxw_;
		Eigen::VectorXd dxc_, udxc_, dxc_over_2_, udxc_over_2_;
		Eigen::VectorXd cenp_, cenc_, cenm_;
		

	private:

		void Allocate();
		void Build();
		void BuildFromCoordinates(const Eigen::VectorXd& coordinates);
	};
}

#include "Grid1D.hpp"

#endif /* OpenSMOKE_Grid1D_H */

