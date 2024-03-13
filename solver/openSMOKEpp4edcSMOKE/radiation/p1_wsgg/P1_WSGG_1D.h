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

#ifndef OpenSMOKE_P1_WSGG_1D_H
#define	OpenSMOKE_P1_WSGG_1D_H

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp> 
#include <boost/math/special_functions/pow.hpp> 
#include "utilities/grids/multi-purpose/Grid1DMultiPurpose.h"

namespace OpenSMOKE
{
	//!  A class to model radiative heat transfer in 1D geometries according to the P1-WSGG model
	/*!
	This class allows to calculate radiative heat transfer in 1D geometries (planar, cylindrical, and spherical)
	The user has to supply a proper 1D grid (planar, cylindrical, or spherical) and the P1_WSGG_1D class solves the
	G equation 
	The P1-WSGG model implemented here is taken from: Modest, Radiative Heat Transfer, 3rd Edition, Chapter 20
	*/

	class P1_WSGG_1D
	{

	public:

		/**
		*@brief Constructor
		*@param grid the 1D grid ((planar, cylindrical, and spherical)
		*@param model the model to calculate weights and absorption coefficients
		*/
		P1_WSGG_1D(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model);

		/**
		*@brief Returns the number of spatial points (according to the grid adopted during the construction phase)
		*/
		unsigned int np() const { return np_; }

		/**
		*@brief Returns the number of spectranl bands
		*/
		unsigned int nbands() const { return nbands_; }

		/**
		*@brief Returns the temperature (in K) of left wall (planar) or internal walls (cylindrical and spherical geometries)
		*/
		double T_wall1() const { return T_wall1_; }

		/**
		*@brief Returns the temperature (in K) of right wall (planar) or external walls (cylindrical and spherical geometries)
		*/
		double T_wall2() const { return T_wall2_; }

		/**
		*@brief Returns the emissivity of left wall (planar) or internal walls (cylindrical and spherical geometries)
		*/
		double epsilon_wall1() const { return epsilon_wall1_; }

		/**
		*@brief Returns the emissivity of right wall (planar) or external walls (cylindrical and spherical geometries)
		*/
		double epsilon_wall2() const { return epsilon_wall2_; }

		/**
		*@brief Returns the net heat flux (in W/m2) originating from the left wall (planar) or internal walls (cylindrical and spherical geometries)
		*/
		double q_from_wall1();

		/**
		*@brief Returns the net heat flux (in W/m2) originating from the right wall (planar) or external walls (cylindrical and spherical geometries)
		*/
		double q_from_wall2();
		
		/**
		*@brief Turn on/off the radiative equilibrium assumption (default is false)
		*@param flag true/false
		*/
		void SetRadiativeEquilibrium(const bool flag);

		/**
		*@brief Set the uniform temperature profile
		*@param T temperature (in K)
		*/
		void SetTemperature(const double T);

		/**
		*@brief Set the uniform partial pressure profile of CO2
		*@param PCO2_Pa pressure (in Pa)
		*/
		void SetPressureCO2(const double PCO2_Pa);

		/**
		*@brief Set the uniform partial pressure profile of H2O
		*@param PH2O_Pa pressure (in Pa)
		*/
		void SetPressureH2O(const double PH2O_Pa);

		/**
		*@brief Set the uniform soot volume fraction
		*@param fv soot volume fraction [-]
		*/
		void SetSootVolumeFraction(const double fv);		

		/**
		*@brief Set the temperature profile
		*@param T temperature (in K)
		*/
		void SetTemperature(const Eigen::VectorXd& T);

		/**
		*@brief Set the pressure partial profile of CO2
		*@param PCO2_Pa temperature (in Pa)
		*/
		void SetPressureCO2(const Eigen::VectorXd& PCO2_Pa);

		/**
		*@brief Set the pressure partial profile of H2O
		*@param PH2O_Pa temperature (in Pa)
		*/
		void SetPressureH2O(const Eigen::VectorXd& PH2O_Pa);

		/**
		*@brief Set the soot volume fraction profile
		*@param fv soot volume fraction [-]
		*/
		void SetSootVolumeFraction(const Eigen::VectorXd& fv);

		/**
		*@brief Set the temperatures of left and right walls (planar) or internal and external walls (cylindrical and spherical geometries)
		*@param T_wall1 left side or internal side temperature (in K)
		*@param T_wall2 right side or external side temperature (in K)
		*/
		void SetWallsTemperatures(const double T_wall1, const double T_wall2);

		/**
		*@brief Set the emissivities of left and right walls (planar) or internal and external walls (cylindrical and spherical geometries)
		*@param epsilon_wall1 left side or internal side emissivity (dimensionless)
		*@param epsilon_wall2 right side or external side emissivity (dimensionless)
		*/
		void SetWallsEmissivities(const double epsilon_wall1, const double epsilon_wall2);

		/**
		*@brief Returns the heat flux (in W/m2)
		*@param i index of mesh point (0-based)
		*/
		double q(const unsigned int i);

		/**
		*@brief Returns the incident radiation (in W/m2)
		*@param i index of mesh point (0-based)
		*/
		double G(const unsigned int i);

		/**
		*@brief Returns the divergence of heat flux (in W/m3)
		*@param i index of mesh point (0-based)
		*/
		double divq(const unsigned int i);

		/**
		*@brief Solves the G equation
		*/
		void SolveSpatialEquations();

		/**
		*@brief Print the solution of G equation on a file
		*@param file_name the path+name of file
		*/
		void PrintFile(const boost::filesystem::path file_name);

	private:

		bool is_scattering_;				//!< scattering on/off
		bool is_anisotropic_scattering_;	//!< anisotropic scattering on/off
		bool radiative_equilibrium_;		//!< radiative equilibrium  on/off

		Eigen::VectorXd T_;					//!< temperature profile (in K)
		Eigen::VectorXd PH2O_Pa_;			//!< pressure profile (in Pa)
		Eigen::VectorXd PCO2_Pa_;			//!< pressure profile (in Pa)
		Eigen::VectorXd fv_;				//!< soot volume fraction
		Eigen::VectorXd dummy_;				//!< dummy variable

		std::vector<Eigen::VectorXd>	kappa_;
		std::vector<Eigen::VectorXd>	gamma_;
		std::vector<Eigen::VectorXd>	a_;
		std::vector<Eigen::VectorXd>	G_;
		std::vector<Eigen::VectorXd>	dG_;
		std::vector<Eigen::VectorXd>	diffG_;



		double T_wall1_;					//!< temperature of left or internal wall (in K)
		double T_wall2_;					//!< temperature of right or external wall (in K)
		double epsilon_wall1_;				//!< emissivity of left or internal wall (dimensionless)
		double epsilon_wall2_;				//!< emissivity of right or external wall (dimensionless)

		Grid1DMultiPurpose& grid_;						//!< 1D grid

		unsigned int np_;								//!< number of space (grid) points
		unsigned int nbands_;							//!< number of space (grid) points
		unsigned int jMin_;								//!< spectral windows

		OpenSMOKE::ExtinctionCoefficientModel wsgg_model_;				//!< model for calculating the weights and the absorption coefficients

	private:

		/**
		*@brief Allocation of memory
		*/
		void MemoryAllocation();	

		/**
		*@brief Updates the physical coefficients (beta, gamma, etc.) before solving the non linear system
		*/
		void UpdateProperties();		

		/**
		*@brief Updates the spectral band properties
		*/
		void UpdateWSGGCoefficients();

		/**
		*@brief Returns the mean Planck absorption coefficient (in 1/m)
		*@param i index of mesh point (0-based)
		*/
		double kappa(const unsigned int i);

		/**
		*@brief Returns the sum of weights of bands (must be equal to 1, just for checking purposes)
		*@param i index of mesh point (0-based)
		*/
		double a(const unsigned int i);

	public:

		/**
		*@brief System of non linear equations corresponding to the G equations
		        The equations are solved in a fully coupled way (needed only for radiative equilibrium problem)
		        [TODO] Should be private, but it is not possible for compatibility issues with the NLS solver
		*@param y unknown (G)
		*@param t time (dummy variable)
		*@param f residuals
		*/
		void SpatialEquationsGlobal(const double* y, const double t, double* f);

		/**
		*@brief System of non linear equations corresponding to the G equations
		        The G equations are solved in a segregated approach
				[TODO] Should be private, but it is not possible for compatibility issues with the NLS solver
		*@param j the band to be solved
		*@param y unknown (G)
		*@param t time (dummy variable)
		*@param f residuals
		*/
		void SpatialEquationsPerBand(const unsigned int j, const double* y, const double t, double* f);

		/**
		*@brief Calls the non linear solver for the fully-coupled equations
		*/
		void SolveSpatialEquationsGlobal();

		/**
		*@brief Calls the non linear solver for the segregated approach
		*/
		void SolveSpatialEquationsPerBand();

	private:

		static const double SIGMA_;					//!< Stefan-Boltzmann constant (in W/m2/K4)
		static const double PI_;
		static const double PI_TIMES_4_;
		static const double PI_TIMES_2_;
		static const double SIGMA_OVER_PI_;
		static const double SIGMA_TIMES_4_;
	};
}

#include "P1_WSGG_1D.hpp"

#endif	/* OpenSMOKE_P1_WSGG_1D_H */
