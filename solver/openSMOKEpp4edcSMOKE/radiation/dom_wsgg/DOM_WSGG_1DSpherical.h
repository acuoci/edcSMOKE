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

#ifndef OpenSMOKE_DOM_WSGG_1DSpherical_H
#define	OpenSMOKE_DOM_WSGG_1DSpherical_H

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp> 
#include <boost/math/special_functions/pow.hpp> 
#include "utilities/grids/multi-purpose/Grid1DMultiPurpose.h"

namespace OpenSMOKE
{
	//!  A class to model radiative heat transfer in 1D geometries according to the discrete ordinates method (DOM)
	/*!
	This class allows to calculate radiative heat transfer in 1D, spherically symmetric geometries
	The user has to supply a proper 1D spherically symmetric grid and the class solves the
	governing equations for the discrete ordinates
	The DOM implemented here is taken from: Modest, Radiative Heat Transfer, 3rd Edition, Chapter 17
	*/

	class DOM_WSGG_1DSpherical
	{

	public:

		/**
		*@brief Constructor
		*@param grid the 1D grid ((planar, cylindrical, and spherical)
		*@param model the model to calculate weights and absorption coefficients
		*/
		DOM_WSGG_1DSpherical(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model);

		/**
		*@brief Returns the number of discrete ordinates
		*/
		unsigned int n() const { return n_; }

		/**
		*@brief Returns half the number of discrete ordinates
		*/
		unsigned int N() const { return N_; }

		/**
		*@brief Returns the number of spatial points (according to the grid adopted during the construction phase)
		*/
		unsigned int np() const { return np_; }

		/**
		*@brief Returns the number of bands (the first one is always the clear window)
		*/
		unsigned int nbands() const { return nbands_; }
		
		/**
		*@brief Returns the direction cosines
		*/
		const double* mu() const { return mu_; }

		/**
		*@brief Returns the direction weights
		*/
		const double* wp() const { return wp_; }

		/**
		*@brief Returns the product of direction cosines and weights
		*/
		const double* mu_times_wp() const { return mu_times_wp_; }

		/**
		*@brief Returns the black-body intensity of radiation of internal wall (in W/m2)
		*/
		double Iblack_wall1() const { return Iblack_wall1_; }

		/**
		*@brief Returns the black-body intensity of radiation of external wall (in W/m2)
		*/
		double Iblack_wall2() const { return Iblack_wall2_; }

		/**
		*@brief Returns the temperature of internal wall (in K)
		*/
		double T_wall1() const { return T_wall1_; }

		/**
		*@brief Returns the temperature of external wall (in K)
		*/
		double T_wall2() const { return T_wall2_; }

		/**
		*@brief Returns the emissivity of internal wall
		*/
		double epsilon_wall1() const { return epsilon_wall1_; }

		/**
		*@brief Returns the emissivity of external wall
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
		*@brief Turn on/off the Legendre formulation for calculation of angular derivative (default is false)
		*@param flag true/false
		*/
		void SetLegendre(const bool flag);

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
		*@brief Set the soot volume fraction profile
		*@param fv soot volume fraction [-]
		*/
		void SetSootVolumeFraction(const Eigen::VectorXd& fv);

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
		*@brief Set the temperatures of internal and external walls (cylindrical and spherical geometries)
		*@param T_wall1 left side or internal side temperature (in K)
		*@param T_wall2 right side or external side temperature (in K)
		*/
		void SetWallsTemperatures(const double T_wall1, const double T_wall2);

		/**
		*@brief Set the emissivities of internal and external walls (cylindrical and spherical geometries)
		*@param epsilon_wall1 internal side emissivity (dimensionless)
		*@param epsilon_wall2 external side emissivity (dimensionless)
		*/
		void SetWallsEmissivities(const double epsilon_wall1, const double epsilon_wall2);

		/**
		*@brief Returns the heat flux (in W/m2)
		*@param j index of mesh point (0-based)
		*/
		double q(const unsigned int j);

		/**
		*@brief Returns the incident radiation (in W/m2)
		*@param j index of mesh point (0-based)
		*/
		double G(const unsigned int j);

		/**
		*@brief Returns the heat flux (in W/m2)
		*@param j index of mesh point (0-based)
		*@param k spectral band (0-based)
		*/
		double q(const unsigned int j, const unsigned int k);

		/**
		*@brief Returns the incident radiation (in W/m2)
		*@param j index of mesh point (0-based)
		*@param k spectral band (0-based)
		*/
		double G(const unsigned int j, const unsigned int k);

		/**
		*@brief Returns the divergence of heat flux (in W/m3)
		*@param j index of mesh point (0-based)
		*/
		double divq(const unsigned int j);

		/**
		*@brief Returns the divergence of heat flux (in W/m3)
		*@param kappa Planck mean absorption coefficient (in 1/m)
		*@param T     temperature (in K)
		*@param I     intensities (in W/m2)
		*/
		double divq(const double kappa, const double T, const Eigen::VectorXd &I);
		
		/**
		*@brief Solves the G equation
		*/
		void SolveSpatialEquations();

		/**
		*@brief Print the solution of G equation on a file
		*@param file_name the path+name of file
		*/
		void PrintFile(const boost::filesystem::path file_name);

		void UpdateWSGGCoefficients();
		void UpdateProperties();


	protected:

		bool is_scattering_;				//!< scattering on/off
		bool is_anisotropic_scattering_;	//!< anisotropic scattering on/off
		bool radiative_equilibrium_;		//!< radiative equilibrium  on/off
		bool legendre_;						//!< legendre formulation on/off

		unsigned int n_;					//!< number of discrete ordinates
		unsigned int N_;					//!< half number of discrete ordinates
		unsigned int np_;					//!< number of space (grid) points
		unsigned int nbands_;				//!< number of spectral bands (the first one is always the clear window)

		double* mu_;						//!< direction cosines
		double* wp_;						//!< direction weights
		double* alfa_;						//!< coefficients for angular derivative
		double* mu_times_wp_;				//!< product of direction cosines and weights

		Eigen::VectorXd T_;							//!< temperature profile (in K)
		Eigen::VectorXd PCO2_Pa_;					//!< CO2 partial pressure profile (in Pa)
		Eigen::VectorXd PH2O_Pa_;					//!< H2O partial pressure profile (in Pa)
		Eigen::VectorXd fv_;				//!< soot volume fraction

		std::vector<Eigen::VectorXd> kappa_;		//!< absorption coefficient (in 1/m)
		std::vector<Eigen::VectorXd> beta_;			//!< extinction coefficient (in 1/m)
		std::vector<Eigen::VectorXd> a_;			//!< spectral weights (-)

	//	Eigen::VectorXd sigma_;						//!< scattering coefficient (in 1/m)
	//	Eigen::VectorXd A1_;						//!< anisotropic scatter coefficient (linear contribution)
	//	Eigen::VectorXd omega_;						//!< albedo (dimensionless)

		double T_wall1_;								//!< temperature of left or internal wall (in K)
		double T_wall2_;								//!< temperature of right or external wall (in K)
		double epsilon_wall1_;							//!< emissivity of left or internal wall (dimensionless)
		double epsilon_wall2_;							//!< emissivity of right or external wall (dimensionless)
		double Iblack_wall1_;
		double Iblack_wall2_;

		std::vector< std::vector< Eigen::VectorXd > > I_;				//!< discretized intesnisities
		std::vector< std::vector< Eigen::VectorXd > > dI_;				//!< spatial derivatives of intensities
		std::vector< std::vector< Eigen::VectorXd > > fI_;				//!< residuals

		Grid1DMultiPurpose& grid_;						//!< 1D grid

		Eigen::VectorXd S_anistropic_scattering_;		//!< source term for anisotroipc scattering
		Eigen::FullPivLU<Eigen::MatrixXd>* ALU_;		//!< linear system matrix for Legendre formulation

		OpenSMOKE::ExtinctionCoefficientModel wsgg_model_;				//!< model for calculating the weights and the absorption coefficients

	protected:

		/**
		*@brief Calculates the black-body intensity of radiation (in W/m2)
		*@param T temperature (in K)
		*/
		double Iblack(const double T);

		/**
		*@brief Performs the final operations after construction
		*/
		void Finalize();

		/**
		*@brief Updates the first-order spatial derivatives of intensisities
		*/
		void FirstOrderDerivatives();

		/**
		*@brief Calculates the spectral radiosity of internal wall (in W/m2)
		*@param k the spectral band
		*/
		double J_wall1(const unsigned int k);

		/**
		*@brief Calculates the spectral radiosity of external wall (in W/m2)
		*@param k the spectral band
		*/
		double J_wall2(const unsigned int k);

		/**
		*@brief Allocation of memory
		*/
		void MemoryAllocation();

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
		*@brief System of non linear equations corresponding to the G equation
		        [TODO] Should be private, but it is not possible for compatibility issues with the NLS solver
		*@param y unknown (G)
		*@param t time (dummy variable)
		*@param f residuals
		*/
		void SpatialEquations(const double* y, const double t, double* f);

	protected:

		static const double SIGMA_;
		static const double PI_;
		static const double PI_TIMES_4_;
		static const double PI_TIMES_2_;
		static const double SIGMA_OVER_PI_;
		static const double SIGMA_TIMES_4_;
	};
}

#include "DOM_WSGG_1DSpherical.hpp"

#endif	/* OpenSMOKE_DOM_WSGG_1DSpherical_H */
