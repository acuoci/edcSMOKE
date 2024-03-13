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

#ifndef OpenSMOKE_SensitivityMap_BlockTridiagonal_SteadyState_H
#define OpenSMOKE_SensitivityMap_BlockTridiagonal_SteadyState_H

#include <Eigen/Dense>

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"
#include "maps/KineticsMap_CHEMKIN.h"
#include "math/OpenSMOKEBandMatrix.h"

#if OPENSMOKE_USE_MKL == 1
#include <Eigen/PardisoSupport>
#endif

#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
#include <Eigen/SuperLUSupport>
#endif

#if OPENSMOKE_USE_UMFPACK == 1
#include <Eigen/UmfPackSupport>
#endif

namespace OpenSMOKE
{
	//!  A class to perform sensitivity analysis
	/*!
	A class to perform sensitivity analysis
	*/

	class SensitivityMap_BlockTridiagonal_SteadyState
	{

	public:

		/**
		*@brief Default constructor
		*@param kineticsMap map containing the kinetic mechanism
		*@param number_of_equations total number of equations of the ODE (or NLS) associated to the sensitivity analysis
		*@param block_dimension number of equations per block
		*/
		SensitivityMap_BlockTridiagonal_SteadyState(KineticsMap_CHEMKIN& kineticMap, const unsigned int number_of_equations, const unsigned int block_dimension);

		/**
		*@brief Calculates the steady state sensitivity coefficients (dense version)
		*@param T current temperature [K]
		*@param P_Pa current pressure [Pa]
		*@param X current mole fractions [-]
		*@param J current Jacobian matrix [NExNE]
		*@param scaling_Jp current scaling vector[NE]
		*/
		void Calculate(const Eigen::VectorXd& T, const double P_Pa, const std::vector<Eigen::VectorXd>& X, OpenSMOKE::OpenSMOKEBandMatrixDouble& J, const std::vector<Eigen::VectorXd>& scaling_Jp);

		/**
		*@brief Saves sensitivity coefficients on xml file for post-processing analyses
		*@param folder_name name of folder where output file will be written
		*/
		void SaveOnXMLFile(const std::string folder_name);

		/**
		*@brief Sets the sensitivity type (default SENSITIVITY_KINETIC_CONSTANT)
		*/
		void SetSensitivityType(const PhysicalConstants::sensitivity_type type);

		/**
		*@brief Sets the type of energy equation (default CONSTANT_VOLUME_SYSTEM)
		*/
		void SetEnergyEquationType(const EnergyEquationType type);

		/**
		*@brief Returns the list of sensitivity parameters, according to the sensitivity_type
		*/
		const OpenSMOKE::OpenSMOKEVectorDouble& parameters() const { return parameters_; }

		/**
		*@brief Returns the total number of parameters
		*/
		unsigned int number_of_parameters() const { return number_of_parameters_; }

		/**
		*@brief Sets the index of temperature equation (0-index based)
		*/
		void SetIndexOfTemperature(const unsigned int index);

		/**
		*@brief Sets the index of mass flow rate equation (0-index based)
		*/
		void SetIndexOfMassFlowRate(const unsigned int index);

		/**
		*@brief Sets indices of species to be written on file (0-index based)
		*/
		void SetIndicesOfSpecies(const Eigen::VectorXi& indices_species);

		/**
		*@brief Returns the cpu time to perform a single factorization
		*/
		double cpuTimeSingleFactorization() const { return cpuTimeSingleFactorization_; }

		/**
		*@brief Returns the cpu time to perform a all the single factorizations (cumulative)
		*/
		double cpuTimeFactorization() const { return cpuTimeFactorization_; }

		/**
		*@brief Returns the cpu time to perform a single solution (NP linear systems)
		*/
		double cpuTimeSingleSolution() const { return cpuTimeSingleSolution_; }

		/**
		*@brief Returns the cpu time to perform the single solutions (NP linear systems) (cumulative)
		*/
		double cpuTimeSolution() const { return cpuTimeSolution_; }

		/**
		*@brief Returns the cpu time to perform a single assembling
		*/
		double cpuTimeSingleAssembling() const { return cpuTimeSingleAssembling_; }

		/**
		*@brief Returns the cpu time to perform the single assembling (cumulative)
		*/
		double cpuTimeAssembling() const { return cpuTimeAssembling_; }

	protected:

		KineticsMap_CHEMKIN& kinetics_;	//!< reference to the kinetic map

		int index_of_temperature_;				//!< index of temperature equation (if available) (0-index based)
		int index_of_mass_flow_rate_;			//!< index of mass flow rate (if available) (0-index based)
		int index_of_species_;					//!< index of first species equation (0-index based)

		unsigned int number_of_species_;		//!< total number of species
		unsigned int number_of_reactions_;		//!< total number of reactions
		unsigned int number_of_parameters_;		//!< total number of parameters
		unsigned int number_of_equations_;		//!< number of equations (for the jacobian construction)
		unsigned int block_dimension_;			//!< dimension of single block
		unsigned int number_of_blocks_;			//!< number of blocks

		Eigen::MatrixXd  sensitivity_temperature_;		//!< matrix containing the sensitivity coefficients
		Eigen::MatrixXd  sensitivity_mass_flow_rate_;	//!< matrix containing the sensitivity coefficients
		Eigen::MatrixXd* sensitivity_species_;			//!< matrix containing the sensitivity coefficients

		Eigen::VectorXi indices_species_;		//!< indices of species tobe written on file (0-based)

		OpenSMOKE::OpenSMOKEVectorDouble parameters_;			//!< vector containing the sensitivity parameters

		EnergyEquationType					energy_type_;		//!< type of energy equation (default CONSTANT_VOLUME_SYSTEM)
		PhysicalConstants::sensitivity_type sensitivity_type_;	//!< type of parameters

		double cpuTimeSingleFactorization_;			//!< cpu time to perform a single factorization
		double cpuTimeFactorization_;				//!< cpu time to perform a all the single factorizations (cumulative)
		double cpuTimeSingleSolution_;				//!< cpu time to perform a single solution (NP linear systems)
		double cpuTimeSolution_;					//!< cpu time to perform the single solutions (NP linear systems) (cumulative)
		double cpuTimeSingleAssembling_;			//!< cpu time to perform a single overall operations (factorization, solution assembling) (NP linear systems)
		double cpuTimeAssembling_;					//!< cpu time to perform the single overall operations (factorization, solution assembling) (cumulative)
	};
}

#include "SensitivityAnalysisMap_BlockTridiagonal_SteadyState.hpp"

#endif // OpenSMOKE_SensitivityMap_BlockTridiagonal_SteadyState_H







