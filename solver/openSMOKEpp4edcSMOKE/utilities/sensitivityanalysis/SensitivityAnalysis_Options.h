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

#ifndef OpenSMOKE_SensitivityAnalysis_Options_H
#define	OpenSMOKE_SensitivityAnalysis_Options_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"
#include "math/PhysicalConstants.h"

namespace OpenSMOKE
{
	class SensitivityAnalysis_Options
	{
	public:

		/**
		*@brief Default constructor
		*/
		SensitivityAnalysis_Options();

		/**
		*@brief Initializes the object from an external dictionary
		*@param dictionary external dictionary
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Transfers all the options to the sensitivity analysis map
		*@param map the sensitivity analysis map receiving the options
		*/
		void TransferOptions(OpenSMOKE::SensitivityMap& map);

		/**
		*@brief Sets the number of substeps to perform the sensitivity analysis (default 2)
		*/
		void SetNumberOfSubSteps(const unsigned int number_of_substeps);
		
		/**
		*@brief Sets the sensitivity type (default SENSITIVITY_KINETIC_CONSTANT)
		*/
		void SetSensitivityType(const PhysicalConstants::sensitivity_type type);

		/**
		*@brief Sets the linear algebra package for solving the linear systems (default SOLVER_DENSE_EIGEN)
		*/
		void SetDenseSolverType(const DenseSolverType type);

		/**
		*@brief Sets the linear algebra package for solving the linear systems (default SOLVER_SPARSE_NONE)
		*/
		void SetSparseSolverType(const SparseSolverType type);

		/**
		*@brief Sets the LU factorization type (default DENSE_DECOMPOSITION_FULL_PIVOTING_LU)
		*/
		void SetDenseDecompositionType(const DenseDecompositionType type);

		/**
		*@brief Sets the sparse preconditioner type (default PRECONDITIONER_SPARSE_ILUT)
		*/
		void SetSparsePreconditionerType(const SparsePreconditionerType type);

		/**
		*@brief Returns the dense solver type
		*/
		DenseSolverType dense_solver_type() const { return dense_solver_type_; }

		/**
		*@brief Returns the sparse solver type
		*/
		SparseSolverType sparse_solver_type() const { return sparse_solver_type_; }

		/**
		*@brief Returns the sparse solver type
		*/
		SparsePreconditionerType sparse_preconditioner_type() const { return sparse_preconditioner_type_; }

		/**
		*@brief Returns the decomposition type
		*/
		DenseDecompositionType dense_decomposition_type() const { return dense_decomposition_type_; }

		/**
		*@brief Returns the sparse preconditioner fill factor
		*/
		int sparse_preconditioner_fill_factor() const { return sparse_preconditioner_fill_factor_;  };

		/**
		*@briefReturns the sparse preconditioner drop tolerance
		*/
		double sparse_preconditioner_drop_tolerance() const { return sparse_preconditioner_drop_tolerance_; };

		/**
		*@brief Returns the number of substeps
		*/
		int number_of_substeps() const { return number_of_substeps_; }

		/**
		*@brief Returns the sensitivity type
		*/
		PhysicalConstants::sensitivity_type sensitivity_type() const { return sensitivity_type_; }

		/**
		*@brief Returns the list of species
		*/
		const std::vector<std::string>& list_of_species() const { return list_of_species_; }

	private:

		PhysicalConstants::sensitivity_type sensitivity_type_;	//!< type of parameters

		int number_of_substeps_;								//!< number of substeps (default 2)
		std::vector<std::string>	list_of_species_;			//!< list of species for which the sensitivity coefficients will be written on output	

		DenseSolverType dense_solver_type_;						//!< linear algebra package for solving the linear systems (default LINEAR_ALGEBRA_PACKAGE_EIGEN)
		DenseDecompositionType dense_decomposition_type_;		//!< LU factorization type (default LINEAR_DECOMPOSITION_FULL_PIVOTING_LU)

		SparseSolverType sparse_solver_type_;
		SparsePreconditionerType sparse_preconditioner_type_;		
		int sparse_preconditioner_fill_factor_;
		double sparse_preconditioner_drop_tolerance_;
	};
}

#include "SensitivityAnalysis_Options.hpp"

#endif	/* OpenSMOKE_SensitivityAnalysis_Options_H */

