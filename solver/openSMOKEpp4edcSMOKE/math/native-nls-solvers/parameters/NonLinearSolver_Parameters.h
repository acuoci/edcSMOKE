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

#ifndef OpenSMOKE_NLS_Parameters_H
#define	OpenSMOKE_NLS_Parameters_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include "math/OpenSMOKEFunctions.h"
#include "math/OpenSMOKEStdInclude.h"
#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace NlsSMOKE
{
	class NonLinearSolver_Parameters
	{
	public:

		enum NLS_SOLVER { NLS_SOLVER_OPENSMOKEPP, NLS_SOLVER_BZZNLS, NLS_SOLVER_KINSOL };
		enum NlsStrategy { NLS_STRATEGY_NEWTON_BASIC, NLS_STRATEGY_NEWTON_GLOBALIZATION, NLS_STRATEGY_FIXED_POINT, NLS_STRATEGY_PICARD };

	public:

		/**
		*@brief Default constructor
		*/
		NonLinearSolver_Parameters();

		/**
		*@brief Reads the options from a file
		*@param dictionary dictionary from which the oprions are extracted
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Sets the type of integrator
		*@param type integrator type: NLS_SOLVER_OPENSMOKEPP, NLS_SOLVER_BZZNLS, NLS_SOLVER_KINSOL
		*/
		void SetType(const NLS_SOLVER type) { type_ = type; }
		
		/**
		*@brief Turns on/off the minimum constraints on the variables/unknowns
		*@param flag true means the minimum constraints are turned on (by default turned on)
		*/
		void SetMinimumConstraints(const bool flag) { minimum_constraints_ = flag; }

		/**
		*@brief Turns on/off the maximum constraints on the variables/unknowns
		*@param flag true means the maximum constraints are turned on (by default turned off)
		*/
		void SetMaximumConstraints(const bool flag) { maximum_constraints_ = flag; }

		/**
		*@brief True if all the unknowns must be non negative
		*@param flag true means all the unknowns must be non negative (by default turned off)
		*/
		void SetNonNegativeUnknowns(const bool flag) { non_negative_unknowns_ = flag; }
	
		/**
		*@brief Sets the verbosity level
		*@param verbosity_level the verbosity level (0 means no output)
		*/
		void SetVerbosityLevel(const int verbosity_level) { verbosity_level_ = verbosity_level; }
		
		/**
		*@brief Returns the solver type: NLS_SOLVER_OPENSMOKEPP, NLS_SOLVER_BZZNLS, NLS_SOLVER_KINSOL
		*/
		NLS_SOLVER type() const { return type_; }

		/**
		*@brief Returns true is the minimum constraints are turned on
		*/
		bool minimum_constraints() const { return minimum_constraints_; }

		/**
		*@brief Returns true is the maximum constraints are turned on
		*/
		bool maximum_constraints() const { return maximum_constraints_; }

		/**
		*@brief Returns true is all the unknowns are non negative
		*/
		bool non_negative_unknowns() const { return non_negative_unknowns_; }

		/**
		*@brief Returns the (TODO)
		*/
		double tolerance_function() const { return tolerance_function_; }

		/**
		*@brief Returns the (TODO)
		*/
		double tolerance_step() const { return tolerance_step_; }

		/**
		*@brief Returns the absolute tolerance
		*/
		double tolerance_absolute() const { return tolerance_absolute_; }

		/**
		*@brief Returns the relative tolerance
		*/
		double tolerance_relative() const { return tolerance_relative_; }

		/**
		*@brief Returns the absolute tolerances
		*/
		const Eigen::VectorXd& tolerances_absolute() const { return tolerances_absolute_; }

		/**
		*@brief Returns the relative tolerances
		*/
		const Eigen::VectorXd& tolerances_relative() const { return tolerances_relative_; }

		/**
		*@brief Returns the (TODO)
		*/
		double relative_error() const { return relative_error_; }

		/**
		*@brief Returns the maximum number of iterations
		*/
		int maximum_number_iterations() const { return maximum_number_iterations_; }

		/**
		*@brief Returns the (TODO)
		*/
		int maximum_setup_calls() const { return maximum_setup_calls_; }

		/**
		*@brief Returns the (TODO)
		*/
		int maximum_sub_setup_calls() const { return maximum_sub_setup_calls_; }

		/**
		*@brief Returns the maximum step to be used during the Newton's correction
		*/
		double maximum_newton_step() const { return maximum_newton_step_; }

		/**
		*@brief Returns the (TODO)
		*/
		int maximum_beta_fails() const { return maximum_beta_fails_; }

		/**
		*@brief Returns the (TODO)
		*/
		int eta_choice() const { return eta_choice_; }

		/**
		*@brief Returns the (TODO)
		*/
		double omega_min() const { return omega_min_; }

		/**
		*@brief Returns the (TODO)
		*/
		double omega_max() const { return omega_max_; }

		/**
		*@brief Returns the verbosity level
		*/
		int verbosity_level() const { return verbosity_level_; }

		/**
		*@brief Returns the strategy adopted for solving the non linear system: 
		*       NLS_STRATEGY_NEWTON_BASIC, NLS_STRATEGY_NEWTON_GLOBALIZATION, NLS_STRATEGY_FIXED_POINT, NLS_STRATEGY_PICARD
		*/
		NlsStrategy strategy() const { return strategy_; }

		/**
		*@brief Returns the solver to be used for solving the linear systems associated to the Jacobian matrix
		*/
		OpenSMOKE::SparseSolverType jacobian_solver() const { return jacobian_solver_; }

		/**
		*@brief Returns the Jacobian structure adopted
		*/
		OpenSMOKE::JacobianStructureType jacobian_structure() const { return jacobian_structure_; }
		
		/**
		*@brief Returns the preconditioner to be used for solving the linear systems associated to the Jacobian matrix
		*/
		OpenSMOKE::SparsePreconditionerType preconditioner() const { return preconditioner_; }
		
		/**
		*@brief Returns true is a sparse linear algebra solver is enabled
		*/
		bool sparse_linear_algebra() const { return sparse_linear_algebra_; }

		/**
		*@brief Returns the (TODO)
		*/
		int scaling_policy() const { return scaling_policy_; }

		/**
		*@brief Sets the relative tolerance
		*@param relative_tolerance the relative tolerance (same for all the variables)
		*/
		void SetRelativeTolerance(const double tolerance_relative) { tolerance_relative_ = tolerance_relative; }

		/**
		*@brief Sets the absolute tolerance
		*@param absolute_tolerance the absolute tolerance (same for all the variables)
		*/
		void SetAbsoluteTolerance(const double tolerance_absolute) { tolerance_absolute_ = tolerance_absolute; }

		/**
		*@brief Sets the relative tolerances
		*@param relative_tolerance the relative tolerances
		*/
		void SetRelativeTolerances(const Eigen::VectorXd& tolerances_relative) { tolerances_relative_ = tolerances_relative; }

		/**
		*@brief Sets the absolute tolerances
		*@param absolute_tolerance the absolute tolerances
		*/
		void SetAbsoluteTolerances(const Eigen::VectorXd& tolerances_absolute) { tolerances_absolute_ = tolerances_absolute; }


	private:

		NLS_SOLVER type_;						//!< NLS solver type

		bool minimum_constraints_;				//!< true if the minimum constraints are turned on (by default turned on)
		bool maximum_constraints_;				//!< true if the maximum constraints are turned on (by default turned off)
		bool non_negative_unknowns_;			//!< true if the all the variables must be non negative (by default turned off)

		double tolerance_function_;				//!< TODO
		double tolerance_step_;					//!< TODO
		double tolerance_absolute_;				//!< absolute tolerance
		double tolerance_relative_;				//!< relative tolerance
		Eigen::VectorXd tolerances_absolute_;	//!< absolute tolerance
		Eigen::VectorXd tolerances_relative_;	//!< relative tolerance
		double relative_error_;					//!< TODO

		int maximum_number_iterations_;			//!< maximum number of iterations
		int maximum_setup_calls_;				//!< TODO
		int maximum_sub_setup_calls_;			//!< TODO
		int maximum_beta_fails_;				//!< TODO
		double maximum_newton_step_;			//!< maximum step for Newton's correction

		int scaling_policy_;					//!< TODO

		double omega_min_;						//!< TODO
		double omega_max_;						//!< TODO
			
		int verbosity_level_;					//!< verbosity level
		int eta_choice_;						//!< TODO

		NlsStrategy strategy_;									//!< strategy to be used for solving the non linear system
		OpenSMOKE::JacobianStructureType jacobian_structure_;			//!< Jacobian structure type
		OpenSMOKE::SparseSolverType jacobian_solver_;			//!< solver to be used for solving linear systems associated to the Jacobian matrix
		OpenSMOKE::SparsePreconditionerType preconditioner_;	//!< preconditioner to be used for solving linear systems associated to the Jacobian matrix
		bool sparse_linear_algebra_;							//!< true if sparse linear solvers have to be used
	};
}

#include "NonLinearSolver_Parameters.hpp"

#endif	/* NonLinearSolver_Parameters_H */

