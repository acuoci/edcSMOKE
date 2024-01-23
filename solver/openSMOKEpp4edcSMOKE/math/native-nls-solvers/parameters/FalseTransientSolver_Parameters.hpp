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

#include "FalseTransientSolver_Parameters_Grammar.h"

namespace NlsSMOKE
{
	FalseTransientSolver_Parameters::FalseTransientSolver_Parameters() 
	{
		// Default values
		type_ = FALSETRANSIENT_SOLVER_BZZNLS;
		
		minimum_constraints_ = true;
		maximum_constraints_ = false;
		non_negative_unknowns_ = false;

		tolerance_function_ = 0.;		// default: uroundoff_1_over_3 = 6.055454E-06;
		tolerance_step_ = 0.;			// default: uroundoff_2_over_3 = 3.666853E-11;
		relative_error_ = 0.;			// default: uroundoff_1_over_2 = 1.490116E-08;
		tolerance_absolute_ = 0.;		// default: 1e-10
		tolerance_relative_ = 0.;		// default: sqrt(MachEpsFloat)

		maximum_number_iterations_ = 200;	// default: 200
		maximum_setup_calls_ = 10;			// default: 10
		maximum_sub_setup_calls_ = 5;		// default: 5
		maximum_newton_step_ = 0.;			// default: 1000 * ||y0||
		maximum_beta_fails_ = 10;			// default: 10

		omega_min_ = 0.;	// default: 1e-5
		omega_max_ = 0.;	// default: 0.90

		eta_choice_ = 0;	// 0=KIN_ETACHOICE1 1=KIN_ETACHOICE2 2=KIN_ETACHONSTANT
		verbosity_level_ = 0;

		strategy_ = NonLinearSolver_Parameters::NLS_STRATEGY_NEWTON_GLOBALIZATION;
		sparse_linear_algebra_ = false;
        
		jacobian_structure_ = OpenSMOKE::JACOBIAN_STRUCTURE_DENSE;
        	jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU;
		preconditioner_ = OpenSMOKE::PRECONDITIONER_SPARSE_ILUT;
 
		scaling_policy_ = 1;

		initial_step_ = 1.e-6;
		minimum_number_steps_ = 100;

		steps_reusing_jacobian_ = 20;
		steps_before_increasing_ = 20;
		increment_factor_ = 10.;
		maximum_step_ = 1.e-1;

		decrement_factor_ = 2.;
		minimum_step_ = 1e-9;
	}

	void FalseTransientSolver_Parameters::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		FalseTransientSolver_Parameters_Grammar grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@FalseTransientSolver") == true)
		{
			std::string name;
			dictionary.ReadString("@FalseTransientSolver", name);
			
			if (name == "OpenSMOKE++")
			{
				type_ = FALSETRANSIENT_SOLVER_OPENSMOKEPP;
			}
			else if (name == "BzzNls")
			{
				#if OPENSMOKE_USE_BZZMATH == 0
					OpenSMOKE::FatalErrorMessage("OpenSMOKE++ was built without the BzzMath support. Please select a different NLS solver (false-transient)");
				#endif
				type_ = FALSETRANSIENT_SOLVER_BZZNLS;
			}
			else if (name == "KinSol")
			{
				#if OPENSMOKE_USE_SUNDIALS == 0
					OpenSMOKE::FatalErrorMessage("OpenSMOKE++ was built without the Sundials (KinSol) support. Please select a different NLS solver (false-transient)");
				#endif
				type_ = FALSETRANSIENT_SOLVER_KINSOL;
			}

			else OpenSMOKE::FatalErrorMessage("Unknown @FalseTransient Solver: " + name);
		}

		// If the OpenSMOKE++ solver is selected, the user can choose sparse linear algebra
		if (type_ == FALSETRANSIENT_SOLVER_OPENSMOKEPP)
		{
			if (dictionary.CheckOption("@Jacobian") == true)
			{
				std::string name;
				dictionary.ReadString("@Jacobian", name);

				if (name == "Dense")			sparse_linear_algebra_ = false;
				else if (name == "Band")		sparse_linear_algebra_ = false;
				else if (name == "TridiagonalBlock")	sparse_linear_algebra_ = false;
				else if (name == "Sparse")		sparse_linear_algebra_ = true;
				else OpenSMOKE::FatalErrorMessage("Unknown @Jacobian. Available options: Dense | Band | TridiagonalBlock | Sparse");

				if (name == "Dense")			jacobian_structure_ = OpenSMOKE::JACOBIAN_STRUCTURE_DENSE;
				else if (name == "Band")		jacobian_structure_ = OpenSMOKE::JACOBIAN_STRUCTURE_BAND;
				else if (name == "TridiagonalBlock")	jacobian_structure_ = OpenSMOKE::JACOBIAN_STRUCTURE_TRIDIAGONAL_BLOCK;
				else if (name == "Sparse")		jacobian_structure_ = OpenSMOKE::JACOBIAN_STRUCTURE_SPARSE;
				else OpenSMOKE::FatalErrorMessage("Unknown @Jacobian. Available options: Dense | Band | TridiagonalBlock | Sparse");
			}

			if (dictionary.CheckOption("@SparseSolver") == true)
			{
				std::string name;
				dictionary.ReadString("@SparseSolver", name);

				#ifdef __APPLE__
				if (name == "EigenSparseLU")		jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU;
				else if (name == "EigenBiCGSTAB")   jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB;
				else if (name == "EigenGMRES")		jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES;
				else if (name == "EigenDGMRES")		jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES;
				else if (name == "Pardiso")			jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO;
				else if (name == "SuperLUSerial")   jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL;
				else if (name == "UMFPack")			jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK;
                #else
				if (name == "EigenSparseLU")		jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_SPARSE_LU;
				else if (name == "EigenBiCGSTAB")   jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_BICGSTAB;
				else if (name == "EigenGMRES")		jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_GMRES;
				else if (name == "EigenDGMRES")		jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_DGMRES;
				else if (name == "Pardiso")			jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_PARDISO;
				else if (name == "SuperLUSerial")   jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL;
				else if (name == "UMFPack")			jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_UMFPACK;
                #endif

				else OpenSMOKE::FatalErrorMessage("Unknown @SparseSolver. Available options: EigenSparseLU (default) | EigenBiCGSTAB | EigenGMRES | EigenDGMRES | Pardiso | SuperLUSerial | UMFPack");
			}

			if (dictionary.CheckOption("@Preconditioner") == true)
			{
				std::string name;
				dictionary.ReadString("@Preconditioner", name);

                #ifdef __APPLE__
                if (name == "ILUT")			   preconditioner_ = OpenSMOKE::PRECONDITIONER_SPARSE_ILUT;
				else if (name == "diagonal")   preconditioner_ = OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL;
                #else
				if (name == "ILUT")			   preconditioner_ = OpenSMOKE::SparsePreconditionerType::PRECONDITIONER_SPARSE_ILUT;
				else if (name == "diagonal")   preconditioner_ = OpenSMOKE::SparsePreconditionerType::PRECONDITIONER_SPARSE_DIAGONAL;
                #endif

				else OpenSMOKE::FatalErrorMessage("Unknown @Preconditioner. Available options: ILUT (default) | diagonal");
			}
		}


		if (dictionary.CheckOption("@FunctionTolerance") == true)
			dictionary.ReadDouble("@FunctionTolerance", tolerance_function_);

		if (dictionary.CheckOption("@StepTolerance") == true)
			dictionary.ReadDouble("@StepTolerance", tolerance_step_);

		if (dictionary.CheckOption("@AbsoluteTolerance") == true)
			dictionary.ReadDouble("@AbsoluteTolerance", tolerance_absolute_);

		if (dictionary.CheckOption("@RelativeTolerance") == true)
			dictionary.ReadDouble("@RelativeTolerance", tolerance_relative_);

		if (dictionary.CheckOption("@RelativeError") == true)
			dictionary.ReadDouble("@RelativeError", relative_error_);

		if (dictionary.CheckOption("@MaximumNumberOfIterations") == true)
			dictionary.ReadInt("@MaximumNumberOfIterations", maximum_number_iterations_);

		if (dictionary.CheckOption("@MaximumSetupCalls") == true)
			dictionary.ReadInt("@MaximumSetupCalls", maximum_setup_calls_);

		if (dictionary.CheckOption("@MaximumSubSetupCalls") == true)
			dictionary.ReadInt("@MaximumSubSetupCalls", maximum_sub_setup_calls_);

		if (dictionary.CheckOption("@Scaling") == true)
			dictionary.ReadInt("@Scaling", scaling_policy_);

		if (dictionary.CheckOption("@Strategy") == true)
		{
			std::string name;
			dictionary.ReadString("@Strategy", name);
			if (name == "NewtonBasic")					strategy_ = NonLinearSolver_Parameters::NLS_STRATEGY_NEWTON_BASIC;
			else if (name == "NewtonGlobalization")		strategy_ = NonLinearSolver_Parameters::NLS_STRATEGY_NEWTON_GLOBALIZATION;
			else OpenSMOKE::FatalErrorMessage("Unknown @Strategy: " + name);
		}

		if (dictionary.CheckOption("@VerbosityLevel") == true)
			dictionary.ReadInt("@VerbosityLevel", verbosity_level_);

		if (dictionary.CheckOption("@StepsReusingJacobian") == true)
			dictionary.ReadInt("@StepsReusingJacobian", steps_reusing_jacobian_);

		if (dictionary.CheckOption("@StepsBeforeIncreasing") == true)
			dictionary.ReadInt("@StepsBeforeIncreasing", steps_before_increasing_);

		if (dictionary.CheckOption("@MinimumNumberSteps") == true)
			dictionary.ReadInt("@MinimumNumberSteps", minimum_number_steps_);

		if (dictionary.CheckOption("@InitialStep") == true)
			dictionary.ReadDouble("@InitialStep", initial_step_);

		if (dictionary.CheckOption("@IncrementFactor") == true)
			dictionary.ReadDouble("@IncrementFactor", increment_factor_);

		if (dictionary.CheckOption("@MaximumStep") == true)
			dictionary.ReadDouble("@MaximumStep", maximum_step_);

		if (dictionary.CheckOption("@DecrementFactor") == true)
			dictionary.ReadDouble("@DecrementFactor", decrement_factor_);

		if (dictionary.CheckOption("@MinimumStep") == true)
			dictionary.ReadDouble("@MinimumStep", minimum_step_);

		if (dictionary.CheckOption("@MinimumConstraints") == true)
			dictionary.ReadBool("@MinimumConstraints", minimum_constraints_);

		if (dictionary.CheckOption("@MaximumConstraints") == true)
			dictionary.ReadBool("@MaximumConstraints", maximum_constraints_);

		if (dictionary.CheckOption("@NonNegativeVariables") == true)
			dictionary.ReadBool("@NonNegativeVariables", non_negative_unknowns_);
	}
}
