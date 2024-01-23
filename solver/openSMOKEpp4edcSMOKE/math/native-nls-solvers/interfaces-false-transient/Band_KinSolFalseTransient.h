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
|   This file is part of OpenSMOKE++ Suite.                               |
|                                                                         |
|   Copyright(C) 2015, 2014, 2013  Alberto Cuoci                          |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_Solve_Band_KinSolFalseTransient_H
#define OpenSMOKE_Solve_Band_KinSolFalseTransient_H

// Sundials
#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_band.h>  /* access to band SUNMatrix        */
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver  */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */
#include <sundials/sundials_math.h>    /* access to SUNRexp               */

// OpenSMOKE++
#include "math/native-nls-solvers/NonLinearSystemSolver"

namespace NlsSMOKE
{
	template<typename Object>
	int loop_KinSol(void* kmem, Object* object, const NlsSMOKE::FalseTransientSolver_Parameters& parameters, void *user_data, int& global_count, N_Vector& y, N_Vector& yInitial, N_Vector& scale);

	template<typename Object>
	int Solve_Band_KinSolFalseTransient(Object* object, const NlsSMOKE::FalseTransientSolver_Parameters& parameters)
	{
		std::cout << "Band False Transient solution (KinSol)..." << std::endl;

		const int neq = object->NumberOfEquations();

		N_Vector y;
		N_Vector yInitial;
		N_Vector scale;
		FalseTransient_UserData data;
		SUNMatrix J;
		SUNLinearSolver LS;
		SUNContext sunctx;

		int flag;
		void *kmem;

		y = NULL;
		yInitial = NULL;
		scale = NULL;
		kmem = NULL;
		sunctx = NULL;
		data = NULL;
		J = NULL;
		LS = NULL;

  		/* Create the SUNDIALS context that all SUNDIALS objects require */
  		int retval = SUNContext_Create(NULL, &sunctx);
  		if (check_retval(&retval, (char*)"SUNContext_Create", 1)) return(1);

		// Memory allocation
		{
			y = N_VNew_Serial(neq, sunctx);
			if (check_flag((void *)y, (char*)"N_VNew_Serial", 0)) return(1);

			yInitial = N_VNew_Serial(neq, sunctx);
			if (check_flag((void *)yInitial, (char*)"N_VNew_Serial", 0)) return(1);

			scale = N_VNew_Serial(neq, sunctx);
			if (check_flag((void *)scale, (char*)"N_VNew_Serial", 0)) return(1);
		}

		// User data
		{
			data = (FalseTransient_UserData)malloc(sizeof *data);
			data->deltat = parameters.initial_step();
			data->yInitial = NV_DATA_S(yInitial);
		}

		// Initialize and allocate memory for KINSOL
		kmem = KINCreate(sunctx);
		if (check_retval((void *)kmem, (char*)"KINCreate", 0)) return(1);

		// Set user data
		flag = KINSetUserData(kmem, data);
		if (check_flag(&flag, (char*)"KINSetUserData", 1)) return(1);

		// Assign the system
		flag = KINInit(kmem, kinsol_equations_false_transient, y);
		if (check_flag(&flag, (char*)"KINInit", 1)) return(1);

		// Set optional input
		{
			// Function tolerance
			if (parameters.tolerance_function() > 0.)
			{
				flag = KINSetFuncNormTol(kmem, parameters.tolerance_function());
				if (check_flag(&flag, (char*)"KINSetFuncNormTol", 1)) return(1);
			}

			// Step tolerance
			if (parameters.tolerance_step() > 0.)
			{
				flag = KINSetScaledStepTol(kmem, parameters.tolerance_step());
				if (check_flag(&flag, (char*)"KINSetScaledStepTol", 1)) return(1);
			}

			// Relative error
			if (parameters.relative_error() > 0.)
			{
				flag = KINSetRelErrFunc(kmem, parameters.relative_error());
				if (check_flag(&flag, (char*)"KINSetRelErrFunc", 1)) return(1);
			}

			// Force a Jacobian re-evaluation every mset iterations (default 10)
			flag = KINSetMaxSetupCalls(kmem, parameters.maximum_setup_calls());
			if (check_flag(&flag, (char*)"KINSetMaxSetupCalls", 1)) return(1);

			// Every msubset iterations, test if a Jacobian evaluation is necessary (default 5)
			flag = KINSetMaxSubSetupCalls(kmem, parameters.maximum_sub_setup_calls());
			if (check_flag(&flag, (char*)"KINSetMaxSubSetupCalls", 1)) return(1);

			// Maximum number of nonlinear iterations allowed
			flag = KINSetNumMaxIters(kmem, parameters.maximum_number_iterations());
			if (check_flag(&flag, (char*)"KINSetNumMaxIters", 1)) return(1);

			// Maximum beta failures in linesearch algorithm
			flag = KINSetMaxBetaFails(kmem, parameters.maximum_beta_fails());
			if (check_flag(&flag, (char*)"KINSetMaxBetaFails", 1)) return(1);

			// Maximum allowed scaled lenght of Newton step
			if (parameters.maximum_newton_step() > 0.)
			{
				flag = KINSetMaxNewtonStep(kmem, parameters.maximum_newton_step());
				if (check_flag(&flag, (char*)"KINSetMaxNewtonStep", 1)) return(1);
			}

			// Parameters omega_min and omega_max
			if (parameters.omega_min() > 0. || parameters.omega_max() > 0.)
			{
				flag = KINSetResMonParams(kmem, parameters.omega_min(), parameters.omega_max());
				if (check_flag(&flag, (char*)"KINSetResMonParams", 1)) return(1);
			}

			// Eta choice
			{
				if (parameters.eta_choice() == 0)	flag = KINSetEtaForm(kmem, KIN_ETACHOICE1);
				if (parameters.eta_choice() == 1)	flag = KINSetEtaForm(kmem, KIN_ETACHOICE2);
				if (parameters.eta_choice() == 2)	flag = KINSetEtaForm(kmem, KIN_ETACONSTANT);
				if (check_flag(&flag, (char*)"KINSetEtaForm", 1)) return(1);
			}

			// Print level
			flag = KINSetPrintLevel(kmem, parameters.verbosity_level());
			if (check_flag(&flag, (char*)"KINSetPrintLevel", 1)) return(1);
		}

		// Set minimum constraints
		if (parameters.non_negative_unknowns() == true)
		{
			N_Vector constraints;
			constraints = NULL;
			constraints = N_VNew_Serial(neq, sunctx);
			if (check_flag((void *)constraints, (char*)"N_VNew_Serial", 0)) return(1);

			// All the unknowns will be constrained to be >= 0
			N_VConst(1., constraints);
			flag = KINSetConstraints(kmem, constraints);
			if (check_flag(&flag, (char*)"KINSetConstraints", 1)) return(1);

			N_VDestroy_Serial(constraints);
		}

		// Set maximum constraints
		if (parameters.maximum_constraints() == true)
		{
			// Sundials KinSol does not support this function
		}

		// Attach band linear solver
		{
			const int mu = object->UpperBand();
			const int ml = object->LowerBand();

			/* In Sundials version 4.1.x no need to specify the fourth argument */
			// J = SUNBandMatrix(neq, mu, ml, mu + ml);
			J = SUNBandMatrix(neq, mu, ml, sunctx);
			if (check_flag((void *)J, (char*)"SUNBandMatrix", 0)) return(1);

			#if OPENSMOKE_USE_MKL == 1
			{
				LS = SUNLinSol_LapackBand(y, J, sunctx);
				if (check_flag((void *)LS, (char*)"SUNLapackBand", 0)) return(1);
			}
			#else
			{
				LS = SUNLinSol_Band(y, J, sunctx);
				if (check_flag((void *)LS, (char*)"SUNDenseLinearSolver", 0)) return(1);
			}
			#endif

			flag = KINSetLinearSolver(kmem, LS, J);
			if (check_flag(&flag, (char*)"KINSetLinearSolver", 1)) return(1);
		}

		// Initialize
		object->UnknownsVector(NV_DATA_S(y));
		object->UnknownsVector(NV_DATA_S(yInitial));

		// Scaling factors (ONE means no scaling)
		N_VConst_Serial(1., scale);
		if (parameters.scaling_policy() != 0)
		{
			realtype *pt_scale = NV_DATA_S(scale);
			realtype *pt_y = NV_DATA_S(y);
			for (int i = 0; i < neq; i++)
				pt_scale[i] = 1. / std::max(1.e-2, pt_y[i]);
		}

		// Call main solver
		int global_count = 0;
		int flag_solver;
		for (;;)
		{
			const int loop_result = loop_KinSol(kmem, object, parameters, data, global_count, y, yInitial, scale);

			// The loop was successfully accomplished
			if (loop_result == 0)
			{
				flag_solver = 0;
				break;
			}
			// The time step was too large
			else
			{
				data->deltat /= parameters.decrement_factor();
				data->deltat = std::max(data->deltat, parameters.minimum_step());

				std::cout << "Decrease time step..." << std::endl;
				std::cout << "dt=" << data->deltat << " alfa=" << parameters.decrement_factor() << " n=" << global_count << std::endl;

				flag_solver = -1;
			}
		}


		// Free memory
		N_VDestroy_Serial(y);
		N_VDestroy_Serial(yInitial);
		N_VDestroy_Serial(scale);
		KINFree(&kmem);
		free(data);
		SUNLinSolFree(LS);
		SUNMatDestroy(J);

		return flag_solver;
	}

	template<typename Object>
	int loop_KinSol(void* kmem, Object* object, const NlsSMOKE::FalseTransientSolver_Parameters& parameters, void *user_data, int& global_count, N_Vector& y, N_Vector& yInitial, N_Vector& scale)
	{
		int flag_solver;
		FalseTransient_UserData data = (FalseTransient_UserData)user_data;

		unsigned int count = 0;
		do
		{
			// Update the Jacobian at first call: true
			if (count == 0)
			{
				int flag = KINSetNoInitSetup(kmem, false);
				if (check_flag(&flag, (char*)"KINSetNoInitSetup", 1)) return(1);
			}

			// The time step is kept constant for a specific number of times
			if (count == parameters.steps_before_increasing())
			{
				data->deltat *= parameters.increment_factor();
				data->deltat = std::min(data->deltat, parameters.maximum_step());
				count = 0;

				std::cout << "dt=" << data->deltat << " beta=" << parameters.increment_factor() << " n=" << global_count << std::endl;
			}

			if (parameters.strategy() == NlsSMOKE::NonLinearSolver_Parameters::NLS_STRATEGY_NEWTON_BASIC)
			{
				flag_solver = KINSol(kmem, y, KIN_NONE, scale, scale);
				if (check_flag(&flag_solver, (char*)"KINSol::KIN_NONE", 0)) return(1);
			}
			else if (parameters.strategy() == NlsSMOKE::NonLinearSolver_Parameters::NLS_STRATEGY_NEWTON_GLOBALIZATION)
			{
				flag_solver = KINSol(kmem, y, KIN_LINESEARCH, scale, scale);
				if (check_flag(&flag_solver, (char*)"KINSol::KIN_LINESEARCH", 0)) return(1);
			}
			else if (parameters.strategy() == NlsSMOKE::NonLinearSolver_Parameters::NLS_STRATEGY_FIXED_POINT)
			{
				flag_solver = KINSol(kmem, y, KIN_FP, scale, scale);
				if (check_flag(&flag_solver, (char*)"KINSol::KIN_FP", 0)) return(1);
			}
			else if (parameters.strategy() == NlsSMOKE::NonLinearSolver_Parameters::NLS_STRATEGY_PICARD)
			{
				flag_solver = KINSol(kmem, y, KIN_PICARD, scale, scale);
				if (check_flag(&flag_solver, (char*)"KINSol::KIN_PICARD", 0)) return(1);
			}

			realtype* pt_y = NV_DATA_S(y);
			realtype* pt_yInitial = NV_DATA_S(yInitial);
			for (int i = 0; i < NV_LENGTH_S(y); i++)
				pt_yInitial[i] = pt_y[i];

			object->CorrectedUnknownsVector(NV_DATA_S(y));
			object->norm();

			// Update the Jacobian at first call: false
			{
				int flag = KINSetNoInitSetup(kmem, true);
				if (check_flag(&flag, (char*)"KINSetNoInitSetup", 1)) return(1);
			}

			count++;
			global_count++;

		} while (global_count <= parameters.minimum_number_steps());

		return 0;
	}
}

#endif // OpenSMOKE_Solve_Band_KinSolFalseTransient_H


