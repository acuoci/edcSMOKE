/*-----------------------------------------------------------------------*\
|                  _       _____ __  __  ____  _  ________                |
|                 | |     / ____|  \/  |/ __ \| |/ /  ____|               |
|          ___  __| | ___| (___ | \  / | |  | | ' /| |__                  |
|         / _ \/ _` |/ __|\___ \| |\/| | |  | |  < |  __|                 |
|        |  __/ (_| | (__ ____) | |  | | |__| | . \| |____                |
|         \___|\__,_|\___|_____/|_|  |_|\____/|_|\_\______|               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, M.R. Malik, A. Parente                             |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|   Contacts: Mohammad Rafi Malik                                         |
|   Aero-Thermo-Mechanical Department                                     |
|   Université Libre de Bruxelles                                         |
|   Avenue F. D. Roosevelt 50, 1050 Bruxelles (Belgium)                   |
|                                                                         |
|   Contacts: Alessandro Parente                                          |
|   email: alessandro.parente@ulb.ac.be                                   |
|   Aero-Thermo-Mechanical Department                                     |
|   Université Libre de Bruxelles                                         |
|   Avenue F. D. Roosevelt 50, 1050 Bruxelles (Belgium)                   |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of edcSMOKE solver.                                 |
|                                                                         |
|	License                                                           |
|                                                                         |
|   Copyright(C) 2016-2014 A. Cuoci, M.R. Malik, A. Parente               |
|   edcSMOKE is free software: you can redistribute it and/or modify      |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   edcSMOKE is distributed in the hope that it will be useful,           |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with edcSMOKE. If not, see <http://www.gnu.org/licenses/>.      |
|                                                                         |
\*-----------------------------------------------------------------------*/

// The ODE solver is always used at the first iteration
if (solverType == SOLVERTYPE_ODE || runTimeStep == 1)
{
	// Set initial conditions
	odeSolverPSR.SetInitialConditions(t0, yStar0);

	// Additional ODE solver options (every time, because the number of important species can change)

	// Set linear algebra options
	odeSolverPSR.SetLinearAlgebraSolver(odeParametersFineStructures.linear_algebra());
	odeSolverPSR.SetFullPivoting(odeParametersFineStructures.full_pivoting());

	// Set relative and absolute tolerances
	odeSolverPSR.SetAbsoluteTolerances(odeParametersFineStructures.absolute_tolerance());
	odeSolverPSR.SetRelativeTolerances(odeParametersFineStructures.relative_tolerance());

	// Set minimum and maximum values
	odeSolverPSR.SetMinimumValues(yMin);
	odeSolverPSR.SetMaximumValues(yMax);

	// Solve
	OdeSMOKE::OdeStatus status = odeSolverPSR.Solve(tf);
	odeSolverPSR.Solution(yStarf);
}
else
{
	// NLS Solvers
	typedef NlsSMOKE::KernelDense<OpenSMOKE::NLS_PSR_OpenSMOKE> denseNlsPSR;

	// Memory allocation
	unsigned int neq = yStar0.size()-2;
	Eigen::VectorXd y0(neq);
	Eigen::VectorXd yf(neq);
	Eigen::VectorXd residualsStar(neq);

	// First guess solution
	{
		unsigned int count = 0;
		for (unsigned int i=0;i<fine_structure_nls_psr.index_max_species();++i)	
			y0(count++) = yStar0(i);
		for (unsigned int i=fine_structure_nls_psr.index_max_species()+1;i<drg->number_important_species();++i)	
			y0(count++) = yStar0(i);
	}

	// Define non linear system solver
	NlsSMOKE::NonLinearSolver<denseNlsPSR> nlsSolverPSR;
	nlsSolverPSR.SetReactor(&fine_structure_nls_psr);

	// Set initial conditions
	nlsSolverPSR.SetFirstGuessSolution(y0);

	// Disable print message
	nlsSolverPSR.SetMaximumNumberOfIterations(5);

	// Set relative and absolute tolerances
	nlsSolverPSR.SetAbsoluteTolerances(nlsParametersFineStructures.absTolerance);
	nlsSolverPSR.SetRelativeTolerances(nlsParametersFineStructures.relTolerance);

	// Set minimum and maximum values
	nlsSolverPSR.SetMinimumValues(0.);
	nlsSolverPSR.SetMaximumValues(1.);

	// Disable warning messages about constraint violations
	nlsSolverPSR.SetPrint(false);

	// Solve
	NlsSMOKE::NlsStatus status = nlsSolverPSR();
	nlsSolverPSR.Solution(yf, residualsStar);

	// Norm
	const double norm2 = residualsStar.norm();

	// Analysis
	const double norm2_max = 0.01;
	if (status >= 0 && norm2 < norm2_max)
	{
		// Reconstruct data
		{
			unsigned int count = 0;
			for (unsigned int i=0;i<fine_structure_nls_psr.index_max_species();++i)	
				yStarf(i) = yf(count++);
			for (unsigned int i=fine_structure_nls_psr.index_max_species()+1;i<drg->number_important_species();++i)	
				yStarf(i) = yf(count++);

			double TStar, omegaMax;
			fine_structure_nls_psr.ReconstructData(yf, TStar, omegaMax);
			yStarf(fine_structure_nls_psr.index_max_species()) = omegaMax;
			yStarf(yStarf.size()-1) = TStar;
		}	

		nls_norm2_cumulative_ += norm2;
		nls_norm2_max_ = (norm2>nls_norm2_max_)? norm2 : nls_norm2_max_;
		nls_number_iterations_cumulative_ += nlsSolverPSR.number_iterations();
	}

	else
	{
		nls_failures_++;

		// Set initial conditions
		odeSolverPSR.SetInitialConditions(t0, yStar0);

		// Additional ODE solver options (every time, because the number of important species can change)

		// Set linear algebra options
		odeSolverPSR.SetLinearAlgebraSolver(odeParametersFineStructures.linear_algebra());
		odeSolverPSR.SetFullPivoting(odeParametersFineStructures.full_pivoting());

		// Set relative and absolute tolerances
		odeSolverPSR.SetAbsoluteTolerances(odeParametersFineStructures.absolute_tolerance());
		odeSolverPSR.SetRelativeTolerances(odeParametersFineStructures.relative_tolerance());

		// Set minimum and maximum values
		odeSolverPSR.SetMinimumValues(yMin);
		odeSolverPSR.SetMaximumValues(yMax);

		// Solve
		OdeSMOKE::OdeStatus status = odeSolverPSR.Solve(tf);
		odeSolverPSR.Solution(yStarf);
	}
}

