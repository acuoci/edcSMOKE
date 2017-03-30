/*-----------------------------------------------------------------------*\
|                  _       _____ __  __  ____  _  ________                |
|                 | |     / ____|  \/  |/ __ \| |/ /  ____|               |
|          ___  __| | ___| (___ | \  / | |  | | ' /| |__                  |
|         / _ \/ _` |/ __|\___ \| |\/| | |  | |  < |  __|                 |
|        |  __/ (_| | (__ ____) | |  | | |__| | . \| |____                |
|         \___|\__,_|\___|_____/|_|  |_|\____/|_|\_\______|               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, M.R. Malik, Z. Li, A. Parente                      |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|   Contacts: Mohammad Rafi Malik, Zhiyi Li, Alessandro Parente           |
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
|   Copyright(C) 2017-2014 A. Cuoci, A. Parente                           |
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

// ISAT Algorithm
{	
	// Scaled variables
	for(unsigned int i=0;i<ne;i++)
		phiISAT_HOM(i) = yStar0(i)*scalingFactors_ISAT(i);

	chemComp *phi0base = NULL;
	
	// Check if retrieval is possible
	if(isatTable->retrieve(phiISAT_HOM, phi0base)) 
	{					
		double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
		
		// makes interpolation
		isatTable->interpol(phiISAT_HOM, RphiISAT_HOM, phi0base);

		//check negative value
		for(unsigned int i=0;i<ne;i++)
			yStarf(i) = std::max(RphiISAT_HOM(i), 0.)/scalingFactors_ISAT(i);

		nRetHOM++;	
		
		double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
		
		cpuTimeRet += (t2-t1);
	} 
	else
	{
		// Direct integration
		{
			double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();

			// Set initial conditions
			odeSolverPSR.SetInitialConditions(t0, yStar0);

			// Additional ODE solver options
			if (nGrowHOM == 0 && nAddHOM == 0)
			{
				// Set linear algebra options
				odeSolverPSR.SetLinearAlgebraSolver(odeParametersFineStructures.linear_algebra());
				odeSolverPSR.SetFullPivoting(odeParametersFineStructures.full_pivoting());

				// Set relative and absolute tolerances
				odeSolverPSR.SetAbsoluteTolerances(odeParametersFineStructures.absolute_tolerance());
				odeSolverPSR.SetRelativeTolerances(odeParametersFineStructures.relative_tolerance());

				// Set minimum and maximum values
				odeSolverPSR.SetMinimumValues(yMin);
				odeSolverPSR.SetMaximumValues(yMax);
			}

			OdeSMOKE::OdeStatus status = odeSolverPSR.Solve(tf);
			odeSolverPSR.Solution(yStarf);

			// Move the solution from DI to ISAT
			for(unsigned int i=0;i<ne;i++)
				RphiISAT_HOM(i) = std::max(yStarf(i), 0.)*scalingFactors_ISAT(i);

			double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();

			cpuTimeDI += (t2-t1);
		}

		// Growth
		if(isatTable->grow(phiISAT_HOM, RphiISAT_HOM, phi0base)) 
		{
			double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();

			phi0base->growEOA(phiISAT_HOM);
			nGrowHOM++;

			double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();

			cpuTimeGrowth += (t2-t1);
		} 
		// Addition
		else
		{
			double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			// compute mapping gradient
			calcMappingGradient(	phiISAT_HOM, RphiISAT_HOM, mapGrad_HOM, scalingFactors_ISAT, 
						luSolver_ISAT, (tf-t0), numberSubSteps_ISAT, &odeSolverPSR);
	
			// add a new leaf 
			bool flag = isatTable->add(phiISAT_HOM, RphiISAT_HOM, mapGrad_HOM, phi0base); 

			if(flag == false)	
				Info << "ISAT Error - Addition process failed..." << endl;
			
			nAddHOM++;
			
			double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			cpuTimeAddition += (t2-t1);
		}

		// check if it is necessary to balance
		isatTable->cleanAndBalance();
	}
}

