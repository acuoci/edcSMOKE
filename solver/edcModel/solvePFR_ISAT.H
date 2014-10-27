// ISAT Algorithm
{	
	// Scaled variables
	for(unsigned int i=0;i<ne;i++)
		phiISAT_HOM(i) = yStar0(i)*scalingFactors_ISAT(i);

	chemComp *phi0base = NULL;
	
	// Check if retrieval is possible
	if(isat_HOM->retrieve(phiISAT_HOM, phi0base)) 
	{					
		double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
		
		// makes interpolation
		isat_HOM->interpol(phiISAT_HOM, RphiISAT_HOM, phi0base);

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
			odeSolverPFR.SetInitialConditions(t0, yStar0);

			// Additional ODE solver options
			if (nGrowHOM == 0 && nAddHOM == 0)
			{
				// Set linear algebra options
				odeSolverPFR.SetLinearAlgebraSolver(odeParametersFineStructures.linear_algebra());
				odeSolverPFR.SetFullPivoting(odeParametersFineStructures.full_pivoting());

				// Set relative and absolute tolerances
				odeSolverPFR.SetAbsoluteTolerances(odeParametersFineStructures.absolute_tolerance());
				odeSolverPFR.SetRelativeTolerances(odeParametersFineStructures.relative_tolerance());

				// Set minimum and maximum values
				odeSolverPFR.SetMinimumValues(yMin);
				odeSolverPFR.SetMaximumValues(yMax);
			}

			OdeSMOKE::OdeStatus status = odeSolverPFR.Solve(tf);
			odeSolverPFR.Solution(yStarf);

			// Move the solution from DI to ISAT
			for(unsigned int i=0;i<ne;i++)
				RphiISAT_HOM(i) = std::max(yStarf(i), 0.)*scalingFactors_ISAT(i);

			double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();

			cpuTimeDI += (t2-t1);
		}

		// Growth
		if(isat_HOM->grow(phiISAT_HOM, RphiISAT_HOM, phi0base)) 
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
						luSolver_ISAT, (tf-t0), numberSubSteps_ISAT, &odeSolverPFR);
	
			// add a new leaf 
			bool flag = isat_HOM->add(phiISAT_HOM, RphiISAT_HOM, mapGrad_HOM, qrSolver_ISAT, phi0base); 

			if(flag == false)	
				Info << "ISAT Error - Addition process failed..." << endl;
			
			nAddHOM++;
			
			double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			cpuTimeAddition += (t2-t1);
		}

		// check if is necessary to balance
		isat_HOM->cleanAndBalance();
	}
}

