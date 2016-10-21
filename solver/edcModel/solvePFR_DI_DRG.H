{	
	// Set initial conditions
	odeSolverPFR.SetInitialConditions(t0, yStar0);

	// Additional ODE solver options
	if (celli == 0)
	{
		// Set linear algebra options
		odeSolverPFR.SetLinearAlgebraSolver(odeParametersFineStructures.linear_algebra());
		odeSolverPFR.SetFullPivoting(odeParametersFineStructures.full_pivoting());
	}

	// Set relative and absolute tolerances
	odeSolverPFR.SetAbsoluteTolerances(odeParametersFineStructures.absolute_tolerance());
	odeSolverPFR.SetRelativeTolerances(odeParametersFineStructures.relative_tolerance());

	// Set minimum and maximum values
	odeSolverPFR.SetMinimumValues(yMin);
	odeSolverPFR.SetMaximumValues(yMax);
						
	// Solve
	OdeSMOKE::OdeStatus status = odeSolverPFR.Solve(tf);
	odeSolverPFR.Solution(yStarf);
}
