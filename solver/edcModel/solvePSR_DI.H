{
	// Set initial conditions
	odeSolverPSR.SetInitialConditions(t0, yStar0);

	// Additional ODE solver options
	if (celli == 0)
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
	
	// Solve
	OdeSMOKE::OdeStatus status = odeSolverPSR.Solve(tf);
	odeSolverPSR.Solution(yStarf);
}
