{
	volScalarField& he = thermo.he();

	fvScalarMatrix EEqnNonReacting
	(
		fvm::ddt(rho, he)
		+ (
			mvConvection.valid()
			? mvConvection->fvmDiv(phi, he)
			: fvm::div(phi, he)
		)
		+ fvc::ddt(rho, K) + fvc::div(phi, K)
		+ (
			he.name() == "e"
			? mvConvection.valid()
			? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
			: fvc::div(fvc::absolute(phi, rho, U), p/rho)
			: -dpdt
		)
		+ thermophysicalTransport.divq(he)
		==
		rho*(U&g)
		+ reactions.Qdot()
		+ fvModels.source(rho, he)
	);

	EEqnNonReacting.relax();

	fvConstraints.constrain(EEqnNonReacting);

	EEqnNonReacting.solve();

	fvConstraints.constrain(he);

	thermo.correct();

	Info << "Min/Max T:" << min(thermo.T()).value() << ' ' << max(thermo.T()).value() << endl;
}
