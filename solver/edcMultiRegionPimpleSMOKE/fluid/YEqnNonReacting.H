if (Y.size())
{
	mvConvection = tmp<fv::convectionScheme<scalar>>
	(
		fv::convectionScheme<scalar>::New
		(
			mesh,
			fields,
			phi,
			mesh.schemes().div("div(phi,Yi_h)")
		)
	);
}

reactions.correct();

forAll(Y, i)
{
	if (composition.solve(i))
	{
		volScalarField& Yi = Y[i];

		fvScalarMatrix YiEqnNonReacting
		(
			fvm::ddt(rho, Yi) 
			+ mvConvection->fvmDiv(phi, Yi)
			+ thermophysicalTransport.divj(Yi)
			==
			reactions.R(Yi) +
			fvModels.source(rho, Yi)
		);

		YiEqnNonReacting.relax();

		fvConstraints.constrain(YiEqnNonReacting);

		YiEqnNonReacting.solve("Yi");

		fvConstraints.constrain(Yi);
	}
}

composition.normalise();

