{
	T.max(280.);
    	T.min(3000.);
	TStar_[0].max(280.);
    	TStar_[0].min(3000.);

	volScalarField& he = thermo.he();

	fvScalarMatrix EEqn
	(
		fvm::div(phi, he)
		+ 	(
				he.name() == "e"
				? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
				: fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
			)
		- fvm::laplacian(turbulence->alphaEff(), he)
		==
		Sh_ +
		fvOptions(rho, he)
	);

	EEqn.relax();

	fvOptions.constrain(EEqn);

	EEqn.solve();

	fvOptions.correct(he);

	thermo.correct();

	Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;


	{
		scalarField& TCells = thermo.T().internalField();
    		if( 	(spark == true) && 
			(runTime.value() >= ignitionTime)&&(runTime.value() <= (ignitionTime+ignitionDuration)) ) 
   		 {
			Info << "Forcing the spark..." << endl;
        		for(int i = 0; i<ignCellsList.size(); i++)
        		{
           			 TCells[ignCellsList[i]] = ignitionTemperature;
           			 TStar_[0].internalField()[ignCellsList[i]] = ignitionTemperature;
        		}
    		}
	}
}
