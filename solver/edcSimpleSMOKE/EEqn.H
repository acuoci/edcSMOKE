{
	T.max(270.);
    	T.min(3500.);

	if (combustionModel == "EDC")
	{
		TStar_[0].max(270.);
    		TStar_[0].min(3500.);
	}
	volScalarField& he = thermo.he();

	// Energy equation
	{
		fvScalarMatrix EEqn
		(
			fvm::div(phi, he)
			+ 	(
					he.name() == "e"
					? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
					: fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
			)
			- fvm::laplacian(turbulence->alphaEff(), he)
           	     - radiation->Sh(thermo)
			==
			Sh_ +
			fvOptions(rho, he)
		);

		EEqn.relax();

		fvOptions.constrain(EEqn);

		EEqn.solve();

		fvOptions.correct(he);

		thermo.correct();
		radiation->correct();
	}

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
				 if (combustionModel == "EDC")
           			 	TStar_[0].internalField()[ignCellsList[i]] = ignitionTemperature;
        		}
    		}
	}
}
