{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
      - radiation->Sh(thermo)
//    - fvm::laplacian(turbulence->muEff(), he)  // unit lewis no.
     ==
        //reaction->Sh()
	Sh_
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();

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
