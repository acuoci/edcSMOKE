if (iLES == true)
{	 
	// Extract the subgrid fields from OpenFOAM
	kSgs.internalField() =  turbulence->k();
	eSgs.internalField() =  turbulence->epsilon() ;
	bound ( eSgs, eSgsMin );   

	// Reynolds number
	ReTauLES.internalField() = kSgs.internalField() * kSgs.internalField() / eSgs.internalField() / turbulence->mu() * thermo.rho() ; 
	bound ( ReTauLES, 0. );   

	// Fine structure properties
	gammaL = pow(  turbulence->mu() / thermo.rho() * eSgs / pow( kSgs, 2.0), 0.25);		
	tauStar = pow( turbulence->mu() / thermo.rho() / eSgs , 0.5);
	bound ( tauStar, tauMin );
	uTau_ = pow(thermo.rho() * eSgs / (turbulence->muEff()+musmall_), 0.5);
}
    

