{
	bool debug = false;
	const double chi = 1.;
	const double relaxationFactorFineStructures = 1.;

	Info << "Fine structure properties..." << endl;
	{
		uTau_ = turbulence->epsilon()/(turbulence->k()+ksmall_);
		scalarField nu = thermo.mu().internalField()/rho;
		scalarField gammaL = Ccsi * pow(nu*turbulence->epsilon()/pow(turbulence->k()+ksmall_,2.),0.25);	
		gammaStar = pow(gammaL, 2.);
		gammaStar = min(maxGammaSquared,gammaStar);
		mDotStar  = 1./Ctau*sqrt(turbulence->epsilon()/nu);

		edc_coefficient = rho*gammaStar*mDotStar*chi/(1.-gammaStar*chi);
	}

	Info << "ODE Integration..." << endl;
	const int unsigned ns = Y.size();
	const int unsigned ne = ns+1;

	Eigen::VectorXd yStar0(ne);
	Eigen::VectorXd yStarf(ne);
	Eigen::VectorXd yMin(ne);
	Eigen::VectorXd yMax(ne);
	
	for(unsigned int i=0;i<ns;i++)	yMin(i) = 0.; yMin(ne-1) = 200.;
	for(unsigned int i=0;i<ns;i++)	yMax(i) = 1.; yMax(ne-1) = 5000.;

	OpenSMOKE::OpenSMOKEVectorDouble yMean(ns);
	OpenSMOKE::OpenSMOKEVectorDouble xMean(ns);
	OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);
	
	double hMean, MWMean;
	const scalarField& pCells = thermo.p().internalField();
	const scalarField& TCells = thermo.T().internalField();
	//const scalarField& hCells = thermo.he().internalField();
	scalarField& ShCells = Sh_.internalField();

	HomogeneousODE fine_structure_ode(*thermodynamicsMapXML, *kineticsMapXML);

	const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	forAll(TCells, celli)
	{	
		// Extract the mean mass fractions
		for(unsigned int i=0;i<ns;i++)
			yMean[i+1] = Y[i].internalField()[celli];
			
		// Convert to mole fractions
		thermodynamicsMapXML->MoleFractions_From_MassFractions(xMean, MWMean, yMean);

		// Extract the fine strucutres mass fractions
		for(unsigned int i=0;i<ns;i++)
			yStar0[i] = YStar_[i].internalField()[celli];
		yStar0[ne-1] = TStar_[0].internalField()[celli];

		// Initialize the fine structure
		thermodynamicsMapXML->SetTemperature(TCells[celli]);
		thermodynamicsMapXML->SetPressure(pCells[celli]);
		thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(hMean, xMean);
		hMean/=MWMean;	
		fine_structure_ode.Set(TCells[celli], pCells[celli], yMean, hMean, gammaStar[celli], mDotStar[celli], chi);

		// Final solution
		const double t0 = 0.;
		const double tf = numberOfResidenceTimes/max(mDotStar[celli], SMALL);
		ODESystem_OpenSMOKE_Homogeneous ode(fine_structure_ode);
		OpenSMOKE::StiffOdeSolverObject_Dense o( yStar0, t0, &ode, "Eigen" );
		
		o.SetMinimumConstraints(yMin);
		o.SetMaximumConstraints(yMax);
		o.SetTollAbs(absTolerance);
		o.SetTollRel(relTolerance);

		//yStarf = o(tf,tf);
		yStarf = o(tf);	
				
		for(unsigned int i=0;i<ns;i++)
			YStar_[i].internalField()[celli] = (1.-relaxationFactorFineStructures)*yStar0[i] + relaxationFactorFineStructures*yStarf[i];
		TStar_[0].internalField()[celli] = (1.-relaxationFactorFineStructures)*yStar0[ne-1] + relaxationFactorFineStructures*yStarf[ne-1];
 
		for(unsigned int i=0;i<ns;i++)
			R_[i].internalField()[celli] = edc_coefficient[celli]*(YStar_[i].internalField()[celli]-Y[i].internalField()[celli]);

		thermodynamicsMapXML->SetTemperature(TCells[celli]);
		thermodynamicsMapXML->SetPressure(pCells[celli]);
		thermodynamicsMapXML->hMolar_Species_From_MoleFractions(h_species);

		double HeatRelease = 0.;
		for(unsigned int i=0;i<ns;i++)
			HeatRelease += -h_species[i+1]/thermodynamicsMapXML->MW()[i+1] * R_[i].internalField()[celli];
		ShCells[celli] = HeatRelease;

		if (debug == true)
		{
			Info << celli << " " << TCells[celli] << " " << gammaStar[celli] << " " << mDotStar[celli] << " " << HeatRelease << endl;
			for(unsigned int i=0;i<ns;i++)
			{
				Info << Y[i].name() << " " << yStar0[i] << " " << yStarf[i] << " " << yMean[i+1] << " " << R_[i].internalField()[celli] << endl;
			}
			getchar();
		}
	}
	const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
	Info << tEnd - tStart << " s " << "(" << (tEnd - tStart)/double(TCells.size())*1e3 << " ms per cell)" << endl;
	
	WriteScalarField(mesh, mDotStar,  "edc::mDotStar");
	WriteScalarField(mesh, gammaStar, "edc::gammaStar");
}
