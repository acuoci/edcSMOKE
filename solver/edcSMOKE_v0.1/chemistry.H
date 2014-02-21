if (homogeneousReactions == true)
{
	bool debug = false;
	const double relaxFineStructures = 1.;

	Info << "Fine structure properties..." << endl;
	{
		scalarField nu = thermo.mu().internalField()/rho;
		scalarField gammaL = 2.1377 * pow(nu*turbulence->epsilon()/pow(turbulence->k(),2.),0.25);	
		gammaStar = pow(gammaL, 2.);
		gammaStar = min(0.7,gammaStar);
		mDotStar  = 2.45*sqrt(turbulence->epsilon()/nu);

		const double chi = 1.;
		edc_coefficient = rho*gammaStar*mDotStar*chi/(1.-gammaStar*chi);
	}

	Info << "ODE Integration..." << endl;
	const int ns = Y.size();

	Eigen::VectorXd yStar0(ns);
	Eigen::VectorXd yStarf(ns);
	OpenSMOKE::OpenSMOKEVectorDouble yMean(ns);
	OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);
	Eigen::VectorXd yMin(ns);
	Eigen::VectorXd yMax(ns);
	
	for(unsigned int i=0;i<ns;i++)	yMin(i) = 0.;
	for(unsigned int i=0;i<ns;i++)	yMax(i) = 1.;

	const scalarField& pCells = thermo.p().internalField();
	const scalarField& TCells = thermo.T().internalField();
	const scalarField& hCells = thermo.he().internalField();
	scalarField& ShCells = Sh_.internalField();

	HomogeneousODE fine_structure_ode(*thermodynamicsMapXML, *kineticsMapXML);

	const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	forAll(TCells, celli)
	{	
		// Extract the mean mass fractions
		for(unsigned int i=0;i<ns;i++)
			yMean[i+1] = Y[i].internalField()[celli];

		// Extract the fine strucutres mass fractions
		for(unsigned int i=0;i<ns;i++)
			yStar0[i] = YStar_[i].internalField()[celli];

		// Initialize the fine structure
		fine_structure_ode.Set(TCells[celli], pCells[celli], yMean, gammaStar[celli], mDotStar[celli]);

		// Final solution
		const double t0 = 0.;
		const double tf = 1000./max(mDotStar[celli], SMALL);
		ODESystem_OpenSMOKE_Homogeneous ode(fine_structure_ode);
		OpenSMOKE::StiffOdeSolverObject_Dense o( yStar0, t0, &ode, "Eigen" );
		
		o.SetMinimumConstraints(yMin);
		o.SetMaximumConstraints(yMax);
		//o.SetTollAbs(speciesAbsoluteTolerance_);
		//o.SetTollRel(speciesRelativeTolerance_);

		yStarf = o(tf,tf);	
				
		for(unsigned int i=0;i<ns;i++)
			YStar_[i].internalField()[celli] = (1.-relaxFineStructures)*yStar0[i] + relaxFineStructures*yStarf[i];
 
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
}
