if (edcCounter == numberIterationsToUpdateChemistry)
{
	edcCounter = 0;

	bool debug = false;
	const double chi = 1.;
	const double relaxationFactorFineStructures = 1.;

	Info << "Fine structure properties..." << endl;
	if (iLES == true)
	{
		//volSymmTensorField D(symm(fvc::grad(U)));
		//uTau_ = sqrt(2*(dev(D) && D));
	}
	else
	{
		if ( turbulence->db().objectRegistry::foundObject<volScalarField>("epsilon") )
		{
			const volScalarField& epsilonConst = turbulence->db().objectRegistry::lookupObject<volScalarField>("epsilon");
			eps_ = epsilonConst;
		}
		else if ( turbulence->db().objectRegistry::foundObject<volScalarField>("omega") )
		{
			const scalar Cmu = 0.09;
			const volScalarField& om = turbulence->db().objectRegistry::lookupObject<volScalarField>("omega");

			eps_ = Cmu*om*turbulence->k();
		}
		else
		{
			FatalError << "No epsilon or omega from turbulence model!" << exit(FatalError);
		}

		// Fine structures
		{
			uTau_ = eps_/(turbulence->k()+ksmall_);
			scalarField nu = thermo.mu().internalField()/rho;
			scalarField gammaL = Ccsi * pow(nu*eps_/pow(turbulence->k()+ksmall_,2.),0.25);	
			gammaStar = pow(gammaL, 2.);
			gammaStar = min(maxGammaSquared,gammaStar);
			mDotStar  = 1./Ctau*sqrt(eps_/nu);

			edc_coefficient = rho*gammaStar*mDotStar*chi/(1.-gammaStar*chi);
		}
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
			yStar0(i) = YStar_[i].internalField()[celli];
		yStar0(ne-1) = TStar_[0].internalField()[celli];

		// Initialize the fine structure
		thermodynamicsMapXML->SetTemperature(TCells[celli]);
		thermodynamicsMapXML->SetPressure(pCells[celli]);
		thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(hMean, xMean);
		hMean/=MWMean;	
		fine_structure_ode.Set(TCells[celli], pCells[celli], yMean, hMean, gammaStar[celli], mDotStar[celli], chi);

		// Final solution
		const double t0 = 0.;
		const double tf = numberOfResidenceTimes/max(mDotStar[celli], SMALL);

		// Solve the ODE system
		if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			// Set initial conditions
			odeSolver.SetInitialConditions(t0, yStar0);

			// Additional ODE solver options
			if (celli == 0)
			{
				// Set linear algebra options
				odeSolver.SetLinearAlgebraSolver(odeParametersFineStructures.linear_algebra());
				odeSolver.SetFullPivoting(odeParametersFineStructures.full_pivoting());

				// Set relative and absolute tolerances
				odeSolver.SetAbsoluteTolerances(odeParametersFineStructures.absolute_tolerance());
				odeSolver.SetRelativeTolerances(odeParametersFineStructures.relative_tolerance());

				// Set minimum and maximum values
				odeSolver.SetMinimumValues(yMin);
				odeSolver.SetMaximumValues(yMax);
			}
						
			// Solve
			OdeSMOKE::OdeStatus status = odeSolver.Solve(tf);
			odeSolver.Solution(yStarf);
		}
		#if OPENSMOKE_USE_DVODE == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
		{
			SolverDVODE->SetInitialValues(t0, yStar0.data());
			SolverDVODE->Solve(tf);
			SolverDVODE->Solution(yStarf.data());
		}
		#endif
		#if OPENSMOKE_USE_ODEPACK == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
		{
			SolverDLSODE->SetInitialValues(t0, yStar0.data());
			SolverDLSODE->Solve(tf);
			SolverDLSODE->Solution(yStarf.data());
		}
		#endif
		#if OPENSMOKE_USE_ODEPACK == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
		{
			SolverDLSODA->SetInitialValues(t0, yStar0.data());
			SolverDLSODA->Solve(tf);
			SolverDLSODA->Solution(yStarf.data());
		}
		#endif
		#if OPENSMOKE_USE_DASPK == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
		{
			SolverDASPK->SetInitialValues(t0, yStar0.data());
			SolverDASPK->Solve(tf);
			SolverDASPK->Solution(yStarf.data());
		}
		#endif
		#if OPENSMOKE_USE_RADAU == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
		{
			SolverRADAU5->SetInitialValues(t0, yStar0.data());
			SolverRADAU5->Solve(tf);
			SolverRADAU5->Solution(yStarf.data());
		}
		#endif
		#if OPENSMOKE_USE_MEBDF == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
		{
			SolverMEBDF->SetInitialValues(t0, yStar0.data());
			SolverMEBDF->Solve(tf);
			SolverMEBDF->Solution(yStarf.data());
		}
		#endif
		#if OPENSMOKE_USE_SUNDIALS == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
		{
			SolverCVODE->SetInitialValues(t0, yStar0.data());
			SolverCVODE->Solve(tf);
			SolverCVODE->Solution(yStarf.data());
		}
		#endif
				
		for(unsigned int i=0;i<ns;i++)
			YStar_[i].internalField()[celli] = (1.-relaxationFactorFineStructures)*yStar0(i) + relaxationFactorFineStructures*yStarf(i);
		TStar_[0].internalField()[celli] = (1.-relaxationFactorFineStructures)*yStar0(ne-1) + relaxationFactorFineStructures*yStarf(ne-1);
 
		for(unsigned int i=0;i<ns;i++)
			R_[i].internalField()[celli] = edc_coefficient[celli]*(YStar_[i].internalField()[celli]-Y[i].internalField()[celli]);

		thermodynamicsMapXML->SetTemperature(TCells[celli]);
		thermodynamicsMapXML->SetPressure(pCells[celli]);
		thermodynamicsMapXML->hMolar_Species(h_species);

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

edcCounter++;

