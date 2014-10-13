/*----------------------------------------------------------------------*\
|                  _       _____ __  __  ____  _  ________                |
|                 | |     / ____|  \/  |/ __ \| |/ /  ____|               |
|          ___  __| | ___| (___ | \  / | |  | | ' /| |__                  |
|         / _ \/ _` |/ __|\___ \| |\/| | |  | |  < |  __|                 |
|        |  __/ (_| | (__ ____) | |  | | |__| | . \| |____                |
|         \___|\__,_|\___|_____/|_|  |_|\____/|_|\_\______|               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, M.R. Malik, A. Parente                             |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|   Contacts: Mohammad Rafi Malik                                         |
|   Aero-Thermo-Mechanical Department                                     |
|   Université Libre de Bruxelles                                         |
|   Avenue F. D. Roosevelt 50, 1050 Bruxelles (Belgium)                   |
|                                                                         |
|   Contacts: Alessandro Parente                                          |
|   email: alessandro.parente@ulb.ac.be                                   |
|   Aero-Thermo-Mechanical Department                                     |
|   Université Libre de Bruxelles                                         |
|   Avenue F. D. Roosevelt 50, 1050 Bruxelles (Belgium)                   |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of edcSMOKE solver.                                 |
|                                                                         |
|	License                                                           |
|                                                                         |
|   Copyright(C) 2014 A. Cuoci, M.R. Malik, A. Parente                    |
|   edcSMOKE is free software: you can redistribute it and/or modify      |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   edcSMOKE is distributed in the hope that it will be useful,           |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with edcSMOKE. If not, see <http://www.gnu.org/licenses/>.      |
|                                                                         |
\*-----------------------------------------------------------------------*/

if (edcCounter == numberIterationsToUpdateChemistry)
{
	edcCounter = 0;

	bool debug = false;
	const double chi = 1.;
	const double relaxationFactorFineStructures = 1.;

	Info << "Fine structure properties..." << endl;
	#if STEADYSTATE == 0
	if (iLES == true)
	{
		const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("LESProperties");

		epsilonSGS = 1.048*pow(turbulence->k(), 1.5)/lesModel.delta();
		scalarField nu = thermo.mu().internalField()/rho;
		scalar Cnusgs = 0.05;
		scalar Cnusgs_025 = 0.472871;
		scalarField gammaL = (Ccsi*Cnusgs_025) * pow(thermo.mu().internalField()/turbulence->mut(), 0.25);	
		gammaStar = pow(gammaL, 2.);
		gammaStar = min(maxGammaSquared,gammaStar);
		mDotStar  = 1./Ctau*sqrt(epsilonSGS/nu);
		edc_coefficient = rho*gammaStar*mDotStar*chi/(1.-gammaStar*chi);

		Info<< "min/max(sgs_epsilon) = " << min(epsilonSGS).value() << ", " << max(epsilonSGS).value() << endl;
		Info<< "min/max(les_epsilon) = " << min(lesModel.epsilon()).value() << ", " << max(lesModel.epsilon()).value() << endl;
		Info<< "min/max(t_k) =         " << min(turbulence->k()).value() << ", " << max(turbulence->k()).value() << endl;
		Info<< "min/max(les_k) =       " << min(lesModel.k()).value() << ", " << max(lesModel.k()).value() << endl;
		Info<< "min/max(delta) =       " << min(lesModel.delta()).value() << ", " << max(lesModel.delta()).value() << endl;
		Info<< "min/max(mu) =          " << min(thermo.mu()).value() << ", " << max(thermo.mu()).value() << endl;
		Info<< "min/max(mut) =         " << min(turbulence->mut()).value() << ", " << max(turbulence->mut()).value() << endl;
		Info<< "min/max(gammaStar) =   " << min(gammaStar) << ", " << max(gammaStar) << endl;
		Info<< "min/max(mDotStar) =    " << min(mDotStar) << ", " << max(mDotStar) << endl;
	}
	else
	#endif
	{
		const volScalarField& epsilonConst = turbulence->epsilon();

		scalarField nu = thermo.mu().internalField()/rho;
		scalarField gammaL = Ccsi * pow(nu*turbulence->epsilon()/pow(turbulence->k()+ksmall_,2.),0.25);	
		gammaStar = pow(gammaL, 2.);
		gammaStar = min(maxGammaSquared,gammaStar);
		mDotStar  = 1./Ctau*sqrt(turbulence->epsilon()/nu);
		edc_coefficient = rho*gammaStar*mDotStar*chi/(1.-gammaStar*chi);

		Info<< "min/max(epsilon) = " << min(epsilonConst).value() << ", " << max(epsilonConst).value() << endl;
		Info<< "min/max(gammaStar) = " << min(gammaStar) << ", " << max(gammaStar) << endl;
		Info<< "min/max(mDotStar) = " << min(mDotStar) << ", " << max(mDotStar) << endl;
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

	unsigned int counter = 0;
	const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	forAll(TCells, celli)
	{	
		// Integration time
		const double t0 = 0.;
		double tf = 0.;

		// Extract the mean mass fractions
		for(unsigned int i=0;i<ns;i++)
			yMean[i+1] = Y[i].internalField()[celli];
			
		// Convert to mole fractions
		thermodynamicsMapXML->MoleFractions_From_MassFractions(xMean, MWMean, yMean);

		// Extract the fine strucutres mass fractions
		if (reactorType == REACTORTYPE_PSR)
		{
			for(unsigned int i=0;i<ns;i++)
				yStar0(i) = YStar_[i].internalField()[celli];
			yStar0(ne-1) = TStar_[0].internalField()[celli];

			// Initialize the fine structure
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(hMean, xMean);
			hMean/=MWMean;	
			
			fine_structure_ode_psr.Set(TCells[celli], pCells[celli], yMean, hMean, gammaStar[celli], mDotStar[celli], chi);

			tf = numberOfResidenceTimes/max(mDotStar[celli], SMALL);
		}
		else if (reactorType == REACTORTYPE_PFR)
		{
			for(unsigned int i=0;i<ns;i++)
				yStar0(i) = yMean[i+1];
			yStar0(ne-1) = TCells[celli];

			fine_structure_ode_pfr.Set(pCells[celli], gammaStar[celli], mDotStar[celli], chi);


			tf = (1.-gammaStar[celli]*chi)/max(mDotStar[celli], SMALL);
		}			

		// Solve the ODE system
		if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			if (reactorType == REACTORTYPE_PSR)
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
			else if (reactorType == REACTORTYPE_PFR)
			{
				// Set initial conditions
				odeSolverPFR.SetInitialConditions(t0, yStar0);

				// Additional ODE solver options
				if (celli == 0)
				{
					// Set linear algebra options
					odeSolverPFR.SetLinearAlgebraSolver(odeParametersFineStructures.linear_algebra());
					odeSolverPFR.SetFullPivoting(odeParametersFineStructures.full_pivoting());

					// Set relative and absolute tolerances
					odeSolverPFR.SetAbsoluteTolerances(odeParametersFineStructures.absolute_tolerance());
					odeSolverPFR.SetRelativeTolerances(odeParametersFineStructures.relative_tolerance());

					// Set minimum and maximum values
					odeSolverPFR.SetMinimumValues(yMin);
					odeSolverPFR.SetMaximumValues(yMax);
				}
						
				// Solve
				OdeSMOKE::OdeStatus status = odeSolverPFR.Solve(tf);
				odeSolverPFR.Solution(yStarf);
			}

		}
		#if OPENSMOKE_USE_DVODE == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
		{
			if (reactorType == REACTORTYPE_PSR)
			{
				SolverDVODEPSR->SetInitialValues(t0, yStar0.data());
				SolverDVODEPSR->Solve(tf);
				SolverDVODEPSR->Solution(yStarf.data());
			}
			else if (reactorType == REACTORTYPE_PFR)
			{
				SolverDVODEPFR->SetInitialValues(t0, yStar0.data());
				SolverDVODEPFR->Solve(tf);
				SolverDVODEPFR->Solution(yStarf.data());
			}
		}
		#endif
		#if OPENSMOKE_USE_ODEPACK == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
		{
			if (reactorType == REACTORTYPE_PSR)
			{
				SolverDLSODEPSR->SetInitialValues(t0, yStar0.data());
				SolverDLSODEPSR->Solve(tf);
				SolverDLSODEPSR->Solution(yStarf.data());
			}
			else if (reactorType == REACTORTYPE_PFR)
			{
				SolverDLSODEPFR->SetInitialValues(t0, yStar0.data());
				SolverDLSODEPFR->Solve(tf);
				SolverDLSODEPFR->Solution(yStarf.data());
			}
		}
		#endif
		#if OPENSMOKE_USE_ODEPACK == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
		{
			if (reactorType == REACTORTYPE_PSR)
			{
				SolverDLSODAPSR->SetInitialValues(t0, yStar0.data());
				SolverDLSODAPSR->Solve(tf);
				SolverDLSODAPSR->Solution(yStarf.data());
			}
			else if (reactorType == REACTORTYPE_PFR)
			{
				SolverDLSODAPFR->SetInitialValues(t0, yStar0.data());
				SolverDLSODAPFR->Solve(tf);
				SolverDLSODAPFR->Solution(yStarf.data());
			}
		}
		#endif
		#if OPENSMOKE_USE_DASPK == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
		{
			if (reactorType == REACTORTYPE_PSR)
			{
				SolverDASPKPSR->SetInitialValues(t0, yStar0.data());
				SolverDASPKPSR->Solve(tf);
				SolverDASPKPSR->Solution(yStarf.data());
			}
			else if (reactorType == REACTORTYPE_PFR)
			{
				SolverDASPKPFR->SetInitialValues(t0, yStar0.data());
				SolverDASPKPFR->Solve(tf);
				SolverDASPKPFR->Solution(yStarf.data());
			}
		}
		#endif
		#if OPENSMOKE_USE_RADAU == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
		{
			if (reactorType == REACTORTYPE_PSR)
			{
				SolverRADAU5PSR->SetInitialValues(t0, yStar0.data());
				SolverRADAU5PSR->Solve(tf);
				SolverRADAU5PSR->Solution(yStarf.data());
			}
			else if (reactorType == REACTORTYPE_PFR)
			{
				SolverRADAU5PFR->SetInitialValues(t0, yStar0.data());
				SolverRADAU5PFR->Solve(tf);
				SolverRADAU5PFR->Solution(yStarf.data());
			}
		}
		#endif
		#if OPENSMOKE_USE_MEBDF == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
		{
			if (reactorType == REACTORTYPE_PSR)
			{
				SolverMEBDFPSR->SetInitialValues(t0, yStar0.data());
				SolverMEBDFPSR->Solve(tf);
				SolverMEBDFPSR->Solution(yStarf.data());
			}
			else if (reactorType == REACTORTYPE_PFR)
			{
				SolverMEBDFPFR->SetInitialValues(t0, yStar0.data());
				SolverMEBDFPFR->Solve(tf);
				SolverMEBDFPFR->Solution(yStarf.data());
			}
		}
		#endif
		#if OPENSMOKE_USE_SUNDIALS == 1
		else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
		{
			if (reactorType == REACTORTYPE_PSR)
			{
				SolverCVODEPSR->SetInitialValues(t0, yStar0.data());
				SolverCVODEPSR->Solve(tf);
				SolverCVODEPSR->Solution(yStarf.data());
			}
			else if (reactorType == REACTORTYPE_PFR)
			{
				SolverCVODEPFR->SetInitialValues(t0, yStar0.data());
				SolverCVODEPFR->Solve(tf);
				SolverCVODEPFR->Solution(yStarf.data());
			}
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

		if (counter%(int(0.20*mesh.nCells())+1) == 0)
			Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << endl;

		counter++;
	}
	const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
	Info << tEnd - tStart << " s " << "(" << (tEnd - tStart)/double(TCells.size())*1e3 << " ms per cell)" << endl;
	
	WriteScalarField(mesh, mDotStar,  "edc::mDotStar");
	WriteScalarField(mesh, gammaStar, "edc::gammaStar");
}

edcCounter++;

