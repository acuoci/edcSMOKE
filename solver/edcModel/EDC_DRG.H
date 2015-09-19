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

		tfStar = (1.-gammaStar*chi)/(mDotStar);
		Info<< "min/max(tfStar) = " << min(tfStar) << ", " << max(tfStar) << endl;
	}

	double hMean, MWMean;
	const scalarField& pCells = thermo.p().internalField();
	const scalarField& TCells = thermo.T().internalField();
	//const scalarField& hCells = thermo.he().internalField();
	scalarField& ShCells = Sh_.internalField();
	scalarField& drgSpeciesCells = drg_important_species.internalField();
	scalarField& drgReactionsCells = drg_important_reactions.internalField();
	scalarField& cpuReactionsCells = cpu_reaction_times.internalField();

	Info << "ODE Integration..." << endl;
	//const int unsigned ns = Y.size();
	//const int unsigned ne = ns+1+1;
	const unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
	      unsigned int NEQ = 0;
	

	Eigen::VectorXd yStar0;
	Eigen::VectorXd yStarf;
	Eigen::VectorXd yMin;
	Eigen::VectorXd yMax;
	Eigen::VectorXd y0;
	Eigen::VectorXd yf;
	OpenSMOKE::OpenSMOKEVectorDouble omegaMean(NC);
	OpenSMOKE::OpenSMOKEVectorDouble omegaStar(NC);
	OpenSMOKE::OpenSMOKEVectorDouble xMean(NC);
	OpenSMOKE::OpenSMOKEVectorDouble h_species(NC);
	OpenSMOKE::OpenSMOKEVectorDouble c_(NC);
	
	// ISAT Statistics
	#if EDCSMOKE_USE_ISAT == 1
		unsigned int nAddHOM   = 0;
		unsigned int nGrowHOM  = 0;
		unsigned int nRetHOM   = 0;
		double cpuTimeRet      = 0.;
		double cpuTimeDI       = 0.;
		double cpuTimeGrowth   = 0.;
		double cpuTimeAddition = 0.;
	#endif
	unsigned int counter = 0;
	unsigned int counter_skipped = 0;
	const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	forAll(TCells, celli)
	{	
		// Integration time
		const double t0 = 0.;
		double tf = 0.;
		//- Solving for celli:	
		double tLocalStart = OpenSMOKE::OpenSMOKEGetCpuTime();

		if (TCells[celli] > drg_minimum_temperature_for_chemistry)
		{
			// Extract the mean mass fractions
			for(unsigned int i=0;i<NC;i++)
				omegaMean[i+1] = Y[i].internalField()[celli];

			// Convert to mole fractions
			thermodynamicsMapXML->MoleFractions_From_MassFractions(xMean, MWMean, omegaMean);

			// Extract the fine strucutres mass fractions
			if (reactorType == REACTORTYPE_PSR)
			{
				for(unsigned int i=0;i<drg_temperature.size();i++)
					if(TCells[celli] < drg_temperature[i])
					{
						drg->SetEpsilon(drg_epsilon[i]);
						break;
					}

				double cTot_ = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
				Product(cTot_, xMean, &c_);
				drg->Analysis(TCells[celli], pCells[celli], c_);

				NEQ = drg->number_important_species()+2;

				yStar0.resize(NEQ);
				yStarf.resize(NEQ);
				yMax.resize(NEQ);	yMax.setConstant(1.); yMax(NEQ-2) = 5000.;	yMax(NEQ-1) = 1.e16; 
				yMin.resize(NEQ);	yMin.setConstant(0.); yMin(NEQ-2) = 280.;	yMin(NEQ-1) = 0.;

				for (unsigned int i=0;i<drg->number_important_species();++i)	
				{
					const unsigned int j = drg->indices_important_species()[i]+1;
					yStar0(i) = omegaMean[j];
					yStar0(NEQ-2) = TCells[celli];
					yStar0(NEQ-1) = 0.;

				}

				// Extract the star mass fractions
				for(unsigned int i=1;i<=NC;i++)
					omegaStar[i] = omegaMean[i];

				// Set DRG
				fine_structure_ode_psr.SetDRG(drg);

				// Initialize the fine structure
				thermodynamicsMapXML->SetTemperature(TCells[celli]);
				thermodynamicsMapXML->SetPressure(pCells[celli]);
				thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(hMean, xMean);
				hMean/=MWMean;	

				fine_structure_ode_psr.Set(TCells[celli], pCells[celli], omegaMean, omegaStar, hMean, gammaStar[celli], mDotStar[celli], chi);

				tf = numberOfResidenceTimes/max(mDotStar[celli], SMALL);
			}
			else if (reactorType == REACTORTYPE_PFR)	// TODO DRG
			{
				tf = (1.-gammaStar[celli]*chi)/max(mDotStar[celli], SMALL);

				for(unsigned int i=0;i<NC;i++)
				yStar0(i) = omegaMean[i+1];
				yStar0(NEQ-2) = TCells[celli];
				yStar0(NEQ-1) = tf;

				fine_structure_ode_pfr.Set(pCells[celli]);
			}			

			// Solve the ODE system
			if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
			{
				if (reactorType == REACTORTYPE_PSR)
				{
					// Set initial conditions
					odeSolverPSR.SetInitialConditions(t0, yStar0);

					// Additional ODE solver options
					//if (celli == 0)
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
					#if EDCSMOKE_USE_ISAT == 1
						if(isatCheck == true)
						{
							#include "solvePFR_ISAT.H"
						}
						else
						{
							#include "solvePFR_DI.H"
						}
					#else
						#include "solvePFR_DI.H"
					#endif
				}

			}
			
			for(unsigned int i=0;i<NEQ-2;i++)
			{
               			const unsigned int j = drg->indices_important_species()[i];
                		YStar_[j].internalField()[celli] = (1.-relaxationFactorFineStructures)*yStar0(i) + relaxationFactorFineStructures*yStarf(i);
			}
				TStar_[0].internalField()[celli] = (1.-relaxationFactorFineStructures)*yStar0(NEQ-2) + relaxationFactorFineStructures*yStarf(NEQ-2);
	 
			for(unsigned int i=0;i<NEQ-2;i++)
				R_[i].internalField()[celli] = edc_coefficient[celli]*(YStar_[i].internalField()[celli]-Y[i].internalField()[celli]);

			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->hMolar_Species(h_species);

			double HeatRelease = 0.;
			for(unsigned int i=0;i<NC;i++)
				HeatRelease += -h_species[i+1]/thermodynamicsMapXML->MW()[i+1] * R_[i].internalField()[celli];
			ShCells[celli] = HeatRelease;

			if (debug == true)
			{
				Info << celli << " " << TCells[celli] << " " << gammaStar[celli] << " " << mDotStar[celli] << " " << HeatRelease << endl;
				for(unsigned int i=0;i<NC;i++)
				{
					Info << Y[i].name() << " " << yStar0[i] << " " << yStarf[i] << " " << omegaMean[i+1] << " " << R_[i].internalField()[celli] << endl;
				}
				getchar();
			}

			// Fill the relevant DRG fields
			drgSpeciesCells[celli]   = drg->number_important_species();
			drgReactionsCells[celli] = drg->number_important_reactions();
		}
		else
		{
			drgSpeciesCells[celli]   = 0;
			drgReactionsCells[celli] = 0;
			counter_skipped++;
		}	
		double tLocalEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuReactionsCells[celli] = (tLocalEnd - tLocalStart)*1.e6;

			if (counter%(int(0.20*mesh.nCells())+1) == 0)
				Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << endl;

			counter++;
	}
	const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
	Info << "CPU time for chemistry (EDC) : " << tEnd - tStart << " s " << "(" << (tEnd - tStart)/double(TCells.size())*1e3 << " ms per cell)" << endl;

	#if EDCSMOKE_USE_ISAT == 1
	if(isatCheck == true) 
	{
		Info << endl;
		Info << " ********* ISAT HOM stats **********" << endl;
		
		Info << "   Direct Integration : " << isatTable->nAdd()+isatTable->nGrow()  << " (" << nAddHOM+nGrowHOM << ")" << " (" << (nAddHOM+nGrowHOM)/double(mesh.nCells())*100. << "%)" << endl;
		Info << "      Add             : " << isatTable->nAdd()  << " (" << nAddHOM  << ")" << " (" << nAddHOM/double(mesh.nCells())*100. << "%)" << endl;
		Info << "      Grow            : " << isatTable->nGrow() << " (" << nGrowHOM << ")" << " (" << nGrowHOM/double(mesh.nCells())*100. << "%)" << endl;
		Info << "   Retrieve           : " << isatTable->nUse()  << " (" << nRetHOM  << ")" << " (" << nRetHOM/double(mesh.nCells())*100. << "%)" << endl;
		Info << endl;				

		const double cpuTimeIntegration = cpuTimeDI + cpuTimeGrowth + cpuTimeAddition;
		Info << "   CPU Integration  : " << cpuTimeIntegration  << " (" << cpuTimeIntegration/(tEnd-tStart)*100. << "%)" << endl;
		Info << "     CPU DI         : " << cpuTimeDI           << " (" << cpuTimeDI/(tEnd-tStart)*100.          << "%)" << endl;				
		Info << "     CPU Growth     : " << cpuTimeGrowth       << " (" << cpuTimeGrowth/(tEnd-tStart)*100.         << "%)" << endl;
		Info << "     CPU Addition   : " << cpuTimeAddition     << " (" << cpuTimeAddition/(tEnd-tStart)*100.    << "%)" << endl;
		Info << "   CPU Retrieve     : " << cpuTimeRet          << " (" << cpuTimeRet/(tEnd-tStart)*100.         << "%)" << endl;
		Info << endl;

		Info << "      BTS  : " << isatTable->nBTS()  << endl;
		Info << "      MRU  : " << isatTable->nMRU()  << endl;
		Info << "      MFU  : " << isatTable->nMFU()  << endl << endl;
		Info << endl;
	}
	#endif

	WriteScalarField(mesh, mDotStar,  "edc::mDotStar");
	WriteScalarField(mesh, gammaStar, "edc::gammaStar");
	WriteScalarField(mesh, tfStar, 	  "edc::tfStar");
}

edcCounter++;

