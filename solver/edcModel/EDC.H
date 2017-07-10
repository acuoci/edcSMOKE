/*-----------------------------------------------------------------------*\
|                  _       _____ __  __  ____  _  ________                |
|                 | |     / ____|  \/  |/ __ \| |/ /  ____|               |
|          ___  __| | ___| (___ | \  / | |  | | ' /| |__                  |
|         / _ \/ _` |/ __|\___ \| |\/| | |  | |  < |  __|                 |
|        |  __/ (_| | (__ ____) | |  | | |__| | . \| |____                |
|         \___|\__,_|\___|_____/|_|  |_|\____/|_|\_\______|               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, M.R. Malik, Z. Li, A. Parente                      |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|   Contacts: Mohammad Rafi Malik, Zhiyi Li, Alessandro Parente           |
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
|   Copyright(C) 2017-2014 A. Cuoci, A. Parente                           |
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

{
	if (counterToUpdateChemistry == numberIterationsToUpdateChemistry)
	{
		counterToUpdateChemistry = 0;

		bool debug = false;

		Info << "Fine structures properties..." << endl;
		#include "propertiesFineStructures.H"

		Info << "Fine structures integration..." << endl;
		const int unsigned ns = Y.size();
		const int unsigned ne = ne_total;
			  unsigned ne_drg = 0;

		Eigen::VectorXd residuals(ne);
		Eigen::VectorXd yStar0(ne);
		Eigen::VectorXd yStarf(ne);
		Eigen::VectorXd yMin(ne);
		Eigen::VectorXd yMax(ne);
	
		for(unsigned int i=0;i<ns;i++)	yMin(i) = 0.; yMin(ne-3) = 200.;	yMin(ne-2) = 0.;	yMin(ne-1) = 0.;
		for(unsigned int i=0;i<ns;i++)	yMax(i) = 1.; yMax(ne-3) = 5000.;	yMax(ne-2) = 1.e16;	yMax(ne-1) = 1.e16;

		OpenSMOKE::OpenSMOKEVectorDouble omegaMean(ns);
		OpenSMOKE::OpenSMOKEVectorDouble omegaStar(ns);
		OpenSMOKE::OpenSMOKEVectorDouble xMean(ns);
		OpenSMOKE::OpenSMOKEVectorDouble cMean(ns);
		OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);
	
		double hMean, MWMean;
		const scalarField& pCells = thermo.p().internalField();
		const scalarField& TCells = thermo.T().internalField();
		#if OPENFOAM_VERSION == 4
		scalarField& QdotCells = Qdot.ref();
		scalarField& dQCells = dQ.ref();
		#else
		scalarField& QdotCells = Qdot.internalField();
		scalarField& dQCells = dQ.internalField();
		#endif

		// DRG Fields
		#if OPENFOAM_VERSION == 4
		scalarField& drgSpeciesCells = drg_important_species.ref();
		scalarField& drgReactionsCells = drg_important_reactions.ref();
		scalarField& cpuReactionsCells = cpu_reaction_times.ref();
		#else
		scalarField& drgSpeciesCells = drg_important_species.internalField();
		scalarField& drgReactionsCells = drg_important_reactions.internalField();
		scalarField& cpuReactionsCells = cpu_reaction_times.internalField();
		#endif


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

		// Set DRG Analysis
		if (drg_analysis == true) 
			drg->ResetCounters();
	
		// Statistics for PSR/NLS
		double nls_norm2_cumulative_ = 0.;
		double nls_norm2_max_ = 0.;
		int nls_failures_ = 0;
		int nls_number_iterations_cumulative_ = 0;

		// Statistics
		unsigned int counter = 0;
		double cpuTimeAllCells_ODE = 0.;

		const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		forAll(TCells, celli)
		{	
			// Integration time
			const double t0 = 0.;
			double tf = 0.;

			// Extract the mean mass fractions
			for(unsigned int i=0;i<ns;i++)
				omegaMean[i+1] = Y[i].internalField()[celli];

			// Convert to mole fractions
			thermodynamicsMapXML->MoleFractions_From_MassFractions(xMean.GetHandle(), MWMean, omegaMean.GetHandle());

			// Extract the fine strucutres mass fractions
			if (reactorType == REACTORTYPE_PSR)
			{
				if (drg_analysis == false)
				{
					for(unsigned int i=0;i<ns;i++)
						yStar0(i) = YStar_[i].internalField()[celli];
					yStar0(ne-3) = TStar_[0].internalField()[celli];
					yStar0(ne-2) = (1.-gammaStarDen[celli])/max(mDotStar[celli], SMALL);
					yStar0(ne-1) = gammaStarNum[celli];

					// Extract the star mass fractions
					for(unsigned int i=1;i<=ns;i++)
						omegaStar[i] = omegaMean[i];

					// Initialize the fine structure
					thermodynamicsMapXML->SetTemperature(TCells[celli]);
					thermodynamicsMapXML->SetPressure(pCells[celli]);
					hMean = thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(xMean.GetHandle());
					hMean/=MWMean;	
			
					fine_structure_ode_psr.Set(TCells[celli], pCells[celli], omegaMean, omegaStar, hMean, gammaStarDen[celli], mDotStar[celli]);
					if(solverType == SOLVERTYPE_NLS)
					{
						int index_max;
						omegaStar.Max(&index_max);
						fine_structure_nls_psr.Set(TCells[celli], pCells[celli], omegaMean, omegaStar, hMean, gammaStarDen[celli], mDotStar[celli], index_max);
					}

					tf = numberOfResidenceTimes/max(mDotStar[celli], SMALL);
				}
				else
				{
					// Set epsilon for DRG
					{
						for(unsigned int i=0;i<drg_temperature.size();i++)
							if(TCells[celli] < drg_temperature[i])
							{
								drg->SetEpsilon(drg_epsilon[i]);
								break;
							}
					}

					// DRG Analysis
					{
						double cTot_ = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
						Product(cTot_, xMean, &cMean);
						drg->Analysis(TCells[celli], pCells[celli], cMean);
					}

					// Number of import species (+1)
					ne_drg = drg->number_important_species()+1;

					// Set minimum and maximum values
					{
						yStar0.resize(ne_drg);
						yStarf.resize(ne_drg);
						yMax.resize(ne_drg);	yMax.setConstant(1.); yMax(ne_drg-1) = 5000.;
						yMin.resize(ne_drg);	yMin.setConstant(0.); yMin(ne_drg-1) = 280.;
					}

					// Set initial condition
					for (unsigned int i=0;i<drg->number_important_species();++i)	
					{
						const unsigned int j = drg->indices_important_species()[i];
						yStar0(i) = YStar_[j].internalField()[celli];
					}
					yStar0(ne_drg-1) = TCells[celli];

					// Extract the star mass fractions
					for(unsigned int i=1;i<=ns;i++)
					{
						omegaStar[i] = YStar_[i-1].internalField()[celli];
					}

					// Set DRG
					fine_structure_ode_psr.SetDRG(drg);
					if(solverType == SOLVERTYPE_NLS)
						fine_structure_nls_psr.SetDRG(drg);

					// Initialize the fine structure
					thermodynamicsMapXML->SetTemperature(TCells[celli]);
					thermodynamicsMapXML->SetPressure(pCells[celli]);
					hMean = thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(xMean.GetHandle());
					hMean/=MWMean;	

					fine_structure_ode_psr.Set(TCells[celli], pCells[celli], omegaMean, omegaStar, hMean, gammaStarDen[celli], mDotStar[celli]);
					if(solverType == SOLVERTYPE_NLS)
					{
						// The species with the highest mass fraction is selected among the important species
						double max_value = 0.;
						unsigned int index_max = 0;
						for (unsigned int i=0;i<drg->number_important_species();++i)	
						{
							const unsigned int j = drg->indices_important_species()[i]+1;
							if (omegaStar[j] > max_value)
							{
								max_value = omegaStar[j];
								index_max = i;
							}
						}

						fine_structure_nls_psr.Set(TCells[celli], pCells[celli], omegaMean, omegaStar, hMean, gammaStarDen[celli], mDotStar[celli], index_max);
					}

					tf = numberOfResidenceTimes/max(mDotStar[celli], SMALL);
				}
			}
			else if (reactorType == REACTORTYPE_PFR)
			{
				if (drg_analysis == false)
				{
					tf = (1.-gammaStarDen[celli])/max(mDotStar[celli], SMALL);

					for(unsigned int i=0;i<ns;i++)
						yStar0(i) = omegaMean[i+1];
					yStar0(ne-3) = TCells[celli];
					yStar0(ne-2) = tf;
					yStar0(ne-1) = 0.;

					fine_structure_ode_pfr.Set(pCells[celli]);
				}
				else
				{
					// Set epsilon for DRG
					{
						for(unsigned int i=0;i<drg_temperature.size();i++)
							if(TCells[celli] < drg_temperature[i])
							{
								drg->SetEpsilon(drg_epsilon[i]);
								break;
							}
					}

					// DRG Analysis
					{
						double cTot_ = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
						Product(cTot_, xMean, &cMean);
						drg->Analysis(TCells[celli], pCells[celli], cMean);
					}

					// Number of import species (+1)
					ne_drg = drg->number_important_species()+1;

					// Set minimum and maximum values
					{
						yStar0.resize(ne_drg);
						yStarf.resize(ne_drg);
						yMax.resize(ne_drg);	yMax.setConstant(1.); yMax(ne_drg-1) = 5000.;
						yMin.resize(ne_drg);	yMin.setConstant(0.); yMin(ne_drg-1) = 280.;
					}

					// Set integration time
					tf = (1.-gammaStarDen[celli])/max(mDotStar[celli], SMALL);

					// Set initial condition
					for (unsigned int i=0;i<drg->number_important_species();++i)	
					{
						const unsigned int j = drg->indices_important_species()[i]+1;
						yStar0(i) = omegaMean[j];
					}
					yStar0(ne_drg-1) = TCells[celli];

					// Extract the star mass fractions
					for(unsigned int i=1;i<=ns;i++)
						omegaStar[i] = omegaMean[i];

					// Set DRG
					fine_structure_ode_pfr.SetDRG(drg);

					// Set the PFR solver
					fine_structure_ode_pfr.Set(pCells[celli], omegaStar);
				}
			}			

			// Solve the ODE system
			{
				const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

				if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
				{
					if (reactorType == REACTORTYPE_PSR)
					{

						#if EDCSMOKE_USE_ISAT == 1

							if(isatCheck == true)
							{
								#include "solvePSR_ISAT.H"
							}
							else
							{
								if (drg_analysis == false)
								{
									#include "solvePSR_DI.H"
								}
								else
								{
									#include "solvePSR_DI_DRG.H"
								}
							}

						#else
							if (drg_analysis == false)
							{
								#include "solvePSR_DI.H"
							}
							else
							{
								#include "solvePSR_DI_DRG.H"
							}

						#endif
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
								if (drg_analysis == false)
								{
									#include "solvePFR_DI.H"
								}
								else
								{
									#include "solvePFR_DI_DRG.H"
								}
							}

						#else
							if (drg_analysis == false)
							{
								#include "solvePFR_DI.H"
							}
							else
							{
								#include "solvePFR_DI_DRG.H"
							}

						#endif
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
			
				const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

				cpuTimeAllCells_ODE += (tEnd-tStart);
			}
	
			if (drg_analysis == false)
			{		
				#if OPENFOAM_VERSION == 4
						for(unsigned int i=0;i<ns;i++)
							YStar_[i].ref()[celli] = yStarf(i);
						TStar_[0].ref()[celli] = yStarf(ne-3);
				 
						for(unsigned int i=0;i<ns;i++)
							R_[i].ref()[celli] = edc_coefficient[celli]*(YStar_[i].internalField()[celli]-Y[i].internalField()[celli]);
				#else
						for(unsigned int i=0;i<ns;i++)
							YStar_[i].internalField()[celli] = yStarf(i);
						TStar_[0].internalField()[celli] = yStarf(ne-3);
				 
						for(unsigned int i=0;i<ns;i++)
							R_[i].internalField()[celli] = edc_coefficient[celli]*(YStar_[i].internalField()[celli]-Y[i].internalField()[celli]);
				#endif
			}
			else
			{
				#if OPENFOAM_VERSION == 4
					for(unsigned int i=0;i<ne_drg-1;i++)
					{
			       			const unsigned int j = drg->indices_important_species()[i];
						YStar_[j].ref()[celli] = yStarf(i);
					}
					TStar_[0].ref()[celli] = yStarf(ne_drg-1);
			 
					// Reset to 0 every formation rate
					for(unsigned int i=0;i<ns;i++)
						R_[i].ref()[celli] = 0.;

					// Update only formation rates for relevant species
					for(unsigned int i=0;i<ne_drg-1;i++)
					{
						const unsigned int j = drg->indices_important_species()[i];
						R_[j].ref()[celli] = edc_coefficient[celli]*(YStar_[j].internalField()[celli]-Y[j].internalField()[celli]);
					}
				#else
					for(unsigned int i=0;i<ne_drg-1;i++)
					{
			       			const unsigned int j = drg->indices_important_species()[i];
						YStar_[j].internalField()[celli] = yStarf(i);
					}
					TStar_[0].internalField()[celli] = yStarf(ne_drg-1);
			 
					// Reset to 0 every formation rate
					for(unsigned int i=0;i<ns;i++)
						R_[i].internalField()[celli] = 0.;

					// Update only formation rates for relevant species
					for(unsigned int i=0;i<ne_drg-1;i++)
					{
						const unsigned int j = drg->indices_important_species()[i];
						R_[j].internalField()[celli] = edc_coefficient[celli]*(YStar_[j].internalField()[celli]-Y[j].internalField()[celli]);
					}
				#endif
				
			}

			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->hMolar_Species(h_species.GetHandle());

			double HeatRelease = 0.;
			for(unsigned int i=0;i<ns;i++)
				HeatRelease += -h_species[i+1]/thermodynamicsMapXML->MW(i) * R_[i].internalField()[celli];
			QdotCells[celli] = HeatRelease;
			dQCells[celli] = HeatRelease * mesh.V()[celli];

			if (laminarDiffusion == true)
			{
				#if OPENFOAM_VERSION == 4
				for(unsigned int i=0;i<ns;i++)
					HeSpecies_[i].ref()[celli] = -h_species[i+1]/thermodynamicsMapXML->MW(i);
				#else
				for(unsigned int i=0;i<ns;i++)
					HeSpecies_[i].internalField()[celli] = -h_species[i+1]/thermodynamicsMapXML->MW(i);
				#endif
			}

			if (debug == true)
			{
				Info << celli << " " << TCells[celli] << " " << gammaStarNum[celli] << " " << mDotStar[celli] << " " << HeatRelease << endl;
				for(unsigned int i=0;i<ns;i++)
				{
					Info << Y[i].name() << " " << yStar0[i] << " " << yStarf[i] << " " << omegaMean[i+1] << " " << R_[i].internalField()[celli] << endl;
				}
				getchar();
			}

			if (counter%(int(0.20*mesh.nCells())+1) == 0)
			{
				if (solverType == SOLVERTYPE_ODE)
					Info << "   Accomplished: " << counter << "/" << mesh.nCells() << endl;
				else
				{
					Info 	<< "  Accomplished: " 		<< counter << "/" << mesh.nCells() << " ";
					Info 	<< "  Norm2(avg): " 		<< nls_norm2_cumulative_/static_cast<double>(counter+1) 	
						<< "  Norm2(max): " 		<< nls_norm2_max_ 
						<< "  Iterations(avg): " 	<< nls_number_iterations_cumulative_/static_cast<double>(counter+1) 
						<< "  Ode: " 			<< nls_failures_ << endl;
				}
			}

			counter++;

			if (drg_analysis == true)
			{
				// Fill the relevant DRG fields
				drgSpeciesCells[celli]   = drg->number_important_species();
				drgReactionsCells[celli] = drg->number_important_reactions();
			}
		}
		const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
		Info << "CPU time for chemistry (EDC): " << tEnd - tStart << " s " << "(" << (tEnd - tStart)/double(TCells.size())*1e3 << " ms per cell)" << endl;
		Info << "CPU time for ODE solution:    " << cpuTimeAllCells_ODE << " s " << "(" << (cpuTimeAllCells_ODE)/double(TCells.size())*1e3 << " ms per cell)" << endl;
		if (drg_analysis == true) 
		{
			Info << "CPU time for DRG analysis:    " << drg->cpuTimeOverall() << " s " << "(" << (drg->cpuTimeOverall())/double(TCells.size())*1e3 << " ms per cell)" << endl;
			drg->WriteCpuTimesOnTheScreen();
		}

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
	
		// Correct boundary conditions
		Qdot.correctBoundaryConditions();
		dQ.correctBoundaryConditions();
	}

	counterToUpdateChemistry++;
}

