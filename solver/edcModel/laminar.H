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

#if STEADYSTATE == 1
{
	const int ns = Y.size();
	const int ncells = mesh.nCells();
	OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);

	if (chemistryLinearModel_().implicitSourceTerm() == true && chemistryLinearModel_().jacobianUpdate() == 1)
	{
		Info<< " * Implicit chemical source terms evaluation (continuous)... " << endl;
	
		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
		{
			const scalarField& pCells = thermo.p().internalField();
			const scalarField& TCells = thermo.T().internalField();
			scalarField& QdotCells = Qdot.ref();
			scalarField& dQCells = dQ.ref();

			Eigen::VectorXd J(ns+1);
			Eigen::VectorXd source(ns+1);	
			Eigen::VectorXd y(ns+1);

			forAll(TCells, celli)
			{
				for(int i=0;i<ns;i++)
					y(i) = Y[i].internalField()[celli];
				y(ns) = TCells[celli];

				chemistryLinearModel_().reactionSourceTerms( y, pCells[celli], source);

				if (chemistryLinearModel_().sparseJacobian() == false)
					chemistryLinearModel_().reactionJacobian( y, pCells[celli], J );
				else
					chemistryLinearModel_().reactionJacobianSparse( y, pCells[celli], J, true, true);

				for(int i=0;i<ns+1;i++)
					sourceImplicit_[i].ref()[celli] = J(i);
		
				for(int i=0;i<ns+1;i++)
					sourceExplicit_[i].ref()[celli] = source(i) - J(i)*y(i);


				// Explicit heat release rate
				thermodynamicsMapXML->SetTemperature(TCells[celli]);
				thermodynamicsMapXML->SetPressure(pCells[celli]);
				thermodynamicsMapXML->hMolar_Species(h_species.GetHandle());

				double HeatRelease = 0.;
				for(unsigned int i=0;i<ns;i++)
					HeatRelease += -h_species[i+1]/thermodynamicsMapXML->MW(i)*source(i);
				QdotCells[celli] = HeatRelease;
				dQCells[celli] = HeatRelease * mesh.V()[celli];
			}
		}
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(ncells)*1000. << " ms per cell)" << endl;
	}

	else if (chemistryLinearModel_().implicitSourceTerm() == true && chemistryLinearModel_().jacobianUpdate() != 1)
	{
		Info<< " * Implicit chemical source terms evaluation (discrete)... " << endl;

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
		{
			const scalarField& pCells = thermo.p().internalField();
			const scalarField& TCells = thermo.T().internalField();
			scalarField& QdotCells = Qdot.ref();
			scalarField& dQCells = dQ.ref();

			Eigen::VectorXd J(ns+1);
			Eigen::VectorXd source(ns+1);	
			Eigen::VectorXd y(ns+1);

			forAll(TCells, celli)
			{
				for(int i=0;i<ns;i++)
					y(i) = Y[i].internalField()[celli];
				y(ns) = TCells[celli];

				chemistryLinearModel_().reactionSourceTerms( y, pCells[celli], source );

				if (chemistryLinearModel_().jacobianCounter() == chemistryLinearModel_().jacobianUpdate())
				{
					if (chemistryLinearModel_().sparseJacobian() == false)
						chemistryLinearModel_().reactionJacobian( y, pCells[celli], J );
					else
						chemistryLinearModel_().reactionJacobianSparse( y, pCells[celli], J, true, true);

					Jstored_[celli] = J;
				}

				for(int i=0;i<ns+1;i++)
					sourceImplicit_[i].ref()[celli] = Jstored_[celli](i);
	
				for(int i=0;i<ns+1;i++)
					sourceExplicit_[i].ref()[celli] = source(i) - Jstored_[celli](i)*y(i);


				// Explicit heat release rate
				thermodynamicsMapXML->SetTemperature(TCells[celli]);
				thermodynamicsMapXML->SetPressure(pCells[celli]);
				thermodynamicsMapXML->hMolar_Species(h_species.GetHandle());

				double HeatRelease = 0.;
				for(unsigned int i=0;i<ns;i++)
					HeatRelease += -h_species[i+1]/thermodynamicsMapXML->MW(i)*source(i);
				QdotCells[celli] = HeatRelease;
				dQCells[celli] = HeatRelease * mesh.V()[celli];
			}
		}
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(ncells)*1000. << " ms per cell)" << endl;

		// Update the counter
		chemistryLinearModel_().updateJacobianCounter();
	}

	else
	{
		Info<< " * Explicit source terms evaluation... " << endl;

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
		{
			const scalarField& pCells = thermo.p().internalField();
			const scalarField& TCells = thermo.T().internalField();
			scalarField& QdotCells = Qdot.ref();
			scalarField& dQCells = dQ.ref();

			Eigen::VectorXd source(ns+1);	
			Eigen::VectorXd y(ns+1);

			forAll(TCells, celli)
			{
				for(int i=0;i<ns;i++)
					y(i) = Y[i].internalField()[celli];
				y(ns) = TCells[celli];

				chemistryLinearModel_().reactionSourceTerms( y, pCells[celli], source );

				for(int i=0;i<ns+1;i++)
					sourceImplicit_[i].ref()[celli] = 0.;
	
				for(int i=0;i<ns+1;i++)
					sourceExplicit_[i].ref()[celli] = source(i);

				// Explicit heat release rate
				thermodynamicsMapXML->SetTemperature(TCells[celli]);
				thermodynamicsMapXML->SetPressure(pCells[celli]);
				thermodynamicsMapXML->hMolar_Species(h_species.GetHandle());

				double HeatRelease = 0.;
				for(unsigned int i=0;i<ns;i++)
					HeatRelease += -h_species[i+1]/thermodynamicsMapXML->MW(i)*source(i);
				QdotCells[celli] = HeatRelease;
				dQCells[celli] = HeatRelease * mesh.V()[celli];
			}
		}
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(ncells)*1000. << " ms per cell)" << endl;
	}

	if (laminarDiffusion == true)
	{
		const scalarField& TCells = thermo.T().internalField();
		const scalarField& pCells = thermo.p().internalField();

		forAll(TCells, celli)
		{	
			OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);

			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->hMolar_Species(h_species.GetHandle());

			for(unsigned int i=0;i<ns;i++)
				HeSpecies_[i].ref()[celli] = -h_species[i+1]/thermodynamicsMapXML->MW(i);
		}
	}

	// Correct boundary conditions
	Qdot.correctBoundaryConditions();
	dQ.correctBoundaryConditions();
}

#else
{
	if (counterToUpdateChemistry == numberIterationsToUpdateChemistry)
	{
		counterToUpdateChemistry = 0;

		bool debug = false;

		Info << "Chemistry integration (laminar combustion model)..." << endl;
		const int unsigned ns = Y.size();
		const int unsigned ne = ns+1;

		Eigen::VectorXd residuals(ne);
		Eigen::VectorXd yStar0(ne);
		Eigen::VectorXd yStarf(ne);
		Eigen::VectorXd yMin(ne);
		Eigen::VectorXd yMax(ne);
	
		for(unsigned int i=0;i<ns;i++)	yMin(i) = 0.; 	yMin(ns) = minT;
		for(unsigned int i=0;i<ns;i++)	yMax(i) = 1.; 	yMax(ns) = maxT;

		OpenSMOKE::OpenSMOKEVectorDouble omegaMean(ns);
		OpenSMOKE::OpenSMOKEVectorDouble omegaStar(ns);
		OpenSMOKE::OpenSMOKEVectorDouble xMean(ns);
		OpenSMOKE::OpenSMOKEVectorDouble cMean(ns);
		OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);
	
		double hMean, MWMean;
		const scalarField& pCells = thermo.p().internalField();
		const scalarField& TCells = thermo.T().internalField();
		const scalarField& rhoCells = rho.internalField();

		scalarField& QdotCells = Qdot.ref();
		scalarField& dQCells = dQ.ref();

		// Statistics
		unsigned int counter = 0;
		double cpuTimeAllCells_ODE = 0.;

		// Integration time
		const double deltat = runTime.deltaT().value()+1e-12;

		const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		forAll(TCells, celli)
		{	
			// Extract the mean mass fractions
			for(unsigned int i=0;i<ns;i++)
				omegaMean[i+1] = std::max(1.e-16, Y[i].internalField()[celli]);

			// Normalize the mass fractions
			double sumOmegaMean = 0.;
			for(unsigned int i=1;i<=ns;i++)
				sumOmegaMean += omegaMean[i];
			for(unsigned int i=1;i<=ns;i++)
				omegaMean[i] /= sumOmegaMean;		

			// Convert to mole fractions
			thermodynamicsMapXML->MoleFractions_From_MassFractions(xMean.GetHandle(), MWMean, omegaMean.GetHandle());

			// Extract the fine strucutres mass fractions
			{
				for(unsigned int i=0;i<ns;i++)
					yStar0(i) = omegaMean[i+1];
				yStar0(ns) = TCells[celli];

				fine_structure_ode_pfr.Set(pCells[celli]);
				fine_structure_ode_pfr.SetCheckMassFractions(true);
			}			

			// Solve the ODE system
			{
				const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

				#include "solvePFR_DI_laminar.H"
			
				const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

				cpuTimeAllCells_ODE += (tEnd-tStart);
			}
	
			{		
				for(unsigned int i=0;i<ns;i++)
					YStar_[i].ref()[celli] = yStarf(i);
				TStar_[0].ref()[celli] = yStarf(ns);
		 
				const double coefficient = rhoCells[celli] / deltat;
				for(unsigned int i=0;i<ns;i++)
					R_[i].ref()[celli] = coefficient*(YStar_[i].internalField()[celli]-Y[i].internalField()[celli]);		
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
				for(unsigned int i=0;i<ns;i++)
					HeSpecies_[i].ref()[celli] = -h_species[i+1]/thermodynamicsMapXML->MW(i);
			}

			if (counter%(int(0.20*mesh.nCells())+1) == 0)
			{
				Info << "   Accomplished: " << counter << "/" << mesh.nCells() << endl;
			}

			counter++;
		}
		const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
		Info << "CPU time for chemistry (laminar): " << tEnd - tStart << " s " << "(" << (tEnd - tStart)/double(TCells.size())*1e3 << " ms per cell)" << endl;
		Info << "CPU time for ODE solution:        " << cpuTimeAllCells_ODE << " s " << "(" << (cpuTimeAllCells_ODE)/double(TCells.size())*1e3 << " ms per cell)" << endl;
	
		// Correct boundary conditions
		Qdot.correctBoundaryConditions();
		dQ.correctBoundaryConditions();
	}

	counterToUpdateChemistry++;
}
#endif


