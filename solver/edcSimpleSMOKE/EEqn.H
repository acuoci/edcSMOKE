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

if (energyEquation == true)
{
	if (combustionModel == "laminar")
	{
		T.max(minT);
	    	T.min(maxT);

		if (laminarDiffusion == true)
		{
			dimensionedScalar zero("zero", dimensionSet(1,-1,-3,0,0,0,0), scalar(0.)); 
			massDiffusionInEnergyEquation = zero;

			forAll(Y, i)
			{ 
				volScalarField& Yi = Y[i];
				volScalarField& Dmixi = Dmix_[i];
				volScalarField& HeSpeciesi = HeSpecies_[i];

				massDiffusionInEnergyEquation += HeSpeciesi* (fvc::laplacian(rho*Dmixi, Yi)); 
			}
		}

		// Energy equation
		{
			volScalarField& he = thermo.he();

			fvScalarMatrix TEqn
			(
				  thermo.Cp()*fvm::div(phi, T)	    
				==
				  fvm::laplacian(rho*thermo.Cp()*turbulence->nut()/Prt,T) 
		                + massDiffusionInEnergyEquation
				+ ( radiation->Ru() - fvm::Sp(radiation->Rp()*pow3(T), T) )
				+ sourceExplicit_[Y.size()] + fvm::Sp(sourceImplicit_[Y.size()],T)
				+ thermo.Cp()*fvModels.source(rho, T)
			);

			TEqn.relax();

			fvConstraints.constrain(TEqn);

			TEqn.solve("T");

			fvConstraints.constrain(T);

			Info<< " * T gas min/max (after transport)   = " << min(T).value() << ", " << max(T).value() << endl;

	    		thermo.correct();
	    		radiation->correct();
		}

		T.max(minT);
	    	T.min(maxT);
	
		Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;

		#include "sparkPatch.H"
	}
	else
	{
		T.max(minT);
	    	T.min(maxT);

		if (combustionModel == "EDC")
		{
			TStar_[0].max(minT);
	    		TStar_[0].min(maxT);
		}
		volScalarField& he = thermo.he();

		if (laminarDiffusion == true)
		{
			dimensionedScalar zero("zero", dimensionSet(1,-1,-3,0,0,0,0), scalar(0.)); 
			massDiffusionInEnergyEquation = zero;

			forAll(Y, i)
			{ 
				volScalarField& Yi = Y[i];
				volScalarField& Dmixi = Dmix_[i];
				volScalarField& HeSpeciesi = HeSpecies_[i];

				massDiffusionInEnergyEquation += HeSpeciesi* (fvc::laplacian(rho*Dmixi, Yi)); 
			}
		}

		// Energy equation
		{
		
			fvScalarMatrix EEqn
			(
				fvm::div(phi, he)
				+ (
				    he.name() == "e"
				  ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
				  : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
				  )
				- fvm::laplacian(rho*turbulence->nut()/Prt, he)
				- radiation->Sh(thermo, he)
		                + massDiffusionInEnergyEquation
				==
				Qdot +
				fvModels.source(rho, he)
			);

	    		EEqn.relax();

	    		fvConstraints.constrain(EEqn);

	   		EEqn.solve();

	    		fvConstraints.constrain(he);

	    		thermo.correct();
	    		radiation->correct();
		}

		T.max(minT);
	    	T.min(maxT);
	
		Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;

		#include "sparkPatch.H"
	}
}
