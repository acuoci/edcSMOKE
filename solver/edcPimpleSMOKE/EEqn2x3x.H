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
	T.max(minT);
	T.min(maxT);

	if (combustionModel == "EDC")
	{
		TStar_[0].max(minT);
		TStar_[0].min(maxT);
	}

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
		- fvm::laplacian(turbulence->mut()/Prt, he)
                + massDiffusionInEnergyEquation
		==
		Qdot
		+ radiation->Sh(thermo)
		+ fvOptions(rho, he)
	);

	EEqn.relax();
	fvOptions.constrain(EEqn);
	EEqn.solve();
	fvOptions.correct(he);

	radiation->correct();
	thermo.correct();

	Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;

	#include "sparkPatch.H" 
}

