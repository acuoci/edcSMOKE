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

	Info<< "Creating reaction model\n" << endl;

	autoPtr<fluidReactionThermo> pThermo
	(
		fluidReactionThermo::New(mesh)
	);
	fluidReactionThermo& thermo = pThermo();
	thermo.validate(args.executable(), "h", "e");

	basicSpecieMixture& composition = thermo.composition();
	PtrList<volScalarField>& Y = composition.Y();

	volScalarField rho
	(
		IOobject
		(
			"rho",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		thermo.rho()
	);

	volScalarField& p = thermo.p();
	const volScalarField& psi = thermo.psi();
	volScalarField& T = thermo.T();

	Info<< "Reading field U\n" << endl;
	volVectorField U
	(
		IOobject
		(
			"U",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);

	#include "compressibleCreatePhi.H"

	pressureReference pressureReference
	(
	    p,
	    simple.dict(),
	    thermo.incompressible()
	);

	mesh.schemes().setFluxRequired(p.name());


	Info << "Creating turbulence model.\n" << nl;
	autoPtr<compressible::momentumTransportModel> turbulence
	(
	    compressible::momentumTransportModel::New
	    (
		rho,
		U,
		phi,
		thermo
	    )
	);


	Info<< "Creating reaction model\n" << endl;
	autoPtr<combustionModel> reaction(combustionModel::New(thermo, turbulence()));

    	dimensionedScalar initialMass = fvc::domainIntegrate(rho);

	#include "createMRF.H"
	#include "createFvModels.H"
	#include "createFvConstraints.H"

    	multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

	forAll(Y, i)
	{
	    fields.add(Y[i]);
	}
	fields.add(thermo.he());


	Info<< "Creating field number of active species\n" << endl;

	volScalarField drg_important_species
	(
	    IOobject
	    (
		"drg::important_species",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("drg::important_species", dimensionSet(0, 0, 0, 0, 0), 0.0)
	);

	volScalarField drg_important_reactions
	(
	    IOobject
	    (
		"drg::important_reactions",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("drg::important_reactions", dimensionSet(0, 0, 0, 0, 0), 0.0)
	);

	volScalarField cpu_reaction_times
	(
	    IOobject
	    (
		"cpu::reaction_times",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("cpu::reaction_times", dimensionSet(0, 0, 0, 0, 0), 0.0)
	);

	volScalarField mix_nu
	(
	    IOobject
	    (
		"mix::nu",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("mix::nu", dimensionSet(0, 2, -1, 0, 0), 0)
	);

	volScalarField mix_mu
	(
	    IOobject
	    (
		"mix::mu",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("mix::mu", dimensionSet(1, -1, -1, 0, 0), 0)
	);


