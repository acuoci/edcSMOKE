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

#include "createRDeltaT.H"

Info<< nl << "Reading thermophysicalProperties" << endl;

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
        mesh
    ),
    thermo.rho()
);

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

volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
      volScalarField& T = thermo.T();

#include "compressibleCreatePhi.H"

pressureReference pressureReference(p, pimple.dict(), false);
mesh.schemes().setFluxRequired(p.name());


dimensionedScalar rhoMax
(
    dimensionedScalar::lookupOrDefault
    (
	"rhoMax",
	pimple.dict(),
	dimDensity,
	GREAT
    )
);

dimensionedScalar rhoMin
(
    dimensionedScalar::lookupOrDefault
    (
	"rhoMin",
	pimple.dict(),
	dimDensity,
	0
    )
);

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

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidReactionThermophysicalTransportModel>
thermophysicalTransport
(
    fluidReactionThermophysicalTransportModel::New
    (
        turbulence(),
        thermo
    )
);


Info<< "Creating reaction model\n" << endl;
autoPtr<combustionModel> reaction(combustionModel::New(thermo, turbulence()));


Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

fileName path = args.rootPath()+"/"+args.globalCaseName();
OFstream Pav(path+"/p.Vessel");


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


