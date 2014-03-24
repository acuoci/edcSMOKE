/*  
  Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<psiThermo> pThermo
    (
        psiThermo::New(mesh)
    );
    psiThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

*/


	Info<< "Creating reaction model\n" << endl;

	autoPtr<combustionModels::psiCombustionModel> reaction
	(
	    combustionModels::psiCombustionModel::New(mesh)
	);

	psiReactionThermo& thermo = reaction->thermo();
	thermo.validate(args.executable(), "h", "e");

	basicMultiComponentMixture& composition = thermo.composition();
	PtrList<volScalarField>& Y = composition.Y();

	word inertSpecie(thermo.lookup("inertSpecie"));



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


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);

    dimensionedScalar rhoMax(simple.dict().lookup("rhoMax"));
    dimensionedScalar rhoMin(simple.dict().lookup("rhoMin"));

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

	// Set the turbulence into the reaction model
	reaction->setTurbulence(turbulence());

    dimensionedScalar initialMass = fvc::domainIntegrate(rho);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());
