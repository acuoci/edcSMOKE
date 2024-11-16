    const wordList fluidNames
    (
        rp.found("fluid") ? rp["fluid"] : wordList(0)
    );

    PtrList<fvMesh> fluidRegions(fluidNames.size());

    forAll(fluidNames, i)
    {
        Info << "Create fluid mesh for region " << fluidNames[i] << " for time = " << runTime.timeName() << nl << endl;

        fluidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    fluidNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }

    if (fluidNames[0] != "reacting")
	FatalErrorIn(args.executable()) << "The first fluid mesh must have reacting name" << exit(FatalError);
