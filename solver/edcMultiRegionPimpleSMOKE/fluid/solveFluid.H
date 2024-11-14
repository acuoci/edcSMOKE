if (!pimple.flow())
{
    if (pimple.models())
    {
        fvModels.correct();
    }

    if (pimple.thermophysics())
    {
        tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);

	#include "properties.H"

        if (Ecorr == 0)
        {
            #include "YEqn.H"
        }

        #include "EEqn.H"
    }
}
else
{
    if (Ecorr == 0)
    {
        if (!mesh.schemes().steady() && pimples.firstPimpleIter())
        {
            #include "rhoEqn.H"
        }

        if (pimple.models())
        {
            fvModels.correct();
        }

        #include "UEqn.H"
    }

    if (pimple.thermophysics())
    {
        tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);

	#include "properties.H"

        if (Ecorr == 0)
        {
            #include "YEqn.H"
        }

        #include "EEqn.H"
    }

    if (Ecorr == nEcorr - 1)
    {
Info << "Here 1" << endl; 
        tmp<fvVectorMatrix>& tUEqn = UEqns[0];
Info << "Here 2" << endl; 
        fvVectorMatrix& UEqn = tUEqn.ref();
Info << "Here 3" << endl; 
        // --- PISO loop
        while (pimple.correct())
        {
            #include "pEqn.H"
        }
Info << "Here 4" << endl; 
        if (pimples.pimpleTurbCorr(i))
        {
             turbulence->correct();
        }

        if (!mesh.schemes().steady() && pimples.finalPimpleIter())
        {
            rho = thermo.rho();
        }
    }
}
