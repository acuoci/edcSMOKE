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

// This is a unsteady simulation
#define STEADYSTATE 0

// This is a multi-region solver
#define MULTIREGIONSOLVER 1

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// OpenSMOKE++ Dictionaries
#include "dictionary/OpenSMOKE_Dictionary"

// ODE solvers
#include "math/native-ode-solvers/MultiValueSolver"
#include "math/external-ode-solvers/ODE_Parameters.h"

// NLS solvers
#include "math/native-nls-solvers/NonLinearSystemSolver"
#include "math/native-nls-solvers/parameters/NonLinearSolver_Parameters.h"

// OpenFOAM
#include "fvCFD.H"
#include "fluidReactionThermo.H"
#include "combustionModel.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidReactionThermophysicalTransportModel.H"
#include "multivariateScheme.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "coordinateSystem.H"
#include "pimpleMultiRegionControl.H"
#include "pressureReference.H"
#include "hydrostaticInitialisation.H"
#include "radiationModel.H"

// Solid
#include "compressibleCourantNo.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"

// Utilities
#include "Utilities.H"

// DRG
#include "DRG.H"

// ODE systems
#include "ODE_PSR.H"
#include "ODE_PSR_Interface.H"
#include "ODE_PFR.H"
#include "ODE_PFR_Interface.H"
#include "ODE_PFR_Laminar.H"
#include "ODE_PFR_Laminar_Interface.H"

// NLS Systems
#include "NLS_PSR.H"
#include "NLS_PSR_Interface.H"

// Characteristic chemical times
#include "CharacteristicChemicalTimes.H"

// ISAT
#if EDCSMOKE_USE_ISAT == 1
    #include "ISAT.h"
    #include "numericalJacobian4ISAT.H"
    #include "mappingGradients/mappingGradient4OpenFOAM.h"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#define NO_CONTROL
	#define CREATE_MESH createMeshesPostProcess.H
	#include "postProcess.H"

	#include "setRootCaseLists.H"
	#include "createTime.H"

	#include "createMeshes.H"
	fvMesh& mesh = fluidRegions[0];

	pimpleMultiRegionControl pimples(fluidRegions, solidRegions);
	
	#include "createFields.H"
        #include "createSolidFields.H"

	#include "readGravitationalAcceleration.H"
	#include "createOpenSMOKEFields.H"
	#include "createRadiationModel.H"
	#include "initContinuityErrs.H"

	// Complete pressure controls (the fluid region is unique)
	pimpleNoLoopControl& pimple = pimples.pimple(0);
	pressureReference pressureReference(p, pimple.dict(), false);
	scalar cumulativeContErr = 0.;

	// This solver does not support moving mesh but it uses the pressure
	// equation of one which does, so we need a dummy face-momentum field
	autoPtr<surfaceVectorField> rhoUf(nullptr);


	#include "createTimeControls.H"
	#include "readSolidTimeControls.H"
	#include "compressibleMultiRegionCourantNo.H"
	#include "solidRegionDiffusionNo.H"
	#include "setInitialMultiRegionDeltaT.H"

	unsigned int runTimeStep = 0;

	while (pimples.run(runTime))
	{
		#include "readTimeControls.H"
		#include "readSolidTimeControls.H"

		#include "compressibleMultiRegionCourantNo.H"
		#include "solidRegionDiffusionNo.H"
		#include "setMultiRegionDeltaT.H"

		runTime++;
		runTimeStep++;

		Info<< "Time = " << runTime.userTimeName() << nl << endl;

		// Optional number of energy correctors
		const int nEcorr = pimples.dict().lookupOrDefault<int>
		(
			"nEcorrectors",
			1
		);

		// --- PIMPLE loop
		while (pimples.loop())
		{
			List<tmp<fvVectorMatrix>> UEqns(fluidRegions.size());

			for(int Ecorr=0; Ecorr<nEcorr; Ecorr++)
			{
				forAll(solidRegions, i)
				{
					Info << "\nSolving for solid region " << solidRegions[i].name() << endl;
					#include "setRegionSolidFields.H"
					#include "solveSolid.H"
				}

				forAll(fluidRegions, i)
				{
					Info << "\nSolving for fluid region " << fluidRegions[i].name() << endl;
					#include "solveFluid.H"
				}
			}
		}

		runTime.write();

		Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;
	}

	Info<< "End\n" << endl;

	return 0;
}


/*
int main(int argc, char *argv[])
{
    unsigned int runTimeStep = 0;

        #define NO_CONTROL
        #define CREATE_MESH createMeshesPostProcess.H
        #include "postProcess.H"

        #include "setRootCaseLists.H"
        #include "createTime.H"
        #include "createMeshes.H"
	pimpleMultiRegionControl pimples(fluidRegions, solidRegions);
	#include "readGravitationalAcceleration.H"
	#include "createDyMControls.H"
	#include "initContinuityErrs.H"
	#include "createFields.H"
	#include "createOpenSMOKEFields.H"
	#include "createRhoUfIfPresent.H"

	#include "createRadiationModel.H"

	turbulence->validate();

	if (!LTS)
	{
		#include "compressibleCourantNo.H"
		#include "setInitialDeltaT.H"
	}

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );
        }

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        fvModels.preUpdateMesh();

        // Store momentum to set rhoUf for introduced faces.
        autoPtr<volVectorField> rhoU;
        if (rhoUf.valid())
        {
            rhoU = new volVectorField("rhoU", rho*U);
        }

        // Update the mesh for topology change, mesh to mesh mapping
        mesh.update();


        runTime++;
	runTimeStep++;
        Info<< "Time = " << runTime.timeName() << nl << endl;


        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (!pimple.flow())
            {
                if (pimple.models())
                {
                    fvModels.correct();
                }

                if (pimple.thermophysics())
                {
		    #include "properties.H"
		    #include "YEqn.H"
		    #include "EEqn.H"
                }
            }
            else
            {
                if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
                {
                    // Move the mesh
                    mesh.move();

                    if (mesh.changing())
                    {
                        MRF.update();

                        if (correctPhi)
                        {
                            #include "correctPhi.H"
                        }

                        if (checkMeshCourantNo)
                        {
                            #include "meshCourantNo.H"
                        }
                    }
                }

                if (pimple.firstPimpleIter() && !pimple.simpleRho())
                {
                    #include "rhoEqn.H"
                }

                if (pimple.models())
                {
                    fvModels.correct();
                }

                #include "UEqn.H"

                if (pimple.thermophysics())
                {
                    #include "properties.H"
		    #include "YEqn.H"
		    #include "EEqn.H"
                }

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }
        }

        rho = thermo.rho();

        runTime.write();

	Pav << runTime.timeName() << "\t" << p.weightedAverage(mesh.V()).value() << endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
*/
// ************************************************************************* //
