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

// This is a steady state solver
#define STEADYSTATE 1

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
#include "simpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "radiationModel.H"
#include "ChemistryLinearModel.H"

// Utilities
#include "Utilities.H"

// ODE system
#include "DRG.H"
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

// Solid
#include "compressibleCourantNo.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    unsigned int runTimeStep = 0;

	#include "postProcess.H"

	#include "setRootCaseLists.H"
	#include "createTime.H"

	// Create multiple meshes
	#include "createMeshes.H"
	fvMesh& mesh = fluidRegions[0];

    #include "readGravitationalAcceleration.H"
//	#include "createControl.H"
    PtrList<simpleControl> simples(fluidRegions.size());
    for(int i=0;i<fluidRegions.size();i++)
        simples.set(i, new simpleControl(fluidRegions[i]));

    simpleControl& simple = simples[0];
	#include "createFields.H"
	#include "createNonReactingFields.H"
    #include "createSolidFields.H"

    #include "createOpenSMOKEFields.H"
	#include "createRadiationModel.H"
	#include "initContinuityErrs.H"

	turbulence->validate();


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simples[0].loop(runTime))
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

		// Check for regions to be solved
		const bool is_fluid_active = runTime.controlDict().lookupOrDefault<Switch>("fluid", true);
		const bool is_solid_active = runTime.controlDict().lookupOrDefault<Switch>("solid", true);

        fvModels.correct();

        // Solid region equations
		if (is_solid_active == true)
		{
			// Solve solid regions
			forAll(solidRegions, i)
			{
				Info << "\nSolving for solid region " << solidRegions[i].name() << endl;

				#include "setRegionSolidFields.H"
				#include "solveSolid.H"
			}
		}

        // Fluid region equations
		if (is_fluid_active == true)
		{
            // Solving for reacting region
            {
                Info << "\nSolving for reacting fluid region " << fluidRegions[0].name() << endl;

                simpleControl& simple = simples[0];

                scalar cumulativeContErr = cumulativeContErrs[0];

	            if (momentumEquations == true)
	            {
		            // Pressure-velocity SIMPLE corrector
		            {
		                #include "UEqn.H"
		                #include "properties.H"
		                #include "YEqn.H"
		                #include "EEqn.H"
		                #include "pEqn.H"
		            }
	            }
	            else
	            {
		            #include "properties.H"
		            #include "YEqn.H"
		            #include "EEqn.H"
	            }

                turbulence->correct();
            }

        	// Non reacting regions
			for(int i=1;i<fluidRegions.size();i++)
			{
				Info << "\nSolving for non-reacting fluid region " << fluidRegions[i].name() << endl;

				#include "setRegionFluidFields.H"

                // Set simple control
                simpleControl& simple = simples[i];
                simple.read();
                simple.storePrevIterFields();

                tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);

                if (momentumEquations == true)
	            {
		            // Pressure-velocity SIMPLE corrector
		            {
		                #include "UEqnNonReacting.H"
		                #include "properties.H"
		                #include "YEqnNonReacting.H"
		                #include "EEqnNonReacting.H"
		                #include "pEqn.H"
		            }
	            }
	            else
	            {
		            #include "properties.H"
		            #include "YEqnNonReacting.H"
		            #include "EEqnNonReacting.H"
	            }

                turbulence.correct();	
                thermophysicalTransport.correct();		
			}
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
