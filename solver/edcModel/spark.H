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

// Spark
Switch spark = false;
labelList ignCellsList;
scalar ignitionTime = 0.;
scalar ignitionDuration = 0.;
scalar ignitionTemperature = 0.;
scalar ignitionAdditionalMassFraction = 1.e-3;
{
	dictionary sparkSubDictionary(combustionOptionsDictionary.subDict("Spark"));

	spark = readBool(sparkSubDictionary.lookup("spark"));

	Foam::vector position;    
	scalar ignitionDiameter = 0.;

	// Read spark properties if detected
	if(spark == true)
	{
		position = (sparkSubDictionary.lookup("position"));
		ignitionTime = readScalar(sparkSubDictionary.lookup("time"));
		ignitionDuration = readScalar(sparkSubDictionary.lookup("duration"));
		ignitionTemperature = readScalar(sparkSubDictionary.lookup("temperature"));
		ignitionDiameter = readScalar(sparkSubDictionary.lookup("diameter"));   
		ignitionAdditionalMassFraction = readScalar(sparkSubDictionary.lookup("additionalMassFraction"));
	}

	label nIgnCells = 1;
	if(spark == true)
	{
		const volVectorField& centres = mesh.C();

		label ignitionCell;
		ignitionCell = mesh.findCell(position);
		
		if (ignitionCell == -1)
		{
			Info << "OpenFOAM WARNING:" << endl;
			Info << "Theorical ignition position does not correspond to any cell." << endl;
			Info << "Please check position\nIn file constant/solverOptions.H\n"    << endl;
		}
		else 
		{
			Info << "Ignition cell center: " << centres[ignitionCell] << ". Original position: " << position <<endl;

			scalar radius = ignitionDiameter/2.0;
			ignCellsList.setSize(1);
			ignCellsList[0] = ignitionCell;

			forAll(centres, celli)
			{
				scalar dist = mag(centres[celli] - position);

				if (dist < radius && celli != ignitionCell)
				{
					ignCellsList.setSize(nIgnCells+1);

					ignCellsList[nIgnCells] = celli;

					nIgnCells++;
				}
			}

			if (ignCellsList.size() > 0) 
				Info << "Added " << nIgnCells << " cells to ignition point.\n" <<	endl;
			else 
				Info << "OpenFOAM wasn't able to find more than one cell in the ignition point.\n" << endl;
		}
	}
}
