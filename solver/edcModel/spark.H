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
