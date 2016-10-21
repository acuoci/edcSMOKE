// Patch temperature field for ignition
{
	if( 	(spark == true) && (runTime.value() >= ignitionTime)&&(runTime.value() <= (ignitionTime+ignitionDuration)) )
	{
		Info << "Forcing the spark..." << endl;

		scalarField& TCells = thermo.T().internalField();
		for(int i = 0; i<ignCellsList.size(); i++)
		{
   			 TCells[ignCellsList[i]] = ignitionTemperature;

			if (combustionModel == "EDC")
			{
				TStar_[0].internalField()[ignCellsList[i]] = ignitionTemperature;
			}
			else
			{
				for(unsigned int j=0;j<Y.size();j++)
					Y[j].internalField()[ignCellsList[i]] += ignitionAdditionalMassFraction;
				
				const double sum = 1.+Y.size()*ignitionAdditionalMassFraction;
				for(unsigned int j=0;j<Y.size();j++)
					Y[j].internalField()[ignCellsList[i]] /= sum;
			}
		}
	}
}
