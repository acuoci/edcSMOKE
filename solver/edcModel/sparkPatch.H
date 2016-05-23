// Patch temperature field for ignition
{
	if( 	(spark == true) && (runTime.value() >= ignitionTime)&&(runTime.value() <= (ignitionTime+ignitionDuration)) )
	{
		Info << "Forcing the spark..." << endl;
		scalarField& TCells = thermo.T().internalField();
		for(int i = 0; i<ignCellsList.size(); i++)
		{
   			 TCells[ignCellsList[i]] = ignitionTemperature;
   			 TStar_[0].internalField()[ignCellsList[i]] = ignitionTemperature;
		}
	}
}
