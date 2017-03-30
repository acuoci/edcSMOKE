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

// Patch temperature field for ignition
#if OPENFOAM_VERSION == 4
{
	if( 	(spark == true) && (runTime.value() >= ignitionTime)&&(runTime.value() <= (ignitionTime+ignitionDuration)) )
	{
		Info << "Forcing the spark..." << endl;

		scalarField& TCells = thermo.T().ref();
		for(int i = 0; i<ignCellsList.size(); i++)
		{
   			 TCells[ignCellsList[i]] = ignitionTemperature;

			if (combustionModel == "EDC")
			{
				TStar_[0].ref()[ignCellsList[i]] = ignitionTemperature;
			}
			else
			{
				for(unsigned int j=0;j<Y.size();j++)
					Y[j].ref()[ignCellsList[i]] += ignitionAdditionalMassFraction;
				
				const double sum = 1.+Y.size()*ignitionAdditionalMassFraction;
				for(unsigned int j=0;j<Y.size();j++)
					Y[j].ref()[ignCellsList[i]] /= sum;
			}
		}
	}
}
#else
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
#endif
