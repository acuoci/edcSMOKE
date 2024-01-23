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

if (laminarDiffusion == true)
{
	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info<< "Properties evaluation... " ;

 	OpenSMOKE::OpenSMOKEVectorDouble Dmixvector(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble moleFractions(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble massFractions(thermodynamicsMapXML->NumberOfSpecies());
	

	// Internal fields
	{
		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField();
		scalarField& lambdaCells = lambda.ref();

		forAll(TCells, celli)
		{
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			
			transportMapXML->SetPressure(pCells[celli]);
			transportMapXML->SetTemperature(TCells[celli]);

			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].internalField()[celli];
				
			double MWmix;
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),MWmix,massFractions.GetHandle());
			
			lambdaCells[celli] = transportMapXML->ThermalConductivity(moleFractions.GetHandle());
			transportMapXML->MassDiffusionCoefficients(Dmixvector.GetHandle(), moleFractions.GetHandle());

			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				Dmix_[i].ref()[celli] = Dmixvector[i+1];
		}
	}

	// Boundaries
	forAll(T.boundaryField(), patchi)
	{
		fvPatchScalarField& pp = p.boundaryFieldRef()[patchi];
		fvPatchScalarField& pT = T.boundaryFieldRef()[patchi];
		fvPatchScalarField& plambda = lambda.boundaryFieldRef()[patchi];
        
		forAll(pT, facei)
		{
			thermodynamicsMapXML->SetPressure(pp[facei]);
			thermodynamicsMapXML->SetTemperature(pT[facei]);

			transportMapXML->SetPressure(pp[facei]);
			transportMapXML->SetTemperature(pT[facei]);
			
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].boundaryField()[patchi][facei];

			double MWmix;
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),MWmix,massFractions.GetHandle());

			plambda[facei] = transportMapXML->ThermalConductivity(moleFractions.GetHandle());
			transportMapXML->MassDiffusionCoefficients(Dmixvector.GetHandle(),moleFractions.GetHandle());

			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				Dmix_[i].boundaryFieldRef()[patchi][facei] = Dmixvector[i+1];
		}
	}
	
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}
