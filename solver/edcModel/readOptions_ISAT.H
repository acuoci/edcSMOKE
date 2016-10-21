/*-----------------------------------------------------------------------*\
|                  _       _____ __  __  ____  _  ________                |
|                 | |     / ____|  \/  |/ __ \| |/ /  ____|               |
|          ___  __| | ___| (___ | \  / | |  | | ' /| |__                  |
|         / _ \/ _` |/ __|\___ \| |\/| | |  | |  < |  __|                 |
|        |  __/ (_| | (__ ____) | |  | | |__| | . \| |____                |
|         \___|\__,_|\___|_____/|_|  |_|\____/|_|\_\______|               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, M.R. Malik, A. Parente                             |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|   Contacts: Mohammad Rafi Malik                                         |
|   Aero-Thermo-Mechanical Department                                     |
|   Université Libre de Bruxelles                                         |
|   Avenue F. D. Roosevelt 50, 1050 Bruxelles (Belgium)                   |
|                                                                         |
|   Contacts: Alessandro Parente                                          |
|   email: alessandro.parente@ulb.ac.be                                   |
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
|   Copyright(C) 2016-2014 A. Cuoci, M.R. Malik, A. Parente               |
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

//- Reading ISAT parameters
Switch isatCheck(combustionOptionsDictionary.subDict("ISAT").lookup("ISAT"));
ISAT *isatTable;
Eigen::VectorXd scalingFactors_ISAT;
int luSolver_ISAT = 1;
label numberSubSteps_ISAT = 1;

if (isatCheck == true && drg_analysis == true)
{
	FatalError << "ISAT and DRG cannot be used together! Please turn off one between ISAT and DRG" << endl; 		
}

if (isatCheck == true)
{
	// mandatory if ISAT ON
	scalar epsilon_ISAT = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("tolerance", 1e-4);
	       numberSubSteps_ISAT = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("numberSubSteps", 1);

	scalar maxSizeMRU = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("maxSizeMRU", 100);
	scalar maxSizeMFU = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("maxSizeMFU", 100);
	scalar maxSizeBT = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("maxSizeBT", 100000);
	scalar maxSearchMRU = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("maxSearchMRU", 10);
	scalar maxSearchMFU = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("maxSearchMFU", 30);
	Switch searchMRU(combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<Switch>("searchMRU", true)); 
	Switch searchMFU(combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<Switch>("searchMFU", true));
	Switch searchBF(combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<Switch>("searchBruteForce", false));
	Switch clearIfFull(combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<Switch>("clearingIfFull", false));
	Switch cleanAndBalance(combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<Switch>("cleanAndBalance", true));
	scalar maxTimeOldCoeff = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("maxTimeOldCoeff", 0.7);
	scalar maxGrowCoeff = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("maxGrowCoeff", 0.5);
	scalar maxHeightCoeff = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("maxHeightCoeff", 20.);
	scalar minUsedCoeff = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("minUsedCoeff", 0.01);
	scalar balanceFactorRetrieve = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("balanceFactorRetrieve", 2.);
	scalar balanceFactorAddition = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("balanceFactorAddition", 0.1);
	word   luFactorization = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<word>("luFactorization","Partial");
	word   qrFactorization = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<word>("qrFactorization","Full");

	if (luFactorization != "Partial" && luFactorization != "Full")
	{
		Info << "Wrong luFactorization options: Partial || Full" << endl;
		abort();
	} 
	else 
	{
		if (luFactorization == "Full") 
			luSolver_ISAT = 0;
		else
			luSolver_ISAT = 1;		
	}

	int qrSolver_ISAT = 0;
	if (qrFactorization != "NoPivoting" && qrFactorization != "Partial" && qrFactorization != "Full")
	{
		Info << "Wrong qrFactorization options: NoPivoting || Partial || Full" << endl;
		abort();
	} 
	else 
	{
		if (qrFactorization == "Full") 
			qrSolver_ISAT = 0;
		else if (qrFactorization == "Partial") 
			qrSolver_ISAT = 1;		
		else if (qrFactorization == "NoPivoting")
			qrSolver_ISAT = 2;
	}


	//- ISAT scale factor 
	scalingFactors_ISAT.setConstant(ne_total,1.);
	if(isatCheck == true) 
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();

		//massFractionsTol = epsilon_ISAT;
	
		dictionary scalingFactors(combustionOptionsDictionary.subDict("ISAT").subDict("scalingFactors"));
		for(unsigned int i=0;i<NC;i++) 
		{	
			if(scalingFactors.found(thermodynamicsMapXML->NamesOfSpecies()[i])) 
				scalingFactors_ISAT(i) = 1./readScalar(scalingFactors.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
			else 
				scalingFactors_ISAT(i) = 1./readScalar(scalingFactors.lookup("others"));
		}
		scalingFactors_ISAT(NC)   = 1./readScalar(scalingFactors.lookup("T"));
		scalingFactors_ISAT(NC+1) = 1./readScalar(scalingFactors.lookup("tf"));
		scalingFactors_ISAT(NC+2) = 1./readScalar(scalingFactors.lookup("gamma"));
	}
	//- ISAT scaling error factors
	Eigen::VectorXd scalingErrors_ISAT;
	scalingErrors_ISAT.setConstant(ne_total,1.);
	if(isatCheck == true) 
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
	
		dictionary scalingErrors(combustionOptionsDictionary.subDict("ISAT").subDict("scalingErrors"));
	
		for(unsigned int i=0;i<NC;i++) 
		{
			if(scalingErrors.found(thermodynamicsMapXML->NamesOfSpecies()[i])) 
				scalingErrors_ISAT(i) = readScalar(scalingErrors.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
			else 
				scalingErrors_ISAT(i) = readScalar(scalingErrors.lookup("others"));
		}
		scalingErrors_ISAT(NC)   = readScalar(scalingErrors.lookup("T"));
		scalingErrors_ISAT(NC+1) = readScalar(scalingErrors.lookup("tf"));
		scalingErrors_ISAT(NC+2) = readScalar(scalingErrors.lookup("gamma"));
	}

	//- ISAT HOM
	isatTable = new ISAT(scalingErrors_ISAT, epsilon_ISAT, ne_total);	

	// - Setting ISAT_HOM param
	isatTable->setMaxSizeBT(maxSizeBT);
	isatTable->setMaxSizeMRU(maxSizeMRU);
	isatTable->setMaxSizeMFU(maxSizeMFU);
	isatTable->setMaxSearchMRU(maxSearchMRU);
	isatTable->setMaxSearchMFU(maxSearchMFU);
	isatTable->setFlagSearchMRU(searchMRU);
	isatTable->setFlagSearchMFU(searchMFU);
	isatTable->setFlagClearingIfFull(clearIfFull);
	isatTable->setMaxGrowCoeff(maxGrowCoeff);
	isatTable->setMaxHeightCoeff(maxHeightCoeff);
	isatTable->setMaxTimeOldCoeff(maxTimeOldCoeff);
	isatTable->setMinUsedCoeff(minUsedCoeff);
	isatTable->setBalanceFactorRet(balanceFactorRetrieve);
	isatTable->setBalanceFactorAdd(balanceFactorAddition);
	isatTable->setQRType(qrSolver_ISAT);
	isatTable->setFlagCleanAndBalance(cleanAndBalance);

	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
	
		Info << endl << "ISAT parameters " << endl;
		Info << "   tolerance           : " << epsilon_ISAT << endl;
		Info << "   luFactorization     : " << luFactorization << endl; 	
		Info << "   qrFactorization     : " << qrFactorization << endl; 	
	
		Info << "   scalingFactors      : " << endl;
		for(unsigned int i=0;i<NC;i++)  
		{
			if(scalingFactors_ISAT(i) != 1.) 
				Info << "      " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t       : " << scalingFactors_ISAT(i) << endl;
		}
		Info << "       others          : " << readScalar(combustionOptionsDictionary.subDict("ISAT").subDict("scalingFactors").lookup("others")) << endl;
		Info << "       T               : " << scalingFactors_ISAT(NC) << endl;
		Info << "       tf              : " << scalingFactors_ISAT(NC+1) << endl;
		Info << "       gamma           : " << scalingFactors_ISAT(NC+2) << endl;


		Info << "   clear if BT full    : " << clearIfFull << endl;

		Info << "   search in MRU       : " << searchMRU << endl;
		Info << "   search in MFU       : " << searchMFU << endl;
		Info << "   search BF           : " << searchBF << endl;

		Info << "   dimension parameters: " << endl;
		Info << "      max size BT      : " << maxSizeBT << endl;	
		if(searchMRU == true) 
		{
			Info << "      max size MRU     : " << maxSizeMRU << endl;
			Info << "      max search MRU   : " << maxSearchMRU << endl;
		}
		if(searchMFU == true) 
		{
			Info << "      max size MFU     : " << maxSizeMRU << endl;
			Info << "      max search MFU   : " << maxSearchMFU << endl;  
		}
	
		if (cleanAndBalance == true)
		{
			Info << "   balance parameters  : " << endl;
			Info << "      balanceFactorRetrieve   : " << balanceFactorRetrieve << endl;	
			Info << "      balanceFactorAddition   : " << balanceFactorAddition << endl;	
			Info << "      maxHeightCoeff          : " << maxHeightCoeff << endl;	
			Info << "      maxGrowCoeff            : " << maxGrowCoeff << endl;
			Info << "      minUsedCoeff            : " << minUsedCoeff << endl;
			Info << "      maxTimeOldCoeff         : " << maxTimeOldCoeff << endl; 
		}

		Info << "   scaling error       : " << endl;
		for(unsigned int i=0;i<NC;i++)  
		{
			if(scalingErrors_ISAT(i) != 1.) 
				Info << "      " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t       : " << scalingErrors_ISAT(i) << endl;
		}
		Info << "       others          : " << readScalar(combustionOptionsDictionary.subDict("ISAT").subDict("scalingErrors").lookup("others")) << endl;
		Info << "       T               : " << scalingErrors_ISAT(NC) << endl;
		Info << "       tf              : " << scalingErrors_ISAT(NC+1) << endl;
		Info << "       gamma           : " << scalingErrors_ISAT(NC+2) << endl;
	
		Info << endl;
	}

	if (	odeParametersFineStructures.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE &&
		isatCheck == true)
	{
		Info << "In the current version of edcSMOKE the ISAT can be used only together with the native OpenSMOKE Integrator." << endl;
		abort();
	}
}



