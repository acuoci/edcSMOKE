//- Reading ISAT parameters
Switch isatCheck(combustionOptionsDictionary.subDict("ISAT").lookup("ISAT"));


// mandatory if ISAT ON
scalar epsilon_ISAT = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("tolerance", 1e-4);
label  numberSubSteps_ISAT = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("numberSubSteps", 1);

scalar maxSizeMRU = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("maxSizeMRU", 100);
scalar maxSizeMFU = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("maxSizeMFU", 100);
scalar maxSizeBT = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("maxSizeBT", 100000);
scalar maxSearchMRU = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("maxSearchMRU", 10);
scalar maxSearchMFU = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<int>("maxSearchMFU", 30);
Switch searchMRU(combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<Switch>("searchMRU", true)); 
Switch searchMFU(combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<Switch>("searchMFU", true));
Switch searchBF(combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<Switch>("searchBruteForce", false));
Switch clearIfFull(combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<Switch>("clearingIfFull", false));
scalar maxTimeOldCoeff = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("maxTimeOldCoeff", 0.7);
scalar maxGrowCoeff = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("maxGrowCoeff", 0.5);
scalar maxHeightCoeff = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("maxHeightCoeff", 20.);
scalar maxUsedCoeff = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<double>("maxUsedCoeff", 20.);
Switch printStats(combustionOptionsDictionary.subDict("ISAT").lookupOrDefault("statistics", false));
word   luFactorization = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<word>("luFactorization","Partial");
word   qrFactorization = combustionOptionsDictionary.subDict("ISAT").lookupOrDefault<word>("qrFactorization","Full");

int luSolver_ISAT = 1;
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
Eigen::VectorXd scalingFactors_ISAT;
scalingFactors_ISAT.setConstant(thermodynamicsMapXML->NumberOfSpecies()+2,1.);
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
}
//- ISAT scaling error factors
Eigen::VectorXd scalingErrors_ISAT;
scalingErrors_ISAT.setConstant(thermodynamicsMapXML->NumberOfSpecies()+2,1.);
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
}

//- ISAT HOM
ISAT *isat_HOM = new ISAT(scalingErrors_ISAT, epsilon_ISAT, thermodynamicsMapXML->NumberOfSpecies()+2);	

// - Setting ISAT_HOM param
isat_HOM->setMaxSizeBT(maxSizeBT);
isat_HOM->setMaxSizeMRU(maxSizeMRU);
isat_HOM->setMaxSizeMFU(maxSizeMFU);
isat_HOM->setMaxSearchMRU(maxSearchMRU);
isat_HOM->setMaxSearchMFU(maxSearchMFU);
isat_HOM->setFlagSearchMRU(searchMRU);
isat_HOM->setFlagSearchMFU(searchMFU);
isat_HOM->setFlagClearingIfFull(clearIfFull);
isat_HOM->setMaxGrowCoeff(maxGrowCoeff);
isat_HOM->setMaxHeightCoeff(maxHeightCoeff);
isat_HOM->setMaxTimeOldCoeff(maxTimeOldCoeff);
isat_HOM->setMaxUsedCoeff(maxUsedCoeff);

if(isatCheck == true) 
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
	
	Info << "   balance parameters  : " << endl;
	Info << "      maxHeightCoeff   : " << maxHeightCoeff << endl;	
	Info << "      maxGrowCoeff     : " << maxGrowCoeff << endl;
	Info << "      maxUsedCoeff     : " << maxUsedCoeff << endl;
	Info << "      maxTimeOldCoeff  : " << maxTimeOldCoeff << endl; 

	Info << "   scaling error       : " << endl;
	for(unsigned int i=0;i<NC;i++)  
	{
		if(scalingErrors_ISAT(i) != 1.) 
			Info << "      " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t       : " << scalingErrors_ISAT(i) << endl;
	}
	Info << "       others          : " << readScalar(combustionOptionsDictionary.subDict("ISAT").subDict("scalingErrors").lookup("others")) << endl;
	Info << "       T               : " << scalingErrors_ISAT(NC) << endl;
	
	Info << endl;
	Info << "   statistics          : " << printStats << endl << endl;
}

if (	odeParametersFineStructures.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE &&
	isatCheck == true)
{
	Info << "In the current version of edcSMOKE the ISAT can be used only together with the native OpenSMOKE Integrator." << endl;
	abort();
}

if (	reactorType == REACTORTYPE_PSR &&
	isatCheck == true)
{
	Info << "In the current version of edcSMOKE the ISAT can be used only together with the PFR model (fine structures)." << endl;
	abort();
}




