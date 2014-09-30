// Read the kinetic scheme in XML format
OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>* thermodynamicsMapXML; 
OpenSMOKE::KineticsMap_CHEMKIN<double>* kineticsMapXML;

{	
	word kinetics(thermo.lookup("kineticMechanism"));
	boost::filesystem::path path_kinetics = kinetics;

	rapidxml::xml_document<> doc;
	std::vector<char> xml_string;
	OpenSMOKE::OpenInputFileXML(doc, xml_string,path_kinetics / "kinetics.xml");

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>(doc); 
	//transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc); 
	kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN<double>(*thermodynamicsMapXML, doc); 					
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
}

const word combustionOptionsDictionaryName("combustionOptions");

Info<< "Reading CombustionOptions dictionary\n" << endl;
IOdictionary combustionOptionsDictionary
(
	IOobject
	(
		combustionOptionsDictionaryName,
		U.time().constant(),
		U.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

// Selecting combustion model
Switch homogeneousReactions(combustionOptionsDictionary.lookup("homogeneousReactions"));

// RAS or LES
Switch iLES = false;
{
	Info<< "Reading turbulenceProperties dictionary\n" << endl;
	IOdictionary turbulencePropertiesDictionary
	(
		IOobject
		(
			"turbulenceProperties",
			U.time().constant(),
			U.db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);


	word simulationType(turbulencePropertiesDictionary.lookup("simulationType"));
	if (simulationType == "RASModel")
	{
		iLES = false;
		Info << "This is a RANS simulation..." << endl;
	}
	else if (simulationType == "LESModel")
	{
		iLES = true;
		Info << "This is a Large Eddy Simulation (LES)..." << endl;
	}
	else
	{
		FatalError << "The ED models require turbulent conditions (i.e. LESModel or RASModel)" << exit(FatalError);
	}
}

// Check if epsilon can be found (only for RAS)
if (iLES == false)
{
	if ( turbulence->db().objectRegistry::foundObject<volScalarField>("epsilon") )
	{
		Info << "Found epsilon from turbulence model!" << endl;
	}
	else
	{
		FatalError << "No epsilon from turbulence model!" << exit(FatalError);
	}
}

// Models
word combustionModel(combustionOptionsDictionary.lookup("combustionModel"));
const scalar maxT = readScalar(combustionOptionsDictionary.lookup("maxT"));
const scalar minT = readScalar(combustionOptionsDictionary.lookup("minT"));

// Numerical parameters for EDC
// scalar CD1 = 0.134;
// scalar CD2 = 0.5;
scalar Ccsi = 2.1377; 	// Ccsi = pow(3.*CD2/4/CD1/CD1, 0.25);
scalar Ctau = 0.40825;	// Ctau = sqrt(CD2/3.);

scalar maxGammaSquared = 0.75487766248;
scalar numberOfResidenceTimes = 100.;
label  numberIterationsToUpdateChemistry = 1;
label  edcCounter = 0;

// Numerical parameters for ED
scalar constA = 4.;
scalar constB = 0.5;
std::vector< std::vector<unsigned int> > ed_reactant_index_;
std::vector< std::vector<unsigned int> > ed_product_index_;
std::vector< std::vector<double> > ed_reactant_nu_;
std::vector< std::vector<double> > ed_product_nu_;
std::vector< std::vector<double> > ed_denominator_reactants_;
std::vector<double> ed_denominator_products_;

// Reading constants from file
{
// ED Constants
if (combustionModel == "ED")
{
	dictionary edSubDictionary(combustionOptionsDictionary.subDict("ED"));

	constA = readScalar(edSubDictionary.lookup("A"));	// default: 4
	constB = readScalar(edSubDictionary.lookup("B"));	// default: 0.5
}
// EDFR Constants
else if (combustionModel == "EDFR")
{
	dictionary edfrSubDictionary(combustionOptionsDictionary.subDict("EDFR"));

	constA = readScalar(edfrSubDictionary.lookup("A"));	// default: 4
	constB = readScalar(edfrSubDictionary.lookup("B"));	// default: 0.5
}
// EDC Constants
else if (combustionModel == "EDC")
{
	dictionary edcSubDictionary(combustionOptionsDictionary.subDict("EDC"));

	Ccsi = readScalar(edcSubDictionary.lookup("Ccsi"));						// default: 2.1377
	Ctau = readScalar(edcSubDictionary.lookup("Ctau"));						// default: 0.40825
	maxGammaSquared = readScalar(edcSubDictionary.lookup("maxGammaSquared"));			// default: 0.75487766248
	numberOfResidenceTimes = readScalar(edcSubDictionary.lookup("numberOfResidenceTimes"));		// default: 100.
	numberIterationsToUpdateChemistry = readLabel(edcSubDictionary.lookup("numberIterationsToUpdateChemistry")); // default 1
	edcCounter = numberIterationsToUpdateChemistry;
}
// Error Message
else
{
	FatalError << "Combustion model available: ED || EDFR || EDC" << endl;
}


if (combustionModel == "ED" || combustionModel == "EDFR")
{
	if (kineticsMapXML->NumberOfReversibleReactions() != 0)
		FatalError << "The ED and EDFR models cannot be used with reversible reactions" << endl;

	ed_denominator_reactants_.resize(kineticsMapXML->NumberOfReactions());
	ed_denominator_products_.resize(kineticsMapXML->NumberOfReactions());
	ed_reactant_index_.resize(kineticsMapXML->NumberOfReactions());
	ed_product_index_.resize(kineticsMapXML->NumberOfReactions());
	ed_reactant_nu_.resize(kineticsMapXML->NumberOfReactions());
	ed_product_nu_.resize(kineticsMapXML->NumberOfReactions());

	for(unsigned int k=0;k<kineticsMapXML->NumberOfReactions();k++)
		ed_denominator_products_[k] = 0.;

	for (int k=0; k<kineticsMapXML->stoichiometry().stoichiometric_matrix_reactants().outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML->stoichiometry().stoichiometric_matrix_reactants(),k); it; ++it)
		{
			ed_reactant_index_[it.row()].push_back(it.col());
			ed_reactant_nu_[it.row()].push_back(it.value());
			ed_denominator_reactants_[it.row()].push_back( it.value()*thermodynamicsMapXML->MW()[it.col()+1] );
		}
	}
	
	for (int k=0; k<kineticsMapXML->stoichiometry().stoichiometric_matrix_products().outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML->stoichiometry().stoichiometric_matrix_products(),k); it; ++it)
		{
			ed_product_index_[it.row()].push_back(it.col());
			ed_product_nu_[it.row()].push_back(it.value());
			ed_denominator_products_[it.row()]+=it.value()*thermodynamicsMapXML->MW()[it.col()+1];
		};
	}
}
}

// Ode parameters
OpenSMOKE::ODE_Parameters odeParametersFineStructures;
if (combustionModel == "EDC")
{
	const dictionary& odeParametersFineStructuresDictionary = combustionOptionsDictionary.subDict("Ode");

	//- Mass fractions tolerance
	scalar relTolerance = readScalar(odeParametersFineStructuresDictionary.lookup("relTolerance"));
	scalar absTolerance = readScalar(odeParametersFineStructuresDictionary.lookup("absTolerance"));
	odeParametersFineStructures.SetRelativeTolerance(relTolerance);
	odeParametersFineStructures.SetAbsoluteTolerance(absTolerance);
	
	//- Full pivoting (only for OpenSMOKE solver)
	Switch fullPivoting(odeParametersFineStructuresDictionary.lookup("fullPivoting"));
	odeParametersFineStructures.SetFullPivoting(fullPivoting);
	
	//- Maximum order of integration (only for OpenSMOKE solver)
	label maximumOrder = readLabel(odeParametersFineStructuresDictionary.lookup("maximumOrder"));
	odeParametersFineStructures.SetMaximumOrder(maximumOrder);
	
	// Type
	word homogeneousODESolverString(odeParametersFineStructuresDictionary.lookup("odeSolver"));
	if (	homogeneousODESolverString != "OpenSMOKE" 	&& homogeneousODESolverString != "DVODE"  && 
		homogeneousODESolverString != "DLSODE" 		&& homogeneousODESolverString != "DLSODA" && 
		homogeneousODESolverString != "CVODE" 		&& homogeneousODESolverString != "DASPK"  &&
		homogeneousODESolverString != "MEBDF" 		&& homogeneousODESolverString != "RADAU5"  
	   )
	{
		Info << "Wrong homogeneous ODE Solver: OpenSMOKE || DVODE || DLSODE || DLSODA || CVODE || DASPK || MEBDF || RADAU5" << endl;
		abort();
	}

	if (homogeneousODESolverString == "OpenSMOKE")	odeParametersFineStructures.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE);
	if (homogeneousODESolverString == "DVODE") 	odeParametersFineStructures.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE);
	if (homogeneousODESolverString == "DLSODE") 	odeParametersFineStructures.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE);
	if (homogeneousODESolverString == "DLSODA") 	odeParametersFineStructures.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA);
	if (homogeneousODESolverString == "CVODE") 	odeParametersFineStructures.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE);
	if (homogeneousODESolverString == "DASPK") 	odeParametersFineStructures.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK);
	if (homogeneousODESolverString == "MEBDF") 	odeParametersFineStructures.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF);
	if (homogeneousODESolverString == "RADAU5") 	odeParametersFineStructures.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5);	
}

// Check type
{
	if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
	{
		#if OPENSMOKE_USE_DVODE != 1
		{
			Info << "The solver was compiled without the DVODE support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
	{
		#if OPENSMOKE_USE_SUNDIALS != 1
		{
			Info << "The solver was compiled without the CVODE support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if ( odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE ||
	          odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA )
	{
		#if OPENSMOKE_USE_ODEPACK != 1
		{
			Info << "The solver was compiled without the ODEPACK support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
	{
		#if OPENSMOKE_USE_RADAU != 1
		{
			Info << "The solver was compiled without the RADAU support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}	
	else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
	{
		#if OPENSMOKE_USE_DASPK != 1
		{
			Info << "The solver was compiled without the DASPK support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
	{
		#if OPENSMOKE_USE_MEBDF != 1
		{
			Info << "The solver was compiled without the MEBDF support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}			
}
		
//- Species mass fractions
PtrList<volScalarField> YStar_;
PtrList<volScalarField> TStar_;
PtrList<volScalarField> R_;
scalarField Tau(mesh.nCells(),0.);
scalarField gammaStar(mesh.nCells(),0.);
scalarField mDotStar(mesh.nCells(),0.);
scalarField edc_coefficient(mesh.nCells(),0.);
{
	// Check the order of species
	if (thermodynamicsMapXML->NumberOfSpecies() != Y.size())
		FatalError << "Kinetic schemes do not match..." << endl;

	for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		if (Y[i].name() != thermodynamicsMapXML->NamesOfSpecies()[i])
		{
			FatalError << "Kinetic schemes do not match..." << endl;
		}

	YStar_.resize(Y.size());
	TStar_.resize(1);
	R_.resize(Y.size());

	unsigned int species_from_file = 0;
	unsigned int species_from_default = 0;

	forAll(Y,i)
	{
		std::cout << " * Processing species "<< i << " : " << Y[i].name() << " : ";
		IOobject header
		(
			"edc::YStar_" + Y[i].name(),
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		// check if field exists and can be read
		if (combustionModel == "EDC")
		{
			if (header.headerOk())
			{
				std::cout << "from file..." << std::endl;
				YStar_.set
				(
					i,
					new volScalarField
					(
						IOobject
						(
							"edc::YStar_"+ Y[i].name(),
							mesh.time().timeName(),
							mesh,
							IOobject::MUST_READ,
							IOobject::AUTO_WRITE
						),
						mesh
					)
				);
				species_from_file++;
			}
			else
			{
				std::cout << "from mean values..." << std::endl;
				YStar_.set
				(
					i,
					new volScalarField
					(
						IOobject
						(
							"edc::YStar_"+ Y[i].name(),
							mesh.time().timeName(),
							mesh,
							IOobject::NO_READ,
							IOobject::AUTO_WRITE
						),
						mesh,
						dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.)
					)
				);
			
				YStar_[i] = Y[i];
				species_from_default++;
			}
		}
		
		R_.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"edc::R_"+ Y[i].name(),
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.)
			)
		);	
	}
}

// TStar
if (combustionModel == "EDC")
{
	std::cout << " * Reading TStar field..." << std::endl;
	{
		IOobject header
		(
			"edc::TStar",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		// check if field exists and can be read
		if (header.headerOk())
		{
			std::cout << "from file..." << std::endl;
			TStar_.set
			(
				0,
				new volScalarField
				(
					IOobject
					(
						"edc::TStar",
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
			std::cout << "from mean values..." << std::endl;
			TStar_.set
			(
				0,
				new volScalarField
				(
					IOobject
					(
						"edc::TStar",
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					mesh,
					dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.)
				)
			);
			
			TStar_[0] = thermo.T();
		}
	}
}

volScalarField Sh_
(
    IOobject
    (
        "Sh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Sh", dimEnergy/dimTime/dimVolume, 0.0)
);

dimensionedScalar ksmall_("ksmall", dimensionSet(0,2,-2,0,0,0,0), scalar(1e-6) );

volScalarField uTau_
(
    IOobject
    (
        "uTau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("uTau", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

volScalarField epsilonSGS
(
    IOobject
    (
        "epsilonSGS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilonSGS", dimensionSet(0,2,-3,0,0,0,0), 0.0)
);

// Normalize the mass fractions (mean fields)
{
	Info << "Normalizing mass fractions (on entry...)" << endl;

	const scalarField& TCells = thermo.T().internalField();
	
	const double maxSumMassFractionsError = 1e-3;
	double sumTotal = 0.;
	double maxError = 0.;
	forAll(TCells, celli)
	{	
		double sum = 0.;

		for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			sum += Y[i].internalField()[celli];

		double e = fabs(sum-1.);

		if ( e > maxSumMassFractionsError)
			FatalError << "Fatal error: the sum of mass fractions in cell " << celli << " is equal to " << sum << " (error " << sum-1. << endl; 
			
		for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			Y[i].internalField()[celli]/=sum;

		sumTotal += sum;
		if (e>maxError) maxError = e;
	}

	Info << "Mean error (mean field): " << fabs( sumTotal/double(TCells.size())-1.) << endl;
	Info << "Max error (mean field):  " << maxError << endl;
}

// Normalize the mass fractions (fine structures)
if (combustionModel == "EDC")
{
	Info << "Normalizing mass fractions (on entry...)" << endl;

	const scalarField& TCells = thermo.T().internalField();
	
	const double maxSumMassFractionsError = 1e-3;
	double sumTotal = 0.;
	double maxError = 0.;
	forAll(TCells, celli)
	{	
		double sum = 0.;

		for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			sum += YStar_[i].internalField()[celli];

		double e = fabs(sum-1.);

		if ( e > maxSumMassFractionsError)
			FatalError << "Fatal error: the sum of mass fractions in cell " << celli << " is equal to " << sum << " (error " << sum-1. << endl; 
			
		for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			YStar_[i].internalField()[celli]/=sum;

		sumTotal += sum;
		if (e>maxError) maxError = e;
	}

	Info << "Mean error (fine structures): " << fabs( sumTotal/double(TCells.size())-1.) << endl;
	Info << "Max error (fine structures):  " << maxError << endl;
}

// Spark (to ignite the mixture)
#include "spark.H"

// Fine Structure (reactor)
HomogeneousODE fine_structure_ode(*thermodynamicsMapXML, *kineticsMapXML);
if (odeParametersFineStructures.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
	fine_structure_ode.SetCheckMassFractions(true);

// ODE Solvers
typedef OdeSMOKE::KernelDense<OpenSMOKE::HomogeneousODE_OpenSMOKE> denseOde;
typedef OdeSMOKE::MethodGear<denseOde> methodGear;
OdeSMOKE::MultiValueSolver<methodGear> odeSolver;
odeSolver.SetReactor(&fine_structure_ode);

// External ODE solvers (requires third-party libraries)
#if OPENSMOKE_USE_DVODE == 1
OpenSMOKE::HomogeneousODE_DVODE *dvodeSolver;
OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::HomogeneousODE_DVODE> *SolverDVODE;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
{
	dvodeSolver = OpenSMOKE::HomogeneousODE_DVODE::GetInstance();
	dvodeSolver->SetReactor(&fine_structure_ode);

	SolverDVODE = new OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::HomogeneousODE_DVODE>(dvodeSolver);
	SolverDVODE->SetDimensions(fine_structure_ode.NumberOfEquations());
	SolverDVODE->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDVODE->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDVODE->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_ODEPACK == 1
OpenSMOKE::HomogeneousODE_DLSODE *dlsodeSolver;
OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::HomogeneousODE_DLSODE> *SolverDLSODE;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodeSolver = OpenSMOKE::HomogeneousODE_DLSODE::GetInstance();
	dlsodeSolver->SetReactor(&fine_structure_ode);

	SolverDLSODE = new OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::HomogeneousODE_DLSODE>(dlsodeSolver);
	SolverDLSODE->SetDimensions(fine_structure_ode.NumberOfEquations());
	SolverDLSODE->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDLSODE->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDLSODE->SetAnalyticalJacobian(false);
}
OpenSMOKE::HomogeneousODE_DLSODA *dlsodaSolver;
OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::HomogeneousODE_DLSODA> *SolverDLSODA;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodaSolver = OpenSMOKE::HomogeneousODE_DLSODA::GetInstance();
	dlsodaSolver->SetReactor(&fine_structure_ode);

	SolverDLSODA = new OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::HomogeneousODE_DLSODA>(dlsodaSolver);
	SolverDLSODA->SetDimensions(fine_structure_ode.NumberOfEquations());
	SolverDLSODA->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDLSODA->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDLSODA->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_SUNDIALS == 1
OpenSMOKE::HomogeneousODE_CVODE *cvodeSolver;
OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::HomogeneousODE_CVODE> *SolverCVODE;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
{
	cvodeSolver = OpenSMOKE::HomogeneousODE_CVODE::GetInstance();
	cvodeSolver->SetReactor(&fine_structure_ode);

	SolverCVODE = new OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::HomogeneousODE_CVODE>(cvodeSolver);
	SolverCVODE->SetDimensions(fine_structure_ode.NumberOfEquations());
	SolverCVODE->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverCVODE->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverCVODE->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_DASPK == 1
OpenSMOKE::HomogeneousODE_DASPK *daspkSolver;
OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::HomogeneousODE_DASPK> *SolverDASPK;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
{
	daspkSolver = OpenSMOKE::HomogeneousODE_DASPK::GetInstance(thermodynamicsMapXML->NumberOfSpecies()+1);
	daspkSolver->SetReactor(&fine_structure_ode);

	SolverDASPK = new OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::HomogeneousODE_DASPK>(daspkSolver);
	SolverDASPK->SetDimensions(fine_structure_ode.NumberOfEquations());
	SolverDASPK->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDASPK->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDASPK->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_RADAU == 1
OpenSMOKE::HomogeneousODE_RADAU5 *radau5Solver;
OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::HomogeneousODE_RADAU5> *SolverRADAU5;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
{
	radau5Solver = OpenSMOKE::HomogeneousODE_RADAU5::GetInstance();
	radau5Solver->SetReactor(&fine_structure_ode);

	SolverRADAU5 = new OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::HomogeneousODE_RADAU5>(radau5Solver);
	SolverRADAU5->SetDimensions(fine_structure_ode.NumberOfEquations());
	SolverRADAU5->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverRADAU5->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverRADAU5->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_MEBDF == 1
OpenSMOKE::HomogeneousODE_MEBDF *mebdfSolver;
OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::HomogeneousODE_MEBDF> *SolverMEBDF;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
{
	mebdfSolver = OpenSMOKE::HomogeneousODE_MEBDF::GetInstance();
	mebdfSolver->SetReactor(&fine_structure_ode);

	SolverMEBDF = new OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::HomogeneousODE_MEBDF>(mebdfSolver);
	SolverMEBDF->SetDimensions(fine_structure_ode.NumberOfEquations());
	SolverMEBDF->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverMEBDF->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverMEBDF->SetAnalyticalJacobian(false);
}
#endif

