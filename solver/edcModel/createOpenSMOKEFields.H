/*----------------------------------------------------------------------*\
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
|   Copyright(C) 2014 A. Cuoci, M.R. Malik, A. Parente                    |
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

// Read the kinetic scheme in XML format
OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML; 
OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML;
OpenSMOKE::TransportPropertiesMap_CHEMKIN* transportMapXML;

{	
	word kinetics(thermo.lookup("kineticMechanism"));
	boost::filesystem::path path_kinetics = kinetics;

	rapidxml::xml_document<> doc;
	std::vector<char> xml_string;
	OpenSMOKE::OpenInputFileXML(doc, xml_string,path_kinetics / "kinetics.xml");

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc); 
	transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc); 
	kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, doc); 					
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
}

// Total number of equations for ODE systems (without DRG)
const int ne_total = thermodynamicsMapXML->NumberOfSpecies()+3;	// species and temperature + 2 dummy variables

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

// DRG (Directed Relation Graph)
OpenSMOKE::DRG* drg;
Switch drg_analysis = false;
double drg_minimum_temperature_for_chemistry = 300.;
List<double>  drg_epsilon;
List<double>  drg_temperature;
const dictionary& drgDictionary = combustionOptionsDictionary.subDict("DRG");
{
	drg_analysis = Switch(drgDictionary.lookup(word("drg")));

	if ( drg_analysis == true)
	{
		drg_minimum_temperature_for_chemistry  = readScalar(drgDictionary.lookup("minTemperature"));

		std::vector<unsigned int> drgListSpecies;
		List<word>  listSpecies(drgDictionary.lookup("species"));
		drgListSpecies.resize(listSpecies.size());
		for (int i=0;i<drgListSpecies.size();i++)
			drgListSpecies[i] = thermodynamicsMapXML->IndexOfSpecies(listSpecies[i])-1;

		drg_epsilon = readList<double>(drgDictionary.lookup("epsilon"));
		drg_temperature = readList<double>(drgDictionary.lookup("temperature"));

		drg = new OpenSMOKE::DRG(thermodynamicsMapXML, kineticsMapXML);
		drg->SetKeySpecies(drgListSpecies);
	}
}

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
		#if STEADYSTATE == 1
		FatalError << "LES is available only for unsteady solvers!" << exit(FatalError);
		#else
		iLES = true;
		Info << "This is a Large Eddy Simulation (LES)..." << endl;
		#endif
	}
	else
	{
		FatalError << "The ED models require turbulent conditions (i.e. LESModel or RASModel)" << exit(FatalError);
	}
}

// Check if epsilon can be found (only for RAS)
#if STEADYSTATE == 0
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
#endif

// Models
std::string combustion_model_tag;
word combustionModel(combustionOptionsDictionary.lookup("combustionModel"));
const scalar maxT = readScalar(combustionOptionsDictionary.lookup("maxT"));
const scalar minT = readScalar(combustionOptionsDictionary.lookup("minT"));
Switch laminarDiffusion(combustionOptionsDictionary.lookup("laminarDiffusion"));
scalar Sct = readScalar(combustionOptionsDictionary.lookup("Sct"));

Info<< "Creating additional terms for EDC-LES\n" << endl;

volScalarField tauStar
(
    IOobject
    (
        "tauStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("tauStar", dimensionSet(0,0,1,0,0,0,0), 0.0)
);

volScalarField epsilonSGS
(
    IOobject
    (
        "epsilonSGS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilonSGS", dimensionSet(0,2,-3,0,0,0,0), 0.0)
);

volScalarField gammaL
(
    IOobject
    (
        "gammaL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("gammaL", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

volScalarField massDiffusionInEnergyEquation
(
    IOobject
    (
        "massDiffusionInEnergyEquation",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("massDiffusionInEnergyEquation", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
);

// Counters
label  counterToUpdateChemistry = 0;

// Numerical parameters for EDC/PaSR
scalar Ccsi = 2.1377;
scalar Ctau = 0.40825;
scalar Cmix = 1.;
word   edc_version = "v2005";
scalar maxGammaL = 0.9;
scalar numberOfResidenceTimes = 100.;
label  numberIterationsToUpdateChemistry = 1;
scalar edc_exp = 2.;

// EDC v2016: Coefficients revised by A. Parente
scalar Ctau_coefficient   = 0.75;
scalar Cgamma_coefficient = std::sqrt(3./2.);
scalar CRe_coefficient    = 0.09;

scalar threshold_temperature_chemical_time = 300.;
scalar threshold_chemical_time = 0.1;
CharacteristicChemicalTimes::CharacteristicChemicalTimesType characteristic_chemical_times_type = CharacteristicChemicalTimes::CHEMICAL_TIMES_FORMATION_RATES;

const double minimum_Da = 1e-10;
const double maximum_Da = 10.;

// Numerical parameters for ED
scalar constA = 4.;
scalar constB = 0.5;
std::vector< std::vector<unsigned int> > ed_reactant_index_;
std::vector< std::vector<unsigned int> > ed_product_index_;
std::vector< std::vector<double> > ed_reactant_nu_;
std::vector< std::vector<double> > ed_product_nu_;
std::vector< std::vector<double> > ed_denominator_reactants_;
std::vector<double> ed_denominator_products_;

// Reactor type: PSR or PFR
enum { REACTORTYPE_PSR, REACTORTYPE_PFR} reactorType;

// Solver type: ODE or NLS
enum { SOLVERTYPE_ODE, SOLVERTYPE_NLS} solverType;

// Reactor type: PSR or PFR
enum { SCALE_KOLMOGOROV, SCALE_MEAN, SCALE_GLOBAL } tauMixType;

// Reading constants from file
{
	// ED Constants
	if (combustionModel == "ED")
	{
		combustion_model_tag = "ed";

		dictionary edSubDictionary(combustionOptionsDictionary.subDict("ED"));

		constA = readScalar(edSubDictionary.lookup("A"));	// default: 4
		constB = readScalar(edSubDictionary.lookup("B"));	// default: 0.5
	}
	// EDFR Constants
	else if (combustionModel == "EDFR")
	{
		combustion_model_tag = "edfr";

		dictionary edfrSubDictionary(combustionOptionsDictionary.subDict("EDFR"));

		constA = readScalar(edfrSubDictionary.lookup("A"));	// default: 4
		constB = readScalar(edfrSubDictionary.lookup("B"));	// default: 0.5
	}
	// EDC Constants
	else if (combustionModel == "EDC")
	{
		combustion_model_tag = "edc";

		dictionary edcSubDictionary(combustionOptionsDictionary.subDict("EDC"));

		edc_version = word(edcSubDictionary.lookup("version"));								// default: v2005
		Ccsi = readScalar(edcSubDictionary.lookup("Ccsi"));								// default: 2.1377
		Ctau = readScalar(edcSubDictionary.lookup("Ctau"));								// default: 0.40825
		maxGammaL = readScalar(edcSubDictionary.lookup("maxGammaL"));							// default: 0.9
		numberOfResidenceTimes = readScalar(edcSubDictionary.lookup("numberOfResidenceTimes"));				// default: 100.
		numberIterationsToUpdateChemistry = readLabel(edcSubDictionary.lookup("numberIterationsToUpdateChemistry")); 	// default 1
		counterToUpdateChemistry = numberIterationsToUpdateChemistry;

		// Version
		if (edc_version == "v2005")
		{
			edc_exp = 2.;
		}
		else if (edc_version == "v1981")
		{
			edc_exp = 3.;
		}
		else if (edc_version == "v2016")
		{
			edc_exp = 2.;
		}
		else
		{
			Info << "Wrong EDC version: v1981 || v2005 || v2016" << endl;
			abort();
		}


		// Reactor type options
		{
			word reactorTypeString(edcSubDictionary.lookup("reactorType"));
	
			if ( reactorTypeString == "PSR")
			     reactorType = REACTORTYPE_PSR;
			else if ( reactorTypeString == "PFR")
				  reactorType = REACTORTYPE_PFR;
			else
			{
				Info << "Wrong reactor type: PSR || PFR" << endl;
				abort();
			}
		}

		// Solver type options
		{
			word solverTypeString(edcSubDictionary.lookup("solverType"));
			if ( solverTypeString == "ODE")
			     solverType = SOLVERTYPE_ODE;
			else if ( solverTypeString == "NLS")
				  solverType = SOLVERTYPE_NLS;
			else
			{
				Info << "Wrong solver type: ODE || NLS" << endl;
				abort();
			}
		}

		// Coefficients (only in case of v2016)
		if (edc_version == "v2016")
		{
			Ctau_coefficient   = readScalar(edcSubDictionary.lookup("CtauCoeff"));		// default: 0.75
			Cgamma_coefficient = readScalar(edcSubDictionary.lookup("CgammaCoeff"));	// default: sqrt(3/2)
			CRe_coefficient    = readScalar(edcSubDictionary.lookup("CReCoeff"));		// default: 0.09

			// The following coeffients are currently implemented in the official OF version
			//  * Ctau_coefficient   = 0.01
			//  * Cgamma_coefficient = 0.5
			//  * CRe_coefficient    = 1.0
		}

		// Characteristic chemical type options (only in case of v2016)
		if (edc_version == "v2016")
		{
			threshold_temperature_chemical_time = readScalar(edcSubDictionary.lookup("ChemicalTimesTemperatureThreshold"));
			threshold_chemical_time = readScalar(edcSubDictionary.lookup("ChemicalTimesThreshold"));
			word chemicalTimesTypeString(edcSubDictionary.lookup("ChemicalTimesType"));

			if (chemicalTimesTypeString == "formationRates")
				characteristic_chemical_times_type = CharacteristicChemicalTimes::CHEMICAL_TIMES_FORMATION_RATES;
			else if (chemicalTimesTypeString == "reactionRates")
				characteristic_chemical_times_type = CharacteristicChemicalTimes::CHEMICAL_TIMES_REACTION_RATES;
			else if (chemicalTimesTypeString == "eigenvalues")
				characteristic_chemical_times_type = CharacteristicChemicalTimes::CHEMICAL_TIMES_EIGENVALUES;
			else
			{
				Info << "Wrong chemical times type: formationRates | reactionRates | eigenvalues" << endl;
				abort();
			}
		}
	}
	// PaSR Constants
	else if (combustionModel == "PaSR")
	{
		combustion_model_tag = "pasr";

		dictionary PaSRSubDictionary(combustionOptionsDictionary.subDict("PaSR"));
		word tauMixTypestring(PaSRSubDictionary.lookup("tauMixType"));
	
		if ( tauMixTypestring == "kolmogorovScale")
			tauMixType = SCALE_KOLMOGOROV;
		else if ( tauMixTypestring == "geometricMeanScale")
			tauMixType = SCALE_MEAN;
		else if ( tauMixTypestring == "globalScale")
			tauMixType = SCALE_GLOBAL;
		else
		{
			Info << "Wrong tau_mix type: kolmogorovScale || geometriMeanScale || globalScale" << endl;
			abort();
		}

		Cmix = readScalar(PaSRSubDictionary.lookup("Cmix"));                        
		numberIterationsToUpdateChemistry = readLabel(PaSRSubDictionary.lookup("numberIterationsToUpdateChemistry")); 
		counterToUpdateChemistry = numberIterationsToUpdateChemistry;

		// Type options
		{
			word reactorTypeString(PaSRSubDictionary.lookup("reactorType"));
			if ( reactorTypeString == "PFR")
			{
				  reactorType = REACTORTYPE_PFR;
			}
			else
			{
				Info << "Wrong reactor type: only PFR" << endl;
				abort();
			}
		}

		// Characteristic chemical type options
		{
			threshold_temperature_chemical_time = readScalar(PaSRSubDictionary.lookup("ChemicalTimesTemperatureThreshold"));
			threshold_chemical_time = readScalar(PaSRSubDictionary.lookup("ChemicalTimesThreshold"));
			word chemicalTimesTypeString(PaSRSubDictionary.lookup("ChemicalTimesType"));

			if (chemicalTimesTypeString == "formationRates")
				characteristic_chemical_times_type = CharacteristicChemicalTimes::CHEMICAL_TIMES_FORMATION_RATES;
			else if (chemicalTimesTypeString == "reactionRates")
				characteristic_chemical_times_type = CharacteristicChemicalTimes::CHEMICAL_TIMES_REACTION_RATES;
			else if (chemicalTimesTypeString == "eigenvalues")
				characteristic_chemical_times_type = CharacteristicChemicalTimes::CHEMICAL_TIMES_EIGENVALUES;
			else
			{
				Info << "Wrong chemical times type: formationRates | reactionRates | eigenvalues" << endl;
				abort();
			}
		}
	}
	// Error Message
	else
	{
		FatalError << "Combustion model available: ED || EDFR || EDC || PaSR" << endl;
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
				ed_denominator_reactants_[it.row()].push_back( it.value()*thermodynamicsMapXML->MW(it.col()) );
			}
		}
	
		for (int k=0; k<kineticsMapXML->stoichiometry().stoichiometric_matrix_products().outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML->stoichiometry().stoichiometric_matrix_products(),k); it; ++it)
			{
				ed_product_index_[it.row()].push_back(it.col());
				ed_product_nu_[it.row()].push_back(it.value());
				ed_denominator_products_[it.row()]+=it.value()*thermodynamicsMapXML->MW(it.col());
			};
		}
	}
}

// Ode parameters
OpenSMOKE::ODE_Parameters odeParametersFineStructures;

if (combustionModel == "EDC" || combustionModel == "PaSR")
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

// Nls parameters
NonLinearSystemParameters nlsParametersFineStructures;

if (combustionModel == "EDC")
{
	const dictionary& nlsParametersFineStructuresDictionary = combustionOptionsDictionary.subDict("Nls");

	//- Mass fractions tolerance
	scalar relTolerance = readScalar(nlsParametersFineStructuresDictionary.lookup("relTolerance"));
	scalar absTolerance = readScalar(nlsParametersFineStructuresDictionary.lookup("absTolerance"));
	scalar thresholdResiduals = readScalar(nlsParametersFineStructuresDictionary.lookup("thresholdResiduals"));
	scalar maxAttempts = readScalar(nlsParametersFineStructuresDictionary.lookup("maxAttempts"));

	nlsParametersFineStructures.relTolerance = relTolerance;
	nlsParametersFineStructures.absTolerance = absTolerance;
	nlsParametersFineStructures.thresholdResiduals = thresholdResiduals;
	nlsParametersFineStructures.maxAttempts = maxAttempts;
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
PtrList<volScalarField> Dmix_;
PtrList<volScalarField> HeSpecies_;

//scalarField Tau(mesh.nCells(),0.);
scalarField gammaStar(mesh.nCells(),0.);
scalarField tfStar(mesh.nCells(),0.);
scalarField mDotStar(mesh.nCells(),0.);
scalarField edc_coefficient(mesh.nCells(),0.);

// Characteristic mixing time [s]
volScalarField tauMix
(
	IOobject
	(
		"combustion::tauMix",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE,
		false
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.), 
	zeroGradientFvPatchScalarField::typeName
);

// Characteristic chemical time [s]
volScalarField tauChem
(
	IOobject
	(
		"combustion::tauChem",
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE,
		false
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.), 
	zeroGradientFvPatchScalarField::typeName
);

// Damkholer number
volScalarField Da
(
	IOobject
	(
		"combustion::Da",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE,
		false
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.), 
	zeroGradientFvPatchScalarField::typeName
);

// Turbulent Reynolds number
volScalarField ReT
(
	IOobject
	(
		"combustion::ReT",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE,
		false
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.), 
	zeroGradientFvPatchScalarField::typeName
);

// PASR coefficient
volScalarField kappa
(
	IOobject
	(
		"combustion::kappa",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE,
		false
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.), 
	zeroGradientFvPatchScalarField::typeName
);

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
	
	if (laminarDiffusion == true)
	{
		HeSpecies_.resize(Y.size());
		Dmix_.resize(Y.size());
	}

	unsigned int species_from_file = 0;
	unsigned int species_from_default = 0;

	forAll(Y,i)
	{
		std::cout << " * Processing species "<< i << " : " << Y[i].name() << " : ";
		IOobject header
		(
			combustion_model_tag + "::YStar_" + Y[i].name(),
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		// check if field exists and can be read
		if (combustionModel == "EDC" || combustionModel == "PaSR")
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
							combustion_model_tag + "::YStar_"+ Y[i].name(),
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
							combustion_model_tag + "::YStar_"+ Y[i].name(),
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
					combustion_model_tag + "::R_"+ Y[i].name(),
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.)
			)
		);	

		if (laminarDiffusion == true)
		{
			Dmix_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"gas::Dmix_"+ Y[i].name(),
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("Dmix", dimensionSet(0, 2, -1, 0, 0), 1e-10)
				)
			);
		}	

		if (laminarDiffusion == true)
		{
			HeSpecies_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"thermo::He_"+ Y[i].name(),
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("He", dimensionSet(0, 2, -2, 0, 0), 0.)
				)
			);
		}	
	}
}

// TStar
if (combustionModel == "EDC" || combustionModel == "PaSR")
{
	std::cout << " * Reading TStar field..." << std::endl;
	{
		IOobject header
		(
			combustion_model_tag + "::TStar",
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
						combustion_model_tag + "::TStar",
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
						combustion_model_tag + "::TStar",
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

volScalarField dQ
(
	IOobject
	(
       		"dQ",
        	runTime.timeName(),
        	mesh,
        	IOobject::NO_READ,
        	IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);

dimensionedScalar ksmall_("ksmall", dimensionSet(0,2,-2,0,0,0,0), SMALL );
dimensionedScalar musmall_("musmall", dimensionSet(1,-1,-1,0,0,0,0), SMALL );
dimensionedScalar epsilonsmall_("epsilonsmall", dimensionSet(0,2,-3,0,0,0,0), SMALL );

volScalarField one_over_tauStar
(
    IOobject
    (
        "one_over_tauStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one_over_tauStar", dimensionSet(0,0,-1,0,0,0,0), 0.0)
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
if (combustionModel == "EDC" || combustionModel == "PaSR")
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

#if EDCSMOKE_USE_ISAT == 1
    #include "readOptions_ISAT.H"
#endif

// Spark (to ignite the mixture)
#include "spark.H"

// --------------------------------------------------------------------------------------------------------
// Fine structure as a PSR
// --------------------------------------------------------------------------------------------------------

// Fine Structure (reactor)
ODE_PSR fine_structure_ode_psr(*thermodynamicsMapXML, *kineticsMapXML);
if (odeParametersFineStructures.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
	fine_structure_ode_psr.SetCheckMassFractions(true);

// Fine Structure (reactor)
NLS_PSR fine_structure_nls_psr(*thermodynamicsMapXML, *kineticsMapXML);

// ODE Solvers
typedef OdeSMOKE::KernelDense<OpenSMOKE::ODE_PSR_OpenSMOKE> denseOdePSR;
typedef OdeSMOKE::MethodGear<denseOdePSR> methodGearPSR;
OdeSMOKE::MultiValueSolver<methodGearPSR> odeSolverPSR;
odeSolverPSR.SetReactor(&fine_structure_ode_psr);

// External ODE solvers (requires third-party libraries)
#if OPENSMOKE_USE_DVODE == 1
OpenSMOKE::ODE_PSR_DVODE *dvodeSolverPSR;
OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::ODE_PSR_DVODE> *SolverDVODEPSR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
{
	dvodeSolverPSR = OpenSMOKE::ODE_PSR_DVODE::GetInstance();
	dvodeSolverPSR->SetReactor(&fine_structure_ode_psr);

	SolverDVODEPSR = new OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::ODE_PSR_DVODE>(dvodeSolverPSR);
	SolverDVODEPSR->SetDimensions(fine_structure_ode_psr.NumberOfEquations());
	SolverDVODEPSR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDVODEPSR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDVODEPSR->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_ODEPACK == 1
OpenSMOKE::ODE_PSR_DLSODE *dlsodeSolverPSR;
OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::ODE_PSR_DLSODE> *SolverDLSODEPSR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodeSolverPSR = OpenSMOKE::ODE_PSR_DLSODE::GetInstance();
	dlsodeSolverPSR->SetReactor(&fine_structure_ode_psr);

	SolverDLSODEPSR = new OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::ODE_PSR_DLSODE>(dlsodeSolverPSR);
	SolverDLSODEPSR->SetDimensions(fine_structure_ode_psr.NumberOfEquations());
	SolverDLSODEPSR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDLSODEPSR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDLSODEPSR->SetAnalyticalJacobian(false);
}
OpenSMOKE::ODE_PSR_DLSODA *dlsodaSolverPSR;
OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::ODE_PSR_DLSODA> *SolverDLSODAPSR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodaSolverPSR = OpenSMOKE::ODE_PSR_DLSODA::GetInstance();
	dlsodaSolverPSR->SetReactor(&fine_structure_ode_psr);

	SolverDLSODAPSR = new OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::ODE_PSR_DLSODA>(dlsodaSolverPSR);
	SolverDLSODAPSR->SetDimensions(fine_structure_ode_psr.NumberOfEquations());
	SolverDLSODAPSR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDLSODAPSR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDLSODAPSR->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_SUNDIALS == 1
OpenSMOKE::ODE_PSR_CVODE *cvodeSolverPSR;
OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::ODE_PSR_CVODE> *SolverCVODEPSR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
{
	cvodeSolverPSR = OpenSMOKE::ODE_PSR_CVODE::GetInstance();
	cvodeSolverPSR->SetReactor(&fine_structure_ode_psr);

	SolverCVODEPSR = new OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::ODE_PSR_CVODE>(cvodeSolverPSR);
	SolverCVODEPSR->SetDimensions(fine_structure_ode_psr.NumberOfEquations());
	SolverCVODEPSR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverCVODEPSR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverCVODEPSR->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_DASPK == 1
OpenSMOKE::ODE_PSR_DASPK *daspkSolverPSR;
OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::ODE_PSR_DASPK> *SolverDASPKPSR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
{
	daspkSolverPSR = OpenSMOKE::ODE_PSR_DASPK::GetInstance(thermodynamicsMapXML->NumberOfSpecies()+1);
	daspkSolverPSR->SetReactor(&fine_structure_ode_psr);

	SolverDASPKPSR = new OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::ODE_PSR_DASPK>(daspkSolverPSR);
	SolverDASPKPSR->SetDimensions(fine_structure_ode_psr.NumberOfEquations());
	SolverDASPKPSR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDASPKPSR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDASPKPSR->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_RADAU == 1
OpenSMOKE::ODE_PSR_RADAU5 *radau5SolverPSR;
OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::ODE_PSR_RADAU5> *SolverRADAU5PSR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
{
	radau5SolverPSR = OpenSMOKE::ODE_PSR_RADAU5::GetInstance();
	radau5SolverPSR->SetReactor(&fine_structure_ode_psr);

	SolverRADAU5PSR = new OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::ODE_PSR_RADAU5>(radau5SolverPSR);
	SolverRADAU5PSR->SetDimensions(fine_structure_ode_psr.NumberOfEquations());
	SolverRADAU5PSR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverRADAU5PSR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverRADAU5PSR->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_MEBDF == 1
OpenSMOKE::ODE_PSR_MEBDF *mebdfSolverPSR;
OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::ODE_PSR_MEBDF> *SolverMEBDFPSR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
{
	mebdfSolverPSR = OpenSMOKE::ODE_PSR_MEBDF::GetInstance();
	mebdfSolverPSR->SetReactor(&fine_structure_ode_psr);

	SolverMEBDFPSR = new OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::ODE_PSR_MEBDF>(mebdfSolverPSR);
	SolverMEBDFPSR->SetDimensions(fine_structure_ode_psr.NumberOfEquations());
	SolverMEBDFPSR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverMEBDFPSR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverMEBDFPSR->SetAnalyticalJacobian(false);
}
#endif


// --------------------------------------------------------------------------------------------------------
// Fine structure as a PFR
// --------------------------------------------------------------------------------------------------------

// Fine Structure (reactor)
ODE_PFR fine_structure_ode_pfr(*thermodynamicsMapXML, *kineticsMapXML);
if (odeParametersFineStructures.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
	fine_structure_ode_pfr.SetCheckMassFractions(true);

// ODE Solvers
typedef OdeSMOKE::KernelDense<OpenSMOKE::ODE_PFR_OpenSMOKE> denseOdePFR;
typedef OdeSMOKE::MethodGear<denseOdePFR> methodGearPFR;
OdeSMOKE::MultiValueSolver<methodGearPFR> odeSolverPFR;
odeSolverPFR.SetReactor(&fine_structure_ode_pfr);

// External ODE solvers (requires third-party libraries)
#if OPENSMOKE_USE_DVODE == 1
OpenSMOKE::ODE_PFR_DVODE *dvodeSolverPFR;
OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::ODE_PFR_DVODE> *SolverDVODEPFR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
{
	dvodeSolverPFR = OpenSMOKE::ODE_PFR_DVODE::GetInstance();
	dvodeSolverPFR->SetReactor(&fine_structure_ode_pfr);

	SolverDVODEPFR = new OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::ODE_PFR_DVODE>(dvodeSolverPFR);
	SolverDVODEPFR->SetDimensions(fine_structure_ode_pfr.NumberOfEquations());
	SolverDVODEPFR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDVODEPFR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDVODEPFR->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_ODEPACK == 1
OpenSMOKE::ODE_PFR_DLSODE *dlsodeSolverPFR;
OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::ODE_PFR_DLSODE> *SolverDLSODEPFR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodeSolverPFR = OpenSMOKE::ODE_PFR_DLSODE::GetInstance();
	dlsodeSolverPFR->SetReactor(&fine_structure_ode_pfr);

	SolverDLSODEPFR = new OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::ODE_PFR_DLSODE>(dlsodeSolverPFR);
	SolverDLSODEPFR->SetDimensions(fine_structure_ode_pfr.NumberOfEquations());
	SolverDLSODEPFR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDLSODEPFR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDLSODEPFR->SetAnalyticalJacobian(false);
}
OpenSMOKE::ODE_PFR_DLSODA *dlsodaSolverPFR;
OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::ODE_PFR_DLSODA> *SolverDLSODAPFR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodaSolverPFR = OpenSMOKE::ODE_PFR_DLSODA::GetInstance();
	dlsodaSolverPFR->SetReactor(&fine_structure_ode_pfr);

	SolverDLSODAPFR = new OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::ODE_PFR_DLSODA>(dlsodaSolverPFR);
	SolverDLSODAPFR->SetDimensions(fine_structure_ode_pfr.NumberOfEquations());
	SolverDLSODAPFR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDLSODAPFR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDLSODAPFR->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_SUNDIALS == 1
OpenSMOKE::ODE_PFR_CVODE *cvodeSolverPFR;
OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::ODE_PFR_CVODE> *SolverCVODEPFR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
{
	cvodeSolverPFR = OpenSMOKE::ODE_PFR_CVODE::GetInstance();
	cvodeSolverPFR->SetReactor(&fine_structure_ode_pfr);

	SolverCVODEPFR = new OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::ODE_PFR_CVODE>(cvodeSolverPFR);
	SolverCVODEPFR->SetDimensions(fine_structure_ode_pfr.NumberOfEquations());
	SolverCVODEPFR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverCVODEPFR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverCVODEPFR->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_DASPK == 1
OpenSMOKE::ODE_PFR_DASPK *daspkSolverPFR;
OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::ODE_PFR_DASPK> *SolverDASPKPFR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
{
	daspkSolverPFR = OpenSMOKE::ODE_PFR_DASPK::GetInstance(thermodynamicsMapXML->NumberOfSpecies()+1);
	daspkSolverPFR->SetReactor(&fine_structure_ode_pfr);

	SolverDASPKPFR = new OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::ODE_PFR_DASPK>(daspkSolverPFR);
	SolverDASPKPFR->SetDimensions(fine_structure_ode_pfr.NumberOfEquations());
	SolverDASPKPFR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverDASPKPFR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverDASPKPFR->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_RADAU == 1
OpenSMOKE::ODE_PFR_RADAU5 *radau5SolverPFR;
OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::ODE_PFR_RADAU5> *SolverRADAU5PFR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
{
	radau5SolverPFR = OpenSMOKE::ODE_PFR_RADAU5::GetInstance();
	radau5SolverPFR->SetReactor(&fine_structure_ode_pfr);

	SolverRADAU5PFR = new OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::ODE_PFR_RADAU5>(radau5SolverPFR);
	SolverRADAU5PFR->SetDimensions(fine_structure_ode_pfr.NumberOfEquations());
	SolverRADAU5PFR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverRADAU5PFR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverRADAU5PFR->SetAnalyticalJacobian(false);
}
#endif

#if OPENSMOKE_USE_MEBDF == 1
OpenSMOKE::ODE_PFR_MEBDF *mebdfSolverPFR;
OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::ODE_PFR_MEBDF> *SolverMEBDFPFR;
if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
{
	mebdfSolverPFR = OpenSMOKE::ODE_PFR_MEBDF::GetInstance();
	mebdfSolverPFR->SetReactor(&fine_structure_ode_pfr);

	SolverMEBDFPFR = new OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::ODE_PFR_MEBDF>(mebdfSolverPFR);
	SolverMEBDFPFR->SetDimensions(fine_structure_ode_pfr.NumberOfEquations());
	SolverMEBDFPFR->SetAbsoluteTolerance(odeParametersFineStructures.absolute_tolerance());
	SolverMEBDFPFR->SetRelativeTolerance(odeParametersFineStructures.relative_tolerance());
	SolverMEBDFPFR->SetAnalyticalJacobian(false);
}
#endif

#if EDCSMOKE_USE_ISAT == 1
    #include "memoryAllocation_ISAT.H"
#endif
