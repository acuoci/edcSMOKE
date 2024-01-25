/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ Suite.                               |
|                                                                         |
|   Copyright(C) 2014, 2013  Alberto Cuoci                                |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// Thermodynamics
#include "kernel/thermo/Species.h"
#include "kernel/thermo/ThermoPolicy_CHEMKIN.h"
#include "kernel/thermo/ThermoReader.h"
#include "kernel/thermo/ThermoReaderPolicy_CHEMKIN.h"

// Transport
#include "kernel/transport/TransportPolicy_CHEMKIN.h"
#include "kernel/transport/TransportReader.h"
#include "kernel/transport/TransportReaderPolicy_CHEMKIN.h"

// Kinetics
#include "kernel/kinetics/ReactionPolicy_CHEMKIN.h"
#include "kernel/kinetics/ReactionPolicy_Surface_CHEMKIN.h"
#include "kernel/kinetics/ReactionPolicy_Solid_CHEMKIN.h"
#include "kernel/kinetics/ReactionPolicy_Liquid_CHEMKIN.h"

// Preprocessing
#include "preprocessing/PreProcessorSpecies.h"
#include "preprocessing/PreProcessorKinetics.h"
#include "preprocessing/PreProcessorKineticsPolicy_CHEMKIN.h"
#include "preprocessing/PreProcessorSpeciesPolicy_CHEMKIN_WithTransport.h"
#include "preprocessing/PreProcessorSurfaceKineticsPolicy_CHEMKIN.h"
#include "preprocessing/PreProcessorSolidKineticsPolicy_CHEMKIN.h"
#include "preprocessing/PreProcessorLiquidKineticsPolicy_CHEMKIN.h"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/TransportPropertiesMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"
#include "maps/ThermodynamicsMap_Surface_CHEMKIN.h"
#include "maps/KineticsMap_Surface_CHEMKIN.h"
#include "maps/ThermodynamicsMap_Solid_CHEMKIN.h"
#include "maps/KineticsMap_Solid_CHEMKIN.h"
#include "maps/ThermodynamicsMap_Liquid_CHEMKIN.h"
#include "maps/KineticsMap_Liquid_CHEMKIN.h"

// Analyzers
#include "analyzers/AnalyzerKineticMechanism.h"

// Grammar rules
#include "Grammar_CHEMKIN_PreProcessor.H"
#include "Grammar_Comments.H"

// License
// #include "licensegenerator/license/OpenSMOKELicenseUtilities.hpp"


int main(int argc, char** argv)
{
	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_CHEMKIN_PreProcessor", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	unsigned int max_number_allowed_species = 100000;
	//OpenSMOKE::OpenSMOKE_CheckLicense(executable_folder, "CHEMKIN_PreProcessor", max_number_allowed_species);

	typedef OpenSMOKE::Species< OpenSMOKE::ThermoPolicy_CHEMKIN, OpenSMOKE::TransportPolicy_CHEMKIN > SpeciesCHEMKIN;
	typedef OpenSMOKE::PreProcessorSpecies< OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<SpeciesCHEMKIN>  > PreProcessorSpecies_CHEMKIN;
	typedef OpenSMOKE::PreProcessorSpecies< OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<SpeciesCHEMKIN>  > PreProcessorSpecies_CHEMKIN_OnlyThermodynamics;
	typedef OpenSMOKE::PreProcessorSpecies< OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<SpeciesCHEMKIN>  > PreProcessorSpecies_CHEMKIN_WithoutTransport;
	typedef OpenSMOKE::PreProcessorKinetics< OpenSMOKE::PreProcessorKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_CHEMKIN> > PreProcessorKinetics_CHEMKIN;
	typedef OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > > ThermoReader_CHEMKIN;
	typedef OpenSMOKE::PreProcessorKinetics< OpenSMOKE::PreProcessorSurfaceKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_Surface_CHEMKIN> >	PreProcessorKinetics_Surface_CHEMKIN;
	typedef OpenSMOKE::PreProcessorKinetics< OpenSMOKE::PreProcessorSolidKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_Solid_CHEMKIN> >		PreProcessorKinetics_Solid_CHEMKIN;
	typedef OpenSMOKE::PreProcessorKinetics< OpenSMOKE::PreProcessorLiquidKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_Liquid_CHEMKIN> >	PreProcessorKinetics_Liquid_CHEMKIN;

	std::string input_file_name_ = "input.dic";
	std::string main_dictionary_name_ = "CHEMKIN_PreProcessor";

	// Program options from command line
	{
		namespace po = boost::program_options;
		po::options_description description("Options for the OpenSMOKEpp_CHEMKIN_PreProcessor");
		description.add_options()
			("help", "print help messages")
			("input", po::value<std::string>(), "name of the file containing the main dictionary (default \"input.dic\")")
			("dictionary", po::value<std::string>(), "name of the main dictionary to be used (default \"CHEMKIN_PreProcessor\")");

		po::variables_map vm;
		try
		{
			po::store(po::parse_command_line(argc, argv, description), vm); // can throw 

			if (vm.count("help"))
			{
				std::cout << "Basic Command Line Parameters" << std::endl;
				std::cout << description << std::endl;
				return OPENSMOKE_SUCCESSFULL_EXIT;
			}

			if (vm.count("input"))
				input_file_name_ = vm["input"].as<std::string>();

			if (vm.count("dictionary"))
				main_dictionary_name_ = vm["dictionary"].as<std::string>();

			po::notify(vm); // throws on error, so do after help in case  there are any problems 
		}
		catch (po::error& e)
		{
			std::cerr << "Fatal error: " << e.what() << std::endl << std::endl;
			std::cerr << description << std::endl;
			return OPENSMOKE_FATAL_ERROR_EXIT;
		}
	}

	// Defines the grammar rules
	GrammarKineticInterpreter grammar;

	// Define the dictionaries
	OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
	dictionaries.ReadDictionariesFromFile(input_file_name_);
	dictionaries(main_dictionary_name_).SetGrammar(grammar);

	// File names
	boost::filesystem::path path_chemkin;
	boost::filesystem::path path_output;
	boost::filesystem::path thermo_file;
	boost::filesystem::path transport_file;
	boost::filesystem::path kinetics_file;
	boost::filesystem::path kinetics_surface_file;
	boost::filesystem::path kinetics_solid_file;
	boost::filesystem::path kinetics_liquid_file;

	// Default values
	bool preprocess_transport_data_ = false;
	bool preprocess_kinetics_ = false;
	bool preprocess_surface_kinetics_ = false;
	bool preprocess_solid_kinetics_ = false;
	bool preprocess_liquid_kinetics_ = false;
	bool write_ascii_thermodynamic_coefficients_ = true;
	bool write_ascii_files_ = false;
	bool write_xml_files_ = true;
	bool write_reaction_strings_ = true;


	// Output folder (compulsory)
	if (dictionaries(main_dictionary_name_).CheckOption("@Output") == true)
		dictionaries(main_dictionary_name_).ReadPath("@Output", path_output);

	// Thermodynamics (compulsory)
	if (dictionaries(main_dictionary_name_).CheckOption("@Thermodynamics") == true)
		dictionaries(main_dictionary_name_).ReadPath("@Thermodynamics", thermo_file);

	// Kinetic mechanism (optional)
	if (dictionaries(main_dictionary_name_).CheckOption("@Kinetics") == true)
	{
		preprocess_kinetics_ = true;
		dictionaries(main_dictionary_name_).ReadPath("@Kinetics", kinetics_file);
	}

	// Surface Kinetic mechanism (optional)
	if (dictionaries(main_dictionary_name_).CheckOption("@Surface") == true)
	{
		preprocess_surface_kinetics_ = true;
		dictionaries(main_dictionary_name_).ReadPath("@Surface", kinetics_surface_file);
	}

	// Solid Kinetic mechanism (optional)
	if (dictionaries(main_dictionary_name_).CheckOption("@Solid") == true)
	{
		preprocess_solid_kinetics_ = true;
		dictionaries(main_dictionary_name_).ReadPath("@Solid", kinetics_solid_file);
	}

	// Liquid Kinetic mechanism (optional)
	if (dictionaries(main_dictionary_name_).CheckOption("@Liquid") == true)
	{
		preprocess_liquid_kinetics_ = true;
		dictionaries(main_dictionary_name_).ReadPath("@Liquid", kinetics_liquid_file);
	}

	// Transport data (optional)
	if (dictionaries(main_dictionary_name_).CheckOption("@Transport") == true)
	{
		preprocess_transport_data_ = true;
		dictionaries(main_dictionary_name_).ReadPath("@Transport", transport_file);
	}

	// Writes the interpreted kinetic scheme in the old format
	bool write_ascii_files_old_style_ = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@OutputOldStyle") == true)
		dictionaries(main_dictionary_name_).ReadBool("@OutputOldStyle", write_ascii_files_old_style_);

	// Re-writes the kinetic mechanism in CHEMKIN format after correction of mass balances
	bool rewrite_chemkin_format_ = true;
	if (dictionaries(main_dictionary_name_).CheckOption("@RewriteCHEMKIN") == true)
		dictionaries(main_dictionary_name_).ReadBool("@RewriteCHEMKIN", rewrite_chemkin_format_);

	// Stores the original CHEMKIN file with line numbers
	bool store_chemkin_file_ = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@StoreCHEMKIN") == true)
		dictionaries(main_dictionary_name_).ReadBool("@StoreCHEMKIN", store_chemkin_file_);

	// Writes the fitting coefficients for the transport properties
	bool write_ascii_fitting_coefficients_ = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@TransportFittingCoefficients") == true)
		dictionaries(main_dictionary_name_).ReadBool("@TransportFittingCoefficients", write_ascii_fitting_coefficients_);

	// Species bundling
	bool species_bundling_ = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@SpeciesBundling") == true)
		dictionaries(main_dictionary_name_).ReadBool("@SpeciesBundling", species_bundling_);

	// Writes the thermodynamic properties in a consistent way
	bool check_thermodynamics_ = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@CheckThermodynamics") == true)
		dictionaries(main_dictionary_name_).ReadBool("@CheckThermodynamics", check_thermodynamics_);

	// Fitting the reverse kinetic constants
	bool write_fitted_kinetic_constants_ = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@ReverseFitting") == true)
		dictionaries(main_dictionary_name_).ReadBool("@ReverseFitting", write_fitted_kinetic_constants_);

	// Fitting the reverse kinetic constants
	bool write_fitted_chebyshev_polynomials_ = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@ChebyshevFitting") == true)
		dictionaries(main_dictionary_name_).ReadBool("@ChebyshevFitting", write_fitted_chebyshev_polynomials_);

	// Fitting the reverse kinetic constants
	bool violation_allowed_chebyshev_polynomials_ = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@ChebyshevViolationAllowed") == true)
		dictionaries(main_dictionary_name_).ReadBool("@ChebyshevViolationAllowed", violation_allowed_chebyshev_polynomials_);

	// Fitting the reverse kinetic constants
	bool write_reaction_tables_ = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@ReactionTables") == true)
		dictionaries(main_dictionary_name_).ReadBool("@ReactionTables", write_reaction_tables_);

	// Fitting the reverse kinetic constants
	std::vector<double> reaction_tables_list_temperatures_;
	if (dictionaries(main_dictionary_name_).CheckOption("@ReactionTablesListOfTemperatures") == true)
	{
		std::vector<std::string> list_temperatures;
		dictionaries(main_dictionary_name_).ReadOption("@ReactionTablesListOfTemperatures", list_temperatures);
		const unsigned int n = list_temperatures.size() - 1;
		const std::string units = list_temperatures[n];

		for (unsigned int i = 0; i < n; i++)
		{
			if (units == "K")	reaction_tables_list_temperatures_.push_back(boost::lexical_cast<double>(list_temperatures[i]));
			else if (units == "C")	reaction_tables_list_temperatures_.push_back(boost::lexical_cast<double>(list_temperatures[i]) + 273.15);
			else
			{
				std::cerr << "Fatal error: " << "Wrong temperature units in @ReactionTablesListOfTemperatures option" << std::endl;
				return OPENSMOKE_FATAL_ERROR_EXIT;
			}
		}
	}

	// Fitting the reverse kinetic constants
	bool write_collision_rate_analysis_ = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@CollisionRateAnalysis") == true)
		dictionaries(main_dictionary_name_).ReadBool("@CollisionRateAnalysis", write_collision_rate_analysis_);

	// Fitting the reverse kinetic constants
	std::vector<double> collision_rate_analysis_list_temperatures_;
	if (dictionaries(main_dictionary_name_).CheckOption("@CollisionRateAnalysisListOfTemperatures") == true)
	{
		std::vector<std::string> list_temperatures;
		dictionaries(main_dictionary_name_).ReadOption("@CollisionRateAnalysisListOfTemperatures", list_temperatures);
		const unsigned int n = list_temperatures.size() - 1;
		const std::string units = list_temperatures[n];

		for (unsigned int i = 0; i < n; i++)
		{
			if (units == "K")		collision_rate_analysis_list_temperatures_.push_back(boost::lexical_cast<double>(list_temperatures[i]));
			else if (units == "C")	collision_rate_analysis_list_temperatures_.push_back(boost::lexical_cast<double>(list_temperatures[i]) + 273.15);
			else
			{
				std::cerr << "Fatal error: " << "Wrong temperature units in @CollisionRateAnalysisListOfTemperatures option" << std::endl;
				return OPENSMOKE_FATAL_ERROR_EXIT;
			}
		}
	}

	// Sparsity pattern analysis
	bool sparsity_pattern_analysis_ = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@SparsityPatternAnalysis") == true)
		dictionaries(main_dictionary_name_).ReadBool("@SparsityPatternAnalysis", sparsity_pattern_analysis_);

	// Reads the comments
	bool write_comments_ = false;
	std::string author_name("undefined");
	std::string place_name("undefined");
	std::string comments("no comments");
	std::string preprocessing_date(OpenSMOKE::GetCurrentDate());
	std::string preprocessing_time(OpenSMOKE::GetCurrentTime());

	if (dictionaries(main_dictionary_name_).CheckOption("@Comments") == true)
	{
		std::string name_of_comments_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@Comments", name_of_comments_subdictionary);

		GrammarComments grammar_comments;
		dictionaries(name_of_comments_subdictionary).SetGrammar(grammar_comments);

		dictionaries(name_of_comments_subdictionary).ReadSequence("@Author", author_name);
		dictionaries(name_of_comments_subdictionary).ReadSequence("@Place", place_name);
		dictionaries(name_of_comments_subdictionary).ReadSequence("@Comments", comments);
	}

	// Create the output directory
	OpenSMOKE::CreateDirectory(path_output);

	// Log file (Setup)
	boost::filesystem::path file_ascii_log_ = path_output / "log";
	std::ofstream fLog;
	fLog.open(std::string(file_ascii_log_.string()).c_str(), std::ios::out);
	OpenSMOKE::CheckIfFileIsOpen(fLog, file_ascii_log_.string());
	fLog.setf(std::ios::scientific);


	// Thermodynamics Analyzer
	if (preprocess_kinetics_ == false)
	{
		if (preprocess_surface_kinetics_ == true)
		{
			std::cerr << "Fatal error: " << "The surface kinetic schemes can be pre-processed only together with a kinetic scheme for the homogeneous phase, to be specified through the @Kinetics option" << std::endl;
			return OPENSMOKE_FATAL_ERROR_EXIT;
		}

		if (preprocess_solid_kinetics_ == true)
		{
			std::cerr << "Fatal error: " << "The solid kinetic schemes can be pre-processed only together with a kinetic scheme for the homogeneous phase, to be specified through the @Kinetics option" << std::endl;
			return OPENSMOKE_FATAL_ERROR_EXIT;
		}

		if (preprocess_liquid_kinetics_ == true)
		{
			std::cerr << "Fatal error: " << "The liquid kinetic schemes can be pre-processed only together with a kinetic scheme for the homogeneous phase, to be specified through the @Kinetics option" << std::endl;
			return OPENSMOKE_FATAL_ERROR_EXIT;
		}

		OpenSMOKE::TransportReader< OpenSMOKE::TransportReaderPolicy_CHEMKIN<OpenSMOKE::TransportPolicy_CHEMKIN > >* transportreader;

		// Preprocessors
		PreProcessorSpecies_CHEMKIN* preprocessor_species_with_transport = nullptr;
		PreProcessorSpecies_CHEMKIN_WithoutTransport* preprocessor_species_without_transport = nullptr;

		// Reading thermodynamic database
		ThermoReader_CHEMKIN thermoreader;
		CheckForFatalError(thermoreader.ReadFromFile(thermo_file.string()));

		if (preprocess_transport_data_ == true)
		{
			// Reading transport database
			transportreader = new OpenSMOKE::TransportReader< OpenSMOKE::TransportReaderPolicy_CHEMKIN<OpenSMOKE::TransportPolicy_CHEMKIN > >;
			CheckForFatalError(transportreader->ReadFromFile(transport_file.string()));

			// Preprocessing the thermodynamic and transport properties
			if (write_ascii_fitting_coefficients_ == false)
			{
				preprocessor_species_with_transport = new PreProcessorSpecies_CHEMKIN(fLog, thermoreader, *transportreader);
				CheckForFatalError(preprocessor_species_with_transport->Setup());
				//	preprocessor_species_with_transport->ThermodynamicConsistency();
			}
			else
			{
				std::vector<std::string> kinetics;
				preprocessor_species_with_transport = new PreProcessorSpecies_CHEMKIN(thermoreader, *transportreader, kinetics, fLog);
				CheckForFatalError(preprocessor_species_with_transport->Setup());
				//	preprocessor_species_with_transport->ThermodynamicConsistency();
			}
		}
		else
		{
			// Preprocessing the thermodynamic properties
			preprocessor_species_without_transport = new PreProcessorSpecies_CHEMKIN_OnlyThermodynamics(thermoreader, fLog);
			CheckForFatalError(preprocessor_species_without_transport->Setup());
			//	preprocessor_species_without_transport->ThermodynamicConsistency();
		}

		if (check_thermodynamics_ == true)
		{
			// Write the Thermodynamic Status on file
			boost::filesystem::path file_thermo_status = path_output / "Thermodynamics_Status.out";
			//boost::filesystem::path file_thermo_matlab = path_output / "Thermodynamics.mlab";
			boost::filesystem::path file_thermo_reformulation = path_output / "thermo.CHEMKIN.CKT";
			boost::filesystem::path file_thermo_reformulation_fixed_T = path_output / "thermo.CHEMKIN.fixedT.CKT";

			if (preprocess_transport_data_ == true)
			{
				preprocessor_species_with_transport->StatusOfThermodynamics(file_thermo_status.string());
				//preprocessor_species_with_transport->WriteThermodynamicsForMATLAB(file_thermo_matlab.string());
				preprocessor_species_with_transport->ReformulationOfThermodynamics(file_thermo_reformulation.string(), thermo_file.string());
				preprocessor_species_with_transport->ReformulationOfThermodynamicsFixedIntermediateTemperature(file_thermo_reformulation_fixed_T.string(), thermo_file.string());
			}
			else
			{
				preprocessor_species_without_transport->StatusOfThermodynamics(file_thermo_status.string());
				//preprocessor_species_without_transport->WriteThermodynamicsForMATLAB(file_thermo_matlab.string());
				preprocessor_species_without_transport->ReformulationOfThermodynamicsFixedIntermediateTemperature(file_thermo_reformulation_fixed_T.string(), thermo_file.string());
				preprocessor_species_without_transport->ReformulationOfThermodynamics(file_thermo_reformulation.string(), thermo_file.string());
			}
		}

		// Preprocessing transport data
		if (preprocess_transport_data_ == true)
			CheckForFatalError(preprocessor_species_with_transport->Fitting());

		// Write thermodynamic coefficients in a readable format
		if (write_ascii_thermodynamic_coefficients_ == true)
		{
			boost::filesystem::path file_ascii_thermodynamic_coefficients_ = path_output / "Thermodynamics_Coefficients.out";
			boost::filesystem::path file_ascii_thermodynamic_tables_ = path_output / "Thermodynamics_Tables.out";

			if (preprocess_transport_data_ == true)
			{
				preprocessor_species_with_transport->WriteThermodynamicCoefficientsOnASCIIFile(file_ascii_thermodynamic_coefficients_.string());
				preprocessor_species_with_transport->WriteThermodynamicTablesOnASCIIFile(file_ascii_thermodynamic_tables_.string());
			}
			else
			{
				preprocessor_species_without_transport->WriteThermodynamicCoefficientsOnASCIIFile(file_ascii_thermodynamic_coefficients_.string());
				preprocessor_species_without_transport->WriteThermodynamicTablesOnASCIIFile(file_ascii_thermodynamic_tables_.string());
			}
		}

		// Write fitting coefficients in ASCII format
		if (write_ascii_fitting_coefficients_ == true && preprocess_transport_data_ == true)
		{
			boost::filesystem::path file_ascii_fitting_coefficients_ = path_output / "TransportProperties_Coefficients.out";
			preprocessor_species_with_transport->WriteFittingCoefficientsOnASCIIFile(file_ascii_fitting_coefficients_.string());
		}

		// Write on XML files
		if (write_xml_files_ == true)
		{
			std::stringstream xml_string;
			xml_string << std::setprecision(8);
			xml_string.setf(std::ios::scientific);

			xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			xml_string << "<opensmoke version=\"0.1a\">" << std::endl;

			xml_string << "<Properties>" << std::endl;
			xml_string << "  <Author>" << author_name << "</Author>" << std::endl;
			xml_string << "  <Place>" << place_name << "</Place>" << std::endl;
			xml_string << "  <Date>" << preprocessing_date << "</Date>" << std::endl;
			xml_string << "  <Time>" << preprocessing_time << "</Time>" << std::endl;
			xml_string << "  <Comments>" << "\n" << OpenSMOKE::SplitStringIntoSeveralLines(comments, 80, "\t\r ") << "\n" << "</Comments>" << std::endl;
			xml_string << "</Properties>" << std::endl;

			// Thermodynamics properties
			if (preprocess_transport_data_ == true)
			{
				preprocessor_species_with_transport->WriteXMLFile(xml_string);
				if (species_bundling_ == true)
				{
					const double bundling_coefficient = 0.1;
					CheckForFatalError(preprocessor_species_with_transport->SpeciesBundling(xml_string, bundling_coefficient));
				}
			}
			else
			{
				preprocessor_species_without_transport->WriteXMLFile(xml_string);
			}

			// Kinetic mechanism (not available)
			xml_string << "<Kinetics type=\"none\" version=\"04-22-2013\">" << std::endl;
			xml_string << "</Kinetics>" << std::endl;

			xml_string << "</opensmoke>" << std::endl;

			// Write file
			boost::filesystem::path kinetics_xml = path_output / "kinetics.xml";
			std::ofstream fOutput;
			fOutput.open(std::string(kinetics_xml.string()).c_str(), std::ios::out);
			fOutput.setf(std::ios::scientific);
			fOutput << xml_string.str();
			fOutput.close();
		}
	}
	else
	{
		bool write_ascii_short_kinetic_summary_ = true;

		// Readers
		ThermoReader_CHEMKIN* thermoreader;
		OpenSMOKE::TransportReader< OpenSMOKE::TransportReaderPolicy_CHEMKIN<OpenSMOKE::TransportPolicy_CHEMKIN > >* transportreader;

		// Reading thermodynamic database
		thermoreader = new OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > >;
		CheckForFatalError(thermoreader->ReadFromFile(thermo_file.string()));

		// Preprocessing the kinetic mechanism
		PreProcessorKinetics_CHEMKIN preprocessor_kinetics(fLog);
		CheckForFatalError(preprocessor_kinetics.ReadFromASCIIFile(kinetics_file.string()));

		// Preprocessors
		PreProcessorSpecies_CHEMKIN* preprocessor_species_with_transport = nullptr;
		PreProcessorSpecies_CHEMKIN_WithoutTransport* preprocessor_species_without_transport = nullptr;

		if (preprocess_transport_data_ == true)	// With transport data
		{
			// Reading transport database
			transportreader = new OpenSMOKE::TransportReader< OpenSMOKE::TransportReaderPolicy_CHEMKIN<OpenSMOKE::TransportPolicy_CHEMKIN > >;
			CheckForFatalError(transportreader->ReadFromFile(transport_file.string()));

			// Preprocessing the thermodynamic and transport properties
			preprocessor_species_with_transport = new PreProcessorSpecies_CHEMKIN(*thermoreader, *transportreader, preprocessor_kinetics, fLog);
		}
		else
		{
			preprocessor_species_without_transport = new PreProcessorSpecies_CHEMKIN_WithoutTransport(*thermoreader, preprocessor_kinetics, fLog);
		}

		// Preprocessing thermodynamics
		if (preprocess_transport_data_ == true)
		{
			CheckForFatalError(preprocessor_species_with_transport->Setup());
			//	preprocessor_species_with_transport->ThermodynamicConsistency();
		}
		else
		{
			CheckForFatalError(preprocessor_species_without_transport->Setup());
			//	preprocessor_species_without_transport->ThermodynamicConsistency();
		}

		if (check_thermodynamics_ == true)
		{
			// Write the Thermodynamic Status on file
			boost::filesystem::path file_thermo_status = path_output / "Thermodynamics_Status.out";
			//boost::filesystem::path file_thermo_matlab = path_output / "Thermodynamics.mlab";
			boost::filesystem::path file_thermo_reformulation = path_output / "thermo.CHEMKIN.CKT";
			boost::filesystem::path file_thermo_reformulation_fixed_T = path_output / "thermo.CHEMKIN.fixedT.CKT";

			if (preprocess_transport_data_ == true)
			{
				preprocessor_species_with_transport->StatusOfThermodynamics(file_thermo_status.string());
				//preprocessor_species_with_transport->WriteThermodynamicsForMATLAB(file_thermo_matlab.string());
				preprocessor_species_with_transport->ReformulationOfThermodynamics(file_thermo_reformulation.string(), thermo_file.string());
				preprocessor_species_with_transport->ReformulationOfThermodynamicsFixedIntermediateTemperature(file_thermo_reformulation_fixed_T.string(), thermo_file.string());
			}
			else
			{
				preprocessor_species_without_transport->StatusOfThermodynamics(file_thermo_status.string());
				//preprocessor_species_without_transport->WriteThermodynamicsForMATLAB(file_thermo_matlab.string());
				preprocessor_species_without_transport->ReformulationOfThermodynamics(file_thermo_reformulation.string(), thermo_file.string());
				preprocessor_species_without_transport->ReformulationOfThermodynamicsFixedIntermediateTemperature(file_thermo_reformulation_fixed_T.string(), thermo_file.string());
			}
		}

		// Preprocessing transport data
		if (preprocess_transport_data_ == true)
			CheckForFatalError(preprocessor_species_with_transport->Fitting());

		// Read kinetics from file
		if (preprocess_transport_data_ == true)
		{
			CheckForFatalError(preprocessor_kinetics.ReadKineticsFromASCIIFile(preprocessor_species_with_transport->AtomicTable()));
		}
		else
		{
			CheckForFatalError(preprocessor_kinetics.ReadKineticsFromASCIIFile(preprocessor_species_without_transport->AtomicTable()));
		}

		// Write interpreted kinetics in old style
		if (preprocess_transport_data_ == true && write_ascii_files_old_style_ == true)
		{
			boost::filesystem::path idealgas_oldstyle = path_output / "idealgas.oldstyle.ascii";
			boost::filesystem::path reactions_oldstyle = path_output / "reactions.oldstyle.ascii";

			preprocessor_species_with_transport->WriteASCIIFileOldStyle(idealgas_oldstyle.string());
			preprocessor_kinetics.WriteASCIIFileOldStyle(reactions_oldstyle.string());
		}

		// Write interpreted kinetics in new style
		if (write_ascii_files_ == true)
		{
			boost::filesystem::path idealgas_newstyle = path_output / "idealgas.newstyle.ascii";	// TOREMOVE
			boost::filesystem::path reactions_newstyle = path_output / "reactions.newstyle.ascii";	// TOREMOVE

			if (preprocess_transport_data_ == true)
			{
				preprocessor_species_with_transport->WriteASCIIFile(idealgas_newstyle.string());	// TOREMOVE
				preprocessor_species_with_transport->WriteASCIIFile(reactions_newstyle.string());	// TOREMOVE
			}
			else
			{
				preprocessor_species_without_transport->WriteASCIIFile(idealgas_newstyle.string());	// TOREMOVE
				preprocessor_species_without_transport->WriteASCIIFile(reactions_newstyle.string());// TOREMOVE
			}
		}

		if (write_xml_files_ == true)
		{
			std::stringstream xml_string;
			xml_string << std::setprecision(8);
			xml_string.setf(std::ios::scientific);

			xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			xml_string << "<opensmoke version=\"0.1a\">" << std::endl;

			xml_string << "<Properties>" << std::endl;
			xml_string << "  <Author>" << author_name << "</Author>" << std::endl;
			xml_string << "  <Place>" << place_name << "</Place>" << std::endl;
			xml_string << "  <Date>" << preprocessing_date << "</Date>" << std::endl;
			xml_string << "  <Time>" << preprocessing_time << "</Time>" << std::endl;
			xml_string << "  <Comments>" << "\n" << OpenSMOKE::SplitStringIntoSeveralLines(comments, 80, "\t\r ") << "\n" << "</Comments>" << std::endl;
			xml_string << "</Properties>" << std::endl;

			// Thermodynamics and transport properties
			if (preprocess_transport_data_ == true)
			{
				preprocessor_species_with_transport->WriteXMLFile(xml_string);
				if (species_bundling_ == true)
				{
					const double bundling_coefficient = 0.1;
					CheckForFatalError(preprocessor_species_with_transport->SpeciesBundling(xml_string, species_bundling_));
				}
			}
			else
			{
				preprocessor_species_without_transport->WriteXMLFile(xml_string);
			}

			// Kinetic mechanism
			preprocessor_kinetics.WriteXMLFile(xml_string);

			// Chebyshev violation
			if (violation_allowed_chebyshev_polynomials_ == true)
				xml_string << "<ViolationChebyshev flag=\"true\"></ViolationChebyshev>" << std::endl;

			xml_string << "</opensmoke>" << std::endl;

			// Write file
			boost::filesystem::path kinetics_xml = path_output / "kinetics.xml";
			std::ofstream fOutput;
			fOutput.open(std::string(kinetics_xml.string()).c_str(), std::ios::out);
			fOutput.setf(std::ios::scientific);
			fOutput << xml_string.str();
			fOutput.close();
		}

		// Write kinetic summary in ASCII format
		if (write_ascii_short_kinetic_summary_ == true)
		{
			boost::filesystem::path file_ascii_short_kinetic_summary_ = path_output / "Kinetics_Summary.out";
			if (preprocess_transport_data_ == true)
				preprocessor_kinetics.WriteShortSummaryOnASCIIFile(file_ascii_short_kinetic_summary_.string(), *preprocessor_species_with_transport);
			else
				preprocessor_kinetics.WriteShortSummaryOnASCIIFile(file_ascii_short_kinetic_summary_.string(), *preprocessor_species_without_transport);
		}

		// Write thermodynamic coefficients in ASCII format
		if (write_ascii_thermodynamic_coefficients_ == true)
		{
			boost::filesystem::path file_ascii_thermodynamic_coefficients_ = path_output / "Thermodynamics_Coefficients.out";
			boost::filesystem::path file_ascii_thermodynamic_tables_ = path_output / "Thermodynamics_Tables.out";

			if (preprocess_transport_data_ == true)
			{
				preprocessor_species_with_transport->WriteThermodynamicCoefficientsOnASCIIFile(file_ascii_thermodynamic_coefficients_.string());
				preprocessor_species_with_transport->WriteThermodynamicTablesOnASCIIFile(file_ascii_thermodynamic_tables_.string());
			}
			else
			{
				preprocessor_species_without_transport->WriteThermodynamicCoefficientsOnASCIIFile(file_ascii_thermodynamic_coefficients_.string());
				preprocessor_species_without_transport->WriteThermodynamicTablesOnASCIIFile(file_ascii_thermodynamic_tables_.string());
			}
		}

		// Write fitting coefficients in ASCII format
		if (write_ascii_fitting_coefficients_ == true && preprocess_transport_data_ == true)
		{
			boost::filesystem::path file_ascii_fitting_coefficients_ = path_output / "TransportProperties_Coefficients.out";
			preprocessor_species_with_transport->WriteFittingCoefficientsOnASCIIFile(file_ascii_fitting_coefficients_.string());
		}

		// Re-write the kinetic mechanism in CHEMKIN format
		if (rewrite_chemkin_format_ == true)
		{
			boost::filesystem::path file_kinetics_CHEMKIN = path_output / "kinetics.CHEMKIN.CKI";

			if (preprocess_transport_data_ == true)
			{
				CheckForFatalError(preprocessor_kinetics.KineticsToCHEMKINFile(preprocessor_species_with_transport->AtomicTable(), file_kinetics_CHEMKIN));
			}
			else
			{
				CheckForFatalError(preprocessor_kinetics.KineticsToCHEMKINFile(preprocessor_species_without_transport->AtomicTable(), file_kinetics_CHEMKIN));
			}
		}

		// Write detailed kinetic mechanism data (requires the kinetics.xml file to be open)
		if (write_fitted_kinetic_constants_ == true ||
			write_reaction_tables_ == true ||
			write_reaction_strings_ == true ||
			sparsity_pattern_analysis_ == true ||
			write_fitted_kinetic_constants_ == true)
		{
			boost::property_tree::ptree ptree;
    			boost::property_tree::read_xml( (path_output / "kinetics.xml").string(), ptree );

			OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML;
			OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML;
			thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
			kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, ptree);

			OpenSMOKE::AnalyzerKineticMechanism<PreProcessorKinetics_CHEMKIN, OpenSMOKE::KineticsMap_CHEMKIN>
				analyzer(preprocessor_kinetics, *kineticsMapXML);

			if (write_reaction_tables_ == true)
			{
				boost::filesystem::path file_ascii_reaction_tables_ = path_output / "Reaction_Tables.out";
				analyzer.WriteReactionTablesOnASCIIFile(file_ascii_reaction_tables_.string(), reaction_tables_list_temperatures_);
			}
			if (write_fitted_kinetic_constants_ == true)
			{
				boost::filesystem::path file_ascii_fitted_kinetics_ = path_output / "Reaction_FittedKinetics.out";
				analyzer.WriteFittedInverseKineticConstantsOnASCIIFile(file_ascii_fitted_kinetics_.string());
			}
			if (write_fitted_chebyshev_polynomials_ == true)
			{
				boost::filesystem::path file_ascii_fitted_kinetics_ = path_output / "Reaction_FittedChebyshev.out";
				analyzer.WriteFittedChebyshevOnASCIIFile(file_ascii_fitted_kinetics_.string());
			}
			if (sparsity_pattern_analysis_ == true)
			{
				boost::filesystem::path file_ascii_fitted_kinetics_ = path_output / "SparsityPattern.out";
				analyzer.SparsityPatternAnalysis(file_ascii_fitted_kinetics_.string());
			}

			if (preprocess_transport_data_ == true && write_collision_rate_analysis_ == true)
			{
				OpenSMOKE::TransportPropertiesMap_CHEMKIN* transportMapXML;
				transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree);

				boost::filesystem::path file_ascii_reaction_tables_ = path_output / "CollisionRates.out";
				analyzer.WriteCollisionRatesOnASCIIFile(*transportMapXML, file_ascii_reaction_tables_.string(), collision_rate_analysis_list_temperatures_);
			}

			if (write_reaction_strings_ == true)
			{
				std::stringstream xml_string;
				xml_string << std::setprecision(8);
				xml_string.setf(std::ios::scientific);

				xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
				xml_string << "<opensmoke version=\"0.1a\">" << std::endl;

				xml_string << "<reaction-number>" << std::endl;
				xml_string << kineticsMapXML->NumberOfReactions() << std::endl;
				xml_string << "</reaction-number>" << std::endl;

				xml_string << "<reaction-names>" << std::endl;
				for (unsigned int k = 0; k < kineticsMapXML->NumberOfReactions(); k++)
				{
					std::string reaction_string;
					preprocessor_kinetics.reactions()[k].GetReactionString(thermodynamicsMapXML->NamesOfSpecies(), reaction_string);
					boost::erase_all(reaction_string, " ");
					xml_string << reaction_string << std::endl;;
				}
				xml_string << "</reaction-names>" << std::endl;

				if (store_chemkin_file_ == true)
					CheckForFatalError(preprocessor_kinetics.KineticsToXMLCHEMKINFile(xml_string));	
				
				xml_string << "</opensmoke>" << std::endl;

				boost::filesystem::path reaction_names_file = path_output / "reaction_names.xml";
				std::ofstream fOutput;
				fOutput.open(std::string(reaction_names_file.string()).c_str(), std::ios::out);
				fOutput.setf(std::ios::scientific);
				fOutput << xml_string.str();
				fOutput.close();
			}
		}

		// Preprocessing the surface kinetic mechanism

		if (preprocess_surface_kinetics_ == true)
		{
			PreProcessorKinetics_Surface_CHEMKIN preprocessor_surface_kinetics(fLog);
			CheckForFatalError(preprocessor_surface_kinetics.ReadFromASCIIFile(kinetics_surface_file.string()));
			PreProcessorSpecies_CHEMKIN_WithoutTransport* surface_preprocessor_species_without_transport;
			surface_preprocessor_species_without_transport = new PreProcessorSpecies_CHEMKIN_WithoutTransport(*thermoreader, preprocessor_kinetics.names_species(), preprocessor_surface_kinetics, fLog);
			CheckForFatalError(surface_preprocessor_species_without_transport->Setup());
			CheckForFatalError(preprocessor_surface_kinetics.ReadKineticsFromASCIIFile(surface_preprocessor_species_without_transport->AtomicTable(), preprocessor_kinetics));

			// Write on XML files
			if (write_xml_files_ == true)
			{
				std::stringstream xml_string;
				xml_string << std::setprecision(8);
				xml_string.setf(std::ios::scientific);

				xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
				xml_string << "<opensmoke version=\"0.1a\">" << std::endl;

				xml_string << "<Properties>" << std::endl;
				xml_string << "  <Author>" << author_name << "</Author>" << std::endl;
				xml_string << "  <Place>" << place_name << "</Place>" << std::endl;
				xml_string << "  <Date>" << preprocessing_date << "</Date>" << std::endl;
				xml_string << "  <Time>" << preprocessing_time << "</Time>" << std::endl;
				xml_string << "  <Comments>" << "\n" << OpenSMOKE::SplitStringIntoSeveralLines(comments, 80, "\t\r ") << "\n" << "  </Comments>" << std::endl;
				xml_string << "</Properties>" << std::endl;

				// Thermodynamics properties
				surface_preprocessor_species_without_transport->WriteXMLFile(xml_string);

				// Kinetic mechanism
				preprocessor_surface_kinetics.WriteXMLFile(xml_string);

				xml_string << "</opensmoke>" << std::endl;

				// Write file
				boost::filesystem::path kinetics_xml = path_output / "kinetics.surface.xml";
				std::ofstream fOutput;
				fOutput.open(std::string(kinetics_xml.string()).c_str(), std::ios::out);
				fOutput.setf(std::ios::scientific);
				fOutput << xml_string.str();
				fOutput.close();
			}

			// Write kinetic summary in ASCII format
			if (write_ascii_short_kinetic_summary_ == true)
			{
				boost::filesystem::path short_summary_kinetics = path_output / "SurfaceSummary.out";
				preprocessor_surface_kinetics.WriteShortSummaryOnASCIIFile(short_summary_kinetics.string(), *surface_preprocessor_species_without_transport);
			}

			if (write_reaction_strings_ == true)
			{

				OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN* thermodynamicsSurfaceMapXML;
				OpenSMOKE::KineticsMap_Surface_CHEMKIN* kineticsSurfaceMapXML;

				{
					boost::property_tree::ptree ptree;
    					boost::property_tree::read_xml( (path_output / "kinetics.surface.xml").string(), ptree );

					thermodynamicsSurfaceMapXML = new OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN(ptree);
					kineticsSurfaceMapXML = new OpenSMOKE::KineticsMap_Surface_CHEMKIN(*thermodynamicsSurfaceMapXML, ptree);
				}

				std::stringstream xml_string;
				xml_string << std::setprecision(8);
				xml_string.setf(std::ios::scientific);

				xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
				xml_string << "<opensmoke version=\"0.1a\">" << std::endl;
				xml_string << "<reaction-number>" << std::endl;
				xml_string << kineticsSurfaceMapXML->NumberOfReactions() << std::endl;
				xml_string << "</reaction-number>" << std::endl;
				xml_string << "<reaction-names>" << std::endl;
				for (unsigned int k = 0; k<kineticsSurfaceMapXML->NumberOfReactions(); k++)
				{
					std::string reaction_string;
					preprocessor_surface_kinetics.reactions()[0][k].GetReactionString(thermodynamicsSurfaceMapXML->NamesOfSpecies(), reaction_string);
					boost::erase_all(reaction_string, " ");
					xml_string << reaction_string << std::endl;;
				}
				xml_string << "</reaction-names>" << std::endl;
				xml_string << "</opensmoke>" << std::endl;

				boost::filesystem::path reaction_names_file = path_output / "surface_reaction_names.xml";
				std::ofstream fOutput;
				fOutput.open(std::string(reaction_names_file.string()).c_str(), std::ios::out);
				fOutput.setf(std::ios::scientific);
				fOutput << xml_string.str();
				fOutput.close();
			}
		}

		// Preprocessing the solid kinetic mechanism
		if (preprocess_solid_kinetics_ == true)
		{
			PreProcessorKinetics_Solid_CHEMKIN preprocessor_solid_kinetics(fLog);
			CheckForFatalError(preprocessor_solid_kinetics.ReadFromASCIIFile(kinetics_solid_file.string()));

			PreProcessorSpecies_CHEMKIN_WithoutTransport* solid_preprocessor_species_without_transport;
			solid_preprocessor_species_without_transport = new PreProcessorSpecies_CHEMKIN_WithoutTransport(*thermoreader, preprocessor_kinetics.names_species(), preprocessor_solid_kinetics, fLog);
			CheckForFatalError(solid_preprocessor_species_without_transport->Setup());
			CheckForFatalError(preprocessor_solid_kinetics.ReadKineticsFromASCIIFile(solid_preprocessor_species_without_transport->AtomicTable(), preprocessor_kinetics));

			// Write on XML files
			if (write_xml_files_ == true)
			{
				std::stringstream xml_string;
				xml_string << std::setprecision(8);
				xml_string.setf(std::ios::scientific);

				xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
				xml_string << "<opensmoke version=\"0.1a\">" << std::endl;

				xml_string << "<Properties>" << std::endl;
				xml_string << "  <Author>" << author_name << "</Author>" << std::endl;
				xml_string << "  <Place>" << place_name << "</Place>" << std::endl;
				xml_string << "  <Date>" << preprocessing_date << "</Date>" << std::endl;
				xml_string << "  <Time>" << preprocessing_time << "</Time>" << std::endl;
				xml_string << "  <Comments>" << "\n" << OpenSMOKE::SplitStringIntoSeveralLines(comments, 80, "\t\r ") << "\n" << "  </Comments>" << std::endl;
				xml_string << "</Properties>" << std::endl;

				// Thermodynamics properties
				solid_preprocessor_species_without_transport->WriteXMLFile(xml_string);

				// Kinetic mechanism
				preprocessor_solid_kinetics.WriteXMLFile(xml_string);

				xml_string << "</opensmoke>" << std::endl;

				// Write file
				boost::filesystem::path kinetics_xml = path_output / "kinetics.solid.xml";
				std::ofstream fOutput;
				fOutput.open(std::string(kinetics_xml.string()).c_str(), std::ios::out);
				fOutput.setf(std::ios::scientific);
				fOutput << xml_string.str();
				fOutput.close();
			}

			// Write kinetic summary in ASCII format
			if (write_ascii_short_kinetic_summary_ == true)
			{
				boost::filesystem::path short_summary_kinetics = path_output / "SolidSummary.out";
				preprocessor_solid_kinetics.WriteShortSummaryOnASCIIFile(short_summary_kinetics.string(), *solid_preprocessor_species_without_transport);
			}
		}

		// Preprocessing the liquid kinetic mechanism
		if (preprocess_liquid_kinetics_ == true)
		{
			PreProcessorKinetics_Liquid_CHEMKIN preprocessor_liquid_kinetics(fLog);
			CheckForFatalError(preprocessor_liquid_kinetics.ReadFromASCIIFile(kinetics_liquid_file.string()));

			PreProcessorSpecies_CHEMKIN_WithoutTransport* liquid_preprocessor_species_without_transport;
			liquid_preprocessor_species_without_transport = new PreProcessorSpecies_CHEMKIN_WithoutTransport(*thermoreader, preprocessor_kinetics.names_species(), preprocessor_liquid_kinetics, fLog);
			CheckForFatalError(liquid_preprocessor_species_without_transport->Setup());
			CheckForFatalError(preprocessor_liquid_kinetics.ReadKineticsFromASCIIFile(liquid_preprocessor_species_without_transport->AtomicTable(), preprocessor_kinetics));

			// Write on XML files
			if (write_xml_files_ == true)
			{
				std::stringstream xml_string;
				xml_string << std::setprecision(8);
				xml_string.setf(std::ios::scientific);

				xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
				xml_string << "<opensmoke version=\"0.1a\">" << std::endl;

				xml_string << "<Properties>" << std::endl;
				xml_string << "  <Author>" << author_name << "</Author>" << std::endl;
				xml_string << "  <Place>" << place_name << "</Place>" << std::endl;
				xml_string << "  <Date>" << preprocessing_date << "</Date>" << std::endl;
				xml_string << "  <Time>" << preprocessing_time << "</Time>" << std::endl;
				xml_string << "  <Comments>" << "\n" << OpenSMOKE::SplitStringIntoSeveralLines(comments, 80, "\t\r ") << "\n" << "  </Comments>" << std::endl;
				xml_string << "</Properties>" << std::endl;

				// Thermodynamics properties
				liquid_preprocessor_species_without_transport->WriteXMLFile(xml_string);

				// Kinetic mechanism
				preprocessor_liquid_kinetics.WriteXMLFile(xml_string);

				xml_string << "</opensmoke>" << std::endl;

				// Write file
				boost::filesystem::path kinetics_xml = path_output / "kinetics.liquid.xml";
				std::ofstream fOutput;
				fOutput.open(std::string(kinetics_xml.string()).c_str(), std::ios::out);
				fOutput.setf(std::ios::scientific);
				fOutput << xml_string.str();
				fOutput.close();
			}

			// Write kinetic summary in ASCII format
			if (write_ascii_short_kinetic_summary_ == true)
			{
				boost::filesystem::path short_summary_kinetics = path_output / "LiquidSummary.out";
				preprocessor_liquid_kinetics.WriteShortSummaryOnASCIIFile(short_summary_kinetics.string(), *liquid_preprocessor_species_without_transport);
			}

			if (write_reaction_strings_ == true)
			{
				OpenSMOKE::ThermodynamicsMap_Liquid_CHEMKIN* thermodynamicsLiquidMapXML;
				OpenSMOKE::KineticsMap_Liquid_CHEMKIN* kineticsLiquidMapXML;

				{
					boost::property_tree::ptree ptree;
    					boost::property_tree::read_xml( (path_output / "kinetics.liquid.xml").string(), ptree );

					thermodynamicsLiquidMapXML = new OpenSMOKE::ThermodynamicsMap_Liquid_CHEMKIN(ptree);
					kineticsLiquidMapXML = new OpenSMOKE::KineticsMap_Liquid_CHEMKIN(*thermodynamicsLiquidMapXML, ptree, 1);
				}

				std::stringstream xml_string;
				xml_string << std::setprecision(8);
				xml_string.setf(std::ios::scientific);

				xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
				xml_string << "<opensmoke version=\"0.1a\">" << std::endl;
				xml_string << "<reaction-number>" << std::endl;
				xml_string << kineticsLiquidMapXML->NumberOfReactions() << std::endl;
				xml_string << "</reaction-number>" << std::endl;
				xml_string << "<reaction-names>" << std::endl;
				for (unsigned int k = 0; k < kineticsLiquidMapXML->NumberOfReactions(); k++)
				{
					std::string reaction_string;
					preprocessor_liquid_kinetics.reactions(0)[k].GetReactionString(thermodynamicsLiquidMapXML->NamesOfSpecies(), reaction_string);
					boost::erase_all(reaction_string, " ");
					xml_string << reaction_string << std::endl;;
				}
				xml_string << "</reaction-names>" << std::endl;
				xml_string << "</opensmoke>" << std::endl;

				boost::filesystem::path reaction_names_file = path_output / "reaction_names.liquid.xml";
				std::ofstream fOutput;
				fOutput.open(std::string(reaction_names_file.string()).c_str(), std::ios::out);
				fOutput.setf(std::ios::scientific);
				fOutput << xml_string.str();
				fOutput.close();
			}
		}
	}

	fLog.close();

	OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_CHEMKIN_PreProcessor", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	return OPENSMOKE_SUCCESSFULL_EXIT;
}

