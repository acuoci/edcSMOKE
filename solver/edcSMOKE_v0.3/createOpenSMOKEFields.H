// Read the kinetic scheme in XML format
OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>* thermodynamicsMapXML; 
OpenSMOKE::KineticsMap_CHEMKIN<double>* kineticsMapXML;

{	
	word kinetics("opensmoke");
	boost::filesystem::path path_kinetics = kinetics;

	rapidxml::xml_document<> doc;
	std::vector<char> xml_string;
	OpenSMOKE::OpenXMLFile(doc, xml_string,path_kinetics / "kinetics.xml");

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>(doc); 
	//transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc); 
	kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN<double>(*thermodynamicsMapXML, doc); 					
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
}
		
//- Species mass fractions
PtrList<volScalarField> YStar_;
PtrList<volScalarField> TStar_;	//TODO
PtrList<volScalarField> R_;
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

// Homogeneous reactions on/off
Switch homogeneousReactions(thermo.lookup("homogeneousReactions"));

// EDC Constants
scalar Ccsi = readScalar(thermo.lookup("Ccsi"));	// default: 2.1377
scalar Ctau = readScalar(thermo.lookup("Ctau"));	// default: 0.40825
scalar maxGammaSquared = readScalar(thermo.lookup("maxGammaSquared"));	// default: 0.75487766248
scalar relTolerance = readScalar(thermo.lookup("relTolerance"));	// default: 1.e-7
scalar absTolerance = readScalar(thermo.lookup("absTolerance"));	// default: 1.e-12
scalar numberOfResidenceTimes = readScalar(thermo.lookup("numberOfResidenceTimes"));	// default: 100.

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

