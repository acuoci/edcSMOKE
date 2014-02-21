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
PtrList<volScalarField> R_;
scalarField gammaStar(mesh.nCells(),0.);
scalarField mDotStar(mesh.nCells(),0.);
scalarField edc_coefficient(mesh.nCells(),0.);
{
	// Check the order of species
	for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		if (Y[i].name() != thermodynamicsMapXML->NamesOfSpecies()[i])
		{
			FatalError << "Kinetic schemes do not match..." << endl;
		}

	YStar_.resize(Y.size());
	R_.resize(Y.size());

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

Switch homogeneousReactions(thermo.lookup("homogeneousReactions"));

// Spark
Switch spark(thermo.lookup("spark"));
labelList ignCellsList;
scalar ignitionTime = 0.;
scalar ignitionDuration = 0.;
scalar ignitionTemperature = 0.;

{
	Foam::vector position;    
	scalar ignitionDiameter = 0.;

	// Read spark properties if detected
	if(spark == true)
	{
		position = (thermo.lookup("position"));
		ignitionTime = readScalar(thermo.lookup("time"));
		ignitionDuration = readScalar(thermo.lookup("duration"));
		ignitionTemperature = readScalar(thermo.lookup("temperature"));
		ignitionDiameter = readScalar(thermo.lookup("diameter"));   
	}

	label nIgnCells = 1;
	if(spark == true)
		#include "spark.H"
}

