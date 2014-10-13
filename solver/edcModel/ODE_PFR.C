
// OpenSMOKE
#include "OpenSMOKE_Definitions.h"
#include <string>
#include <iostream>
#include <numeric>
#include <Eigen/Dense>

// Base classes
#include "thermo/ThermoPolicy_CHEMKIN.h"
#include "kinetics/ReactionPolicy_CHEMKIN.h"
#include "math/PhysicalConstants.h"
#include "math/OpenSMOKEUtilities.h"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"



ODE_PFR::ODE_PFR(
	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMapXML, 
	OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMapXML) :
	thermodynamicsMapXML_(thermodynamicsMapXML),
	kineticsMapXML_(kineticsMapXML)
{
	number_of_gas_species_ = thermodynamicsMapXML_.NumberOfSpecies();
	number_of_equations_ = number_of_gas_species_ + 1;

	ChangeDimensions(number_of_gas_species_, &omegaStar_, true);
	ChangeDimensions(number_of_gas_species_, &xStar_, true);
	ChangeDimensions(number_of_gas_species_, &cStar_, true);
	ChangeDimensions(number_of_gas_species_, &RStar_, true);

	checkMassFractions_ = false;	
}

int ODE_PFR::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	// Recover mass fractions
	if (checkMassFractions_ == true)
	{	for(unsigned int i=1;i<=number_of_gas_species_;++i)
			omegaStar_[i] = max(y[i], 0.);
	}
	else
	{
		for(unsigned int i=1;i<=number_of_gas_species_;++i)
			omegaStar_[i] = y[i];
	}
	// Recover temperature
	const double TStar_ = y[number_of_equations_];
	
	// Calculates the pressure and the concentrations of species
	thermodynamicsMapXML_.MoleFractions_From_MassFractions(xStar_, MWStar_, omegaStar_);
	cTotStar_ = P_Pa_/(PhysicalConstants::R_J_kmol * TStar_);
	rhoStar_ = cTotStar_*MWStar_;
	Product(cTotStar_, xStar_, &cStar_);

	// Calculates thermodynamic properties
	thermodynamicsMapXML_.SetTemperature(TStar_);
	thermodynamicsMapXML_.SetPressure(P_Pa_);
	thermodynamicsMapXML_.cpMolar_Mixture_From_MoleFractions(cpStar_, xStar_);
	cpStar_/=MWStar_;
	
	// Calculates kinetics
	kineticsMapXML_.SetTemperature(TStar_);
	kineticsMapXML_.SetPressure(P_Pa_);
	kineticsMapXML_.ReactionEnthalpiesAndEntropies();
	kineticsMapXML_.KineticConstants();
	kineticsMapXML_.ReactionRates(cStar_);
	kineticsMapXML_.FormationRates(&RStar_);
	const double QRStar_ = kineticsMapXML_.HeatRelease(RStar_);

	// Recovering residuals
	for (unsigned int i=1;i<=number_of_gas_species_;++i)	
		dy[i] = thermodynamicsMapXML_.MW()[i]*RStar_[i]/rhoStar_;
	
	const double Q = 0.; // radiation contribution
	dy[number_of_gas_species_+1] = (QRStar_ - Q)/(rhoStar_*cpStar_);
	
	return 0;
}

int ODE_PFR::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
	//std::cout << t << std::endl;
	return 0;
}

