
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



HomogeneousODE::HomogeneousODE(
	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMapXML, 
	OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMapXML) :
	thermodynamicsMapXML_(thermodynamicsMapXML),
	kineticsMapXML_(kineticsMapXML)
{
	number_of_gas_species_ = thermodynamicsMapXML_.NumberOfSpecies();
	number_of_equations_ = number_of_gas_species_;

	ChangeDimensions(number_of_gas_species_, &omegaSurr_, true);
	ChangeDimensions(number_of_gas_species_, &omegaStar_, true);
	ChangeDimensions(number_of_gas_species_, &xStar_, true);
	ChangeDimensions(number_of_gas_species_, &cStar_, true);
	ChangeDimensions(number_of_gas_species_, &RStar_, true);	
}

int HomogeneousODE::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	unsigned int k=1;
	for(unsigned int i=1;i<=number_of_gas_species_;i++)
		omegaStar_[i] = y[k++];

	for(unsigned int i=1;i<=number_of_gas_species_;i++)
		omegaSurr_[i] = (omegaMean_[i] - omegaStar_[i]*gammaStar_)/(1.-gammaStar_);
	
	// Calculates the pressure and the concentrations of species
	thermodynamicsMapXML_.MoleFractions_From_MassFractions(xStar_, MWStar_, omegaStar_);
	cTotStar_ = P_Pa_/(PhysicalConstants::R_J_kmol * T_);
	rhoStar_ = cTotStar_*MWStar_;
	Product(cTotStar_, xStar_, &cStar_);

	// Calculates thermodynamic properties
	thermodynamicsMapXML_.SetTemperature(T_);
	thermodynamicsMapXML_.SetPressure(P_Pa_);
	
	// Calculates kinetics
	kineticsMapXML_.SetTemperature(T_);
	kineticsMapXML_.SetPressure(P_Pa_);
	kineticsMapXML_.ReactionEnthalpiesAndEntropies();
	kineticsMapXML_.ArrheniusKineticConstants();
	kineticsMapXML_.ReactionRates(cStar_);
	kineticsMapXML_.FormationRates(&RStar_);
	//HeatReleaseGas_ = kineticsMapXML_.HeatRelease(RStar_);

	// Recovering residuals
	for (unsigned int i=1;i<=number_of_gas_species_;++i)	
		dy[i] = thermodynamicsMapXML_.MW()[i]*RStar_[i]/rhoStar_ + mDotStar_*(omegaSurr_[i]-omegaStar_[i]);

	return 0;
}

int HomogeneousODE::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
	std::cout << t << std::endl;
	return 0;
}

