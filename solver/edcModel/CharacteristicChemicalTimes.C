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

// OpenSMOKE
#include "OpenSMOKE_Definitions.h"
#include <string>
#include <iostream>
#include <numeric>
#include <Eigen/Dense>

// Base classes
#include "kernel/thermo/ThermoPolicy_CHEMKIN.h"
#include "kernel/kinetics/ReactionPolicy_CHEMKIN.h"
#include "math/PhysicalConstants.h"
#include "math/OpenSMOKEUtilities.h"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"

#if OPENSMOKE_USE_MKL == 1
#include "mkl.h"
#include "mkl_lapacke.h"
#endif

CharacteristicChemicalTimes::CharacteristicChemicalTimes(
	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXML, 
	OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapXML) :
	thermodynamicsMapXML_(thermodynamicsMapXML),
	kineticsMapXML_(kineticsMapXML)
{
	tauChem_threshold_ = 0.10;	
	type_ = CHEMICAL_TIMES_FORMATION_RATES;
	lapack_mode_ = false;

	ns_ = thermodynamicsMapXML_.NumberOfSpecies();
	nr_ = kineticsMapXML_.NumberOfReactions();
	ne_ = thermodynamicsMapXML_.elements().size();
	n_  = ns_- ne_;

	// Memory allocation
	ChangeDimensions(ns_, &x_, true);
	ChangeDimensions(ns_, &c_, true);
	ChangeDimensions(ns_, &R_, true);
	ChangeDimensions(nr_, &rf_, true);

	// Local variables
	lambda_real_.resize(ns_);
	lambda_imag_.resize(ns_);
	lambda_mod_.resize(ns_);
	lambda_real_cleaned_.resize(ns_);
	vr_.resize(ns_*ns_);
	vl_.resize(ns_*ns_);

	// Memory allocation: Jacobian matrices
	Jc_.resize(ns_, ns_);

	//
	kineticsMapXML_.stoichiometry().GetSumOfStoichiometricCoefficientsOfProducts(sum_nu_);
}

void CharacteristicChemicalTimes::SetThresholdChemicalTime(const double threshold)
{
	tauChem_threshold_ = threshold;
}

void CharacteristicChemicalTimes::SetThresholdTemperatureChemicalTime(const double threshold)
{
	tauChem_temperature_threshold_ = threshold;
}

void CharacteristicChemicalTimes::SetType(const CharacteristicChemicalTimesType type)
{
	type_ = type;
}

void CharacteristicChemicalTimes::SetLapackMode(const bool flag)
{
	lapack_mode_ = flag;

	#if OPENSMOKE_USE_MKL != 1
	lapack_mode_ = false;
	#endif
}

double CharacteristicChemicalTimes::CalculateCharacteristicChemicalTime(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega)
{
	// If the temperature is below the threshold, no calculation is performed
	if (T <= tauChem_temperature_threshold_)
	{
		return tauChem_threshold_;
	}
	// Estimation of chemical time is carried out only if the temperature is above the threshold
	else
	{
		if (type_ == CHEMICAL_TIMES_FORMATION_RATES)
			return FromFormationRates(T, P, omega);
		else if (type_ == CHEMICAL_TIMES_REACTION_RATES)
			return FromReactionRates(T, P, omega);
		else if (type_ == CHEMICAL_TIMES_EIGENVALUES)
			return FromEigenValueAnalysis(T, P, omega);

		return 0;
	}
}

double CharacteristicChemicalTimes::FromFormationRates(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega)
{
	const double dC_over_dt_threshold = 1e-16;

	double MW;
	thermodynamicsMapXML_.MoleFractions_From_MassFractions(x_.GetHandle(), MW, omega.GetHandle());

	const double cTot = P/(PhysicalConstants::R_J_kmol * T);
	const double rho = cTot*MW;
	Product(cTot, x_, &c_);

	// Calculates reaction rates
	thermodynamicsMapXML_.SetTemperature(T);
	thermodynamicsMapXML_.SetPressure(P);
	kineticsMapXML_.SetTemperature(T);
	kineticsMapXML_.SetPressure(P);
	kineticsMapXML_.KineticConstants();
	kineticsMapXML_.ReactionRates(c_.GetHandle());
	kineticsMapXML_.FormationRates(R_.GetHandle());

	double tau = 0.;
	for(unsigned int j=0;j<ns_;j++)
	{
		// Characteristic time of species
		const double dC_over_dt = std::fabs(R_[j+1]);
		const double tau_local = (dC_over_dt > dC_over_dt_threshold) ? c_[j+1]/dC_over_dt : 0.;

		// Selection of maximum characteristic time
		tau = (tau_local > tau) ? tau_local : tau;
	}

	// The characteristic chemical time is always constrained between 0<tau<tau_threshold
	if (tau == 0.)	tau = tauChem_threshold_;
	tau = min(tau, tauChem_threshold_);

	return tau;
}

double CharacteristicChemicalTimes::FromReactionRates(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega)
{
	const double dC_over_dt_threshold = 1e-16;

	double MW;
	thermodynamicsMapXML_.MoleFractions_From_MassFractions(x_.GetHandle(), MW, omega.GetHandle());

	const double cTot = P/(PhysicalConstants::R_J_kmol * T);
	const double rho = cTot*MW;
	Product(cTot, x_, &c_);

	// Calculates reaction rates
	thermodynamicsMapXML_.SetTemperature(T);
	thermodynamicsMapXML_.SetPressure(P);
	kineticsMapXML_.SetTemperature(T);
	kineticsMapXML_.SetPressure(P);
	kineticsMapXML_.KineticConstants();
	kineticsMapXML_.ReactionRates(c_.GetHandle());
	kineticsMapXML_.GetForwardReactionRates(rf_.GetHandle());

	double sum = 0.;
	for(unsigned int j=0;j<nr_;j++)
		sum += rf_[j+1]*sum_nu_(j);

	// The characteristic chemical time is always constrained between 0<tau<tau_threshold
	double tau = tauChem_threshold_;	
	if (sum > 0.)	
		tau = nr_*cTot/sum;
	tau = min(tau, tauChem_threshold_);

	return tau;
}

double CharacteristicChemicalTimes::FromEigenValueAnalysis(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega)
{
	T_ = T;
	P_ = P;

	double MW;
	thermodynamicsMapXML_.MoleFractions_From_MassFractions(x_.GetHandle(), MW, omega.GetHandle());

	const double cTot = P/(PhysicalConstants::R_J_kmol * T);
	const double rho = cTot*MW;
	Product(cTot, x_, &c_);

	// Populate Jacobian matrix
	NumericalJacobian(c_, Jc_);

	// Calculate eigenvectors/eigenvalues
	if (lapack_mode_ == false)
	{
		Eigen::EigenSolver<Eigen::MatrixXd> eigenJc(Jc_);
		
		for (unsigned int i = 0; i < ns_; i++)
		{
			lambda_real_[i] = eigenJc.eigenvalues()[i].real();
			lambda_imag_[i] = eigenJc.eigenvalues()[i].imag();
			lambda_mod_[i] = std::sqrt(lambda_real_[i] * lambda_real_[i] + lambda_imag_[i] * lambda_imag_[i]);
		}
	}
	#if OPENSMOKE_USE_MKL == 1
	else
	{
		int info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'V', 'V', ns_, Jc_.data(), ns_, lambda_real_.data(), lambda_imag_.data(), vl_.data(), ns_, vr_.data(), ns_);

		if (info != 0)
			OpenSMOKE::FatalErrorMessage("OnTheFlyCEMA::Calculate: LAPACKE_dgeev failure");

		for (unsigned int i = 0; i < ns_; i++)
			lambda_mod_[i] = std::sqrt(lambda_real_[i] * lambda_real_[i] + lambda_imag_[i] * lambda_imag_[i]);
	}
	#endif

	// Search conservative modes
	SearchConservativeModes(lambda_real_, lambda_mod_, lambda_real_cleaned_);

	// Search for minimum eigenvalue (modulus)
	double minimum_lambda = 1e16;
	for(unsigned int i=0;i<n_;i++)
		minimum_lambda = ( lambda_real_cleaned_[i] < minimum_lambda) ? lambda_real_cleaned_[i] : minimum_lambda;

	double tau = tauChem_threshold_;
	if (minimum_lambda > 0.)
		tau = 1./minimum_lambda;
	tau = min(tau, tauChem_threshold_);

	return tau;
}

void CharacteristicChemicalTimes::Equations(const OpenSMOKE::OpenSMOKEVectorDouble& c, OpenSMOKE::OpenSMOKEVectorDouble& dc_over_dt)
{
	thermodynamicsMapXML_.SetTemperature(T_);
	thermodynamicsMapXML_.SetPressure(P_);
	kineticsMapXML_.SetTemperature(T_);
	kineticsMapXML_.SetPressure(P_);
	kineticsMapXML_.KineticConstants();
	kineticsMapXML_.ReactionRates(c.GetHandle());
	kineticsMapXML_.FormationRates(dc_over_dt.GetHandle());
}

void CharacteristicChemicalTimes::NumericalJacobian(const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::MatrixXd& J)
{
	// Calculated as suggested by Buzzi (private communication)

	const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
	const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);			
	const double TOLR = 100. * OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT;
	const double TOLA = 1.e-10;

	OpenSMOKE::OpenSMOKEVectorDouble y_plus = y;
	OpenSMOKE::OpenSMOKEVectorDouble dy_original(y.Size());
	OpenSMOKE::OpenSMOKEVectorDouble dy_plus(y.Size());

	Equations(y, dy_original);

	// Derivatives with respect to y[kd]
	for(int kd=1;kd<=y.Size();kd++)
	{
		double hf = 1.e0;
		double error_weight = 1./(TOLA+TOLR*std::fabs(y[kd]));
		double hJ = ETA2 * std::fabs(std::max(y[kd], 1./error_weight));
		double hJf = hf/error_weight;
		hJ = std::max(hJ, hJf);
		hJ = std::max(hJ, ZERO_DER);

		// This is what is done by Buzzi
		double dy = std::min(hJ, 1.e-3 + 1e-3*std::fabs(y[kd]));
		double udy = 1. / dy;
		y_plus[kd] += dy;
		Equations(y_plus, dy_plus);

		for(int j=1;j<=y.Size();j++)
			J(j-1,kd-1) = (dy_plus[j]-dy_original[j]) * udy;

		y_plus[kd] = y[kd];
	}
}

void CharacteristicChemicalTimes::SearchConservativeModes(const std::vector<double>& lambda_real, const std::vector<double>& lambda_mod, std::vector<double>& lambda_real_cleaned)
{
	// Select elements to be removed
	std::vector<size_t> sorted_indices(ns_);
	sorted_indices = OpenSMOKE::SortAndTrackIndicesIncreasing(lambda_mod);
	std::vector<size_t> elements_to_remove = std::vector<size_t>(sorted_indices.begin(), sorted_indices.begin() + ne_);
	std::sort(elements_to_remove.begin(), elements_to_remove.end());
	std::reverse(elements_to_remove.begin(), elements_to_remove.end());

	// Removing elements corresponding to conservation constraints
	lambda_real_cleaned = lambda_real;
	for (unsigned int i = 0; i<elements_to_remove.size(); i++)
		lambda_real_cleaned.erase(lambda_real_cleaned.begin() + elements_to_remove[i]);

	// Absolute values
	for (unsigned int i = 0; i<n_; i++)
		lambda_real_cleaned[i] = std::fabs(lambda_real_cleaned[i]);
}

