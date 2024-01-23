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
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_TransportPropertiesMap_CHEMKIN_H
#define OpenSMOKE_TransportPropertiesMap_CHEMKIN_H

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "TransportPropertiesMap.h"

// Grammar utilities
#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	//!  A class to efficiently evaluate the transport properties
	/*!
	This class provides the tools to calculate in a very efficient way the transport properties.
	In order to ensure a good efficiency a map is created to store all the data
	depending on the temperature. In this way they are recalculated only if strictly needed, i.e. only
	if the temperature changes.
	*/

	class TransportPropertiesMap_CHEMKIN : public TransportPropertiesMap
	{

	public:
		 
		enum StefanMaxwell_Skip_Policy { STEFANMAXWELL_SKIP_POLICY_MOST_ABUNDANT, STEFANMAXWELL_SKIP_POLICY_POLICY_FIXED };


	public:

		/**
		*@brief Creates a transport map for the evaluation of transport properties (obsolete TOREMOVE)
		*@param nSpecies number of species
		*/
		TransportPropertiesMap_CHEMKIN(const unsigned int nSpecies);

		/**
		*@brief Creates a transport map for the evaluation of transport properties
		*@param ptree file in XML format
		*/
		TransportPropertiesMap_CHEMKIN(boost::property_tree::ptree& ptree);

		/**
		*@brief Copy constructor
		*@param rhs the object to be copied in the current object
		*/
		TransportPropertiesMap_CHEMKIN(const TransportPropertiesMap_CHEMKIN& rhs);

		/**
		*@brief Default destructor
		*/
		~TransportPropertiesMap_CHEMKIN();

		/**
		*@brief Set the temperature at which the properties have to be evaluated
		*@param T the temperature value in K
		*/
		virtual void SetTemperature(const double& T);

		/**
		*@brief Set the pressure at which the properties have to be evaluated
		*@param P the pressure value in Pa
		*/
		virtual void SetPressure(const double& P);

		/**
		*@brief Sets the transport coefficient for the requested species
		*/
		virtual void SetCoefficients(const unsigned int i, const double* coefficients);

		/**
		*@brief Import the coefficients from a file in ASCII format (obsolete, TOREMOVE)
		*/
		virtual void ImportCoefficientsFromASCIIFile(std::istream& fInput);

		/**
		*@brief Import the coefficients from a file in ASCII format (obsolete, TOREMOVE)
		*@param fInput input stream
		*/
		virtual void ImportLennardJonesCoefficientsFromASCIIFile(std::istream& fInput);

		/**
		*@brief Import the raw transport coefficients from a file in ASCII format (obsolete, TOREMOVE)
		*@param fInput input stream
		*/
		virtual void ImportRawTransportParametersFromASCIIFile(std::istream& fInput);

		/**
		*@brief Import the coefficients from a file in XML format
		*/
		virtual void ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree);

		/**
		*@brief Import the species bundling coefficients from a file in XML format
		*/
		virtual void ImportSpeciesBundlingFromXMLFile(boost::property_tree::ptree& ptree, const double epsilon);

		/**
		*@brief Import the species from a file in XML format
		*/
		virtual void ImportSpeciesFromXMLFile(boost::property_tree::ptree& ptree);

		/**
		*@brief Import the viscosity model from a file in XML format
		*/
		virtual void ImportViscosityModelFromXMLFile(boost::property_tree::ptree& ptree);

		/**
		*@brief Returns the reduced collision integral Omega11* according to
		        the regression fit reported in Chen et al., Comb. and Flame 186, p. 208 (2017)
		*@param TStar reduced temperature (dimensionless) calculated as TStar=kb*T/eps
		*/
		double Omega11(const double TStar);

		/**
		*@brief Returns the reduced collision integral Omega11* according to
		        a 2D interpolation table
		*@param TStar reduced temperature (dimensionless) calculated as TStar=kb*T/eps
		*@param DStar reduced dipole moment (see CHEMKIN Theory Guide)
		*/
		double Omega11(const double TStar, const double DStar);

		/**
		*@brief Returns the mass (in kg) of species (0-based index)
		*/
		const std::vector<double>& mu() const { return mu_; }

		/**
		*@brief Returns the collision diameters (in m) of species (0-based index)
		*/
		const std::vector<double>& sigma() const { return sigma_; }

		/**
		*@brief Returns the scaled well depths (in K) of species (0-based index)
		*/
		const std::vector<double>& epsilon_over_kb() const { return epsilon_over_kb_; }

		/**
		*@brief Returns the raw transport parameters
		*/
		const std::vector<std::vector<double>>& raw_transport_parameters() const { return raw_transport_parameters_; }

	protected:

		/**
		*@brief TODO (This is just a test subroutine!)
		*/
		void Test(const int nLoops, const double& T, int* index);

	public:

		/**
		*@brief [Bertrand Naud] Solves the thermodiffusion problem (1+M model)
		*/
		virtual void solve_thermodiffusion_onePlusM(const double* x, const double* cp, const unsigned int Mmain, const std::vector<unsigned int>& kreorder);

		/**
		*@brief Combines the species thermal conductivities to calculate the mixture thermal conductivity
		*/
		virtual double lambdaMix(const double* x);

		/**
		*@brief [Bertrand Naud] Calculates the mixture thermal conductivity (1+M model)
		*/
		virtual double lambdaMult(const double* x, const unsigned int Mmain, const std::vector<unsigned int>& kreorder);

		/**
		*@brief Combines the species dynamic viscosities to calculate the mixture dynamic viscosity
		*/
		virtual double etaMix(const double* x);

		/**
		*@brief Combines the species dynamic viscosities to calculate the mixture dynamic viscosity (1+M model)
		*/
		virtual double etaMix(const double* x, const unsigned int Mmain, const std::vector<unsigned int>& kreorder);

		/**
		*@brief Combines the species mass diffusion coefficients to calculate the mixture mass diffusion coefficients
		*/
		virtual void gammaMix(double* gammamix, const double* x);

		/**
		*@brief [Bertrand Naud] Calculates the mixture mass diffusion coefficients (1+M model)
		*/
		virtual void gammaMult(double* gammamult, const double* x, const unsigned int Mmain, const std::vector<unsigned int>& kreorder);

		/**
		*@brief Combines the species mass diffusion coefficients to calculate the mixture mass diffusion coefficients using the bundling algorithm
		*/
		virtual void bundling_gammaMix(double* gammamix, const double* x);

		/**
		*@brief Combines the species thermal diffusion coefficients to calculate the mixture thermal diffusion coefficients
		*/
		virtual void tetaMix(double* tetamix, const double* x);

		/**
		*@brief [Bertrand Naud] Calculates the mixture thermal diffusion coefficients (1+M model)
		*/
		virtual void tetaMult(double* tetamult, const double* x, const unsigned int Mmain, const std::vector<unsigned int>& kreorder);

		/**
		*@brief Calculates the thermal diffusion coefficients according to the Kuo's model
		*@param GammaSoret thermal diffusion coefficients (kg/m/s)
		*@param T temperature
		*@param X vector of mole fractions
		*@param X vector of mass fractions
		*/
		void ThermalDiffusionCoefficients(double* GammaSoret, const double T, const double* X, const double* Y);

		/**
		*@brief Combines the planck mean absorption coefficients of relevant species, according to their mole fractions
		Returns the Planck mean absorption coefficient of the mixture in [1/m]
		*/
		virtual double kPlanckMix(const double* moleFractions);

		/**
		*@brief Calculates the thermal conductivities for all the species
		*/
		inline virtual void lambda();

		/**
		*@brief Calculates the dynamic viscosities for all the species
		*/
		inline virtual void eta();

		/**
		*@brief Calculates the dynamic viscosities for 1+M species
		*/
		inline virtual void eta(const unsigned int Mmain, const std::vector<unsigned int>& kreorder);

		/**
		*@brief Calculates the mass diffusion coefficients for all the species
		*/
		inline virtual void gamma();

		/**
		*@brief Calculates the mass diffusion coefficients for all the species using the bundling algorithm
		*/
		inline virtual void bundling_gamma();

		/**
		*@brief Calculates the thermal diffusion coefficients for all the species
		*/
		inline virtual void teta();

		/**
		*@brief Return the vector of species (1-index based) for which the Soret effect is active
		*/
		const std::vector<unsigned int>& iThermalDiffusionRatios() const { return iThermalDiffusionRatios_; }

		/**
		*@brief [Bertrand Naud] Reset the vector of species (1-index based) for which the Soret effect is active to all species (1+M model)
		*/
		inline virtual void reset_iThermalDiffusionRatios();

		/**
		*@brief Calculates the collision rate constant for bimolecular reactions [m3/kmol/s]
		*/
		double kCollision(const unsigned int i, const unsigned int k, const double T);

		/**
		*@brief Returns true is raw transport parameters are available (strictly needed for applying multicomponent diffusion models)
		*/
		bool is_raw_transport_parameters_available() const { return is_raw_transport_parameters_available_; }
	
	public:


		/**
		*@brief Calculates the ordinary diffusion coefficients according to the Stefan-Maxwell theory
		*@param GammaBar matrix of ordinary diffusion coefficients (including zero row and column corresponding to species to skip)
		*@param x vector of mole fractions
		*/
		void StefanMaxwellGammaBar(Eigen::MatrixXd& GammaBar, const double* x);


		/**
		*@brief Calculates the Fick diffusion fluxes according to the Stefan-Maxwell theory
		*@param jsm fick diffusion fluxes (kg/m2/s)
		*@param GammaBarL matrix of ordinary diffusion coefficients, left side (m2/s)
		*@param rhoL density, left side (kg/m3)
		*@param GammaBarR matrix of ordinary diffusion coefficients, right side (m2/s)
		*@param rhoR density, right side (kg/m3)
		*@param nablaY gradients of mass fractions (1/m)
		*/
		void StefanMaxwellFickFluxes(double* jsm, const Eigen::MatrixXd& GammaBarL, const double rhoL, const Eigen::MatrixXd& GammaBarR, const double rhoR, const double* nablaY);


		/**
		*@brief Calculates the Soret diffusion fluxes according to the Stefan-Maxwell theory
		*@param jsm Soret diffusion fluxes (kg/m2/s)
		*@param rhoL density, left side (kg/m3)
		*@param rhoR density, right side (kg/m3)
		*@param GammaTL thermal diffusion coefficient, left side (kg/m/s)
		*@param GammaTR thermal diffusion coefficient, right side (kg/m/s)
		*@param TL temperature, left side (K)
		*@param TR temperature, right side (K)
		*@param nablaT_over_T temperature gradient divided by temperature (1/m)
		*/
		void StefanMaxwellSoretFluxes(double* jsm,
			const double rhoL, const double rhoR,
			const Eigen::VectorXd& GammaTL, const Eigen::VectorXd& GammaTR, 
			const double TL, const double TR, const double nablaT);


		/**
		*@brief Calculates the mass diffusion fluxes according to the Stefan-Maxwell theory (Soret effect is included)
		*@param jsm mass diffusion fluxes (kg/m2/s)
		*@param GammaBarL matrix of ordinary diffusion coefficients, left side (m2/s)
		*@param rhoL density, left side (kg/m3)
		*@param GammaBarL matrix of ordinary diffusion coefficients, right side (m2/s)
		*@param rhoL density, right side (kg/m3)
		*@param nablaY gradients of mass fractions (1/m)
		*@param GammaTL thermal diffusion coefficient, left side (kg/m/s)
		*@param GammaTR thermal diffusion coefficient, right side (kg/m/s)
		*@param TL temperature, left side (K)
		*@param TR temperature, right side (K)
		*@param nablaT_over_T temperature gradient divided by temperature (1/m)
		*/
		void StefanMaxwellFluxes(double* jsm,
			const Eigen::MatrixXd& GammaBarL, const double rhoL,
			const Eigen::MatrixXd& GammaBarR, const double rhoR, const double* nablaY,
			const Eigen::VectorXd& GammaTL, const Eigen::VectorXd& GammaTR, const double TL, const double TR, const double nablaT);


		/**
		*@brief Sets the policy to choose the species to skip in Stefan-Maxwell calculations
		*@param policy policy to be adopted: most-abundant | last | name
		*@param names of species
		*/
		void StefanMaxwellSkipPolicy(const std::string policy, const std::vector<std::string>& names);

		/**
		*@brief Sets the Stefan-Maxwell options from dictionary
		*@param dictionary dictionary
		*@param names of species
		*/
		void SetupStefanMaxwellFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, const std::vector<std::string>& names);


	private:

		/**
		*@brief Allocates the memory
		*/
		void MemoryAllocation();

		/**
		*@brief Pre-calculates useful data
		*/
		void CompleteInitialization();

		/**
		*@brief Copies the data from another transport map (used by copy constructors)
		*/
		void CopyFromMap(const TransportPropertiesMap_CHEMKIN& rhs);

	private:

		PhysicalConstants::OpenSMOKE_GasMixture_Viscosity_Model viscosity_model;

		double* M;				//!< molecular weights
		double* fittingLambda;			//!< fitting coefficients for the thermal conductivities
		double* fittingEta;			//!< fitting coefficients for the dynamic viscosities
		double* fittingTeta;			//!< fitting coefficients for the thermal diffusion coefficients
		double* fittingGamma;			//!< fitting coefficients for the mass diffusion coefficients
		double* fittingGammaSelfDiffusion;	//!< fitting coefficients for the mass diffusion coefficients

		double* MWRatio1over4;		//!< auxiliary vector for the dynamic viscosity calculation
		double* phi_eta_sup;			//!< auxiliary vector for the dynamic viscosity calculation
		double* phi_eta_inf;			//!< auxiliary vector for the dynamic viscosity calculation

		std::vector<std::vector<double>>	mMWRatio1over4;	//!< auxiliary vector for the dynamic viscosity calculation (for 1+M model, it can be improved)
		std::vector<std::vector<double>>	mphi_eta_sup;		//!< auxiliary vector for the dynamic viscosity calculation (for 1+M model, it can be improved)

		double* sqrtEta;			//!< auxiliary vector for the dynamic viscosity calculation
		double* usqrtEta;			//!< auxiliary vector for the dynamic viscosity calculation
		Eigen::VectorXd sumK;			//!< auxiliary vector for the dynamic viscosity calculation

		double* sqrtMWRatio_inf;		//!< auxiliary vector for the dynamic viscosity calculation
		double* sqrtMWRatio_sup;		//!< auxiliary vector for the dynamic viscosity calculation

		double* sqrtMW;				//!< auxiliary vector for the dynamic viscosity calculation

		double* sum_diffusion_coefficients;	//!< auxiliary vector used for the calculation of the mass diffusion coefficients
		double* x_corrected;			//!< auxiliary vector used for the calculation of the mass diffusion coefficients

		unsigned int count_species_thermal_diffusion_ratios_;   //!< number of species for which the thermal diffusion coefficients are evaluated 
		std::vector<unsigned int> iThermalDiffusionRatios_;	//!< indices of species tracked for the thermal diffusion coefficients

		static const double threshold_;
		double sum_threshold_;

		bool temperature_lambda_must_be_recalculated_;
		bool temperature_eta_must_be_recalculated_;
		bool temperature_gamma_must_be_recalculated_;
		bool temperature_teta_must_be_recalculated_;
		bool pressure_gamma_must_be_recalculated_;
		bool Mijkl_must_be_recalculated_;

		// Indices of relevant species (1-based)
		unsigned int index_H2O_;
		unsigned int index_CO_;
		unsigned int index_CH4_;
		unsigned int index_CO2_;
		unsigned int index_NH3_;
		unsigned int index_NO_;
		unsigned int index_N2O_;

		// Bundling
		unsigned int bundling_number_groups_;
		std::vector<unsigned int> bundling_reference_species_;
		std::vector< std::vector<unsigned int> > bundling_groups_;
		std::vector<unsigned int> bundling_species_group_;
		std::vector<double> bundling_sum_diffusion_coefficients_;
		std::vector<double> bundling_sum_x_groups_;

		double* bundling_fittingGammaSelfDiffusion_;
		double* bundling_fittingGamma_;
		double* bundling_gammaSpecies_;
		double* bundling_gammaSpeciesSelfDiffusion_;

		// Lennard-Jones parameters
		bool is_lennard_jones_available_;			//!< true if the Lennard-Jones coefficients are available
		std::vector<double> mu_;				//!< species masses (in kg)
		std::vector<double> sigma_;				//!< species collision diameters (in m)
		std::vector<double> epsilon_over_kb_;		//!< species scaled well depths (in K)

		// Raw transport coefficients
		bool is_raw_transport_parameters_available_;				//!< true if the raw transport parameters are available
		std::vector<std::vector<double>> raw_transport_parameters_;		//!< raw transport coefficients

		// Stefan-Maxwell theory
		unsigned int stefanmaxwell_skip_;
		StefanMaxwell_Skip_Policy stefanmaxwell_skip_policy_;
		double* M_0511_;
		double* M_0489_;

		std::vector<std::vector<double>>	epsij_over_kb_;	//!< auxiliary vector for the thermodiffusion matrix definition (1+M model)
		std::vector<double>			fParker_298_;		//!< auxiliary vector for the thermodiffusion matrix definition (1+M model)
	};

	class Grammar_StefanMaxwell : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{

	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SwapPolicy",
				OpenSMOKE::SINGLE_STRING,
				"Type of swap policy: most-abundant (default) | last | species name",
				true));
		}
	};
}

#include "TransportPropertiesMap_CHEMKIN.hpp"

#endif // OpenSMOKE_TransportPropertiesMap_CHEMKIN_H
