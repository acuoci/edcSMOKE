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
|	License                                                           |
|                                                                         |
|   Copyright(C) 2019  Alberto Cuoci                                      |
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

#ifndef OpenSMOKE_KineticsMap_Liquid_CHEMKIN_CHEMKIN_H
#define OpenSMOKE_KineticsMap_Liquid_CHEMKIN_CHEMKIN_H

#include "KineticsMap.h"
#include "StoichiometricMap.h"



namespace OpenSMOKE
{
	enum TYPE_OF_LIQUID_KINETICS { TYPE_OF_LIQUID_KINETICS_CHEMKIN_CONVENTIONAL };

	//!  A class to efficiently evaluate the reaction and formation rates, to be used in production codes
	/*!
		 This class provides the tools to calculate in a very efficient way the reaction rates and the
		 formation rates. In order to ensure a good efficiency a map is created to store all the data
		 depending on the temperature. In this way they are recalculated only if strictly needed, i.e. only
		 if the temperature changes
	*/

	class KineticsMap_Liquid_CHEMKIN : public KineticsMap
	{

	public:	// Rate of Production Analysis (ROPA) utilities

		void RateOfProductionAnalysis(const bool iNormalize) const;
		void RateOfProductionAnalysis(std::ostream& fout) const;
		void RateOfProductionAnalysis(ROPA_Data& ropa) const;
		void RateOfProductionAnalysis(ROPA_Data& ropa, const double* rf, const double* rb) const;

	public:

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param thermo the thermodynamic map
		*@param doc xml file  
		*@param target liquid material index (from 1 to the number of materials)
		*/
		KineticsMap_Liquid_CHEMKIN(ThermodynamicsMap_Liquid_CHEMKIN& thermo, boost::property_tree::ptree& ptree, const unsigned int target);

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param thermo the thermodynamic map
		*@param doc xml file  
		*@param target liquid material name
		*/
		KineticsMap_Liquid_CHEMKIN(ThermodynamicsMap_Liquid_CHEMKIN& thermo, boost::property_tree::ptree& ptree, const std::string target);

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
		*@brief Returns the names of the species
		        Please remember that the returned vector is 0-index based
		*/
		const std::vector<std::string>& NamesOfSpecies() const { return thermodynamics_.names(); }

		/**
		*@brief Imports the kinetic schemes from a file in XML format
		*/
		virtual void ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree);

		/**
		*@brief Imports the kinetic schemes from a file in XML format
		*/
		void ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree, const std::string target_material_name);

		/**
		*@brief Imports the kinetic schemes from a file in XML format
		*/
		void ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree, const unsigned int target_material_index);

		/**
		*@brief Imports the list of species from a file in XML format
		*/
		virtual void ImportSpeciesFromXMLFile(boost::property_tree::ptree& ptree);

		/**
		*@brief Calculates the kinetic constants of the reverse reactions
		*/
		void FittedReverseKineticConstants(const double* x_bath, const unsigned int nparameters, Eigen::MatrixXd& fittedKineticParameters);

		/**
		*@brief Calculates the kinetic constants of the reverse reactions
		*/
		void FittedReverseKineticConstants(const unsigned int k, std::ostream& fOut, Eigen::MatrixXd& fittedKineticParameters);

		/*/**
		*@brief Write the data for the reaction tables
		*/
		void WriteKineticData(std::ostream& fOut, const unsigned int k, const double* c_bath, const double conversion_forward=1., const double conversion_backward=1.);

		/** 
		*@brief Write the data for the reaction tables as in gas-phase kinetics
		*/ 
		void WriteKineticData(std::ostream& fOut, const unsigned int k, const double* c_bath, const std::vector<double> list_of_temperatures, const double conversion_forward = 1., const double conversion_backward = 1.);
		
		/**
		*@brief Write the data for the reaction tables
		*/
		void WriteKineticData(std::ostream& fOut, const unsigned int k);

		/**
		*@brief Calculates the formation rates for all the species in the kinetic mechanism
		*@param Rgas the formation rates of all the gaseous species in [kmol/m3/s]
		*@param Rliquid the formation rates of all the liquid-phase spacies in [kmol/m3/s]
		*/
		void FormationRates(double* Rgas, double* Rliquid);

		/**
		*@brief Calculates the heat release
		*@param Rgas the formation rates of all the gaseous species in [kmol/m3/s] (as returned by the FormationRates function)
		*@param Rliquid the formation rates of all the liquid-phase spacies in [kmol/m3/s] (as returned by the FormationRates function)
		*@return the heat release due to the gas/liquid reactions in [J/m3/s]
		*/
		double HeatRelease(const double* Rgas, const double* RLiquid);

		/**
		*@brief Returns the forward reaction rates for all the reactions in the kinetic scheme
		*@param r the reaction rates of all the gaseous species in [kmol/m3/s]
		*/
		void GetForwardReactionRates(double* r);

		/**
		*@brief Returns the backward reaction rates for all the reactions in the kinetic scheme
		        If a reaction is irreversible, it returns zero
		*@param r the reaction rates of all the gaseous species in [kmol/m3/s]
		*/
		void GetBackwardReactionRates(double* r);

		/**
		*@brief Returns the net reaction rates in [kmol/m3/s]
		*/
		const std::vector<double>& GiveMeReactionRates();

		/**
		*@brief Returns the net reaction rates in [kmol/m3/s]
		*@param r the net reaction rates of all the reactions in [kmol/m3/s]
		*/
		void GiveMeReactionRates(double* r);

		/**
		*@brief Multiplies the required reaction rate (net) by a given correction coefficient
		*@param j index (1-based) of reaction rate to be corrected
		*@param Cc correction coefficient
		*/
		void CorrectReactionRate(const unsigned int j, const double Cc);

		/**
		*@brief Calculates the production and the destruction rates for all the species in the kinetic mechanism
		*/
		void ProductionAndDestructionRates(double* P, double* D);

		/**
		*@brief Calculates the reaction rates for all the reactions in the kinetic scheme
		*@param cGas the concentrations of all the gas-phase species in [kmol/m3/s]
		*@param cLiquid the concentrations of all the liquid-phase species in [kmol/m3/s]
		*/
		void ReactionRates(const double* cGas, const double* cLiquid);

		/**
		*@brief Returns the indices of the reversible reactions
		*/
		const std::vector<unsigned int>& IndicesOfReversibleReactions() const { return indices_of_reversible_reactions__; }

		/**
		*@brief Returns the indices of USRPROG reactions
		*/
		const std::vector<unsigned int>& IndicesOfUsrProgReactions() const { return indices_of_usrprog_reactions__; }

		/**
		*@brief Returns the indices of USRPROG reactions
		*/
		const std::vector<std::string>& LabelsOfUsrProgReactions() const { return labels_of_usrprog_reactions__; }

		/**
		*@brief Returns the stoichiometric map
		*@brief the stoichiometric map
		*/
		StoichiometricMap& stoichiometry() { return *stoichiometry_; }
                
        /**
		*@brief Calculates the reaction enthalpies and entropies (to be used for the kinetic constants)
		*/
		void ReactionEnthalpiesAndEntropies();

		/**
		*@brief Calculates the kinetic constants
		*/
		void KineticConstants();
        
        /**
		*@brief Return the frequency factor of a single reaction [kmol, m, s]
		*@param j index of reaction (starting from zero)
		*/
		double A(const unsigned int j) const { return std::exp(lnA__[j]); }

		/**
		*@brief Return the temperature exponent a single reaction
		*@param j index of reaction (starting from zero)
		*/
		double Beta(const unsigned int j) const { return Beta__[j]; }

		/**
		*@brief Return the activation temperature of a single reaction [K]
		*@param j index of reaction (starting from zero)
		*/
		double E_over_R(const unsigned int j) const { return E_over_R__[j]; }


	private:

		/**
		*@brief Calculates the modified Arrhenius constants
		*/
		const std::vector<double>& KArrheniusModified() const { return kArrheniusModified__; }

		/**
		*@brief Calculates the Arrhenius constants
		*/
		const std::vector<double>& KArrhenius() const { return kArrhenius__; }

	private:

		ThermodynamicsMap_Liquid_CHEMKIN& thermodynamics_;		//!< reference to the thermodynamics
		
		std::vector<double> c__;
		std::vector<double> aux_vector__;

		std::vector<unsigned int> indices_of_irreversible_reactions__;				//!< indices of irreversible reactions
		std::vector<unsigned int> indices_of_reversible_reactions__;					//!< indices of reversible reactions
		std::vector<unsigned int> indices_of_thermodynamic_reversible_reactions__;	//!< indices of reversible (thermodynamic) reactions
		std::vector<unsigned int> indices_of_explicitly_reversible_reactions__;		//!< indices of reversible (explicit) reactions

		unsigned int number_of_irreversible_reactions_;
		unsigned int number_of_reversible_reactions_;
		unsigned int number_of_thermodynamic_reversible_reactions_;
		unsigned int number_of_explicitly_reversible_reactions_;

		std::vector<unsigned int> indices_of_usrprog_reactions__;
		std::vector<std::string> labels_of_usrprog_reactions__;
		unsigned int number_of_usrprog_reactions_;

		std::vector<double> lnA__;								//!< frequency factors (log)
		std::vector<double> Beta__;								//!< temperature exponents
		std::vector<double> E_over_R__;							//!< activation temperatures
		std::vector<double> forward_kinetic_order__;			//!< global kinetic order for forward reactions

		std::vector<unsigned int> indices_of_reactions_needing_conversion_;

		std::vector<double> lnA_reversible__;					//!< frequency factors (log) for explicitly reversible reactions
		std::vector<double> Beta_reversible__;					//!< temperature exponents for explicitly reversible reactions
		std::vector<double> E_over_R_reversible__;				//!< activation temperatures for explicitly reversible reactions

		std::vector<double> changeOfMoles__;		//!< list of change of moles

		StoichiometricMap* stoichiometry_;			//!< pointer to the stoichiometry

		bool arrhenius_kinetic_constants_must_be_recalculated_;
		bool nonconventional_kinetic_constants_must_be_recalculated_;
		bool reaction_h_and_s_must_be_recalculated_;

		std::vector<double> reaction_s_over_R__;
		std::vector<double> reaction_h_over_RT__;
		std::vector<double> kArrheniusModified__;
		std::vector<double> kArrhenius__;
		std::vector<double> kArrhenius_reversible__;
		std::vector<double> uKeq__;

		std::vector<double> forwardReactionRates__;
		std::vector<double> reverseReactionRates__;
		std::vector<double> netReactionRates__;

		double Patm_over_RT_;
		double log_Patm_over_RT_;

		std::vector<unsigned int> isThermodynamicallyReversible__;		//!< vector containing the local index of thermodynamically reversible reactions
		std::vector<unsigned int> isExplicitlyReversible__;				//!< vector containing the local index of explicitly reversible reactions

		VectorReactionTags  type_of_reaction__;
		std::vector<unsigned int> local_family_index__;

		// Thermodynamic reversible reactions
		std::vector<double> delta_nu_gas_;

		// Type of kinetics
		TYPE_OF_LIQUID_KINETICS type_of_liquid_kinetics_;
	};
}

#include "KineticsMap_Liquid_CHEMKIN.hpp"

#endif /* OpenSMOKE_KineticsMap_Liquid_CHEMKIN_H */
