/*----------------------------------------------------------------------*\
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
|	License                                                               |
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

#ifndef OpenSMOKE_ThermodynamicsMap_Solid_CHEMKIN_H
#define OpenSMOKE_ThermodynamicsMap_Solid_CHEMKIN_H

#include "ThermodynamicsMap.h"

namespace OpenSMOKE
{
	//!  A class to efficiently evaluate the thermodynamic properties of gas/solid mixtures
	/*!
		 This class provides the tools to calculate in a very efficient way the thermodynamic properties.
		 In order to ensure a good efficiency a map is created to store all the data
		 depending on the temperature. In this way they are recalculated only if strictly needed, i.e. only
		 if the temperature changes.
	*/

	class ThermodynamicsMap_Solid_CHEMKIN : public ThermodynamicsMap
	{

	public:

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param doc file in XML format
		*/
		ThermodynamicsMap_Solid_CHEMKIN(boost::property_tree::ptree& ptree);

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
		*@brief Returns the molecular weight of the mixture from the mole fractions
		*@param x mole fractions of species
		*@param MW the molecular weight in kg/kmol
		*/
		virtual inline double MolecularWeight_From_MoleFractions(const double* x);

		/**
		*@brief Returns the molecular weight of the mixture from the mass fractions
		*@param y mass fractions of species
		*@param MW the molecular weight in kg/kmol
		*/
		virtual inline double MolecularWeight_From_MassFractions(const double* y);

		/**
		*@brief Calculates the mass fractions and the molecular weight from the mole fractions
		*/
		virtual inline void MassFractions_From_MoleFractions(double* y, double& MW, const double* x);

		/**
		*@brief Calculates the mole fractions and the molecular weight from the mass fractions
		*/
		virtual inline void MoleFractions_From_MassFractions(double* x, double& MW, const double* y);

		/**
		*@brief Returns the molecular weight of the solid mixture from the solid mole fractions
		*@param x mole fractions of species
		*@param MW the molecular weight of solid mixture in kg/kmol
		*/
		virtual inline double SolidMolecularWeight_From_SolidMoleFractions(const double* x);

		/**
		*@brief Returns the molecular weight of the solid mixture from the solid mass fractions
		*@param y mass fractions of solid species
		*@param MW the molecular weight of solid mixture in kg/kmol
		*/
		virtual inline double SolidMolecularWeight_From_SolidMassFractions(const double* y);

		/**
		*@brief Calculates the solid mass fractions and the solid molecular weight from the solid mole fractions
		*/
		virtual inline void SolidMassFractions_From_SolidMoleFractions(double* y, double& MW, const double* x);
		
		/**
		*@brief Calculates the solid mole fractions and the molecular weight from the solid mass fractions
		*/	
		virtual inline void SolidMoleFractions_From_SolidMassFractions(double* x, double& MW, const double* y);
                                
		/**
		*@brief Calculates the mixture averaged specific heat (equivalent to CKCPBL in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param cpmix the mixture-averaged specific heat in J/kmol/K
		*/
		virtual double cpMolar_Mixture_From_MoleFractions(const double* x);

		/**
		*@brief Calculates the mixture averaged enthalpy (equivalent to CKHBML in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param hmix the mixture-averaged enthalpy in J/kmol
		*/
		virtual double hMolar_Mixture_From_MoleFractions(const double* x);

		/**
		*@brief Calculates the mixture averaged entropy (equivalent to CKSBML in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param smix the mixture-averaged entropy in J/kmol/K
		*/
		virtual double sMolar_Mixture_From_MoleFractions(const double* x);

		/**
		*@brief Calculates the mixture averaged internal energy (equivalent to CKUBML in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param umix mixture-averaged internal energy in J/kmol
		*/
		virtual double uMolar_Mixture_From_MoleFractions(const double* x);

		/**
		*@brief Calculates the mixture averaged Gibb's free energy (equivalent to CKGBML in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param gmix the mixture-averaged Gibb's free energy in J/kmol
		*/
		virtual double gMolar_Mixture_From_MoleFractions(const double* x);

		/**
		*@brief Calculates the mixture averaged Helmholtz free energy (equivalent to CKABML in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param amix the mixture-averaged Helmholtz free energy in J/kmol
		*/
		virtual double aMolar_Mixture_From_MoleFractions(const double* x);

		/**
		*@brief Calculates the standard specific heats of species at constant pressure (equivalent to CKCPOR in CHEMKIN software)
		*@return the specific heats at constant pressure in J/kmol/K
		*/
		virtual void cpMolar_Species(double* cp_species);

		/**
		*@brief Calculates the standard enthalpies of species (equivalent to CKHORT in CHEMKIN software)
		*@return the standard enthalpies in J/kmol
		*/
		virtual void hMolar_Species(double* h_species);

		/**
		*@brief Calculates the standard entropies of species (equivalent to CKSOR in CHEMKIN software)
		*@return the standard entropies in J/kmol/K
		*/
		virtual void sMolar_Species(double* s_species);
                
        /**
		*@brief Calculates the standard specific heats of species at constant pressure (only for the solid species)
		*@return the specific heats at constant pressure in J/kmol/K
		*/
		virtual void cpMolar_SolidSpecies(double* cp_solidspecies);

		/**
		*@brief Calculates the standard enthalpies of species (only for the solid species)
		*@return the standard enthalpies in J/kmol
		*/
		virtual void hMolar_SolidSpecies(double* h_solidspecies);

		/**
		*@brief Calculates the standard entropies of species (only for the solid species)
		*@return the standard entropies in J/kmol/K
		*/
		virtual void sMolar_SolidSpecies(double* s_solidspecies);

		/**
		*@brief Calculates the standard internal energies of species (equivalent to CKUML in CHEMKIN software)
		*@return the standard entropies in J/kmol
		*/
		virtual void uMolar_Species(double* u_species);

		/**
		*@brief Calculates the standard Gibb's free energies of species (equivalent to CKGML in CHEMKIN software)
		*@return the standard Gibb's free energies in J/kmol
		*/
		virtual void gMolar_Species(double* g_species);

		/**
		*@brief Calculates the standard Helmholtz's free energies of species (equivalent to CKAML in CHEMKIN software)
		*@return the standard Helmholtz's free energies in J/kmol
		*/
		virtual void aMolar_Species(double* a_species);
                
		/**
		*@brief Calculates the standard internal energies of species (only for solid species)
		*@return the standard entropies in J/kmol
		*/
		virtual void uMolar_SolidSpecies(double* u_species);

		/**
		*@brief Calculates the standard Gibb's free energies of species (only for solid species)
		*@return the standard Gibb's free energies in J/kmol
		*/
		virtual void gMolar_SolidSpecies(double* g_species);

		/**
		*@brief Calculates the standard Helmholtz's free energies of species (only for solid species)
		*@return the standard Helmholtz's free energies in J/kmol
		*/
		virtual void aMolar_SolidSpecies(double* a_species);

		/**
		*@brief Calculates the mixture averaged entropies of species (no equivalent in CHEMKIN software)
		*       The mixture averaged species entropies are not the same as the standard state values
		*       and must account for the appropriate pressure abd entropy-of-mixing terms
		*@return the mixture averaged entropies in J/kmol/K
		*/
		virtual void sMolar_Species_MixtureAveraged_From_MoleFractions(double* s_species, const double* x);

		/**
		*@brief Calculates the mixture averaged Gibb's free energies of species (no equivalent in CHEMKIN software)
		*       The mixture averaged Gibb's free energies are not the same as the standard state values
		*       and must account for the appropriate pressure abd entropy-of-mixing terms
		*@return the mixture averaged Gibb's free energies in J/kmol
		*/
		virtual void gMolar_Species_MixtureAveraged_From_MoleFractions(double* g_species, const double* x);

		/**
		*@brief Calculates the mixture averaged Helmholtz's free energies of species (no equivalent in CHEMKIN software)
		*       The mixture averaged Helmholtz's free energies are not the same as the standard state values
		*       and must account for the appropriate pressure abd entropy-of-mixing terms
		*@return the mixture averaged Helmholtz's free energies in J/kmol
		*/
		virtual void aMolar_Species_MixtureAveraged_From_MoleFractions(double* a_species, const double* x);

		/**
		*@brief Returns the normalized species specific heats at constant pressure
		*/
		virtual const std::vector<double>& Species_Cp_over_R() { cp_over_R(); return species_cp_over_R__; }

		/**
		*@brief Returns the normalized species enthalpies
		*/
		virtual const std::vector<double>& Species_H_over_RT() { h_over_RT(); return species_h_over_RT__; }
		
		/**
		*@brief Returns the normalized species entropies
		*/
		virtual const std::vector<double>& Species_S_over_R()  { s_over_R(); return species_s_over_R__; }

		/**
		*@brief Import the coefficients from a file in ASCII format (obsolete, TOREMOVE)
		*/
		virtual void ImportCoefficientsFromASCIIFile(std::ifstream& fInput);

		/**
		*@brief Import the coefficients from a file in XML format
		*/
		virtual void ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree);

		/**
		*@brief Import the species from a file in XML format
		*/
		virtual void ImportSpeciesFromXMLFile(boost::property_tree::ptree& ptree);

		/**
		*@brief TODO (This is just a test subroutine!)
		*/
		void Test(const int nLoops, const double& T, int* index);

		/**
		*@brief Calculates the normalized specific heats for the species
		*/
		inline void cp_over_R();

		/**
		*@brief Returns the number of materials
		*/
		inline unsigned int number_of_materials() const { return number_of_materials_; }

		/**
		*@brief Returns the total number of gas species
		*/
		inline unsigned int number_of_gas_species() const { return number_of_gas_species_; }

		/**
		*@brief Returns the total number of solid species (i.e. all the materials are accounted for)
		*/
		inline unsigned int number_of_solid_species() const { return number_of_solid_species_; }

		/**
		*@brief Returns the names of all the solid species (i.e. all the materials are accounted for)
		        The returned vector is 0-index based
		*/
		inline const std::vector<std::string>& vector_names_solid_species() const { return vector_names_solid_species_; }

		/**
		*@brief Returns the names of the solid species for each material
		        Example: matrix_names_solid_species()[k][j] returns the j-th species of the k-th material
				Be careful: both the indices runs from 0, i.e. the first material is the material number 0
 		*/
		inline const std::vector< std::vector<std::string> >&	matrix_names_solid_species()	const {return matrix_names_solid_species_;}
		
		/**
		*@brief Returns the solid indices (running from 1 to the total number of solid species) of the solid species for each material
		        Example: matrix_indices_solid_species()[k][j] returns the j-th species of the k-th material
				Be careful: both the indices runs from 0, i.e. the first material is the material number 0 and the first local species is 
				the number 0, BUT the returned index is 1-based
 		*/
		inline const std::vector< std::vector<unsigned int>  >&	matrix_indices_solid_species()		const {return matrix_indices_solid_species_;}

		/**
		*@brief Returns the atomic composition matrix
 		*/
		const Eigen::MatrixXd& atomic_composition() const { return this->atomic_composition_; }

		/**
		*@brief Returns the list of names of the elements 
 		*/
		const std::vector<std::string>& elements() const { return this->elements_; }

	protected:
		
		/**
		*@brief Allocates the memory
		*/
		void MemoryAllocation();

		/**
		*@brief Sets the thermodynamic coefficients for each species in the map
		*/
		void SetCoefficients(const unsigned int i, const double* coefficients);

		/**
		*@brief Calculates the normalized specific enthalpies for the species
		*/
		inline void h_over_RT();

		/**
		*@brief Calculates the normalized entropies for the species
		*/
		inline void s_over_R();

		/**
		*@brief Calculates the free Gibbs' energies for the species
		*/
		inline void g_over_RT();

		/**
		*@brief Returns the derivatives of the concentrations with respect to the mass fractions
		*/
		void DerivativesOfConcentrationsWithRespectToMassFractions(const double cTot, const double MW, const double* omega, Eigen::MatrixXd* dc_over_omega);


	private:

		std::vector<double> species_cp_over_R__;
		std::vector<double> species_h_over_RT__;
		std::vector<double> species_g_over_RT__;
		std::vector<double> species_s_over_R__;

		double*	Cp_LT;	/**< correlation coefficients, specific heat, low temperature region */
		double*	Cp_HT;	/**< correlation coefficients, specific heat, high temperature region */
		double*	DH_LT;	/**< correlation coefficients, enthalpy, low temperature region */
		double*	DH_HT;	/**< correlation coefficients, enthalpy, high temperature region */
		double*	DS_LT;	/**< correlation coefficients, entropy, low temperature region */
		double*	DS_HT;	/**< correlation coefficients, entropy, high temperature region */

		double* TL;			/**< low  temperature limit */
		double* TM;			/**< mean temperature limit */
		double* TH;			/**< high temperature limit */

		bool cp_must_be_recalculated_;
		bool h_must_be_recalculated_;
		bool s_must_be_recalculated_;
                
        std::vector<double> aux_vector_total_number_species_;

	private:

		unsigned int number_of_materials_; 
		unsigned int number_of_solid_species_;
		unsigned int number_of_gas_species_;

		std::vector<std::string>					vector_names_solid_species_;
		std::vector< std::vector<std::string> >		matrix_names_solid_species_;
		std::vector< std::vector<unsigned int> >	matrix_indices_solid_species_;

	};
}

#include "ThermodynamicsMap_Solid_CHEMKIN.hpp"

#endif /* OpenSMOKE_ThermodynamicsMap_Solid_CHEMKIN_H */
