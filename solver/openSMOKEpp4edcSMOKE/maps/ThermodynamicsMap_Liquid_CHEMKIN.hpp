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

#include "math/OpenSMOKEUtilities.h"

namespace OpenSMOKE
{ 
	ThermodynamicsMap_Liquid_CHEMKIN::ThermodynamicsMap_Liquid_CHEMKIN(boost::property_tree::ptree& ptree)
	{
		ImportSpeciesFromXMLFile(ptree);
		this->ImportElementsFromXMLFile(ptree);
		MemoryAllocation();
		ImportCoefficientsFromXMLFile(ptree);
	}

	void ThermodynamicsMap_Liquid_CHEMKIN::MemoryAllocation()
	{
		Cp_LT = new double[5*this->nspecies_];	
		Cp_HT = new double[5*this->nspecies_];
		DH_LT = new double[6*this->nspecies_];
		DH_HT = new double[6*this->nspecies_];
		DS_LT = new double[6*this->nspecies_];
		DS_HT = new double[6*this->nspecies_];
		TL = new double[this->nspecies_];	
		TH = new double[this->nspecies_];
		TM = new double[this->nspecies_];

		this->MW__.resize(this->nspecies_);

		species_cp_over_R__.resize(this->nspecies_);
		std::fill(species_cp_over_R__.begin(), species_cp_over_R__.end(), 0.);
		species_h_over_RT__.resize(this->nspecies_);
		std::fill(species_h_over_RT__.begin(), species_h_over_RT__.end(), 0.);
		species_g_over_RT__.resize(this->nspecies_);
		std::fill(species_g_over_RT__.begin(), species_g_over_RT__.end(), 0.);
		species_s_over_R__.resize(this->nspecies_);
		std::fill(species_s_over_R__.begin(), species_s_over_R__.end(), 0.);       
		aux_vector_total_number_species_.resize(this->nspecies_);
		std::fill(aux_vector_total_number_species_.begin(), aux_vector_total_number_species_.end(), 0.);

		cp_must_be_recalculated_ = true;
		h_must_be_recalculated_ = true;
		s_must_be_recalculated_ = true;
	}
 
	void ThermodynamicsMap_Liquid_CHEMKIN::SetTemperature(const double& T)
	{
		this->T_ = T;
		cp_must_be_recalculated_ = true;
		h_must_be_recalculated_ = true;
		s_must_be_recalculated_ = true;
	}

	void ThermodynamicsMap_Liquid_CHEMKIN::SetPressure(const double& P)
	{
		this->P_ = P;
	}
 
	void ThermodynamicsMap_Liquid_CHEMKIN::SetCoefficients(const unsigned k, const double* coefficients)
	{
		const double one_third = 1./3.;

		// Specific heat: high temperature
		{
			unsigned int i = k*5;
			Cp_HT[i++] = coefficients[0];
			Cp_HT[i++] = coefficients[1];
			Cp_HT[i++] = coefficients[2];
			Cp_HT[i++] = coefficients[3];
			Cp_HT[i++] = coefficients[4];
		}

		// Specific heat: low temperature
		{
			unsigned int i = k*5;
			Cp_LT[i++] = coefficients[7];
			Cp_LT[i++] = coefficients[8];
			Cp_LT[i++] = coefficients[9];
			Cp_LT[i++] = coefficients[10];
			Cp_LT[i++] = coefficients[11];
		}

		// Enthalpy: high temperature
		{
			unsigned int j = k*5;
			unsigned int i = k*6;
			DH_HT[i++] = Cp_HT[j++];
			DH_HT[i++] = 0.50 *Cp_HT[j++];
			DH_HT[i++] = one_third*Cp_HT[j++];
			DH_HT[i++] = 0.25 *Cp_HT[j++];
			DH_HT[i++] = 0.20 *Cp_HT[j++];
			DH_HT[i++] = coefficients[5];
		}
	
		// Enthalpy: low temperature
		{
			unsigned int j = k*5;
			unsigned int i = k*6;
			DH_LT[i++] = Cp_LT[j++];
			DH_LT[i++] = 0.50 *Cp_LT[j++];
			DH_LT[i++] = one_third*Cp_LT[j++];
			DH_LT[i++] = 0.25 *Cp_LT[j++];
			DH_LT[i++] = 0.20 *Cp_LT[j++];
			DH_LT[i++] = coefficients[12];
		}

		// Entropy: high temperature
		{
			unsigned int j = k*5;
			unsigned int i = k*6;
			DS_HT[i++] = Cp_HT[j++];
			DS_HT[i++] = Cp_HT[j++];
			DS_HT[i++] = 0.50 *Cp_HT[j++];
			DS_HT[i++] = one_third*Cp_HT[j++];
			DS_HT[i++] = 0.25 *Cp_HT[j++];
			DS_HT[i++] = coefficients[6];
		}

		// Entropy: low temperature
		{
			unsigned int j = k*5;
			unsigned int i = k*6;
			DS_LT[i++] = Cp_LT[j++];
			DS_LT[i++] = Cp_LT[j++];
			DS_LT[i++] = 0.50 *Cp_LT[j++];
			DS_LT[i++] = one_third*Cp_LT[j++];
			DS_LT[i++] = 0.25 *Cp_LT[j++];
			DS_LT[i++] = coefficients[13];
		}

		// Temperature limits
		{
			TL[k] = coefficients[14];
			TH[k] = coefficients[15];
			TM[k] = coefficients[16];
		}

		this->MW__[k] = coefficients[17];
	}

	void ThermodynamicsMap_Liquid_CHEMKIN::ImportCoefficientsFromASCIIFile(std::ifstream& fInput)
	{
		std::cout << " * Reading thermodynamic coefficients of species..." << std::endl;
		double coefficients[18];
		for(unsigned int i=0;i<this->nspecies_;i++)
		{
			for(unsigned int j=0;j<18;j++)
				fInput >> coefficients[j];
			SetCoefficients(i, coefficients);
		}
	}

	void ThermodynamicsMap_Liquid_CHEMKIN::ImportSpeciesFromXMLFile(boost::property_tree::ptree& ptree)
	{
		// Names of species
		{
			try
			{
				this->nspecies_ = ptree.get<unsigned int>("opensmoke.NumberOfSpecies");  
				this->names_.resize(this->nspecies_);

				std::stringstream stream;
				stream.str( ptree.get< std::string >("opensmoke.NamesOfSpecies") );  
				for(unsigned int i=0;i<this->nspecies_;i++)
					stream >> this->names_[i];

				number_of_materials_    = ptree.get<unsigned int>("opensmoke.NumberOfMaterials");
				number_of_liquid_species_ = ptree.get<unsigned int>("opensmoke.NumberOfLiquidSpecies");
				number_of_gas_species_  = this->nspecies_ - number_of_liquid_species_;

				// Liquid species
				{
					std::stringstream stream;
					stream.str( ptree.get< std::string >("opensmoke.LiquidSpecies") ); 
					vector_names_liquid_species_.resize(number_of_liquid_species_);
					for(unsigned int i=0;i<number_of_liquid_species_;i++)
						stream >> vector_names_liquid_species_[i]; 	
				}
			}
			catch(...)
			{
				ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::ImportSpeciesFromXMLFile", "Error in reading the list of species.");
			}
		}

		// Reading all the available materials
		{
			matrix_names_liquid_species_.resize(number_of_materials_);
			matrix_indices_liquid_species_.resize(number_of_materials_);

			unsigned int k=0;
			BOOST_FOREACH( boost::property_tree::ptree::value_type const& node, ptree.get_child( "opensmoke" ) ) 
			{
				boost::property_tree::ptree subtree = node.second;  

				if( node.first == "MaterialDescription" ) 
				{
					std::stringstream stream;
					stream.str( subtree.get< std::string >("") ); 

					unsigned int n;

					stream >> n;
					matrix_names_liquid_species_[k].resize(n);
					matrix_indices_liquid_species_[k].resize(n);

					for(unsigned int i=0;i<n;i++)
					{
						stream >> matrix_names_liquid_species_[k][i];
						stream >> matrix_indices_liquid_species_[k][i];
					}

					k++;
				}
			}
		}
	}
 
	void ThermodynamicsMap_Liquid_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)
	{
		std::cout << " * Reading thermodynamic coefficients of species from XML file..." << std::endl;

		std::stringstream stream;
		stream.str( ptree.get< std::string >("opensmoke.Thermodynamics.NASA-coefficients") );

		double coefficients[18];
		for(unsigned int i=0;i<this->nspecies_;i++)
		{
			for(unsigned int j=0;j<18;j++)
				stream >> coefficients[j];
			SetCoefficients(i, coefficients);
		}
	}
 
	inline double ThermodynamicsMap_Liquid_CHEMKIN::MolecularWeight_From_MoleFractions(const double* x)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::MolecularWeight_From_MoleFractions", "Function not available for surface thermodynamics!");
		return 0;
	}

	inline double ThermodynamicsMap_Liquid_CHEMKIN::MolecularWeight_From_MassFractions(const double* y)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::MolecularWeight_From_MassFractions", "Function not available for surface thermodynamics!");
		return 0;
	}

	inline void ThermodynamicsMap_Liquid_CHEMKIN::MassFractions_From_MoleFractions(double* y, double& MW, const double* x)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::MassFractions_From_MoleFractions", "Function not available for surface thermodynamics!");
	}

	inline void ThermodynamicsMap_Liquid_CHEMKIN::MoleFractions_From_MassFractions(double* x, double& MW, const double* y)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::MoleFractions_From_MassFractions", "Function not available for surface thermodynamics!");
	}

	inline double ThermodynamicsMap_Liquid_CHEMKIN::LiquidMolecularWeight_From_LiquidMoleFractions(const double* x)
	{	
		double MW = 0.;
		for (unsigned int j=0;j<number_of_liquid_species_;j++)
		    MW += x[j] * this->MW__[number_of_gas_species_ + j];
			return MW;
	}

	inline double ThermodynamicsMap_Liquid_CHEMKIN::LiquidMolecularWeight_From_LiquidMassFractions(const double* y)
	{
		double  MW = 0.;
		for (unsigned int j=0;j<number_of_liquid_species_;j++)
		    MW += y[j] / this->MW__[number_of_gas_species_ + j];

		MW = 1./MW;
			return MW;
	}

	inline void ThermodynamicsMap_Liquid_CHEMKIN::LiquidMassFractions_From_LiquidMoleFractions(double* y, double& MW, const double* x)
	{
		MW = 0.;
		for (unsigned int j=0;j<number_of_liquid_species_;j++)
		    MW += x[j] * this->MW__[number_of_gas_species_ + j];
		for (unsigned int j=0;j<number_of_liquid_species_;j++)
		    y[j] = x[j] * this->MW__[number_of_gas_species_ + j] / MW;
	}
	
	inline void ThermodynamicsMap_Liquid_CHEMKIN::LiquidMoleFractions_From_LiquidMassFractions(double* x, double& MW, const double* y)
	{
		MW = 0.;
		for (unsigned int j=0;j<number_of_liquid_species_;j++)
		    MW += y[j] / this->MW__[number_of_gas_species_ + j];
		MW = 1./MW;
		for (unsigned int j=0;j<number_of_liquid_species_;j++)
		    x[j] = y[j] / this->MW__[number_of_gas_species_ + j] * MW;
	}
 
	inline void ThermodynamicsMap_Liquid_CHEMKIN::cp_over_R()
	{
		if (cp_must_be_recalculated_ == true)
		{
			const double T2 = this->T_*this->T_;
			const double T3 = T2*this->T_;
			const double T4 = T3*this->T_;

			unsigned int j = 0;
			for (unsigned int k=0;k<this->nspecies_;k++)
			{
				species_cp_over_R__[k] = (this->T_>TM[k]) ?		Cp_HT[j] + this->T_*Cp_HT[j+1] + T2*Cp_HT[j+2] + T3*Cp_HT[j+3] + T4*Cp_HT[j+4] :
																Cp_LT[j] + this->T_*Cp_LT[j+1] + T2*Cp_LT[j+2] + T3*Cp_LT[j+3] + T4*Cp_LT[j+4] ;
				j += 5;
			}

			cp_must_be_recalculated_ = false;
		}
	}

	inline void ThermodynamicsMap_Liquid_CHEMKIN::h_over_RT()
	{
		if (h_must_be_recalculated_ == true)
		{
			const double T2 = this->T_*this->T_;
			const double T3 = T2*this->T_;
			const double T4 = T3*this->T_;
			const double uT = 1./this->T_;

			int j = 0;
			for (unsigned int k=0;k<this->nspecies_;k++)
			{
				species_h_over_RT__[k] = (this->T_>TM[k]) ?	DH_HT[j] + this->T_*DH_HT[j+1] + T2*DH_HT[j+2] + T3*DH_HT[j+3] + T4*DH_HT[j+4] + uT*DH_HT[j+5]:
															DH_LT[j] + this->T_*DH_LT[j+1] + T2*DH_LT[j+2] + T3*DH_LT[j+3] + T4*DH_LT[j+4] + uT*DH_LT[j+5];
				j += 6;
			}

			h_must_be_recalculated_ = false;
		}
	}
	
	inline void ThermodynamicsMap_Liquid_CHEMKIN::s_over_R()
	{
		if (s_must_be_recalculated_ == true)
		{
			const double logT = log(this->T_);
			const double T2 = this->T_*this->T_;
			const double T3 = T2*this->T_;
			const double T4 = T3*this->T_;

			unsigned int j = 0;
			for (unsigned int k=0;k<this->nspecies_;k++)
			{
				species_s_over_R__[k] = (this->T_>TM[k]) ?	DS_HT[j]*logT + this->T_*DS_HT[j+1] + T2*DS_HT[j+2] + T3*DS_HT[j+3] + T4*DS_HT[j+4] + DS_HT[j+5] :
															DS_LT[j]*logT + this->T_*DS_LT[j+1] + T2*DS_LT[j+2] + T3*DS_LT[j+3] + T4*DS_LT[j+4] + DS_LT[j+5] ;
				j += 6;
			}

			s_must_be_recalculated_ = false;
		}
	}

	inline void ThermodynamicsMap_Liquid_CHEMKIN::g_over_RT()
	{
		h_over_RT();
		s_over_R();

		Difference(this->nspecies_, species_h_over_RT__.data(), species_s_over_R__.data(), species_g_over_RT__.data());
	}

	double ThermodynamicsMap_Liquid_CHEMKIN::cpMolar_Mixture_From_MoleFractions(const double* x)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::cpMolar_Mixture_From_MoleFractions", "Function not available for Liquid thermodynamics!");
		return 0.;
	}
	 
	double ThermodynamicsMap_Liquid_CHEMKIN::hMolar_Mixture_From_MoleFractions(const double* x)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::hMolar_Mixture_From_MoleFractions", "Function not available for Liquid thermodynamics!");
		return 0.;
	}

	double ThermodynamicsMap_Liquid_CHEMKIN::sMolar_Mixture_From_MoleFractions(const double* x)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::sMolar_Mixture_From_MoleFractions", "Function not available for Liquid thermodynamics!");
		return 0.;
	}

	double ThermodynamicsMap_Liquid_CHEMKIN::uMolar_Mixture_From_MoleFractions(const double* x)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::uMolar_Mixture_From_MoleFractions", "Function not available for Liquid thermodynamics!");
		return 0.;
	}
 
	double ThermodynamicsMap_Liquid_CHEMKIN::gMolar_Mixture_From_MoleFractions(const double* x)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::gMolar_Mixture_From_MoleFractions", "Function not available for Liquid thermodynamics!");
		return 0.;
	}
 
	double ThermodynamicsMap_Liquid_CHEMKIN::aMolar_Mixture_From_MoleFractions(const double* x)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::aMolar_Mixture_From_MoleFractions", "Function not available for Liquid thermodynamics!");
		return 0.;
	}

	void ThermodynamicsMap_Liquid_CHEMKIN::cpMolar_Species(double* cp_species)
	{
		cp_over_R();
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol, species_cp_over_R__.data(), cp_species);
	}

	void ThermodynamicsMap_Liquid_CHEMKIN::cpMolar_LiquidSpecies(double* cp_liquidspecies)
	{
		cpMolar_Species(aux_vector_total_number_species_.data());
		for (unsigned int j = 0; j < number_of_liquid_species_; j++)
		    cp_liquidspecies[j] = aux_vector_total_number_species_[number_of_gas_species_+j];
	}
        
	void ThermodynamicsMap_Liquid_CHEMKIN::hMolar_Species(double* h_species)
	{
		h_over_RT();
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol*this->T_, species_h_over_RT__.data(), h_species);
	}
        
	void ThermodynamicsMap_Liquid_CHEMKIN::hMolar_LiquidSpecies(double* h_liquidspecies)
	{
		hMolar_Species(aux_vector_total_number_species_.data());
		for (unsigned int j = 0; j < number_of_liquid_species_; j++)
			h_liquidspecies[j] = aux_vector_total_number_species_[number_of_gas_species_+j];
	}        

	void ThermodynamicsMap_Liquid_CHEMKIN::sMolar_Species(double* s_species)
	{
		s_over_R();
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol, species_s_over_R__.data(), s_species);
	}
        
	void ThermodynamicsMap_Liquid_CHEMKIN::sMolar_LiquidSpecies(double* s_liquidspecies)
	{
		sMolar_Species(aux_vector_total_number_species_.data());
		for (unsigned int j = 0; j < number_of_liquid_species_; j++)
		    s_liquidspecies[j] = aux_vector_total_number_species_[number_of_gas_species_+j];
	}
         
	void ThermodynamicsMap_Liquid_CHEMKIN::sMolar_Species_MixtureAveraged_From_MoleFractions(double* s_species, const double* x)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::sMolar_Species_MixtureAveraged_From_MoleFractions", "Function not available for Liquid thermodynamics!");
	}
 
	void ThermodynamicsMap_Liquid_CHEMKIN::uMolar_Species(double* u_species)
	{
		h_over_RT();
		Sum(this->nspecies_, species_h_over_RT__.data(), -1., u_species);
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol*this->T_, u_species);
	}
         
	void ThermodynamicsMap_Liquid_CHEMKIN::uMolar_LiquidSpecies(double* u_liquidspecies)
	{
		uMolar_Species(aux_vector_total_number_species_.data());
		for (unsigned int j = 0; j < number_of_liquid_species_; j++)
			u_liquidspecies[j] = aux_vector_total_number_species_[number_of_gas_species_+j];
	}        

	void ThermodynamicsMap_Liquid_CHEMKIN::gMolar_Species(double* g_species)
	{
		h_over_RT();
		s_over_R();
		Difference(this->nspecies_, species_h_over_RT__.data(), species_s_over_R__.data(), g_species);
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol*this->T_, g_species);
	}
        
	void ThermodynamicsMap_Liquid_CHEMKIN::gMolar_LiquidSpecies(double* g_liquidspecies)
	{
		gMolar_Species(aux_vector_total_number_species_.data());
		for (unsigned int j = 0; j < number_of_liquid_species_; j++)
		    g_liquidspecies[j] = aux_vector_total_number_species_[number_of_gas_species_+j];
	}

	void ThermodynamicsMap_Liquid_CHEMKIN::gMolar_Species_MixtureAveraged_From_MoleFractions(double* g_species, const double* x)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::gMolar_Species_MixtureAveraged_From_MoleFractions", "Function not available for Liquid thermodynamics!");
	}

	void ThermodynamicsMap_Liquid_CHEMKIN::aMolar_Species(double* a_species)
	{
		h_over_RT();
		s_over_R();
		Difference(this->nspecies_, species_h_over_RT__.data(), species_s_over_R__.data(), a_species);
		Sum(this->nspecies_, -1., a_species);
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol*this->T_, a_species);
	}
        
	void ThermodynamicsMap_Liquid_CHEMKIN::aMolar_LiquidSpecies(double* a_liquidspecies)
	{
		aMolar_Species(aux_vector_total_number_species_.data());
		for (unsigned int j = 0; j < number_of_liquid_species_; j++)
			a_liquidspecies[j] = aux_vector_total_number_species_[number_of_gas_species_+j];
	}        

	void ThermodynamicsMap_Liquid_CHEMKIN::aMolar_Species_MixtureAveraged_From_MoleFractions(double* a_species, const double* x)
	{
		gMolar_Species_MixtureAveraged_From_MoleFractions(a_species, x);
		Sum(this->nspecies_, -PhysicalConstants::R_J_kmol*this->T_, a_species);
	}

	void ThermodynamicsMap_Liquid_CHEMKIN::DerivativesOfConcentrationsWithRespectToMassFractions(const double cTot, const double MW, const double* omega, Eigen::MatrixXd* dc_over_omega)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::DerivativesOfConcentrationsWithRespectToMassFractions", "Function not available for Liquid thermodynamics!");
	}

	void ThermodynamicsMap_Liquid_CHEMKIN::Test(const int nLoops, const double& T, int* index)
	{
		ErrorMessage("ThermodynamicsMap_Liquid_CHEMKIN::Test", "Function not available for Liquid thermodynamics!");
	}

}
