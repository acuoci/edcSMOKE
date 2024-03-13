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

#include "math/OpenSMOKEUtilities.h"
#include "ThermodynamicsMap.h"

namespace OpenSMOKE
{ 
	KineticsMap_Liquid_CHEMKIN::KineticsMap_Liquid_CHEMKIN(ThermodynamicsMap_Liquid_CHEMKIN& thermo, boost::property_tree::ptree& ptree, const unsigned int target) :
	thermodynamics_(thermo)
	{
		ImportSpeciesFromXMLFile(ptree);
		ImportCoefficientsFromXMLFile(ptree, target);
		this->T_ = this->P_ = 0.;
	}

	KineticsMap_Liquid_CHEMKIN::KineticsMap_Liquid_CHEMKIN(ThermodynamicsMap_Liquid_CHEMKIN& thermo, boost::property_tree::ptree& ptree, const std::string target) :
	thermodynamics_(thermo)
	{
		ImportSpeciesFromXMLFile(ptree);
		ImportCoefficientsFromXMLFile(ptree, target);
		this->T_ = this->P_ = 0.;
	}

	void KineticsMap_Liquid_CHEMKIN::SetTemperature(const double& T)
	{
		this->T_old_ = this->T_;
		this->T_ = T;
		this->uT_ = 1./this->T_;
		this->logT_ = log(this->T_);
		Patm_over_RT_ = 101325./PhysicalConstants::R_J_kmol/this->T_;
		log_Patm_over_RT_ = log(Patm_over_RT_);

		if (std::fabs(this->T_-this->T_old_)/this->T_>1.e-14)
		{
			arrhenius_kinetic_constants_must_be_recalculated_ = true;
			nonconventional_kinetic_constants_must_be_recalculated_ = true;
			reaction_h_and_s_must_be_recalculated_ = true;
		}
	}

	void KineticsMap_Liquid_CHEMKIN::SetPressure(const double& P)
	{
		this->P_old_ = this->P_;
		this->P_ = P;
	//	if (std::fabs(this->P_-this->P_old_)/this->P_>1.e-14)
		{
			nonconventional_kinetic_constants_must_be_recalculated_ = true;
		}
	}

	void KineticsMap_Liquid_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)
	{
		ErrorMessage("void KineticsMap_Liquid_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)", "The liquid kinetic map require the user specifies tha material name.");
	}

	void KineticsMap_Liquid_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree, const std::string target_material_name)
	{
		std::string kinetics_type =  ptree.get<std::string>("opensmoke.Kinetics.<xmlattr>.type");
		std::string kinetics_version = ptree.get<std::string>("opensmoke.Kinetics.<xmlattr>.version");
		if (kinetics_type != "OpenSMOKE" || kinetics_version != "01-02-2014")
			ErrorMessage("void KineticsMap_Liquid_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)", "The current kinetic scheme is not supported.");

		unsigned int target_material_index = 0;
		BOOST_FOREACH( boost::property_tree::ptree::value_type const& node, ptree.get_child( "opensmoke.Kinetics" ) )
		{
			boost::property_tree::ptree subtree = node.second;  

			if( node.first == "MaterialKinetics" ) 
			{
				std::string material_name =  subtree.get<std::string>("<xmlattr>.name");
				if ( target_material_name == material_name )
				{
					target_material_index = subtree.get<unsigned int>("<xmlattr>.index");
					break;
				}
			}
		}

		if (target_material_index == 0)
			ErrorMessage("void KineticsMap_Liquid_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)", "The requested liquid material is not available. Please check the name.");
		else
			ImportCoefficientsFromXMLFile(ptree, target_material_index);
	}
 
	void KineticsMap_Liquid_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree, const unsigned int target_material_index)
	{
		if (target_material_index <=0 || target_material_index > thermodynamics_.number_of_materials())
			ErrorMessage("void KineticsMap_Liquid_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)", "The requested liquid material is not available. Please check the name.");

		std::string kinetics_type =  ptree.get<std::string>("opensmoke.Kinetics.<xmlattr>.type");
		std::string kinetics_version = ptree.get<std::string>("opensmoke.Kinetics.<xmlattr>.version");
		if (kinetics_type != "OpenSMOKE" || kinetics_version != "01-02-2014")
			ErrorMessage("void KineticsMap_Liquid_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)", "The current kinetic scheme is not supported.");


		BOOST_FOREACH( boost::property_tree::ptree::value_type const& node, ptree.get_child( "opensmoke.Kinetics" ) )
		{
			boost::property_tree::ptree subtree = node.second;  

			if( node.first == "MaterialKinetics" ) 
			{
				if (target_material_index == subtree.get<unsigned int>("<xmlattr>.index") )
				{
					this->number_of_species_ = subtree.get<unsigned int>("NumberOfSpecies"); 
					this->number_of_reactions_ = subtree.get<unsigned int>("NumberOfReactions");

					// Irreversible reactions
					{
						std::cout << "Reading irreversible..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("Irreversible") );  
						Load(indices_of_irreversible_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_irreversible_reactions_ = static_cast<unsigned int>(indices_of_irreversible_reactions__.size());
					}
		
					// Reversible reactions
					{
						std::cout << "Reading reversible..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("Reversible") );  
						Load(indices_of_reversible_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_reversible_reactions_ = static_cast<unsigned int>(indices_of_reversible_reactions__.size());
					}

					// Thermodynamic Reversible reactions
					{
						std::cout << "Reading reversible thermodynamics..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("Reversible-Thermodynamics") );  
						Load(indices_of_thermodynamic_reversible_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_thermodynamic_reversible_reactions_ = static_cast<unsigned int>(indices_of_thermodynamic_reversible_reactions__.size());
					}

					// Explicit Reversible reactions
					{
						std::cout << "Reading reversible explicit..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("Reversible-Explicit") );  
						Load(indices_of_explicitly_reversible_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_explicitly_reversible_reactions_ = static_cast<unsigned int>(indices_of_explicitly_reversible_reactions__.size());
					}

					// USRPROG reactions
					{
						number_of_usrprog_reactions_ = 0;

						std::cout << "Reading USRPROG reaction..." << std::endl;

						std::stringstream stream;
						stream.str( subtree.get< std::string >("UsrProg") );

						//if (current_node != 0)
						{
							Load(indices_of_usrprog_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
							number_of_usrprog_reactions_ = static_cast<unsigned int>(indices_of_usrprog_reactions__.size());

							if (number_of_usrprog_reactions_ != 0)
							{
								std::stringstream substream;
								substream.str( subtree.get< std::string >("UsrProg-Labels") );
								Load(labels_of_usrprog_reactions__, substream, OPENSMOKE_FORMATTED_FILE);
							}
						}
					}

					// Reading if the kinetic scheme is conventional
					{
						std::cout << "Reading type of kinetic scheme..." << std::endl;
						std::string dummy = subtree.get<std::string>("TypeOfKinetics"); 
						dummy.erase(std::remove(dummy.begin(), dummy.end(), '\n'), dummy.end());

						if (dummy == "chemkin_conventional")	type_of_liquid_kinetics_ = TYPE_OF_LIQUID_KINETICS_CHEMKIN_CONVENTIONAL;
						else
						{
							std::cout << "This type of kinetic mechanism is not available: " << dummy << std::endl;
							std::cout << "Available types: chemkin_conventional" << std::endl;
							std::cout << "Press enter to exit..." << std::endl;						
							getchar();
							exit(OPENSMOKE_FATAL_ERROR_EXIT);
						}
					}

					std::cout << " * Reading kinetic parameters of reactions..." << std::endl;	
					{
						// Direct side
						{
							// lnA
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Direct.lnA") );  
			 					Load(lnA__, stream, OPENSMOKE_FORMATTED_FILE);
							}

							// Beta
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Direct.Beta") );  
			 					Load(Beta__, stream, OPENSMOKE_FORMATTED_FILE);

							}

							// E_over_R
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Direct.E_over_R") );  
			 					Load(E_over_R__, stream, OPENSMOKE_FORMATTED_FILE);
							}

							// Global kinetic order of forward reaction
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Direct.ForwardKineticOrder") );  
			 					Load(forward_kinetic_order__, stream, OPENSMOKE_FORMATTED_FILE);
							}
						}

						// Reverse side
						if (number_of_explicitly_reversible_reactions_ != 0)
						{
							// lnA
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Reverse.lnA") );  
			 					Load(lnA_reversible__, stream, OPENSMOKE_FORMATTED_FILE);
							}

							// Beta
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Reverse.Beta") );  
			 					Load(Beta_reversible__, stream, OPENSMOKE_FORMATTED_FILE);

							}

							// E_over_R
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Reverse.E_over_R") );  
			 					Load(E_over_R_reversible__, stream, OPENSMOKE_FORMATTED_FILE);
							}
						}
					}

					// Reactions needing conversion
					{			
						std::stringstream stream;
						stream.str( subtree.get< std::string >("ReactionsNeedingConversion") ); 

						std::string dummy;

						stream >> dummy;
						unsigned int number_of_reactions_needing_conversion_ = boost::lexical_cast<unsigned int>(dummy);
						indices_of_reactions_needing_conversion_.resize(number_of_reactions_needing_conversion_);

						for (unsigned int j=0;j<number_of_reactions_needing_conversion_;j++)
						{		
							stream >> dummy;
							indices_of_reactions_needing_conversion_[j] = boost::lexical_cast<unsigned int>(dummy);
						}
					}

					// Thermodynamic kinetic constants
					{
						std::stringstream stream;
						stream.str( subtree.get< std::string >("ThermodynamicReversibleReactions") ); 
						
						std::string dummy;

						stream >> dummy;
						unsigned int n = boost::lexical_cast<unsigned int>(dummy);

						delta_nu_gas_.resize(n);

						for (unsigned int j=0;j<n;j++)
						{
							stream >> dummy;
							const unsigned int index_phase = boost::lexical_cast<unsigned int>(dummy);
							stream >> dummy;
							const double delta_sigma = boost::lexical_cast<double>(dummy);
							stream >> dummy;
							delta_nu_gas_[j] = boost::lexical_cast<double>(dummy);
						}				
					}

					// Stoichiometry
					{
						std::string stoichiometry_type =  subtree.get<std::string>("Stoichiometry.<xmlattr>.type");
						std::string stoichiometry_version = subtree.get<std::string>("Stoichiometry.<xmlattr>.version");
						if (stoichiometry_type != "OpenSMOKE" || stoichiometry_version != "01-02-2014")
							ErrorMessage("void KineticsMap_Liquid_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)", "The current stoichiometric data are not supported.");

						std::stringstream stream;
						stream.str( subtree.get< std::string >("Stoichiometry") ); 

						stoichiometry_ = new StoichiometricMap(this->number_of_species_, this->number_of_reactions_);
						stoichiometry_->ReadFromASCIIFile(stream);
						changeOfMoles__ = stoichiometry_->ChangeOfMoles();

						{
							OpenSMOKE::OpenSMOKEVectorBool tmp(this->number_of_reactions_);
							tmp = false;
							for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
								tmp[indices_of_thermodynamic_reversible_reactions__[k-1]] = true;
							stoichiometry_->CompleteChangeOfMoles(tmp.GetHandle());
	 					}
			

						// Memory allocation			
						aux_vector__.resize(this->number_of_species_);
						std::fill(aux_vector__.begin(), aux_vector__.end(), 0.);

						c__.resize(this->number_of_species_);
						std::fill(c__.begin(), c__.end(), 0.);

						reaction_s_over_R__.resize(this->number_of_reactions_);
						std::fill(reaction_s_over_R__.begin(), reaction_s_over_R__.end(), 0.);

						reaction_h_over_RT__.resize(this->number_of_reactions_);
						std::fill(reaction_h_over_RT__.begin(), reaction_h_over_RT__.end(), 0.);

						kArrhenius__.resize(this->number_of_reactions_);
						std::fill(kArrhenius__.begin(), kArrhenius__.end(), 0.);

						kArrheniusModified__.resize(this->number_of_reactions_);
						std::fill(kArrheniusModified__.begin(), kArrheniusModified__.end(), 0.);

						uKeq__.resize(number_of_thermodynamic_reversible_reactions_);
						std::fill(uKeq__.begin(), uKeq__.end(), 0.);

						kArrhenius_reversible__.resize(number_of_explicitly_reversible_reactions_);
						std::fill(kArrhenius_reversible__.begin(), kArrhenius_reversible__.end(), 0.);

						forwardReactionRates__.resize(this->number_of_reactions_);
						std::fill(forwardReactionRates__.begin(), forwardReactionRates__.end(), 0.);

						reverseReactionRates__.resize(this->number_of_reactions_);
						std::fill(reverseReactionRates__.begin(), reverseReactionRates__.end(), 0.);

						netReactionRates__.resize(this->number_of_reactions_);
						std::fill(netReactionRates__.begin(), netReactionRates__.end(), 0.);

						isThermodynamicallyReversible__.resize(this->number_of_reactions_);
						std::fill(isThermodynamicallyReversible__.begin(), isThermodynamicallyReversible__.end(), 0);

						isExplicitlyReversible__.resize(this->number_of_reactions_);
						std::fill(isExplicitlyReversible__.begin(), isExplicitlyReversible__.end(), 0);

						for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
							isThermodynamicallyReversible__[indices_of_thermodynamic_reversible_reactions__[k-1]-1] = k;
						for(unsigned int k=1;k<=number_of_explicitly_reversible_reactions_;k++)
							isExplicitlyReversible__[indices_of_explicitly_reversible_reactions__[k-1]-1] = k;

			
						// Additional indices for sensitivity analysis
						{
							type_of_reaction__.resize(this->number_of_reactions_);
							std::fill(type_of_reaction__.begin(), type_of_reaction__.end(), PhysicalConstants::REACTION_SIMPLE);
						
							local_family_index__.resize(this->number_of_reactions_);
							std::fill(local_family_index__.begin(), local_family_index__.end(), 0);
						}

						std::cout << std::endl;
						std::cout << "----------------------------------------------------------------------------" << std::endl;
						std::cout << " Kinetic Mechanism Summary"<< std::endl;
						std::cout << "----------------------------------------------------------------------------" << std::endl;
						std::cout << " Total number of species:          " << this->number_of_species_ << std::endl;
						std::cout << " Total number of reactions:        " << this->number_of_reactions_ << std::endl;
						std::cout << "   Reversible reactions:           " << number_of_reversible_reactions_ << " (" << number_of_reversible_reactions_/std::max(1.,double(this->number_of_reactions_))*100. << "%)" << std::endl;
						std::cout << "    * by thermodynamics:           " << number_of_thermodynamic_reversible_reactions_ << " (" << number_of_thermodynamic_reversible_reactions_/std::max(1.,double(number_of_reversible_reactions_))*100. << "%)" << std::endl;
						std::cout << "    * by Arrhenius' law:           " << number_of_explicitly_reversible_reactions_ << " (" << number_of_explicitly_reversible_reactions_/std::max(1.,double(number_of_reversible_reactions_))*100. << "%)" << std::endl;
						std::cout << " USRPRG reactions:                 " << number_of_usrprog_reactions_ << std::endl;
						std::cout << std::endl;

						stoichiometry_->Summary(std::cout);
					}

					break;
				}
			}
		}
	}

	void KineticsMap_Liquid_CHEMKIN::ImportSpeciesFromXMLFile(boost::property_tree::ptree& ptree)
	{
		try
		{
			this->number_of_species_ = ptree.get<unsigned int>("opensmoke.NumberOfSpecies"); 				
		}
		catch(...)
		{
			ErrorMessage("KineticsMap_Surface_CHEMKIN::ImportSpeciesFromXMLFile", "Error in reading the number of species.");
		}
	}

	void KineticsMap_Liquid_CHEMKIN::ReactionEnthalpiesAndEntropies()
	{
		if (reaction_h_and_s_must_be_recalculated_ == true)
		{
			stoichiometry_->ReactionEnthalpyAndEntropy(	reaction_h_over_RT__, reaction_s_over_R__,  thermodynamics_.Species_H_over_RT(), thermodynamics_.Species_S_over_R() );

			reaction_h_and_s_must_be_recalculated_ = false;
		}
	}

	void KineticsMap_Liquid_CHEMKIN::KineticConstants()
	{
                ReactionEnthalpiesAndEntropies();
                        
		if (arrhenius_kinetic_constants_must_be_recalculated_ == true)
		{
			// Forward kinetic constants (Arrhenius' Law)
			{
				double *pt_lnA = lnA__.data();
				double *pt_Beta = Beta__.data();
				double *pt_E_over_R = E_over_R__.data();
				double *pt_kArrheniusT = kArrhenius__.data();
			
				for(unsigned int j=0;j<this->number_of_reactions_;j++)
					*pt_kArrheniusT++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;

				Exp(kArrhenius__, &kArrhenius__);
			}

			// Equilibrium constants (inverse value)
			{
				for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
				{
					unsigned int j = indices_of_thermodynamic_reversible_reactions__[k-1];
					uKeq__[k-1] = -reaction_s_over_R__[j-1] + reaction_h_over_RT__[j-1] - log_Patm_over_RT_ * delta_nu_gas_[k-1];
				}
				Exp(uKeq__, &uKeq__);
			}

			// Explicit reverse Arrhenius constants
			if (number_of_explicitly_reversible_reactions_ != 0)
			{
				double *pt_lnA = lnA_reversible__.data();
				double *pt_Beta = Beta_reversible__.data();
				double *pt_E_over_R = E_over_R_reversible__.data();
				double *pt_kArrhenius = kArrhenius_reversible__.data();

				for(unsigned int k=0;k<number_of_explicitly_reversible_reactions_;k++)
					*pt_kArrhenius++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;
			
				Exp(kArrhenius_reversible__, &kArrhenius_reversible__);
			}

			arrhenius_kinetic_constants_must_be_recalculated_ = false;
		}

		//if (nonconventional_kinetic_constants_must_be_recalculated_ == true)
		{
			nonconventional_kinetic_constants_must_be_recalculated_ = false;
		}

		// Conversions
		if (indices_of_reactions_needing_conversion_.size() > 0)
		{
			const double R_times_T = PhysicalConstants::R_J_kmol*this->T_;
			for(unsigned int k=0;k<indices_of_reactions_needing_conversion_.size();k++)
			{
				unsigned int j = indices_of_reactions_needing_conversion_[k];
				kArrhenius__[j-1] *= std::pow(R_times_T, forward_kinetic_order__[j-1]);
			}
		}

		kArrheniusModified__ = kArrhenius__;
	}

	void KineticsMap_Liquid_CHEMKIN::ReactionRates(const double* cGas, const double* cLiquid)
	{
		double cTot = 0.;
		for (unsigned int j = 0; j<thermodynamics_.number_of_gas_species(); j++)
			cTot += cGas[j];
		
		unsigned int count = 0;
		for(unsigned int j=0;j<thermodynamics_.number_of_gas_species();j++)
			c__[count++] = cGas[j];
		
		for(unsigned int j=0;j<thermodynamics_.number_of_liquid_species();j++)
			c__[count++] = cLiquid[j];
		
		if (type_of_liquid_kinetics_ == TYPE_OF_LIQUID_KINETICS_CHEMKIN_CONVENTIONAL)
		{
			// 1. Kinetic constants
			KineticConstants();
		}

		// Calculates the product of conenctrations (for forward and reverse reactions)
		// Be careful: the reverseReactionRates_ vector is defined for all the reactions
		// in the kinetic scheme, not only for the reversible reactions. After calling the
		// function reported below the value of reverseReactionRates_ vector for non reversible
		// reactions is put equal to 1.
		stoichiometry_->ProductOfConcentrations(forwardReactionRates__, reverseReactionRates__, c__.data());

		// Corrects the product of concentrations for reverse reaction by the 
		// thermodynamic equilibrium constant
		for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_thermodynamic_reversible_reactions__[k-1];
			reverseReactionRates__[j-1] *= uKeq__[k-1];
		}

		// Corrects the product of concentrations for reverse reaction by the 
		// explicit Arrhenius kinetic parameters (if provided)
		for(unsigned int k=1;k<=number_of_explicitly_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_explicitly_reversible_reactions__[k-1];
//			reverseReactionRates_[j] *= kArrhenius_reversible_[k-1]/kArrheniusModified__[j-1];
			reverseReactionRates__[j-1] *= kArrhenius_reversible__[k-1]/kArrhenius__[j-1];
		}

		// Calculates the net reaction rate
		// Be careful: the netReactionRates_ vector must be multiplied by the effective 
		// forward kinetic constant, to obtain the real reaction rates in [kmol/m3/s]
		netReactionRates__ = forwardReactionRates__;
		for(unsigned int k=1;k<=number_of_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_reversible_reactions__[k-1]-1;
			netReactionRates__[j] -= reverseReactionRates__[j];
		}

		// Multiplies the net reaction rate by the effective kinetic constant (accounting for 
		// third-body effects, fall-off, etc.). At the end of this function the netReactionRates_
		// vector contains the net reaction rates of all the reactions in [kmol/m3/s]
		ElementByElementProduct(static_cast<int>(netReactionRates__.size()), netReactionRates__.data(), kArrheniusModified__.data(), netReactionRates__.data());
	}
 
	void KineticsMap_Liquid_CHEMKIN::FormationRates(double* Rgas, double* Rliquid)
	{
		std::vector<double> R(thermodynamics_.NumberOfSpecies());
		stoichiometry_->FormationRatesFromReactionRates(R.data(), netReactionRates__.data());
		
		unsigned int count = 0;
		
		for(unsigned int j=0;j<thermodynamics_.number_of_gas_species();j++)
			Rgas[j] = R[count++];
		
		for(unsigned int j=0;j<thermodynamics_.number_of_liquid_species();j++)
			Rliquid[j] = R[count++];
	}
 
	double KineticsMap_Liquid_CHEMKIN::HeatRelease(const double* Rgas, const double* RLiquid)
	{
		unsigned int k = 0;
		for (unsigned int j = 0; j<thermodynamics_.number_of_gas_species(); j++)
			aux_vector__[k++] = Rgas[j];
		for (unsigned int j = 0; j < thermodynamics_.number_of_liquid_species(); j++)
			aux_vector__[k++] = RLiquid[j];

		return -Dot(static_cast<int>(aux_vector__.size()), aux_vector__.data(), thermodynamics_.Species_H_over_RT().data()) * PhysicalConstants::R_J_kmol * this->T_;
	}

	void KineticsMap_Liquid_CHEMKIN::ProductionAndDestructionRates(double* P, double* D)
	{
		stoichiometry_->ProductionAndDestructionRatesFromReactionRates(P, D, netReactionRates__.data());
	}
	
	const std::vector<double>& KineticsMap_Liquid_CHEMKIN::GiveMeReactionRates()
	{
		return netReactionRates__;
	}

	void KineticsMap_Liquid_CHEMKIN::CorrectReactionRate(const unsigned int j, const double Cc)
	{
		netReactionRates__[j-1] *= Cc;
		forwardReactionRates__[j - 1] *= Cc;
		reverseReactionRates__[j - 1] *= Cc;
	}

	void KineticsMap_Liquid_CHEMKIN::GiveMeReactionRates(double* r)
	{
		r = netReactionRates__.data();
	}

	void KineticsMap_Liquid_CHEMKIN::GetForwardReactionRates(double* r)
	{
		ElementByElementProduct(static_cast<int>(forwardReactionRates__.size()), forwardReactionRates__.data(), kArrheniusModified__.data(), r);
	}
 
	void KineticsMap_Liquid_CHEMKIN::GetBackwardReactionRates(double* r)
	{
		for (unsigned int j = 0; j < this->number_of_reactions_; j++)
			r[j] = 0.;

		for(unsigned int k=1;k<=number_of_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_reversible_reactions__[k-1]-1;
			r[j] = reverseReactionRates__[j]*kArrheniusModified__[j];
		}
		for(unsigned int k=1;k<=number_of_explicitly_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_explicitly_reversible_reactions__[k-1]-1;
			r[j] = reverseReactionRates__[j]*kArrheniusModified__[j];
		}                
	}

	void KineticsMap_Liquid_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k)
	{			
		thermodynamics_.SetPressure(101325.);
		thermodynamics_.SetTemperature(298.15);
		SetTemperature(298.15);
		SetPressure(101325.);
		ReactionEnthalpiesAndEntropies();
		
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << (reaction_h_over_RT__[k-1]-reaction_s_over_R__[k-1])*PhysicalConstants::R_kcal_mol*this->T_;	// [kcal/mol]
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_h_over_RT__[k-1] *PhysicalConstants::R_kcal_mol*this->T_;;						// [kcal/mol]
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_s_over_R__[k-1] *PhysicalConstants::R_cal_mol;;							// [cal/mol/K]
	}
 
	void KineticsMap_Liquid_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k, const double* c_bath, const double conversion_forward, const double conversion_backward)
	{			
		// TODO
	}

	void KineticsMap_Liquid_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k, const double* c_bath, const std::vector<double> list_of_temperatures, const double conversion_forward, const double conversion_backward)
	{
		const double patm = 101325.;
		SetPressure(patm);
		thermodynamics_.SetPressure(patm);

		int n_gas, n_liq;		
		n_gas = thermodynamics_.number_of_gas_species();
		n_liq = thermodynamics_.number_of_liquid_species();

		OpenSMOKE::OpenSMOKEVectorDouble Cgas(n_gas), Cliq(n_liq);

		for (int i = 1; i <= n_gas; i++)
			Cgas = c_bath[i];
		for (int i = 1; i <= n_liq; i++)
			Cliq= c_bath[i+n_gas];

		OpenSMOKEVectorDouble temperatures;

		if (list_of_temperatures.size() == 0)
		{
			OpenSMOKE::ChangeDimensions(5, &temperatures, true);
			temperatures[1] = 300.;
			temperatures[2] = 600.;
			temperatures[3] = 900.;
			temperatures[4] = 1200.;
			temperatures[5] = 1500.;
		}
		else
		{
			OpenSMOKE::ChangeDimensions(static_cast<int>(list_of_temperatures.size()), &temperatures, true);
			for (int i = 1; i <= temperatures.Size(); i++)
				temperatures[i] = list_of_temperatures[i - 1];
		}

		fOut << " -------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOut << "    Temperature   kF            Keq           kR            DG            DH            DS            kF            kR" << std::endl;
		fOut << "    [K]           [kmol,m3,s]   [-]           [kmol,m3,s]   [kcal/mol]    [kcal/mol]    [cal/mol/K]   [mol,cm3,s]   [mol,cm3,s]" << std::endl;
		fOut << " -------------------------------------------------------------------------------------------------------------------------------" << std::endl;

		for (int i = 1; i <= temperatures.Size(); i++)
		{
			SetTemperature(temperatures[i]);
			thermodynamics_.SetTemperature(temperatures[i]);

			ReactionEnthalpiesAndEntropies();
			KineticConstants();
			ReactionRates(Cgas.GetHandle(),Cliq.GetHandle());

			// Temperature
			fOut << "    " << std::setw(14) << std::left << std::setprecision(0) << std::fixed << temperatures[i];

			// Forward kinetic constant [kmol, m3, s]
			fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified__[k - 1];

			// Equilibrium and Backward kinetic constant [kmol, m3, s]
			if (isThermodynamicallyReversible__[k - 1] != 0)
			{
				const unsigned int j = isThermodynamicallyReversible__[k - 1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << 1. / uKeq__[j - 1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified__[k - 1] * uKeq__[j - 1];
			}
			else if (isExplicitlyReversible__[k - 1] != 0)
			{
				const unsigned int j = isExplicitlyReversible__[k - 1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified__[k - 1] / kArrhenius_reversible__[j - 1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrhenius_reversible__[j - 1];
			}
			else
			{
				fOut << std::setw(14) << std::left << std::setprecision(0) << std::fixed << 0.;
				fOut << std::setw(14) << std::left << std::setprecision(0) << std::fixed << 0.;
			}

			// Thermodynamic data
			fOut << std::setw(14) << std::left << std::setprecision(3) << std::fixed << (reaction_h_over_RT__[k - 1] - reaction_s_over_R__[k - 1]) * PhysicalConstants::R_kcal_mol * this->T_;	// [kcal/mol]
			fOut << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_h_over_RT__[k - 1] * PhysicalConstants::R_kcal_mol * this->T_;;						// [kcal/mol]
			fOut << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_s_over_R__[k - 1] * PhysicalConstants::R_cal_mol;;							// [cal/mol/K]

			// Forward kinetic constant [mol, cm3, s]
			fOut << std::setw(14) << std::left << std::setprecision(4) << std::scientific << kArrheniusModified__[k - 1] * conversion_forward;

			// Equilibrium and Backward kinetic constant [mol, cm3, s]
			if (isThermodynamicallyReversible__[k - 1] != 0)
			{
				const unsigned int j = isThermodynamicallyReversible__[k - 1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified__[k - 1] * uKeq__[j - 1] * conversion_backward;
			}
			else if (isExplicitlyReversible__[k - 1] != 0)
			{
				const unsigned int j = isExplicitlyReversible__[k - 1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrhenius_reversible__[j - 1] * conversion_backward;
			}
			else
			{
				fOut << std::setw(14) << std::left << std::setprecision(0) << std::fixed << 0.;
			}

			fOut << std::endl;
		}
		fOut << " -------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOut << std::endl;
	}

	void KineticsMap_Liquid_CHEMKIN::FittedReverseKineticConstants(const double* x_bath, const unsigned int nparameters, Eigen::MatrixXd& fittedKineticParameters)
	{			
		// TODO
	}

	void KineticsMap_Liquid_CHEMKIN::FittedReverseKineticConstants(const unsigned int k, std::ostream& fOut, Eigen::MatrixXd& fittedKineticParameters)
	{
		if (isThermodynamicallyReversible__[k-1] != 0)
		{
			const unsigned int j = isThermodynamicallyReversible__[k-1];

			fOut << std::setw(18) << std::right << std::scientific << std::setprecision(4) << std::exp(fittedKineticParameters(0, j-1));
			if (fittedKineticParameters.rows() == 2)
				fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << 0.;
			else
				fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << fittedKineticParameters(2, j - 1);
			fOut << std::setw(16) << std::right << std::fixed << std::setprecision(2) << fittedKineticParameters(1,j-1)/Conversions::J_from_kcal;
			fOut << std::setw(5)  << "";
		}
		else if (isExplicitlyReversible__[k-1] != 0)
		{
			const unsigned int j = isExplicitlyReversible__[k-1];
				
			fOut << std::setw(18) << std::right << std::scientific << std::setprecision(4) << std::exp(lnA_reversible__[j-1]);
			fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << Beta_reversible__[j-1];
			fOut << std::setw(16) << std::right << std::fixed << std::setprecision(2) << E_over_R_reversible__[j-1];
			fOut << std::setw(5)  << "";
		}
	}

	void KineticsMap_Liquid_CHEMKIN::RateOfProductionAnalysis(const bool iNormalize) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates__.data(), iNormalize);
	}

	void KineticsMap_Liquid_CHEMKIN::RateOfProductionAnalysis(std::ostream& fout) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates__.data(), false);
		stoichiometry_->WriteRateOfProductionAnalysis(fout);
	}

	void KineticsMap_Liquid_CHEMKIN::RateOfProductionAnalysis(ROPA_Data& ropa) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates__.data(), false);
		stoichiometry_->WriteRateOfProductionAnalysis(ropa);
	}

	void KineticsMap_Liquid_CHEMKIN::RateOfProductionAnalysis(ROPA_Data& ropa, const double* rf, const double* rb) const
	{
		stoichiometry_->RateOfProductionAnalysis(rf, rb);
		stoichiometry_->WriteRateOfProductionAnalysis(ropa);
	}
}

