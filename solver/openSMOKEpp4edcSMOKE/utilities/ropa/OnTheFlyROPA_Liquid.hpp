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

#include "Grammar_OnTheFlyROPA.h"

namespace OpenSMOKE
{
	void ReorderROPACoefficients(const std::vector<unsigned int>& positive_indices, const std::vector<unsigned int>& negative_indices,
		const std::vector<double>& positive_coefficients, const std::vector<double>& negative_coefficients,
		std::vector<int>& reordered_positive_indices, std::vector<double>& reordered_positive_coefficients,
		std::vector<int>& reordered_negative_indices, std::vector<double>& reordered_negative_coefficients);

		OnTheFlyROPA_Liquid::OnTheFlyROPA_Liquid(OpenSMOKE::ThermodynamicsMap_Liquid_CHEMKIN& thermodynamicsMap,
			OpenSMOKE::KineticsMap_Liquid_CHEMKIN& kineticsMap) :
		
		thermodynamicsMap_(thermodynamicsMap),
		kineticsMap_(kineticsMap)
	
	{
		writing_policy_ = FIXED_TIME_STEPS_;
		number_steps_ = 25;
		merge_reaction_rates_ = false;
		total_liquid_conversion = false;
		reference_species_ = "none";
		is_active_ = false;
		am_i_analyzing = false;
		threshold_contributions_ = 0.03;
		compact_mode_ = true;
		threshold_formation_rate_ = 1.e-32;

		next_step_ = number_steps_;
		next_index_conversion_ = 0;
		next_index_time_ = 0;
	}

	void OnTheFlyROPA_Liquid::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, boost::filesystem::path path_kinetics_output)
	{
		is_active_ = true;

		Grammar_OnTheFlyROPA grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@MergeForwardAndBackwardReactions") == true)
			dictionary.ReadBool("@MergeForwardAndBackwardReactions", merge_reaction_rates_);

		if (dictionary.CheckOption("@CompactOutput") == true)
			dictionary.ReadBool("@CompactOutput", compact_mode_);

		if (dictionary.CheckOption("@NumberOfSteps") == true)
		{
			dictionary.ReadInt("@NumberOfSteps", number_steps_);
			next_step_ = number_steps_;
			writing_policy_ = FIXED_TIME_STEPS_;
		}

		if (dictionary.CheckOption("@Conversions") == true)
		{
			dictionary.ReadOption("@Conversions", list_conversions_);
			list_conversions_.push_back(1.e64);
			writing_policy_ = LIST_CONVERSIONS;
		}

		if (dictionary.CheckOption("@Times") == true)
		{
			std::vector<std::string> list_;
			dictionary.ReadOption("@Times", list_);

			list_times_.resize(list_.size() - 1);
			const std::string units = list_[list_.size() - 1];
			for (unsigned int i = 0; i < list_.size() - 1; i++)
			{
				if (units == "s")
					list_times_[i] = boost::lexical_cast<double>(list_[i]);
				else if (units == "ms")
					list_times_[i] = boost::lexical_cast<double>(list_[i]) / 1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown time units.");
			}
			list_times_.push_back(1.e64);
			writing_policy_ = LIST_TIMES;
		}

		if (dictionary.CheckOption("@Threshold") == true)
			dictionary.ReadDouble("@Threshold", threshold_contributions_);

		if (dictionary.CheckOption("@ReferenceSpecies") == true)
			dictionary.ReadString("@ReferenceSpecies", reference_species_);

		if (dictionary.CheckOption("@Species") == true)
		{
			dictionary.ReadOption("@Species", list_species_);

			if (list_species_[0] == "ALL" || list_species_[0] == "all")
			{
				list_species_.resize(thermodynamicsMap_.number_of_liquid_species());
				list_species_ = thermodynamicsMap_.vector_names_liquid_species();
			}
			else
			{	// check if the species is available in the mechanism
				int index_to_check;
				for (unsigned int j = 0; j < list_species_.size(); j++)
				{
					index_to_check = thermodynamicsMap_.IndexOfSpeciesWithoutError(list_species_[j]) - 1;
					if (index_to_check <0)
					{
						std::cout << "\nSpecies " << list_species_[j] << " is not available in the mechanism" << std::endl;
						OpenSMOKE::FatalErrorMessage("Check the ROPA dictionary ");
						
					}
				}
			}
		}
		
		if (dictionary.CheckOption("@ReferenceLiquidMass") == true)
			total_liquid_conversion = true;
		
		index_ref_species_ = thermodynamicsMap_.IndexOfSpecies(reference_species_ + "(L)") - thermodynamicsMap_.number_of_gas_species();

		Setup(path_kinetics_output);
	}

	void OnTheFlyROPA_Liquid::Setup(const boost::filesystem::path& path_kinetics_output)
	{
		boost::filesystem::path file_name = path_kinetics_output / "reaction_names.liquid.xml";
		OpenSMOKE::ImportReactionNames(file_name, kineticsMap_.NumberOfReactions(), reaction_names_);

		OpenSMOKE::ChangeDimensions(thermodynamicsMap_.number_of_gas_species(), &Rg_, true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMap_.number_of_liquid_species(), &Rl_, true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &P_, true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &D_, true);
		OpenSMOKE::ChangeDimensions(kineticsMap_.NumberOfReactions(), &rf_, true);
		OpenSMOKE::ChangeDimensions(kineticsMap_.NumberOfReactions(), &rb_, true);
	}

	void OnTheFlyROPA_Liquid::WriteHead(std::ofstream& fOut, const std::string reactor_type)
	{
		OpenSMOKE::OpenSMOKE_logo(fOut, "Rate of Production Analysis");

		fOut << "Frequency factors:      [kmol,m,s]" << std::endl;
		fOut << "Activation energies:    [cal/mol]" << std::endl;
		fOut << "Formation rates:        [kmol/m3/s]" << std::endl;
		fOut << "Reaction rates:         [kmol/m3/s]" << std::endl;
		fOut << "Evaporation rates:      [kmol/m3/s]" << std::endl;
		fOut << "Compact mode:           " << compact_mode_ << std::endl;
		fOut << "Merge forward/backward: " << merge_reaction_rates_ << std::endl;
		fOut << "Threshold:              " << threshold_contributions_ * 100 << "%" << std::endl;
		fOut << std::endl;

		fOut << "Reactor: " << reactor_type << std::endl;
		fOut << "Date:    " << OpenSMOKE::GetCurrentDate() << std::endl;
		fOut << "Time:    " << OpenSMOKE::GetCurrentTime() << std::endl;
		fOut << std::endl;
		
		fOut << "Number of total species:	 " << thermodynamicsMap_.NumberOfSpecies() << std::endl;
		fOut << "Number of liquid species:   " << thermodynamicsMap_.number_of_liquid_species() << std::endl;
		fOut << "Number of liquid reactions: " << kineticsMap_.NumberOfReactions() << std::endl;
		fOut << std::endl;
	}

	bool OnTheFlyROPA_Liquid::CheckForROPA(const int current_iteration, const double time, OpenSMOKE::OpenSMOKEVectorDouble& conversions)
	{
		if (current_iteration == 0)
			return true;

		if (writing_policy_ == FIXED_TIME_STEPS_)
		{
			if (current_iteration >= next_step_)
			{
				next_step_ += number_steps_;
				return true;
			}
		}
		else if (writing_policy_ == LIST_CONVERSIONS)
		{
			double X_to_check = 0.;
			if (total_liquid_conversion == true)
				X_to_check = conversions[conversions.Size()];
			else
				X_to_check = conversions[index_ref_species_];
			if (X_to_check >= list_conversions_[next_index_conversion_])
			{
				next_index_conversion_++;
				return true;
			}
		}
		else if (writing_policy_ == LIST_TIMES)
		{
			if (time >= list_times_[next_index_time_])
			{
				next_index_time_++;
				return true;
			}
		}

		return false;
	}	

	void OnTheFlyROPA_Liquid::Analyze(std::ofstream& fOut, const int iteration, const double t, const double T, const double P,
		const OpenSMOKE::OpenSMOKEVectorDouble& cGas, const OpenSMOKE::OpenSMOKEVectorDouble& cLiq,
		const OpenSMOKE::OpenSMOKEVectorDouble& massL, const OpenSMOKE::OpenSMOKEVectorDouble& omega0,
		const OpenSMOKE::OpenSMOKEVectorDouble& mL_ev)
	{
		OpenSMOKE::OpenSMOKEVectorDouble conversions(thermodynamicsMap_.number_of_liquid_species() + 1);
		for (unsigned int i = 1; i <= thermodynamicsMap_.number_of_liquid_species(); i++)
			conversions[i] = 1. - massL[i] / (omega0[i] + 1.e-32);

		conversions[thermodynamicsMap_.number_of_liquid_species() + 1] = 1. - massL.SumAbsElements() / (omega0.SumAbsElements() + 1.e-32);

		am_i_analyzing = false;

		if (CheckForROPA(iteration, t, conversions) == true)
		{
			am_i_analyzing = true;
			int NSL_ = thermodynamicsMap_.number_of_liquid_species();
			int NSG_ = thermodynamicsMap_.number_of_gas_species();

			// Gas-phase
			// Concentration [kmol/m3]
			const double cGTot = cGas.SumElements();

			// Molecular weight [kg/kmol]
			double MW = 0.;
			for (unsigned int i = 1; i <= NSG_; i++)
				MW += cGas[i] * thermodynamicsMap_.MW(i - 1);
			MW /= cGTot;

			// Liquid-phase
			// Concentration [kmol/m3]
			const double cLTot = cLiq.SumElements();

			// Molecular weight [kg/kmol]
			double MWL = 0.;
			for (unsigned int i = 1; i <= NSL_; i++)
				MWL += cLiq[index_Lmixture_Lkin_(i - 1) + 1] * thermodynamicsMap_.MW(index_of_fuel_species_(i - 1) - 1);

			MWL /= cLTot;

			//const double rho = P*MW / PhysicalConstants::R_J_kmol / T;

			// Calculates thermodynamic properties
			thermodynamicsMap_.SetTemperature(T);
			thermodynamicsMap_.SetPressure(P);

			// Calculates kinetics
			kineticsMap_.SetTemperature(T);
			kineticsMap_.SetPressure(P);

			// Reaction rates
			kineticsMap_.KineticConstants();
			kineticsMap_.ReactionRates(cGas.GetHandle(), cLiq.GetHandle());
			kineticsMap_.ProductionAndDestructionRates(P_.GetHandle(), D_.GetHandle());

			// Rate of Production Analysis
			OpenSMOKE::ROPA_Data ropa;
			if (merge_reaction_rates_ == true)
			{
				kineticsMap_.FormationRates(Rg_.GetHandle(), Rl_.GetHandle());
				kineticsMap_.RateOfProductionAnalysis(ropa);
			}
			else
			{
				kineticsMap_.GetForwardReactionRates(rf_.GetHandle());
				kineticsMap_.GetBackwardReactionRates(rb_.GetHandle());
				kineticsMap_.RateOfProductionAnalysis(ropa, rf_.GetHandle(), rb_.GetHandle());
			}

			fOut << "********************************************************************************************************************" << std::endl;
			fOut << std::setw(14) << std::left << "Time[s]";
			fOut << std::setw(14) << std::left << "T[K]";
			fOut << std::setw(14) << std::left << "P[atm]";			
			fOut << std::setw(14) << std::left << "MW[kg/kmol]";
			fOut << std::setw(14) << std::left << "Conversion[%]";
			fOut << std::endl;
			fOut << std::setw(14) << std::left << std::scientific << std::setprecision(2) << t;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(2) << T;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(2) << P / 101325.;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(2) << MWL;
			if (total_liquid_conversion == true)
				fOut << std::setw(14) << std::left << std::fixed << std::setprecision(3) << conversions[conversions.Size()] * 100.;
			else
				fOut << std::setw(14) << std::left << std::fixed << std::setprecision(3) << conversions[index_ref_species_] * 100.;
			fOut << std::endl;
			fOut << "********************************************************************************************************************" << std::endl;
			fOut << std::endl;

			fOut << std::setw(30) << std::left << "Species";
			fOut << std::setw(16) << std::left << "Conc.[kmol/m3]";
			fOut << std::setw(16) << std::left << "Mole fract.";
			fOut << std::setw(16) << std::left << "Mass [kg/kg0]";
			fOut << std::endl;
			fOut << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

			for (unsigned int i = 1; i <= NSL_; i++)
			{				
				if (compact_mode_ == false || massL[i] >= 1.e-16)
				{
					const unsigned int index_of_species = index_of_fuel_species_(i - 1) - 1;
					fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[index_of_species];
					fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << cLiq[i];
					fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << cLiq[i] / cLTot;
					fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << massL[i];
					fOut << std::endl;
				}
			}
			fOut << std::endl;
			fOut << std::endl;

			int index_of_species = 0;
			// Loop over all the species
			for (unsigned int i = 0; i < list_species_.size(); i++)
			{
				
				//index_of_species = index_of_liquid_species(thermodynamicsMap_.IndexOfSpecies(list_species_[i]) - 1);
				index_of_species = thermodynamicsMap_.IndexOfSpeciesWithoutError(list_species_[i] + "(L)") - 1;
				if (index_of_species == -1)
					index_of_species = thermodynamicsMap_.IndexOfSpeciesWithoutError(list_species_[i]) - 1;
				//const unsigned int index_of_species = thermodynamicsMap_.IndexOfSpecies(list_species_[i] + "(L)") - 1; 
				
				// Reorder species
				std::vector<int> reordered_positive_indices;
				std::vector<double> reorder_positive_coefficients;
				std::vector<int> reordered_negative_indices;
				std::vector<double> reorder_negative_coefficients;
				ReorderROPACoefficients(ropa.production_reaction_indices[index_of_species],
					ropa.destruction_reaction_indices[index_of_species],
					ropa.production_coefficients[index_of_species], ropa.destruction_coefficients[index_of_species],
					reordered_positive_indices, reorder_positive_coefficients,
					reordered_negative_indices, reorder_negative_coefficients);

				if (P_[index_of_species + 1] >= threshold_formation_rate_ || D_[index_of_species + 1] >= threshold_formation_rate_)
				{
					// Write 
					fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[index_of_species];
					fOut << "Formation: " << std::setw(12) << std::left << std::setprecision(2) << std::scientific << P_[index_of_species + 1];
					fOut << "Consumption: " << std::setw(12) << std::left << std::setprecision(2) << std::scientific << -D_[index_of_species + 1];
					fOut << "Net: " << std::setw(12) << std::left << std::setprecision(2) << std::scientific << P_[index_of_species + 1] - D_[index_of_species + 1];
					if (thermodynamicsMap_.IndexOfSpeciesWithoutError(list_species_[i] + "(L)") > 0)
						fOut << "Evaporation: " << std::setw(12) << std::left << std::setprecision(2) << std::scientific << mL_ev[index_of_species - NSG_ + 1];
					fOut << std::endl;
					fOut << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

					// Production contribution
					const double sum_positive_coefficients = std::accumulate(reorder_positive_coefficients.begin(), reorder_positive_coefficients.end(), 0.);
					for (unsigned int i = 0; i < reordered_positive_indices.size(); i++)
					{
						const double percentage = reorder_positive_coefficients[i] / sum_positive_coefficients;
						if (percentage > threshold_contributions_)
						{
							unsigned int j = reordered_positive_indices[i];
							fOut << std::setw(6) << std::left << j + 1;
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << reorder_positive_coefficients[i];
							fOut << std::setw(10) << std::right << std::setprecision(1) << std::fixed << percentage * 100. << "%";
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << kineticsMap_.A(j);
							fOut << std::setw(8) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.Beta(j);
							fOut << std::setw(12) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.E_over_R(j) * PhysicalConstants::R_cal_mol;
							fOut << std::left << "   " << reaction_names_[j] << std::endl;
						}
					}

					// Consumption contributions
					const double sum_negative_coefficients = std::accumulate(reorder_negative_coefficients.begin(), reorder_negative_coefficients.end(), 0.);
					for (unsigned int i = 0; i < reordered_negative_indices.size(); i++)
					{
						const double percentage = reorder_negative_coefficients[i] / sum_negative_coefficients;
						if (percentage > threshold_contributions_)
						{
							unsigned int j = reordered_negative_indices[i];
							fOut << std::setw(6) << std::left << j + 1;
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << -reorder_negative_coefficients[i];
							fOut << std::setw(10) << std::right << std::setprecision(1) << std::fixed << -percentage * 100. << "%";
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << kineticsMap_.A(j);
							fOut << std::setw(8) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.Beta(j);
							fOut << std::setw(12) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.E_over_R(j) * PhysicalConstants::R_cal_mol;
							fOut << std::left << "   " << reaction_names_[j] << std::endl;
						}
					}

					fOut << std::endl;					
				}
			}
					
		}
	}

	void OnTheFlyROPA_Liquid::Print_Evaporation(std::ofstream& fOut, const OpenSMOKE::OpenSMOKEVectorDouble& mL_ev)
	{		
		
		{			
			int NSL_ = thermodynamicsMap_.number_of_liquid_species();
			int NSG_ = thermodynamicsMap_.number_of_gas_species();						
			
			std::vector<std::vector<unsigned int>>	all_indices;
			std::vector<std::string>				all_names;
			all_names = thermodynamicsMap_.NamesOfSpecies();

			std::vector<unsigned int>				dummy_prod_indices(NSL_,0);
			std::vector<double>						dummy_prod_coefficients(NSL_, -1.);
			std::vector<double>						percentaged_mevL(NSL_);
			
			for (int k = 0; k < NSL_; k++)
				percentaged_mevL[k] = mL_ev[k + 1];
			

			all_indices = thermodynamicsMap_.matrix_indices_liquid_species();
			if (all_indices.size() > 1)
			{
				std::cout << "Evaporation ROPA not yet implemented for more than one Liquid-phase" << std::endl;
				OpenSMOKE::FatalErrorMessage("Not yet implemented");
			}
			for (int k = 0;k< thermodynamicsMap_.number_of_materials();k++)
			{				

				// Reorder species
				std::vector<int> reordered_positive_indices;
				std::vector<double> reorder_positive_coefficients;
				std::vector<int> reordered_negative_indices;
				std::vector<double> reorder_negative_coefficients;
				ReorderROPACoefficients(dummy_prod_indices, all_indices[k],
					dummy_prod_coefficients, percentaged_mevL,
					reordered_positive_indices, reorder_positive_coefficients,
					reordered_negative_indices, reorder_negative_coefficients);
				
				if (mL_ev.SumAbsElements() > threshold_formation_rate_)
				{
					// Write 
					fOut << std::setw(30) << std::left << "Total Liquid Mass variation";
					fOut << "From evaporation model [kg/s/m3]: " << std::setw(12) << std::left << std::setprecision(2) << std::scientific << -mL_ev.SumAbsElements();
					fOut << std::endl;
					fOut << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

					// Consumption contributions
					const double sum_negative_coefficients = std::accumulate(reorder_negative_coefficients.begin(), reorder_negative_coefficients.end(), 0.);
					for (unsigned int i = 0; i < reordered_negative_indices.size(); i++)
					{
						const double percentage = reorder_negative_coefficients[i] / sum_negative_coefficients;
						if (percentage > threshold_contributions_)
						{
							unsigned int j = reordered_negative_indices[i];
							fOut << std::setw(6) << std::left << j + 1;
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << -reorder_negative_coefficients[i];
							fOut << std::setw(10) << std::right << std::setprecision(1) << std::fixed << -percentage * 100. << "%";
							fOut << std::left << "   " << all_names[index_of_fuel_species_(j - 1) - 1] << std::endl;
						}
					}

					fOut << std::endl;
				}
			}
		}
	}

	void OnTheFlyROPA_Liquid::Print_Evaporation(std::ofstream& fOut, const OpenSMOKE::OpenSMOKEVectorDouble& mL_ev, const OpenSMOKE::OpenSMOKEVectorDouble& Omega_Gas_from_Liquid)
	{

		{
			int NSL_ = thermodynamicsMap_.number_of_liquid_species();
			int NSG_ = thermodynamicsMap_.number_of_gas_species();

			std::vector<unsigned int>				all_indices(NSG_);
			std::vector<std::string>				all_names;
			all_names = thermodynamicsMap_.NamesOfSpecies();

			std::vector<unsigned int>				dummy_prod_indices(NSG_, 0);
			std::vector<double>						dummy_prod_coefficients(NSG_, -1.);
			std::vector<double>						percentaged_mevL(NSG_, 0.);

			
			for (unsigned int k = 0; k < NSG_; k++)
			{
				percentaged_mevL[k] = -Omega_Gas_from_Liquid[k + 1];
				percentaged_mevL[k] += mL_ev[k + 1];
				all_indices[k] = k + 1;
			}
			
			if (thermodynamicsMap_.number_of_materials() > 1)
			{
				std::cout << "Evaporation ROPA not really implemented for more than one Liquid-phase" << std::endl;
				OpenSMOKE::FatalErrorMessage("Not yet implemented");
			}
			for (int k = 0; k < thermodynamicsMap_.number_of_materials(); k++)
			{

				// Reorder species
				std::vector<int> reordered_positive_indices;
				std::vector<double> reorder_positive_coefficients;
				std::vector<int> reordered_negative_indices;
				std::vector<double> reorder_negative_coefficients;
				ReorderROPACoefficients(dummy_prod_indices, all_indices,
					dummy_prod_coefficients, percentaged_mevL,
					reordered_positive_indices, reorder_positive_coefficients,
					reordered_negative_indices, reorder_negative_coefficients);

				if (mL_ev.SumAbsElements() > threshold_formation_rate_)
				{
					// Write 
					fOut << std::setw(30) << std::left << "Total Liquid Mass variation";
					fOut << "From evaporation model [kg/s/m3]: " << std::setw(12) << std::left << std::setprecision(2) << std::scientific << -mL_ev.SumAbsElements();
					fOut << "From reaction rates [kg/s/m3]: " << std::setw(12) << std::left << std::setprecision(2) << std::scientific << -Omega_Gas_from_Liquid.SumAbsElements();
					fOut << std::endl;
					fOut << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

					// Consumption contributions
					const double sum_negative_coefficients = std::accumulate(reorder_negative_coefficients.begin(), reorder_negative_coefficients.end(), 0.);
					for (unsigned int i = 0; i < reordered_negative_indices.size(); i++)
					{
						const double percentage = reorder_negative_coefficients[i] / sum_negative_coefficients;
						if (percentage > threshold_contributions_)
						{
							unsigned int j = reordered_negative_indices[i];
							fOut << std::setw(6) << std::left << j + 1;
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << -reorder_negative_coefficients[i];
							fOut << std::setw(10) << std::right << std::setprecision(1) << std::fixed << -percentage * 100. << "%";
							fOut << std::left << "   " << all_names[j - 1] << std::endl;
						}
					}

					fOut << std::endl;
				}
			}
		}
	}

	
	void OnTheFlyROPA_Liquid::SetIndexLiquidSpecies(const Eigen::VectorXi& index_fuel_species, const Eigen::VectorXi& index_Lmixture_Lkin)
	{
		index_of_fuel_species_ = index_fuel_species;
		index_Lmixture_Lkin_ = index_Lmixture_Lkin;
	}

}

