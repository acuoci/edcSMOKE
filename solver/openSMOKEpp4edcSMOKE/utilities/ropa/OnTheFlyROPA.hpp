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
	template<class T>
	std::string FormatXMLSingleVariable(const std::string label, const T value)
	{
		return "<" + label + ">" + boost::lexical_cast<std::string>(value) + "</" + label + ">\n";
	}

	void ReorderROPACoefficients(const std::vector<unsigned int>& positive_indices, const std::vector<unsigned int>& negative_indices,
		const std::vector<double>& positive_coefficients, const std::vector<double>& negative_coefficients,
		std::vector<int>& reordered_positive_indices, std::vector<double>& reordered_positive_coefficients,
		std::vector<int>& reordered_negative_indices, std::vector<double>& reordered_negative_coefficients);

		OnTheFlyROPA::OnTheFlyROPA(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
									OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap) :
		
		thermodynamicsMap_(thermodynamicsMap),
		kineticsMap_(kineticsMap)
	
	{
		writing_policy_ = FIXED_TIME_STEPS_;
		number_steps_ = 25;
		merge_reaction_rates_ = false;
		reference_species_ = "none";
		is_active_ = false;
		do_the_analysis_ = false;
		threshold_contributions_ = 0.03;
		compact_mode_ = true;
		threshold_formation_rate_ = 1.e-32;

		next_step_ = number_steps_;
		next_index_conversion_ = 0;
		next_index_time_ = 0;
		current_step_ = 0;
		told_ = 0.;

		is_write_xml_ = false;
	}

	void OnTheFlyROPA::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, boost::filesystem::path path_kinetics_output)
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
				list_species_.resize(thermodynamicsMap_.NumberOfSpecies());
				list_species_ = thermodynamicsMap_.NamesOfSpecies();
			}
			else
			{	// check if the species is available in the mechanism
				int index_to_check;
				for (unsigned int j = 0; j < list_species_.size(); j++)
				{
					index_to_check = thermodynamicsMap_.IndexOfSpeciesWithoutError(list_species_[j]) - 1;
					if (index_to_check < 0)
					{
						std::cout << "\nSpecies " << list_species_[j] << " is not available in the mechanism" << std::endl;
						OpenSMOKE::FatalErrorMessage("Check the ROPA dictionary ");

					}
				}
			}
		}

		if (dictionary.CheckOption("@WriteXML") == true)
			dictionary.ReadBool("@WriteXML", is_write_xml_);

		Setup(path_kinetics_output);
	}

	void OnTheFlyROPA::Setup(const boost::filesystem::path& path_kinetics_output)
	{
		boost::filesystem::path file_name = path_kinetics_output / "reaction_names.xml";
		OpenSMOKE::ImportReactionNames(file_name, kineticsMap_.NumberOfReactions(), reaction_names_);

		OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &R_, true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &P_, true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &D_, true);
		OpenSMOKE::ChangeDimensions(kineticsMap_.NumberOfReactions(), &rf_, true);
		OpenSMOKE::ChangeDimensions(kineticsMap_.NumberOfReactions(), &rb_, true);
	}

	void OnTheFlyROPA::WriteHead(std::ofstream& fOut, const std::string reactor_type)
	{
		OpenSMOKE::OpenSMOKE_logo(fOut, "Rate of Production Analysis");

		fOut << "Frequency factors:      [kmol,m,s]" << std::endl;
		fOut << "Activation energies:    [cal/mol]" << std::endl;
		fOut << "Formation rates:        [kmol/m3/s]" << std::endl;
		fOut << "Reaction rates:         [kmol/m3/s]" << std::endl;
		fOut << "Compact mode:           " << compact_mode_ << std::endl;
		fOut << "Merge forward/backward: " << merge_reaction_rates_ << std::endl;
		fOut << "Threshold:              " << threshold_contributions_ * 100 << "%" << std::endl;
		fOut << std::endl;

		fOut << "Reactor: " << reactor_type << std::endl;
		fOut << "Date:    " << OpenSMOKE::GetCurrentDate() << std::endl;
		fOut << "Time:    " << OpenSMOKE::GetCurrentTime() << std::endl;
		fOut << std::endl;

		//	fOut << "Kinetic mechanism:   " << kineticsMap_.NameOfKineticMechanism() << std::endl;
		fOut << "Number of species:   " << thermodynamicsMap_.NumberOfSpecies() << std::endl;
		fOut << "Number of reactions: " << kineticsMap_.NumberOfReactions() << std::endl;
		fOut << std::endl;

		// Write XML file
		if (is_write_xml_ == true)
		{
			xml_string_.open("ROPA.xml", std::ios::out);
			WriteHeadXML(reactor_type);
		}
	}

	bool OnTheFlyROPA::CheckForROPA(const int current_iteration, const double time, OpenSMOKE::OpenSMOKEVectorDouble& conversions)
	{
		if (current_iteration == 0)
			return true;

		if (do_the_analysis_ == true)
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
			if (conversions[thermodynamicsMap_.IndexOfSpecies(reference_species_)] >= list_conversions_[next_index_conversion_])
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

	void OnTheFlyROPA::Analyze(std::ofstream& fOut, const int iteration, const double t, const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEVectorDouble& omega, const OpenSMOKE::OpenSMOKEVectorDouble& omega0)
	{
		// Calculate conversions
		OpenSMOKE::OpenSMOKEVectorDouble conversions(thermodynamicsMap_.NumberOfSpecies());
		for (unsigned int i = 1; i <= thermodynamicsMap_.NumberOfSpecies(); i++)
			conversions[i] = 1. - omega[i] / (omega0[i] + 1.e-32);

		if (CheckForROPA(iteration, t, conversions) == true)
		{
			// Update current step and previous time
			current_step_++;
			if (current_step_ == 1)	told_ = t;

			// Concentration [kmol/m3]
			const double cTot = c.SumElements();

			// Molecular weight [kg/kmol]
			double MW = 0.;
			for (unsigned int i = 1; i <= thermodynamicsMap_.NumberOfSpecies(); i++)
				MW += c[i] * thermodynamicsMap_.MW(i-1);
			MW /= cTot;

			// Density [kg/m3]
			const double rho = P*MW / PhysicalConstants::R_J_kmol / T;

			// Calculates thermodynamic properties
			thermodynamicsMap_.SetTemperature(T);
			thermodynamicsMap_.SetPressure(P);

			// Calculates kinetics
			kineticsMap_.SetTemperature(T);
			kineticsMap_.SetPressure(P);

			// Reaction rates
			kineticsMap_.KineticConstants();
			kineticsMap_.ReactionRates(c.GetHandle());
			kineticsMap_.ProductionAndDestructionRates(P_.GetHandle(), D_.GetHandle());

			// Rate of Production Analysis
			OpenSMOKE::ROPA_Data ropa;
			if (merge_reaction_rates_ == true)
			{
				kineticsMap_.FormationRates(R_.GetHandle());
				kineticsMap_.RateOfProductionAnalysis(ropa);
			}
			else
			{
				kineticsMap_.GetForwardReactionRates(rf_.GetHandle());
				kineticsMap_.GetBackwardReactionRates(rb_.GetHandle());
				kineticsMap_.RateOfProductionAnalysis(ropa, rf_.GetHandle(), rb_.GetHandle());
			}

			fOut << "***********************************************************************************************************" << std::endl;
			fOut << std::setw(14) << std::left << "Time[s]";
			fOut << std::setw(14) << std::left << "T[K]";
			fOut << std::setw(14) << std::left << "P[atm]";
			fOut << std::setw(14) << std::left << "rho[kg/m3]";
			fOut << std::setw(14) << std::left << "MW[kg/kmol]";
			fOut << std::setw(14) << std::left << "Conversion[%]";
			fOut << std::endl;
			fOut << std::setw(14) << std::left << std::scientific << std::setprecision(2) << t;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(2) << T;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(2) << P / 101325.;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(4) << rho;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(2) << MW;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(3) << conversions[thermodynamicsMap_.IndexOfSpecies(reference_species_)] * 100.;
			fOut << std::endl;
			fOut << "***********************************************************************************************************" << std::endl;
			fOut << std::endl;

			fOut << std::setw(30) << std::left << "Species";
			fOut << std::setw(16) << std::left << "Conc.[kmol/m3]";
			fOut << std::setw(16) << std::left << "Mole fract.";
			fOut << std::setw(16) << std::left << "Mass fract.";
			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------------------------------------" << std::endl;

			for (unsigned int i = 1; i <= thermodynamicsMap_.NumberOfSpecies(); i++)
			{
				if (compact_mode_ == false)
				{
					fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[i - 1];
					fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << c[i];
					fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << c[i] / cTot;
					fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << omega[i];
					fOut << std::endl;
				}
				else
				{
					if (omega[i] >= 1.e-16)
					{
						fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[i - 1];
						fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << c[i];
						fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << c[i] / cTot;
						fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << omega[i];
						fOut << std::endl;
					}
				}
			}
			fOut << std::endl;
			fOut << std::endl;

			// Write xml
			if (is_write_xml_ == true)
			{
				xml_string_ << "<ropa " << "time=" << "\"" << t << "\" "
										<< "weight=" << "\"" << t-told_ << "\" "
										<< "id=" << "\"" << current_step_ << "\" "
										<< "T=" << "\"" << T << "\" "
										<< "P=" << "\"" << P/101325. << "\">"
										<< std::endl;
			}

			// Loop over all the species
			for (unsigned int i = 0; i < list_species_.size(); i++)
			{
				const unsigned int index_of_species = thermodynamicsMap_.IndexOfSpecies(list_species_[i]) - 1;

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
					fOut << std::endl;
					fOut << "-----------------------------------------------------------------------------------------------------------" << std::endl;

					// Production contribution
					const double sum_positive_coefficients = std::accumulate(reorder_positive_coefficients.begin(), reorder_positive_coefficients.end(), 0.);
					for (unsigned int i = 0; i < reordered_positive_indices.size(); i++)
					{
						const double percentage = reorder_positive_coefficients[i] / sum_positive_coefficients;
						if (percentage > threshold_contributions_)
						{
							unsigned int j = reordered_positive_indices[i];
							fOut << std::setw(6) << std::left << j;
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << reorder_positive_coefficients[i];
							fOut << std::setw(10) << std::right << std::setprecision(1) << std::fixed << percentage*100. << "%";
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << kineticsMap_.A(j);
							fOut << std::setw(8) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.Beta(j);
							fOut << std::setw(12) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.E_over_R(j)* PhysicalConstants::R_cal_mol;
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
							fOut << std::setw(6) << std::left << j;
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << -reorder_negative_coefficients[i];
							fOut << std::setw(10) << std::right << std::setprecision(1) << std::fixed << -percentage*100. << "%";
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << kineticsMap_.A(j);
							fOut << std::setw(8) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.Beta(j);
							fOut << std::setw(12) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.E_over_R(j)* PhysicalConstants::R_cal_mol;
							fOut << std::left << "   " << reaction_names_[j] << std::endl;
						}
					}

					fOut << std::endl;
				}

				// Write XML
				if (is_write_xml_ == true)
				{
					// Main tag
					xml_string_ << "<species "	<< "name=" << "\"" << thermodynamicsMap_.NamesOfSpecies()[index_of_species] << "\" " 
												<< "id=" << "\"" << index_of_species+1 << "\">" << std::endl;

					// Production contribution
					const double sum_positive_coefficients = std::accumulate(reorder_positive_coefficients.begin(), reorder_positive_coefficients.end(), 0.);
					xml_string_ << FormatXMLSingleVariable("tot-form", sum_positive_coefficients);
					for (unsigned int i = 0; i < reordered_positive_indices.size(); i++)
					{
						const double percentage = reorder_positive_coefficients[i] / sum_positive_coefficients;
						if (percentage > threshold_contributions_)
						{
							unsigned int j = reordered_positive_indices[i];

							xml_string_ << "<form "
										<< "id=" << "\"" << j << "\" "
										<< "react=" << "\"" << reaction_names_[j] << "\">"
										<< reorder_positive_coefficients[i]
										<< "</form>" << std::endl;
						}
					}

					// Consumption contributions
					const double sum_negative_coefficients = std::accumulate(reorder_negative_coefficients.begin(), reorder_negative_coefficients.end(), 0.);
					xml_string_ << FormatXMLSingleVariable("tot-consumption", sum_negative_coefficients);
					for (unsigned int i = 0; i < reordered_negative_indices.size(); i++)
					{
						const double percentage = reorder_negative_coefficients[i] / sum_negative_coefficients;
						if (percentage > threshold_contributions_)
						{
							unsigned int j = reordered_negative_indices[i];

							xml_string_ << "<cons "
										<< "id=" << "\"" << j << "\" "
										<< "react=" << "\"" << reaction_names_[j] << "\">"
										<< reorder_negative_coefficients[i]
										<< "</cons>" << std::endl;
						}
					}

					xml_string_ << "</species>" << std::endl;
				}
			}

			if (is_write_xml_ == true)
			{
				xml_string_ << "</ropa>" << std::endl;
			}

			// Store previous time
			told_ = t;
		}
	}

	void ReorderROPACoefficients(	const std::vector<unsigned int>& positive_indices, const std::vector<unsigned int>& negative_indices,
									const std::vector<double>& positive_coefficients, const std::vector<double>& negative_coefficients,
									std::vector<int>& reordered_positive_indices, std::vector<double>& reordered_positive_coefficients,
									std::vector<int>& reordered_negative_indices, std::vector<double>& reordered_negative_coefficients )
	{
		// Positive coefficients
		reordered_positive_indices.resize(positive_indices.size());
		reordered_positive_coefficients.resize(positive_indices.size());
		for (unsigned int i = 0; i<positive_coefficients.size(); i++)
		{
			reordered_positive_indices[i] = positive_indices[i];
			reordered_positive_coefficients[i] = positive_coefficients[i];
		}
		OpenSMOKE_Utilities::ReorderPairsOfVectors(reordered_positive_coefficients, reordered_positive_indices);
		std::reverse(reordered_positive_indices.begin(), reordered_positive_indices.end());
		std::reverse(reordered_positive_coefficients.begin(), reordered_positive_coefficients.end());

		// Negative coefficients
		reordered_negative_indices.resize(negative_indices.size());
		reordered_negative_coefficients.resize(negative_indices.size());
		for (unsigned int i = 0; i<negative_coefficients.size(); i++)
		{
			reordered_negative_indices[i] = negative_indices[i];
			reordered_negative_coefficients[i] = -negative_coefficients[i];
		}
		OpenSMOKE_Utilities::ReorderPairsOfVectors(reordered_negative_coefficients, reordered_negative_indices);
		//	std::reverse(reordered_negative_indices.begin(), reordered_negative_indices.end());
		//	std::reverse(reordered_negative_coefficients.begin(), reordered_negative_coefficients.end());
	}

	void OnTheFlyROPA::WriteHeadXML(const std::string reactor_type)
	{
		xml_string_ << std::setprecision(8);
		xml_string_.setf(std::ios::scientific);
		xml_string_ << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
		xml_string_ << "<opensmoke version=\"0.1a\">" << std::endl;

		xml_string_ << FormatXMLSingleVariable("AUnits", "[kmol, m, s]");
		xml_string_ << FormatXMLSingleVariable("EUnits", "[cal/mol]");
		xml_string_ << FormatXMLSingleVariable("omegaUnits", "[kmol/m3/s]");
		xml_string_ << FormatXMLSingleVariable("rUnits", "[kmol/m3/s]");
		xml_string_ << FormatXMLSingleVariable("TUnits", "[K]");
		xml_string_ << FormatXMLSingleVariable("PUnits", "[atm]");
		xml_string_ << FormatXMLSingleVariable("Compact", compact_mode_);
		xml_string_ << FormatXMLSingleVariable("Merging", merge_reaction_rates_);
		xml_string_ << FormatXMLSingleVariable("Threshold", threshold_contributions_);

		xml_string_ << FormatXMLSingleVariable("Reactor", reactor_type);
		xml_string_ << FormatXMLSingleVariable("Date", OpenSMOKE::GetCurrentDate());
		xml_string_ << FormatXMLSingleVariable("Time", OpenSMOKE::GetCurrentTime());

		xml_string_ << FormatXMLSingleVariable("NumberOfSpecies", thermodynamicsMap_.NumberOfSpecies());
		xml_string_ << FormatXMLSingleVariable("NumberOfReactions", kineticsMap_.NumberOfReactions());
	}

	void OnTheFlyROPA::Close()
	{
		if (is_write_xml_ == true)
		{
			xml_string_ << "</opensmoke>" << std::endl;
			xml_string_.close();
		}
	}
}

