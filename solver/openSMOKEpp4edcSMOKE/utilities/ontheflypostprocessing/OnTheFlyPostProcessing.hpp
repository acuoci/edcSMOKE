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

#include "Grammar_OnTheFlyPostProcessing.h"
#include <boost/algorithm/string.hpp>

namespace OpenSMOKE
{

	OnTheFlyPostProcessing::OnTheFlyPostProcessing(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
													OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
													const boost::filesystem::path path_output) :

		thermodynamicsMap_(thermodynamicsMap),
		kineticsMap_(kineticsMap)
	{
		path_output_ = path_output;
		is_active_ = false;
		print_formation_rates_ = false;
		print_reaction_rates_ = false;
		formation_rates_moles_ = true;
		class_grouping_ = false;
	}

	void OnTheFlyPostProcessing::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_OnTheFlyPostProcessing grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@FormationRates") == true)
		{
			std::vector<std::string> formation_rates_species;
			dictionary.ReadOption("@FormationRates", formation_rates_species);

			if (formation_rates_species[0] == "ALL" || formation_rates_species[0] == "all")
			{
				print_formation_rates_ = true;
			}
			else if (formation_rates_species[0] == "NONE" || formation_rates_species[0] == "none")
			{
				print_formation_rates_ = false;
			}
			else
			{
				print_formation_rates_ = true;
				indices_of_formation_rates_species_.resize(formation_rates_species.size());
				for (unsigned int i = 0; i < formation_rates_species.size(); i++)
					indices_of_formation_rates_species_[i] = thermodynamicsMap_.IndexOfSpecies(formation_rates_species[i]);
			}
		}

		if (dictionary.CheckOption("@ReactionRates") == true)
		{
			std::vector<std::string> output_reaction_rates;
			dictionary.ReadOption("@ReactionRates", output_reaction_rates);

			if (output_reaction_rates[0] == "ALL" || output_reaction_rates[0] == "all")
			{
				print_reaction_rates_ = true;
			}
			else if (output_reaction_rates[0] == "NONE" || output_reaction_rates[0] == "none")
			{
				print_reaction_rates_ = false;
			}
			else
			{
				print_reaction_rates_ = true;
				indices_of_reaction_rates_.resize(output_reaction_rates.size());
				for (unsigned int i = 0; i < output_reaction_rates.size(); i++)
					indices_of_reaction_rates_[i] = boost::lexical_cast<unsigned int>(output_reaction_rates[i]);
			}
		}

		if (dictionary.CheckOption("@FormationRatesUnits") == true)
		{
			std::string formation_rates_units;
			dictionary.ReadString("@FormationRatesUnits", formation_rates_units);

			if (formation_rates_units == "mass")
				formation_rates_moles_ = false;
			else if (formation_rates_units == "moles")
				formation_rates_moles_ = true;
			else
				OpenSMOKE::FatalErrorMessage("@FormationRatesUnits: wrong definition. Allowed options are: mass | moles (default)");
		}

		if (dictionary.CheckOption("@ClassGrouping") == true)
		{
			boost::filesystem::path input_file;
			dictionary.ReadPath("@ClassGrouping", input_file);

			if (boost::filesystem::exists(input_file))
			{
				std::ifstream fInput(input_file.c_str(), std::ios::in);

				std::vector<std::string> list_families;
				std::vector<std::string> list_species;

				std::string line;
				unsigned int iline = 0;
				while (!fInput.eof())
				{
					iline++;

					getline(fInput, line);
					std::stringstream ss(line);

					if (boost::trim_copy(line) != "")
					{
						std::string family = "";
						std::string species = "";
						ss >> species;
						ss >> family;

						if (family != "" && species != "")
						{
							for (int i = 0; i < list_species.size(); i++)
								if (boost::to_upper_copy(species) == list_species[i])
								{
									std::string message = "@ClassGrouping list: the following species is specified more than once: " + species;
									OpenSMOKE::FatalErrorMessage(message);
								}

							list_species.push_back(species);
							list_families.push_back(boost::to_upper_copy(family));
						}
						else
						{
							std::stringstream label; label << iline;
							std::string message = "@ClassGrouping list: error in line " + label.str();
							OpenSMOKE::FatalErrorMessage(message);
						}
					}
				}

				fInput.close();

				
				class_grouping_names_families_ = list_families;
				std::sort(class_grouping_names_families_.begin(), class_grouping_names_families_.end());
				class_grouping_names_families_.erase(std::unique(class_grouping_names_families_.begin(), class_grouping_names_families_.end()), class_grouping_names_families_.end());

				for (unsigned int j = 0; j < list_families.size(); j++)
				{
					for (unsigned int i = 0; i < class_grouping_names_families_.size(); i++)
					{
						if (list_families[j] == class_grouping_names_families_[i])
						{
							class_grouping_indices_species_.push_back(thermodynamicsMap_.IndexOfSpecies(list_species[j]) - 1);
							class_grouping_indices_families_.push_back(i);
							break;
						}
					}
				}

				std::cout << std::endl;
				std::cout << "------------------------------------------------------" << std::endl;
				std::cout << "Summary: Class grouping" << std::endl;
				std::cout << "------------------------------------------------------" << std::endl;
				for (unsigned int i = 0; i < class_grouping_names_families_.size(); i++)
				{
					std::cout << "Group: " << class_grouping_names_families_[i] << std::endl;
					for (unsigned int j = 0; j < class_grouping_indices_families_.size(); j++)
						if (class_grouping_indices_families_[j] == i)
							std::cout << " * " << thermodynamicsMap_.NamesOfSpecies()[class_grouping_indices_species_[j]] << std::endl;
				}
				std::cout << "------------------------------------------------------" << std::endl;
				std::cout << std::endl;

				class_grouping_ = true;
			}
			else
			{
				OpenSMOKE::FatalErrorMessage("@ClassGrouping: input file does not exist.");
			}
		}

		if (print_formation_rates_ == true || print_reaction_rates_ == true || class_grouping_ == true)
			is_active_ = true;

		MemoryAllocation();
	}

	void OnTheFlyPostProcessing::MemoryAllocation()
	{
		// Basic variables
		NS_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();

		if (is_active_ == true)
		{
			OpenSMOKE::ChangeDimensions(NS_, &x_, true);
			OpenSMOKE::ChangeDimensions(NS_, &c_, true);
			OpenSMOKE::ChangeDimensions(NS_, &R_, true);
		}
	}

	void OnTheFlyPostProcessing::CloseOutputFiles()
	{
		fFormationRates_.Close();
		fReactionRates_.Close();

		if (class_grouping_ == true)
			fClassGrouping_.Close();
	}

	void OnTheFlyPostProcessing::PrepareOutputFiles(const std::vector<std::string>& additional)
	{
		// Output file
		if (print_formation_rates_ == true)
		{
			const boost::filesystem::path file_name = path_output_ / "FormationRates.out";
			fFormationRates_.Open(file_name);

			const unsigned int precision = 12;
			fFormationRates_.AddColumn("t[s]", precision);
			fFormationRates_.AddColumn("x[m]", precision);
			fFormationRates_.AddColumn("y[m]", precision);
			fFormationRates_.AddColumn("z[m]", precision);
			fFormationRates_.AddColumn("T[K]", precision);
			fFormationRates_.AddColumn("P[Pa]", precision);
			fFormationRates_.AddColumn("QR[W/m3]", precision);

			// Additional variables
			for (unsigned int i = 0; i<additional.size(); i++)
				fFormationRates_.AddColumn(additional[i], precision);

			// Species
			if (indices_of_formation_rates_species_.size() != 0)
			{
				for (unsigned int i = 0; i < indices_of_formation_rates_species_.size(); i++)
					fFormationRates_.AddColumn(thermodynamicsMap_.NamesOfSpecies()[indices_of_formation_rates_species_[i] - 1], precision);
			}
			else
			{
				for (unsigned int i = 0; i < NS_; i++)
					fFormationRates_.AddColumn(thermodynamicsMap_.NamesOfSpecies()[i], precision);
			}

			fFormationRates_.Complete();
		}

		if (print_reaction_rates_ == true)
		{
			const boost::filesystem::path file_name = path_output_ / "ReactionRates.out";
			fReactionRates_.Open(file_name);

			const unsigned int precision = 12;
			fReactionRates_.AddColumn("t[s]", precision);
			fReactionRates_.AddColumn("x[m]", precision);
			fReactionRates_.AddColumn("y[m]", precision);
			fReactionRates_.AddColumn("z[m]", precision);
			fReactionRates_.AddColumn("T[K]", precision);
			fReactionRates_.AddColumn("P[Pa]", precision);
			fReactionRates_.AddColumn("QR[W/m3]", precision);

			// Additional variables
			for (unsigned int i = 0; i<additional.size(); i++)
				fReactionRates_.AddColumn(additional[i], precision);

			// Reactions
			if (indices_of_reaction_rates_.size() != 0)
			{
				for (unsigned int i = 0; i < indices_of_reaction_rates_.size(); i++)
					fReactionRates_.AddColumn("r" + boost::lexical_cast<std::string>(indices_of_reaction_rates_[i]), precision);
			}
			else
			{
				for (unsigned int i = 0; i < NR_; i++)
					fReactionRates_.AddColumn("r" + boost::lexical_cast<std::string>(i + 1), precision);
			}

			fReactionRates_.Complete();
		}

		if (class_grouping_ == true)
		{
			const boost::filesystem::path file_name = path_output_ / "ClassGrouping.out";
			fClassGrouping_.Open(file_name);

			const unsigned int precision = 12;
			fClassGrouping_.AddColumn("t[s]", precision);
			fClassGrouping_.AddColumn("x[m]", precision);
			fClassGrouping_.AddColumn("y[m]", precision);
			fClassGrouping_.AddColumn("z[m]", precision);
			fClassGrouping_.AddColumn("T[K]", precision);
			fClassGrouping_.AddColumn("P[Pa]", precision);
			fClassGrouping_.AddColumn("QR[W/m3]", precision);

			// Additional variables
			for (unsigned int i = 0; i < additional.size(); i++)
				fClassGrouping_.AddColumn(additional[i], precision);

			// Groups
			for (unsigned int i = 0; i < class_grouping_names_families_.size(); i++)
				fClassGrouping_.AddColumn(class_grouping_names_families_[i] + "_w", precision);

			for (unsigned int i = 0; i < class_grouping_names_families_.size(); i++)
				fClassGrouping_.AddColumn(class_grouping_names_families_[i] + "_x", precision);

			fClassGrouping_.Complete();
		}
	}

	void OnTheFlyPostProcessing::PrepareOutputFiles()
	{
		std::vector<std::string> additional(0);
		PrepareOutputFiles(additional);
	}

	void OnTheFlyPostProcessing::WriteOnFile(const double t, const double x, const double y, const double z, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omega, const std::vector<double>& additional)
	{
		double QR = 0;

		// Kinetic data
		{
			thermodynamicsMap_.SetTemperature(T);
			thermodynamicsMap_.SetPressure(P_Pa);

			double MW;
			thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW, omega.GetHandle());
			const double cTot = P_Pa / (PhysicalConstants::R_J_kmol * T);
			Product(cTot, x_, &c_);
			
			kineticsMap_.SetTemperature(T);
			kineticsMap_.SetPressure(P_Pa);
			kineticsMap_.ReactionRates(c_.GetHandle());
			kineticsMap_.FormationRates(R_.GetHandle());
			QR = kineticsMap_.HeatRelease(R_.GetHandle());
		}

		// Print formation rates
		if (print_formation_rates_ == true)
		{
			// Basic variables
			fFormationRates_ << t;
			fFormationRates_ << x;
			fFormationRates_ << y;
			fFormationRates_ << z;
			fFormationRates_ << T;
			fFormationRates_ << P_Pa;

			// Reaction heat [W/m3]
			fFormationRates_ << QR;

			// Additional variables
			for (unsigned int i=0;i<additional.size();i++)
				fFormationRates_ << additional[i];

			// Formation Rates [kmol/m3/s]
			if (formation_rates_moles_ == true)
			{
				if (indices_of_formation_rates_species_.size() != 0)
				{
					for (unsigned int i = 0; i < indices_of_formation_rates_species_.size(); i++)
						fFormationRates_ << R_[indices_of_formation_rates_species_[i]];
				}
				else
				{
					for (unsigned int i = 0; i < NS_; i++)
						fFormationRates_ << R_[i+1];
				}
			}
			else
			{
				if (indices_of_formation_rates_species_.size() != 0)
				{
					for (unsigned int i = 0; i < indices_of_formation_rates_species_.size(); i++)
						fFormationRates_ << thermodynamicsMap_.MW(indices_of_formation_rates_species_[i]-1) * R_[indices_of_formation_rates_species_[i]];
				}
				else
				{
					for (unsigned int i = 0; i < NS_; i++)
						fFormationRates_ << thermodynamicsMap_.MW(i) * R_[i+1];
				}
			}

			fFormationRates_.NewRow();
		}

		// Print reaction rates
		if (print_reaction_rates_ == true)
		{
			// Reaction rates
			OpenSMOKE::OpenSMOKEVectorDouble r(NR_);
			kineticsMap_.GiveMeReactionRates(r.GetHandle());

			// Basic variables
			fReactionRates_ << t;
			fReactionRates_ << x;
			fReactionRates_ << y;
			fReactionRates_ << z;
			fReactionRates_ << T;
			fReactionRates_ << P_Pa;

			// Reaction heat [W/m3]
			fReactionRates_ << QR;

			// Additional variables
			for (unsigned int i = 0; i<additional.size(); i++)
				fReactionRates_ << additional[i];

			if (indices_of_reaction_rates_.size() != 0)
			{
				for (unsigned int i = 0; i < indices_of_reaction_rates_.size(); i++)
					fReactionRates_ << r[indices_of_reaction_rates_[i]];
			}
			else
			{
				for (unsigned int i = 1; i <= NR_; i++)
					fReactionRates_ << r[i];
			}

			fReactionRates_.NewRow();
		}

		if (class_grouping_ == true)
		{
			std::vector<double> sum_x(class_grouping_names_families_.size());
			std::vector<double> sum_omega(class_grouping_names_families_.size());
			std::fill(sum_x.begin(), sum_x.end(), 0.);
			std::fill(sum_omega.begin(), sum_omega.end(), 0.);
			for (unsigned int j = 0; j < class_grouping_names_families_.size(); j++)
			{
				sum_x[class_grouping_indices_families_[j]] += x_[class_grouping_indices_species_[j] + 1];
				sum_omega[class_grouping_indices_families_[j]] += omega[class_grouping_indices_species_[j] + 1];
			}

			// Basic variables
			fClassGrouping_ << t;
			fClassGrouping_ << x;
			fClassGrouping_ << y;
			fClassGrouping_ << z;
			fClassGrouping_ << T;
			fClassGrouping_ << P_Pa;

			// Reaction heat [W/m3]
			fClassGrouping_ << QR;

			// Additional variables
			for (unsigned int i = 0; i < additional.size(); i++)
				fClassGrouping_ << additional[i];

			// Groups (mass fractions)
			for (unsigned int i = 0; i < class_grouping_names_families_.size(); i++)
				fClassGrouping_ << sum_omega[i];

			// Groups (mole fractions)
			for (unsigned int i = 0; i < class_grouping_names_families_.size(); i++)
				fClassGrouping_ << sum_x[i];

			fClassGrouping_.NewRow();
		}
	}

	void OnTheFlyPostProcessing::WriteOnFile(const double t, const double x, const double y, const double z, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omega)
	{
		std::vector<double> additional(0);
		WriteOnFile(t, x, y, z, T, P_Pa, omega, additional);
	}
}

