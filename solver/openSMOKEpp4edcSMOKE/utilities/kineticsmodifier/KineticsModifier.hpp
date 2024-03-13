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
|	License                                                               |
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

#include "Grammar_KineticsModifier.h"

namespace OpenSMOKE
{
	KineticsModifier::KineticsModifier()
	{
		is_active_ = false;
		is_verbose_ = false;
	}

	void KineticsModifier::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		is_active_ = false;

		Grammar_KineticsModifier grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@A") == true)
		{
			std::vector<std::string> list_values;
			dictionary.ReadOption("@A", list_values);

			if (list_values.size() % 2 != 0)
				OpenSMOKE::FatalErrorMessage("@A option has a wrong number of elements.");

			const int n = static_cast<int>(list_values.size() / 2);
			for (int i = 0; i < n; i++)
			{
				list_A_index_.push_back(boost::lexical_cast<unsigned int>(list_values[i * 2]));
				list_A_value_.push_back(boost::lexical_cast<double>(list_values[i * 2 + 1]));
			}

			is_active_ = true;
		}

		if (dictionary.CheckOption("@n") == true)
		{
			std::vector<std::string> list_values;
			dictionary.ReadOption("@n", list_values);

			if (list_values.size() % 2 != 0)
				OpenSMOKE::FatalErrorMessage("@n option has a wrong number of elements.");

			const int n = static_cast<int>(list_values.size() / 2);
			for (int i = 0; i < n; i++)
			{
				list_n_index_.push_back(boost::lexical_cast<unsigned int>(list_values[i * 2]));
				list_n_value_.push_back(boost::lexical_cast<double>(list_values[i * 2 + 1]));
			}

			is_active_ = true;
		}

		if (dictionary.CheckOption("@EoverR") == true)
		{
			std::vector<std::string> list_values;
			dictionary.ReadOption("@EoverR", list_values);

			if (list_values.size() % 2 != 0)
				OpenSMOKE::FatalErrorMessage("@EoverR option has a wrong number of elements.");

			const int n = static_cast<int>(list_values.size() / 2);
			for (int i = 0; i < n; i++)
			{
				list_EoverR_index_.push_back(boost::lexical_cast<unsigned int>(list_values[i * 2]));
				list_EoverR_value_.push_back(boost::lexical_cast<double>(list_values[i * 2 + 1]));
			}

			is_active_ = true;
		}

		if (dictionary.CheckOption("@ThirdBody") == true)
		{
			std::vector<std::string> list_values;
			dictionary.ReadOption("@ThirdBody", list_values);

			if (list_values.size() % 3 != 0)
				OpenSMOKE::FatalErrorMessage("@ThirdBody option has a wrong number of elements.");

			const int n = static_cast<int>(list_values.size() / 3);
			for (int i = 0; i < n; i++)
			{
				list_ThirdBody_index_.push_back(boost::lexical_cast<unsigned int>(list_values[i * 3]));
				list_ThirdBody_name_.push_back(list_values[i * 3 + 1]);
				list_ThirdBody_value_.push_back(boost::lexical_cast<double>(list_values[i * 3 + 2]));
			}

			is_active_ = true;
		}

		if (dictionary.CheckOption("@FrozenSpecies") == true)
		{
			dictionary.ReadOption("@FrozenSpecies", list_frozen_species_name_);
			is_active_ = true;
		}

	}

	void KineticsModifier::SetupFromIndices(const std::vector<unsigned int>& indices, const std::vector<double>& coefficients)
	{
		list_A_multiplier_index_ = indices;
		std::for_each(std::begin(list_A_multiplier_index_), std::end(list_A_multiplier_index_), [](unsigned int& x) { x += 1; });
		list_A_multiplier_value_ = coefficients;
		is_active_ = true;
	}

	void KineticsModifier::Setup(ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap)
	{
		is_verbose_ = true;
		
		if (is_active_ == true && is_verbose_ == true)
		{
			if (list_A_index_.size() != 0)
			{
				std::cout << std::endl;
				std::cout << "---------------------------------------------------------------------------" << std::endl;
				std::cout << "Frequency factors [kmol, m, s]" << std::endl;
				std::cout << "---------------------------------------------------------------------------" << std::endl;
				for (unsigned int j = 0; j < list_A_index_.size(); j++)
				{
					std::cout << "  Reaction: ";
					std::cout << std::setw(6) << std::left << list_A_index_[j];
					std::cout << std::setw(12) << std::right << kineticsMap.A(list_A_index_[j] - 1);
					std::cout << " -> ";
					std::cout << std::setw(12) << std::right << list_A_value_[j];
					std::cout << std::endl;

					kineticsMap.Set_A(list_A_index_[j] - 1, list_A_value_[j]);
				}
			}

			if (list_n_index_.size() != 0)
			{
				std::cout << std::endl;
				std::cout << "---------------------------------------------------------------------------" << std::endl;
				std::cout << "Temperature exponents [-]" << std::endl;
				std::cout << "---------------------------------------------------------------------------" << std::endl;
				for (unsigned int j = 0; j < list_n_index_.size(); j++)
				{
					std::cout << "  Reaction: ";
					std::cout << std::setw(6) << std::left << list_n_index_[j];
					std::cout << std::setw(12) << std::right << kineticsMap.Beta(list_n_index_[j] - 1);
					std::cout << " -> ";
					std::cout << std::setw(12) << std::right << list_n_value_[j];
					std::cout << std::endl;

					kineticsMap.Set_Beta(list_n_index_[j] - 1, list_n_value_[j]);
				}
			}

			if (list_EoverR_index_.size() != 0)
			{
				std::cout << std::endl;
				std::cout << "---------------------------------------------------------------------------" << std::endl;
				std::cout << "Activation temperatures [K]" << std::endl;
				std::cout << "---------------------------------------------------------------------------" << std::endl;
				for (unsigned int j = 0; j < list_EoverR_index_.size(); j++)
				{
					std::cout << "  Reaction: ";
					std::cout << std::setw(6) << std::left << list_EoverR_index_[j];
					std::cout << std::setw(12) << std::right << kineticsMap.E_over_R(list_EoverR_index_[j] - 1);
					std::cout << " -> ";
					std::cout << std::setw(12) << std::right << list_EoverR_value_[j];
					std::cout << std::endl;
					kineticsMap.Set_E_over_R(list_EoverR_index_[j] - 1, list_EoverR_value_[j]);
				}
			}

			if (list_ThirdBody_index_.size() != 0)
			{
				std::cout << std::endl;
				std::cout << "---------------------------------------------------------------------------" << std::endl;
				std::cout << "Third-Body efficiencies" << std::endl;
				std::cout << "---------------------------------------------------------------------------" << std::endl;
				for (unsigned int j = 0; j < list_ThirdBody_index_.size(); j++)
				{
					unsigned int species_index = thermodynamicsMap.IndexOfSpecies(list_ThirdBody_name_[j]);
					
					std::cout << "  Reaction: ";
					std::cout << std::setw(6) << std::left << list_ThirdBody_index_[j];
					std::cout << std::setw(12) << std::left << list_ThirdBody_name_[j];
					std::cout << std::setw(6) << std::right << kineticsMap.ThirdBody(list_ThirdBody_index_[j]-1, species_index-1);
					std::cout << " -> ";
					std::cout << std::setw(6) << std::right << list_ThirdBody_value_[j];
					std::cout << std::endl;
					
					kineticsMap.Set_ThirdBody(list_ThirdBody_index_[j] - 1, species_index-1, list_ThirdBody_value_[j]);
				}
			}

			if (list_A_multiplier_index_.size() != 0)
			{
				std::cout << std::endl;
				std::cout << "---------------------------------------------------------------------------" << std::endl;
				std::cout << "Frequency factors [kmol, m, s]" << std::endl;
				std::cout << "---------------------------------------------------------------------------" << std::endl;
				for (unsigned int j = 0; j < list_A_multiplier_index_.size(); j++)
				{
					std::cout << "  Reaction: ";
					std::cout << std::setw(6) << std::left << list_A_multiplier_index_[j];
					std::cout << std::setw(12) << std::right << kineticsMap.A(list_A_multiplier_index_[j]-1);
					std::cout << " -> ";
					std::cout << std::setw(12) << std::right << kineticsMap.A(list_A_multiplier_index_[j]-1)*list_A_multiplier_value_[j];
					std::cout << std::endl;

					kineticsMap.Set_A(list_A_multiplier_index_[j]-1, kineticsMap.A(list_A_multiplier_index_[j]-1)*list_A_multiplier_value_[j] );
				}
			}

			if (list_frozen_species_name_.size() != 0)
			{
				if (list_frozen_species_name_[0] == "ALL" || list_frozen_species_name_[0] == "all")
				{
					list_frozen_species_index_.resize(thermodynamicsMap.NumberOfSpecies());
					for (unsigned int i = 0; i < thermodynamicsMap.NumberOfSpecies(); i++)
						list_frozen_species_index_[i] = i;
				}
				else if (list_frozen_species_name_[0] == "NONE" || list_frozen_species_name_[0] == "none")
				{
					list_frozen_species_name_.resize(0);
				}
				else
				{
					list_frozen_species_index_.resize(list_frozen_species_name_.size());
					for (unsigned int i = 0; i < list_frozen_species_name_.size(); i++)
						list_frozen_species_index_[i] = thermodynamicsMap.IndexOfSpecies(list_frozen_species_name_[i])-1;
				}

				kineticsMap.stoichiometry().SetFrozenSpecies(list_frozen_species_index_);
			}
		}

		if (is_active_ == true && is_verbose_ == false)
		{

		}
	}
}

