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
|   Copyright(C) 2020  Alberto Cuoci                                      |
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

#include <boost/algorithm/string.hpp>

namespace OpenSMOKE
{
    SpeciesClasses::SpeciesClasses()
    {
        is_active_ = false;
    }
    	
    SpeciesClasses::SpeciesClasses(const std::vector<std::string>& species_class_names, const std::vector<int>& species_class_lines_abs)
	{
		is_active_ = false;
		if (species_class_names.size() != 0)
		{
			is_active_ = true;
			species_class_names_ = species_class_names;
			species_class_lines_abs_ = species_class_lines_abs;
			species_indices_.resize(species_class_names_.size());
			n_classes_ = static_cast<unsigned int>(species_class_names_.size());
		}
	}

    SpeciesClasses::SpeciesClasses(const boost::filesystem::path& file_name)
	{
		is_active_ = false;
		ReadXMLFile(file_name);
	}

    bool SpeciesClasses::Checking()
	{
		return true;
	}

    void SpeciesClasses::Summary(std::ostream& out)
	{
		if (is_active_ == true)
		{
			out << std::endl;
			out << "-----------------------------------------------------------------------------------------" << std::endl;
			out << "                                 Species Classes                                   " << std::endl;
			out << "-----------------------------------------------------------------------------------------" << std::endl;
			out << std::left << std::setw(15) << "Index"
			    << std::left << std::setw(30) << "Name"
				<< std::left << std::setw(15) << "#species"
				<< std::left << std::setw(15) << "min-species"
				<< std::left << std::setw(15) << "max-species" << std::endl;
			out << "-----------------------------------------------------------------------------------------" << std::endl;

			unsigned int total_number_species = 0;
			for (unsigned int k = 0; k < n_classes_; k++)
			{
				// out << " ";
				out << std::left << std::setw(15) << k; // indice della classe
				out << std::left << std::setw(30) << species_class_names_[k]; // nome della classe definito nel CHEMKIN
				out << std::left << std::setw(15) << species_indices_[k].back() - species_indices_[k].front() + 1;
				out << std::left << std::setw(15) << *std::min_element(species_indices_[k].begin(), species_indices_[k].end()) + 1;
				out << std::left << std::setw(15) << *std::max_element(species_indices_[k].begin(), species_indices_[k].end()) + 1;
				out << std::endl;

				total_number_species += static_cast<unsigned int>(species_indices_[k].back() - species_indices_[k].front() + 1);
			}
			out << "-----------------------------------------------------------------------------------------" << std::endl;
			out << std::left << std::setw(45) << " " << total_number_species << std::endl;
			out << "-----------------------------------------------------------------------------------------" << std::endl;
		}
    }

	unsigned int SpeciesClasses::index_from_species_name(const std::string name)
	{
		if (is_active_ == true)
		{
			for (unsigned int k = 0; k < n_classes_; k++)
				if (species_class_names_[k] == name)
					return k;

			OpenSMOKE::FatalErrorMessage("The specified species class is not available: " + name);
		}
		else
		{
			OpenSMOKE::FatalErrorMessage("No classes have been defined for species");
		}

		return 0;
	}
	
	bool SpeciesClasses::filter_species(const unsigned int j, const unsigned int i)
	{

		if (i <= species_class_lines_abs_[0])
			return false;

		for (unsigned int k=n_classes_-1;k>=0;k--)
			if (i >= species_class_lines_abs_[k])
			{
				species_indices_[k].push_back(j);
				return true;
			}

		OpenSMOKE::FatalErrorMessage("SpeciesClasses: Error in filtering the species indices");
		return false;
	}

	bool SpeciesClasses::filter_species_nobili(const std::vector<std::string>& names_species,
		const std::vector<std::string>& firstspecies_names_)
	{
		for (unsigned int i = 0; i < firstspecies_names_.size(); i++) {
			std::vector<int> firstspecies_indices;
			auto it = std::find (names_species.begin(), 
				names_species.end(), 
				firstspecies_names_[i]); //same size than n_classes

			int pos = std::distance(names_species.begin(), it);

			firstspecies_indices.push_back(pos);
			species_indices_[i].push_back(pos); 
			if (i>0) species_indices_[i - 1].push_back(pos - 1); //index of the last species in each speciesclass
		}
		species_indices_[n_classes_ - 1].push_back(names_species.size() - 1); //index of the last species in the last speciesclass
		return true;

		//Come gestire il check?
		//OpenSMOKE::FatalErrorMessage("SpeciesClasses: Error in filtering the species indices");
		//return false;
	}

	void SpeciesClasses::WriteXMLFile(std::stringstream& fOutput) const
	{
		std::cout << " * Writing the species classes in XML format..." << std::endl;

		fOutput << "<SpeciesClasses>" << std::endl;
		{
			fOutput << "<NumberOfClasses>" << std::endl;
			fOutput << n_classes_ << std::endl;
			fOutput << "</NumberOfClasses>" << std::endl;

			fOutput << "<ClassNames>" << std::endl;
			for (unsigned int k = 0; k < n_classes_; k++)
				fOutput << species_class_names_[k] << std::endl;
			fOutput << "</ClassNames>" << std::endl;

			fOutput << "<ClassSizes>" << std::endl;
			for (unsigned int k = 0; k < n_classes_; k++)
				fOutput << species_indices_[k].back() - species_indices_[k].front() + 1  << std::endl;
			fOutput << "</ClassSizes>" << std::endl;


			for (unsigned int k = 0; k < n_classes_; k++)
			{
				fOutput << "<ClassSpecies name=\"" << species_class_names_[k] << "\" n=\"" << species_indices_[k].back() - species_indices_[k].front() + 1 << "\">" << std::endl;
				for (unsigned int j = species_indices_[k].front(); j < species_indices_[k].back() + 1; j++)
					fOutput << j << std::endl;
				fOutput << "</ClassSpecies>" << std::endl;
			}
		}
		fOutput << "</SpeciesClasses>" << std::endl;
	}

	void SpeciesClasses::ReadXMLFile(const boost::filesystem::path& file_name)
	{
		// ATTENZIONE DA RIVEDERE BENE
		std::cout << " * Reading the species classes in XML format..." << std::endl;

		boost::property_tree::ptree ptree;
		boost::property_tree::read_xml((file_name).string(), ptree);

		auto child = ptree.get_child_optional("opensmoke.Kinetics.SpeciesClasses");

		if (!child)
		{
			std::cout << "No species classes have been found in the kinetics.xml file" << std::endl;
		}
		else
		{
			// Is active
			is_active_ = true;

			// Number of classes
			n_classes_ = ptree.get<unsigned int>("opensmoke.Kinetics.SpeciesClasses.NumberOfClasses");
			species_indices_.resize(n_classes_);
			species_class_names_.resize(n_classes_);

			// Class sizes
			{
				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.Kinetics.SpeciesClasses.ClassSizes"));
				for (unsigned int k = 0; k < n_classes_; k++)
				{
					unsigned int dummy;
					stream >> dummy;
					species_indices_[k].resize(dummy);
				}
			}

			// Class names
			{
				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.Kinetics.SpeciesClasses.ClassNames"));
				for (unsigned int k = 0; k < n_classes_; k++)
				{
					stream >> species_class_names_[k];
				}
			}

			{
				// Species indices
				unsigned int k = 0;
				BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree.get_child("opensmoke.Kinetics.SpeciesClasses"))
				{
					boost::property_tree::ptree subtree = node.second;

					if (node.first == "ClassSpecies")
					{
						std::string dummy = subtree.get_value< std::string >();
						std::stringstream stream(dummy);
						for (unsigned int j = 0; j < species_indices_[k].back() - species_indices_[k].front() + 1; j++)
						{
							stream >> species_indices_[k][j];
						}

						k++;
					}
				}
			}

			// Summary on the screen
			Summary(std::cout);
		}
	}
} // namespace OpenSMOKE
