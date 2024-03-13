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
	PolimiSootClasses::PolimiSootClasses()
	{
		is_active_ = false;
	}

	PolimiSootClasses::PolimiSootClasses(const std::vector<std::string>& soot_class_names, const std::vector<int>& soot_class_lines_abs)
	{
		is_active_ = false;
		if (soot_class_names.size() != 0)
		{
			is_active_ = true;
			soot_class_names_ = soot_class_names;
			soot_class_lines_abs_ = soot_class_lines_abs;
			reaction_indices_.resize(soot_class_names_.size());
			n_classes_ = static_cast<unsigned int>(soot_class_names_.size());
		}
	}

	PolimiSootClasses::PolimiSootClasses(const boost::filesystem::path& file_name)
	{
		is_active_ = false;
		ReadXMLFile(file_name);
	}

	bool PolimiSootClasses::Checking()
	{
		return true;
	}

	void PolimiSootClasses::Summary(std::ostream& out)
	{
		if (is_active_ == true)
		{
			out << std::endl;
			out << "---------------------------------------------------------------------------------" << std::endl;
			out << "                              POLIMI Soot Classes                                " << std::endl;
			out << "---------------------------------------------------------------------------------" << std::endl;
			out << " Index  Name                 #reactions    min-reaction  max-reaction            " << std::endl;
			out << "---------------------------------------------------------------------------------" << std::endl;

			unsigned int total_number_reactions = 0;
			for (unsigned int k = 0; k < n_classes_; k++)
			{
				out << " ";
				out << std::left << std::setw(7) << k;
				out << std::left << std::setw(21) << soot_class_names_[k]; 
				out << std::left << std::setw(14) << reaction_indices_[k].size();
				out << std::left << std::setw(14) << *std::min_element(reaction_indices_[k].begin(), reaction_indices_[k].end()) + 1;
				out << std::left << std::setw(14) << *std::max_element(reaction_indices_[k].begin(), reaction_indices_[k].end()) + 1;
				out << std::endl;

				total_number_reactions += static_cast<unsigned int>(reaction_indices_[k].size());
			}
			out << "---------------------------------------------------------------------------------" << std::endl;
			out << std::left << std::setw(29) << " " << total_number_reactions << std::endl;
			out << "---------------------------------------------------------------------------------" << std::endl;
		}
	}

	unsigned int PolimiSootClasses::index_from_class_name(const std::string name)
	{
		if (is_active_ == true)
		{
			for (unsigned int k = 0; k < n_classes_; k++)
				if (soot_class_names_[k] == name)
					return k;

			OpenSMOKE::FatalErrorMessage("The specified soot class is not available: " + name);
		}
		else
		{
			OpenSMOKE::FatalErrorMessage("No classes have been defined for soot");
		}

		return 0;
	}
	
	bool PolimiSootClasses::filter_reaction(const unsigned int j, const unsigned int i)
	{
		if (i <= soot_class_lines_abs_[0])
			return false;

		for (int k=n_classes_-1;k>=0;k--)
			if (i >= soot_class_lines_abs_[k])
			{
				reaction_indices_[k].push_back(j);
				return true;
			}

		OpenSMOKE::FatalErrorMessage("PolimiSootClasses: Error in filtering the reaction indices");
		return false;
	}

	void PolimiSootClasses::WriteXMLFile(std::stringstream& fOutput) const
	{
		std::cout << " * Writing the soot classes in XML format..." << std::endl;

		fOutput << "<PolimiSootClasses>" << std::endl;
		{
			fOutput << "<NumberOfClasses>" << std::endl;
			fOutput << n_classes_ << std::endl;
			fOutput << "</NumberOfClasses>" << std::endl;

			fOutput << "<ClassNames>" << std::endl;
			for (unsigned int k = 0; k < n_classes_; k++)
				fOutput << soot_class_names_[k] << std::endl;
			fOutput << "</ClassNames>" << std::endl;

			fOutput << "<ClassSizes>" << std::endl;
			for (unsigned int k = 0; k < n_classes_; k++)
				fOutput << reaction_indices_[k].size() << std::endl;
			fOutput << "</ClassSizes>" << std::endl;


			for (unsigned int k = 0; k < n_classes_; k++)
			{
				fOutput << "<ClassReactions name=\"" << soot_class_names_[k] << "\" n=\"" << reaction_indices_[k].size() << "\">" << std::endl;
				for (unsigned int j = 0; j < reaction_indices_[k].size(); j++)
					fOutput << reaction_indices_[k][j] << std::endl;
				fOutput << "</ClassReactions>" << std::endl;
			}
		}
		fOutput << "</PolimiSootClasses>" << std::endl;
	}

	void PolimiSootClasses::ReadXMLFile(const boost::filesystem::path& file_name)
	{
		std::cout << " * Reading the soot classes in XML format..." << std::endl;

		boost::property_tree::ptree ptree;
		boost::property_tree::read_xml((file_name).string(), ptree);

		auto child = ptree.get_child_optional("opensmoke.Kinetics.PolimiSootClasses");

		if (!child)
		{
			std::cout << "No soot classes have been found in the kinetics.xml file" << std::endl;
		}
		else
		{
			// Is active
			is_active_ = true;

			// Number of classes
			n_classes_ = ptree.get<unsigned int>("opensmoke.Kinetics.PolimiSootClasses.NumberOfClasses");
			reaction_indices_.resize(n_classes_);
			soot_class_names_.resize(n_classes_);

			// Class sizes
			{
				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.Kinetics.PolimiSootClasses.ClassSizes"));
				for (unsigned int k = 0; k < n_classes_; k++)
				{
					unsigned int dummy;
					stream >> dummy;
					reaction_indices_[k].resize(dummy);
				}
			}

			// Class names
			{
				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.Kinetics.PolimiSootClasses.ClassNames"));
				for (unsigned int k = 0; k < n_classes_; k++)
				{
					stream >> soot_class_names_[k];
				}
			}

			{
				// Reaction indices
				unsigned int k = 0;
				BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree.get_child("opensmoke.Kinetics.PolimiSootClasses"))
				{
					boost::property_tree::ptree subtree = node.second;

					if (node.first == "ClassReactions")
					{
						std::string dummy = subtree.get_value< std::string >();
						std::stringstream stream(dummy);
						for (unsigned int j = 0; j < reaction_indices_[k].size(); j++)
						{
							stream >> reaction_indices_[k][j];
						}

						k++;
					}
				}
			}

			// Summary on the screen
			Summary(std::cout);
		}
	}
}

