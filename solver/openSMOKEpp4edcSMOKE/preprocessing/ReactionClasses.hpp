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
|   Copyright(C) 2022  Alberto Cuoci                                      |
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
	ReactionClasses::ReactionClasses()
	{
		is_active_ = false;
	}

	ReactionClasses::ReactionClasses(	const std::vector<std::string>& reaction_mainclass_names, const std::vector<int>& reaction_mainclass_lines_abs,
										const std::vector<std::string>& endreaction_mainclass_names, const std::vector<int>& endreaction_mainclass_lines_abs,
										const std::vector<std::string>& reaction_subclass_names, const std::vector<int>& reaction_subclass_lines_abs,
										const std::vector<std::string>& endreaction_subclass_names, const std::vector<int>& endreaction_subclass_lines_abs )
	{
		is_active_ = false;
		if (reaction_mainclass_names.size() != 0)
		{
			is_active_ = true;
			reaction_mainclass_names_ = reaction_mainclass_names;
			reaction_mainclass_lines_abs_ = reaction_mainclass_lines_abs;
			endreaction_mainclass_names_ = endreaction_mainclass_names;
			endreaction_mainclass_lines_abs_ = endreaction_mainclass_lines_abs;

			reaction_subclass_names_ = reaction_subclass_names;
			reaction_subclass_lines_abs_ = reaction_subclass_lines_abs;
			endreaction_subclass_names_ = endreaction_subclass_names;
			endreaction_subclass_lines_abs_ = endreaction_subclass_lines_abs;

			// Checking main classes
			for (unsigned int j=0;j< reaction_mainclass_names_.size();j++)
				if (endreaction_mainclass_names_.size() >= j + 1)
				{
					if (endreaction_mainclass_names_[j] != reaction_mainclass_names_[j])
						OpenSMOKE::FatalErrorMessage("Missing [ENDREACTIONCLASS] for the following [REACTIONCLASS]: " + reaction_mainclass_names_[j]);
				}
				else
				{
					OpenSMOKE::FatalErrorMessage("Missing [ENDREACTIONCLASS] for the following [REACTIONCLASS]: " + reaction_mainclass_names_[j]);
				}

			// Checking main classes
			for (unsigned int j = 0; j < endreaction_mainclass_names_.size(); j++)
				if (reaction_mainclass_names_.size() >= j + 1)
				{
					if (reaction_mainclass_names_[j] != endreaction_mainclass_names_[j])
						OpenSMOKE::FatalErrorMessage("Missing [REACTIONCLASS] for the following [ENDREACTIONCLASS]: " + endreaction_mainclass_names_[j]);
				}
				else
				{
					OpenSMOKE::FatalErrorMessage("Missing [REACTIONCLASS] for the following [ENDREACTIONCLASS]: " + endreaction_mainclass_names_[j]);
				}

			// Checking sub classes
			for (unsigned int j = 0; j < reaction_subclass_names_.size(); j++)
				if (endreaction_subclass_names_.size() >= j + 1)
				{
					if (endreaction_subclass_names_[j] != reaction_subclass_names_[j])
						OpenSMOKE::FatalErrorMessage("Missing [ENDREACTIONCLASS] for the following [REACTIONCLASS]: " + reaction_subclass_names_[j]);
				}
				else
				{
					OpenSMOKE::FatalErrorMessage("Missing [ENDREACTIONCLASS] for the following [REACTIONCLASS]: " + reaction_subclass_names_[j]);
				}

			// Checking sub classes
			for (unsigned int j = 0; j < endreaction_subclass_names_.size(); j++)
				if (reaction_subclass_names_.size() >= j + 1)
				{
					if (reaction_subclass_names_[j] != endreaction_subclass_names_[j])
						OpenSMOKE::FatalErrorMessage("Missing [REACTIONCLASS] for the following [ENDREACTIONCLASS]: " + endreaction_subclass_names_[j]);
				}
				else
				{
					OpenSMOKE::FatalErrorMessage("Missing [REACTIONCLASS] for the following [ENDREACTIONCLASS]: " + endreaction_subclass_names_[j]);
				}

			// Analysis of main classes: selection of unique values
			{
				std::set<std::string> s{};
				std::vector<int> uniqueidx{};

				for (size_t i = 0; i < reaction_mainclass_names_.size(); i++)
					if (s.insert(reaction_mainclass_names_[i]).second)
						uniqueidx.push_back(i);

				for (const auto i : uniqueidx)
					mainclass_names_.push_back(reaction_mainclass_names_[i]);
			}

			n_mainclasses_ = static_cast<unsigned int>(mainclass_names_.size());
			subclass_names_.resize(n_mainclasses_);
			n_subclasses_.resize(n_mainclasses_);
			reaction_indices_.resize(n_mainclasses_);
			subclass_names_lines_abs_.resize(n_mainclasses_);
			endsubclass_names_lines_abs_.resize(n_mainclasses_);

			// Identification of subclasses
			for (unsigned int k = 0; k < n_mainclasses_; k++)
			{
				for (unsigned int j = 0; j < reaction_subclass_names_.size(); j++)
					if (reaction_mainclass_names_[j] == mainclass_names_[k])
					{
						subclass_names_[k].push_back(reaction_subclass_names_[j]);
						subclass_names_lines_abs_[k].push_back(reaction_subclass_lines_abs_[j]);
						endsubclass_names_lines_abs_[k].push_back(endreaction_subclass_lines_abs_[j]);
					}

				n_subclasses_[k] = subclass_names_[k].size();
				reaction_indices_[k].resize(n_subclasses_[k]);
			}
			

		//	Summary(std::cout);
		}
	}

	ReactionClasses::ReactionClasses(const boost::filesystem::path& file_name)
	{
		is_active_ = false;
		ReadXMLFile(file_name);
	}

	bool ReactionClasses::Checking()
	{
		return true;
	}

	void ReactionClasses::Summary(std::ostream& out)
	{
		if (is_active_ == true)
		{
			out << std::endl;
			out << "-------------------------------------------------------------------------------------------------------------" << std::endl;
			out << "                                              Reaction Classes                                               " << std::endl;
			out << "-------------------------------------------------------------------------------------------------------------" << std::endl;
			out << " Index  Index  Main                 Sub                            #reactions    min-reaction  max-reaction  " << std::endl;
			out << "-------------------------------------------------------------------------------------------------------------" << std::endl;

			unsigned int total_number_reactions = 0;
			for (unsigned int k = 0; k < n_mainclasses_; k++)
			{
				for (unsigned int j = 0; j < n_subclasses_[k]; j++)
				{
					out << " ";
					out << std::left << std::setw(7) << k;
					out << std::left << std::setw(7) << j;
					out << std::left << std::setw(21) << mainclass_names_[k];
					out << std::left << std::setw(31) << subclass_names_[k][j];

					out << std::left << std::setw(14) << reaction_indices_[k][j].size();
					if (reaction_indices_[k][j].size() != 0)
					{
						out << std::left << std::setw(14) << *std::min_element(reaction_indices_[k][j].begin(), reaction_indices_[k][j].end()) + 1;
						out << std::left << std::setw(14) << *std::max_element(reaction_indices_[k][j].begin(), reaction_indices_[k][j].end()) + 1;
					}	
					else
					{
						out << std::left << std::setw(14) << 0;
						out << std::left << std::setw(14) << 0;
					}								
					out << std::endl;

					total_number_reactions += static_cast<unsigned int>(reaction_indices_[k][j].size());
				}
			}
			out << "-------------------------------------------------------------------------------------------------------------" << std::endl;
			out << std::left << std::setw(81) << " " << total_number_reactions << std::endl;
			out << "-------------------------------------------------------------------------------------------------------------" << std::endl;
		}
	}

	unsigned int ReactionClasses::index_from_mainclass_name(const std::string name)
	{
		if (is_active_ == true)
		{
			for (unsigned int k = 0; k < n_mainclasses_; k++)
				if (mainclass_names_[k] == name)
					return k;

			OpenSMOKE::FatalErrorMessage("The specified reaction mainclass is not available: " + name);
		}
		else
		{
			OpenSMOKE::FatalErrorMessage("No mainclasses have been defined for reactions");
		}

		return 0;
	}
	
	bool ReactionClasses::filter_reaction(const unsigned int j, const unsigned int i)
	{
		for (int k=n_mainclasses_-1;k>=0;k--)
			for (int kk = n_subclasses_[k] - 1; kk >= 0; kk--)
				if (i >= subclass_names_lines_abs_[k][kk] && i < endsubclass_names_lines_abs_[k][kk])
					reaction_indices_[k][kk].push_back(j);
				
		return true;
	}

	void ReactionClasses::WriteXMLFile(std::stringstream& fOutput, const std::vector<std::string>& names_species, const Eigen::MatrixXd& nu) const
	{
		std::cout << " * Writing the reaction classes in XML format..." << std::endl;

		fOutput << "<ReactionClasses>" << std::endl;
		{
			fOutput << "<NumberOfMainClasses>" << std::endl;
			fOutput << n_mainclasses_ << std::endl;
			fOutput << "</NumberOfMainClasses>" << std::endl;

			fOutput << "<MainClassNames>" << std::endl;
			for (unsigned int k = 0; k < n_mainclasses_; k++)
				fOutput << mainclass_names_[k] << std::endl;
			fOutput << "</MainClassNames>" << std::endl;

			fOutput << "<MainClassSizes>" << std::endl;
			for (unsigned int k = 0; k < n_mainclasses_; k++)
				fOutput << n_subclasses_[k] << std::endl;
			fOutput << "</MainClassSizes>" << std::endl;

			// Sub classes
			fOutput << "<!--ReactionIndices starts from 0-->" << std::endl;
			fOutput << "<!--Species reports the stoichiometric coefficients for each reaction in ReactionIndices-->" << std::endl;
			for (unsigned int k = 0; k < n_mainclasses_; k++)
			{
				fOutput << "<MainClass name=\"" << mainclass_names_[k] <<"\" index=\"" << k << "\" size=\"" << n_subclasses_[k] << "\">" << std::endl;

				for (unsigned int j = 0; j < n_subclasses_[k]; j++)
				{
					fOutput << "<SubClass name=\"" << subclass_names_[k][j] << "\" index=\"" << j << "\" size=\"" << reaction_indices_[k][j].size() << "\">" << std::endl;
					
					fOutput << "<ReactionIndices>" << std::endl;
					for (unsigned int i=0;i< reaction_indices_[k][j].size();i++)
						fOutput << reaction_indices_[k][j][i] << " ";
					fOutput << std::endl;
					fOutput << "</ReactionIndices>" << std::endl;

					const unsigned int ns = nu.rows();
					const unsigned int nr = nu.cols();
					for (unsigned int kk = 0; kk < ns; kk++)
					{
						// Check if the species is involved in reaction classes
						bool is_present = false;
						for (unsigned int i = 0; i < reaction_indices_[k][j].size(); i++)
							if (nu(kk, reaction_indices_[k][j][i]) != 0.)
							{
								is_present = true;
								break;
							}

						// Write on XML file
						if (is_present == true)
						{
							fOutput << "<Species index=\"" << kk << "\" name = \"" << names_species[kk] << "\">" << std::endl;
							for (unsigned int i = 0; i < reaction_indices_[k][j].size(); i++)
								fOutput << nu(kk, reaction_indices_[k][j][i]) << " ";
							fOutput << std::endl;
							fOutput << "</Species>" << std::endl;
						}
					}

					fOutput << "</SubClass>" << std::endl;
				}

				fOutput << "</MainClass>" << std::endl;
			}
		}
		fOutput << "</ReactionClasses>" << std::endl;
	}

	void ReactionClasses::ReadXMLFile(const boost::filesystem::path& file_name)
	{
		/*
		std::cout << " * Reading the reaction classes in XML format..." << std::endl;

		boost::property_tree::ptree ptree;
		boost::property_tree::read_xml((file_name).string(), ptree);

		auto child = ptree.get_child_optional("opensmoke.Kinetics.ReactionClasses");

		if (!child)
		{
			std::cout << "No reaction classes have been found in the kinetics.xml file" << std::endl;
		}
		else
		{
			// Is active
			is_active_ = true;

			// Number of classes
			n_mainclasses_ = ptree.get<unsigned int>("opensmoke.Kinetics.ReactionClasses.NumberOfMainClasses");
			reaction_indices_.resize(n_mainclasses_);
			mainclass_names_.resize(n_mainclasses_);

			// Main Class sizes
			{
				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.Kinetics.ReactionClasses.MainClassSizes"));
				for (unsigned int k = 0; k < n_mainclasses_; k++)
				{
					unsigned int dummy;
					stream >> dummy;
					reaction_indices_[k].resize(dummy);
				}
			}

			// Class names
			{
				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.Kinetics.ReactionClasses.MainClassNames"));
				for (unsigned int k = 0; k < n_mainclasses_; k++)
				{
					stream >> reaction_mainclass_names_[k];
				}
			}

			{
				// Reaction indices
				unsigned int k = 0;
				BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree.get_child("opensmoke.Kinetics.ReactionClasses"))
				{
					boost::property_tree::ptree subtree = node.second;

					if (node.first == "MainClass")
					{
					//	std::string dummy = subtree.get_value< std::string >();
					//	std::stringstream stream(dummy);
					//	for (unsigned int j = 0; j < reaction_indices_[k].size(); j++)
					//	{
					//		stream >> reaction_indices_[k][j];
					//	}

						std::cout << k << std::endl;
						k++;
					}
				}
			}

			// Summary on the screen
			Summary(std::cout);
		}
		*/
	}
}

