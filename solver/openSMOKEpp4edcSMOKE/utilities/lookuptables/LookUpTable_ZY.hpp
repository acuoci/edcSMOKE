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
|   Copyright(C) 2021  Alberto Cuoci                                      |
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

#include "Grammar_LookUpTable_ZY.h"

namespace OpenSMOKE
{
	LookUpTable_ZY::LookUpTable_ZY(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap) :
		thermodynamicsMap_(thermodynamicsMap)
	{
		is_active_ = false;

		boost::filesystem::path path_output_;

		type_ = LookUpTable_ZY::UNIFORM_Z_UNIFORM_Y;
		is_apriori_ = true;
		is_pah_consumption_ = false;
		thermophoretic_model_ = 0;
		Le_ = 1.0;
		initialization_ = Initialization_Types::FROM_SCRATCH;

		nz_ = 0;
		nc_ = 0;
		dz_ = 0.;
		dc_ = 0.;
		dctilde_ = 0.;
		min_z_ = 1.e16;
		max_z_ = -1.e16;
		min_overall_c_ = 1.e16;
		max_overall_c_ = -1.e16;

		query_z_ = -1.e16;
		query_c_ = -1.e16;

		jC_ = thermodynamicsMap_.IndexOfElementWithoutError("C") - 1;
		jO_ = thermodynamicsMap_.IndexOfElementWithoutError("O") - 1;
		jH_ = thermodynamicsMap_.IndexOfElementWithoutError("H") - 1;

		WC_ = OpenSMOKE::AtomicWeights["C"];
		WO_ = OpenSMOKE::AtomicWeights["O"];
		WH_ = OpenSMOKE::AtomicWeights["H"];

		Y_def_PAH12_ = false;
		Y_def_PAH34_ = false;
		Y_def_PAHLP_ = false;
		Y_def_SP_ = false;
		Y_def_AGG_ = false;

		Y_def_alpha_PAH12_ = 0.;
		Y_def_alpha_PAH34_ = 0.;
		Y_def_alpha_PAHLP_ = 0.;
		Y_def_alpha_SP_ = 0.;
		Y_def_alpha_AGG_ = 0.;
	}

	void LookUpTable_ZY::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_LookUpTable_ZY grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@Path") == true)
			dictionary.ReadPath("@Path", path_output_);

		if (dictionary.CheckOption("@Apriori") == true)
			dictionary.ReadBool("@Apriori", is_apriori_);

		if (dictionary.CheckOption("@Fields") == true)
			dictionary.ReadOption("@Fields", list_fields_);

		std::string initialization = "";
		if (dictionary.CheckOption("@Initialization") == true)
			dictionary.ReadString("@Initialization", initialization);
		if (initialization == "FromScratch")
			initialization_ = Initialization_Types::FROM_SCRATCH;
		else if (initialization == "FromBackup")
			initialization_ = Initialization_Types::FROM_BACKUP;
		else if (initialization == "FromReconstruction")
			initialization_ = Initialization_Types::FROM_RECONSTRUCTION;
		else if (initialization == "FromOnTheFlySolution")
			initialization_ = Initialization_Types::FROM_ONTHEFLY_SOLUTION;
		else
		{
			OpenSMOKE::FatalErrorMessage("@Initialization options: FromScratch | FromBackup | FromReconstruction | FromOnTheFlySolution");
		}

		path_output_ = path_output_ / "lookuptable.main.xml";

		ReadTableFromXMLFile();

		Summary();

		CheckInput();
	}

	void LookUpTable_ZY::Setup(const boost::filesystem::path& path_to_xml, const std::vector<std::string>& list_fields, const bool is_apriori)
	{
		path_output_ = path_to_xml;
		list_fields_ = list_fields;
		is_apriori_ = is_apriori;

		ReadTableFromXMLFile();

		Summary();

		CheckInput();
	}

	void LookUpTable_ZY::ReadTableFromXMLFile()
	{
		{
			std::cout << "Open XML file: " << path_output_.string() << std::endl;

			boost::property_tree::ptree ptree;
			boost::property_tree::read_xml((path_output_).string(), ptree);

			// Reading fuel side
			{
				std::cout << "Reading fuel side composition..." << std::endl;

				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.YFuel"));

				unsigned int n = 0;
				stream >> n;

				Y_fuel_.resize(thermodynamicsMap_.NumberOfSpecies());
				Y_fuel_.setZero();
				for (unsigned int i = 0; i < n; i++)
				{
					std::string name_species;
					double Y;
					stream >> name_species;
					stream >> Y;

					const unsigned int j = thermodynamicsMap_.IndexOfSpecies(name_species) - 1;
					Y_fuel_[j] = Y;
				}

				// Mixture fraction
				{
					double mw;
					Eigen::VectorXd X(thermodynamicsMap_.NumberOfSpecies());
					thermodynamicsMap_.MoleFractions_From_MassFractions(X.data(), mw, Y_fuel_.data());

					// Reconstruct mixture fractions
					const double ZC = X.dot(thermodynamicsMap_.atomic_composition().col(jC_)) * WC_ / mw;
					const double ZH = X.dot(thermodynamicsMap_.atomic_composition().col(jH_)) * WH_ / mw;
					const double ZO = X.dot(thermodynamicsMap_.atomic_composition().col(jO_)) * WO_ / mw;

					ZstarFuel_ = 2. * ZC / WC_ + 0.50 * ZH / WH_ - ZO / WO_;
				}
			}

			// Reading oxidizer side
			{
				std::cout << "Reading oxidizer side composition..." << std::endl;

				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.YOx"));

				unsigned int n = 0;
				stream >> n;

				Y_ox_.resize(thermodynamicsMap_.NumberOfSpecies());
				Y_ox_.setZero();
				for (unsigned int i = 0; i < n; i++)
				{
					std::string name_species;
					double Y;
					stream >> name_species;
					stream >> Y;

					const unsigned int j = thermodynamicsMap_.IndexOfSpecies(name_species) - 1;
					Y_ox_[j] = Y;
				}

				// Mixture fraction 
				{
					double mw;
					Eigen::VectorXd X(thermodynamicsMap_.NumberOfSpecies());
					thermodynamicsMap_.MoleFractions_From_MassFractions(X.data(), mw, Y_ox_.data());

					// Reconstruct mixture fractions
					const double ZC = X.dot(thermodynamicsMap_.atomic_composition().col(jC_)) * WC_ / mw;
					const double ZH = X.dot(thermodynamicsMap_.atomic_composition().col(jH_)) * WH_ / mw;
					const double ZO = X.dot(thermodynamicsMap_.atomic_composition().col(jO_)) * WO_ / mw;

					ZstarOx_ = 2. * ZC / WC_ + 0.50 * ZH / WH_ - ZO / WO_;
				}
			}

			// Progress variable definition
			{
				std::cout << "Reading progress variable definition..." << std::endl;

				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.ProgressVariableDefinition"));

				unsigned int n = 0;
				stream >> n;

				Y_def_indices_.resize(0);
				Y_def_weights_.resize(0);
				for (unsigned int i = 0; i < n; i++)
				{
					std::string name_species;

					double alpha;
					stream >> name_species;
					stream >> alpha;

					if (name_species == "PAH12")
					{
						Y_def_PAH12_ = true;
						Y_def_alpha_PAH12_ = alpha;
					}
					else if (name_species == "PAH34")
					{
						Y_def_PAH34_ = true;
						Y_def_alpha_PAH34_ = alpha;
					}
					else if (name_species == "PAHLP")
					{
						Y_def_PAHLP_ = true;
						Y_def_alpha_PAHLP_ = alpha;
					}
					else if (name_species == "SP")
					{
						Y_def_SP_ = true;
						Y_def_alpha_SP_ = alpha;
					}
					else if (name_species == "AGG")
					{
						Y_def_AGG_ = true;
						Y_def_alpha_AGG_ = alpha;
					}
					else
					{
						Y_def_indices_.push_back(thermodynamicsMap_.IndexOfSpecies(name_species) - 1);
						Y_def_weights_.push_back(alpha / thermodynamicsMap_.MWs()[Y_def_indices_[i]]);
					}
				}

				// List of PAHs (if available)
				std::vector<std::string> PAH12_list = { "C6H6", "C7H8", "INDENE", "C10H8", "C12H8", "BIPHENYL", "FLUORENE", "C6H5C2H",
									"C6H5C2H3", "C6H5C2H5", "C10H7CH3",
									"C6H5", "C7H7", "C10H7", "C12H7", "C12H9", "CH3C6H4", "C6H4C2H", "C6H5C2H2",
									"C10H7CH2", "C10H6CH3", "XYLENE", "RXYLENE", "INDENYL" };

				std::vector<std::string> PAH34_list = { "C14H10", "C16H10", "C18H10", "C18H14", "C14H9", "C16H9", "C18H9", "C6H5C2H4C6H5", "C6H5CH2C6H5" };

				std::vector<std::string> PAHLP_list = { "BIN1A", "BIN1B", "BIN1C", "BIN1AJ", "BIN1BJ", "BIN1CJ",
									"BIN2A", "BIN2B", "BIN2C", "BIN2AJ", "BIN2BJ", "BIN2CJ",
									"BIN3A", "BIN3B", "BIN3C", "BIN3AJ", "BIN3BJ", "BIN3CJ",
									"BIN4A", "BIN4B", "BIN4C", "BIN4AJ", "BIN4BJ", "BIN4CJ" };

				std::vector<std::string> SP_list = { "BIN5A", "BIN5B", "BIN5C", "BIN5AJ", "BIN5BJ", "BIN5CJ",
									"BIN6A", "BIN6B", "BIN6C", "BIN6AJ", "BIN6BJ", "BIN6CJ",
									"BIN7A", "BIN7B", "BIN7C", "BIN7AJ", "BIN7BJ", "BIN7CJ",
									"BIN8A", "BIN8B", "BIN8C", "BIN8AJ", "BIN8BJ", "BIN8CJ",
																		"BIN9A", "BIN9B", "BIN9C", "BIN9AJ", "BIN9BJ", "BIN9CJ",
									"BIN10A", "BIN10B", "BIN10C", "BIN10AJ", "BIN10BJ", "BIN10CJ",
									"BIN11A", "BIN11B", "BIN11C", "BIN11AJ", "BIN11BJ", "BIN11CJ",
									"BIN12A", "BIN12B", "BIN12C", "BIN12AJ", "BIN12BJ", "BIN12CJ" };

				std::vector<std::string> AGG_list = { "BIN13A", "BIN13B", "BIN13C", "BIN13AJ", "BIN13BJ", "BIN13CJ",
									"BIN14A", "BIN14B", "BIN14C", "BIN14AJ", "BIN14BJ", "BIN14CJ",
									"BIN15A", "BIN15B", "BIN15C", "BIN15AJ", "BIN15BJ", "BIN15CJ",
									"BIN16A", "BIN16B", "BIN16C", "BIN16AJ", "BIN16BJ", "BIN16CJ",
									"BIN17A", "BIN17B", "BIN17C", "BIN17AJ", "BIN17BJ", "BIN17CJ",
									"BIN18A", "BIN18B", "BIN18C", "BIN18AJ", "BIN18BJ", "BIN18CJ",
																		"BIN19A", "BIN19B", "BIN19C", "BIN19AJ", "BIN19BJ", "BIN19CJ",
									"BIN20A", "BIN20B", "BIN20C", "BIN20AJ", "BIN20BJ", "BIN20CJ",
									"BIN21A", "BIN21B", "BIN21C", "BIN21AJ", "BIN21BJ", "BIN21CJ",
									"BIN22A", "BIN22B", "BIN22C", "BIN22AJ", "BIN22BJ", "BIN22CJ",
									"BIN23A", "BIN23B", "BIN23C", "BIN23AJ", "BIN23BJ", "BIN23CJ",
									"BIN24A", "BIN24B", "BIN24C", "BIN24AJ", "BIN24BJ", "BIN24CJ",
									"BIN25A", "BIN25B", "BIN25C", "BIN25AJ", "BIN25BJ", "BIN25CJ" };



				indices_pah12_.resize(0);
				for (unsigned int i = 0; i < PAH12_list.size(); i++)
					if (thermodynamicsMap_.IndexOfSpeciesWithoutError(PAH12_list[i]) > 0)
					{
						std::cout << "PAH12: " << PAH12_list[i] << std::endl;
						indices_pah12_.push_back(thermodynamicsMap_.IndexOfSpeciesWithoutError(PAH12_list[i]) - 1);
					}

				indices_pah34_.resize(0);
				for (unsigned int i = 0; i < PAH34_list.size(); i++)
					if (thermodynamicsMap_.IndexOfSpeciesWithoutError(PAH34_list[i]) > 0)
					{
						std::cout << "PAH34: " << PAH34_list[i] << std::endl;
						indices_pah34_.push_back(thermodynamicsMap_.IndexOfSpeciesWithoutError(PAH34_list[i]) - 1);
					}

				indices_pahlp_.resize(0);
				for (unsigned int i = 0; i < PAHLP_list.size(); i++)
					if (thermodynamicsMap_.IndexOfSpeciesWithoutError(PAHLP_list[i]) > 0)
					{
						std::cout << "PAHLP: " << PAHLP_list[i] << std::endl;
						indices_pahlp_.push_back(thermodynamicsMap_.IndexOfSpeciesWithoutError(PAHLP_list[i]) - 1);
					}

				indices_sp_.resize(0);
				for (unsigned int i = 0; i < SP_list.size(); i++)
					if (thermodynamicsMap_.IndexOfSpeciesWithoutError(SP_list[i]) > 0)
					{
						std::cout << "SP: " << SP_list[i] << std::endl;
						indices_sp_.push_back(thermodynamicsMap_.IndexOfSpeciesWithoutError(SP_list[i]) - 1);
					}

				indices_agg_.resize(0);
				for (unsigned int i = 0; i < AGG_list.size(); i++)
					if (thermodynamicsMap_.IndexOfSpeciesWithoutError(AGG_list[i]) > 0)
					{
						std::cout << "AGG: " << AGG_list[i] << std::endl;
						indices_agg_.push_back(thermodynamicsMap_.IndexOfSpeciesWithoutError(AGG_list[i]) - 1);
					}
			}

			// Mixture fraction points
			{
				std::cout << "Reading mixture fraction..." << std::endl;

				nz_ = ptree.get<unsigned int>("opensmoke.Z-points");
				z_.resize(nz_);

				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.Z-coordinates"));
				for (unsigned int i = 0; i < nz_; i++)
					stream >> z_[i];

				dz_ = z_[1] - z_[0];
				min_z_ = z_[0];
				max_z_ = z_[nz_ - 1];
			}

			// Progress Variable
			{
				std::cout << "Reading progress variable..." << std::endl;

				nc_ = ptree.get<unsigned int>("opensmoke.C-points");

				{
					std::stringstream stream;
					stream.str(ptree.get< std::string >("opensmoke.C-coordinates"));

					c_.resize(nc_);
					for (unsigned int i = 0; i < nc_; i++)
						stream >> c_[i];
					dc_ = c_[1] - c_[0];
				}

				{

					std::stringstream stream;
					stream.str(ptree.get< std::string >("opensmoke.Ctilde-coordinates"));

					ctilde_.resize(nc_);
					for (unsigned int i = 0; i < nc_; i++)
						stream >> ctilde_[i];
					dctilde_ = ctilde_[1] - ctilde_[0];
				}
			}

			// Min and max values
			{
				{
					std::stringstream stream;
					stream.str(ptree.get< std::string >("opensmoke.Cmin-values"));

					c_min_.resize(nz_);
					for (unsigned int i = 0; i < nz_; i++)
						stream >> c_min_[i];

					min_overall_c_ = *std::min_element(c_min_.begin(), c_min_.end());
				}

				{
					std::stringstream stream;
					stream.str(ptree.get< std::string >("opensmoke.Cmax-values"));

					c_max_.resize(nz_);
					for (unsigned int i = 0; i < nz_; i++)
						stream >> c_max_[i];

					max_overall_c_ = *std::max_element(c_max_.begin(), c_max_.end());
				}

				// Minimum/Maximum local progress variable
				std::cout << "Checking Cmin and Cmax ..." << std::endl;
				for (unsigned int i = 0; i < nz_; i++)
				{
					const double Cmin = MinimumLocalProgressVariable(z_[i]);
					const double Cmax = MaximumLocalProgressVariable(z_[i]);

					if (Cmax - Cmin <= 0.)
						std::cout << "WARNING: Cmax<=Cmin at point: z: " << z_[i] << " Cmin: " << Cmin << " Cmax: " << Cmax << std::endl;
				}
			}
		}

		// Mass fractions (-)
		tables_Y_.resize(list_fields_.size());
		min_fields_Y_.resize(list_fields_.size());
		max_fields_Y_.resize(list_fields_.size());
		mean_fields_Y_.resize(list_fields_.size());
		interpolated_Y_.resize(list_fields_.size());
		std::fill(interpolated_Y_.begin(), interpolated_Y_.end(), 0.);

		// Production terms (kg/m3/s)
		tables_omegap_.resize(list_fields_.size());
		min_fields_omegap_.resize(list_fields_.size());
		max_fields_omegap_.resize(list_fields_.size());
		mean_fields_omegap_.resize(list_fields_.size());
		interpolated_omegap_.resize(list_fields_.size());
		std::fill(interpolated_omegap_.begin(), interpolated_omegap_.end(), 0.);

		// Consumption terms (kg/m3/s)
		tables_omegad_.resize(list_fields_.size());
		min_fields_omegad_.resize(list_fields_.size());
		max_fields_omegad_.resize(list_fields_.size());
		mean_fields_omegad_.resize(list_fields_.size());
		interpolated_omegad_.resize(list_fields_.size());
		std::fill(interpolated_omegad_.begin(), interpolated_omegad_.end(), 0.);


		// Fields for which the production/consumption source terms are available
		available_omega_.resize(list_fields_.size());
		list_available_omega_.resize(0);
		std::fill(available_omega_.begin(), available_omega_.end(), false);

		// Read individual lookup tables
		for (unsigned int k = 0; k < list_fields_.size(); k++)
		{
			std::string path_lookup_table = (path_output_).string();
			boost::replace_all(path_lookup_table, ".main.xml", ("." + list_fields_[k] + ".xml"));

			boost::property_tree::ptree ptree;
			boost::property_tree::read_xml(path_lookup_table, ptree);

			std::cout << "Reading " << list_fields_[k] << " table..." << std::endl;

			// Mass fractions (-)
			{
				tables_Y_[k].resize(nc_);
				for (unsigned int i = 0; i < nc_; i++)
					tables_Y_[k][i].resize(nz_);

				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.Table_" + list_fields_[k]));

				for (unsigned int i = 0; i < nc_; i++)
					for (unsigned int j = 0; j < nz_; j++)
						stream >> tables_Y_[k][i][j];
			}

			// Production terms (kg/m3/s)
			if (ptree.get_child_optional("opensmoke.Table_OmegaP_" + list_fields_[k]))
			{
				tables_omegap_[k].resize(nc_);
				for (unsigned int i = 0; i < nc_; i++)
					tables_omegap_[k][i].resize(nz_);

				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.Table_OmegaP_" + list_fields_[k]));

				for (unsigned int i = 0; i < nc_; i++)
					for (unsigned int j = 0; j < nz_; j++)
						stream >> tables_omegap_[k][i][j];
			}

			// Consumption rates (kg/m3/s)
			if (ptree.get_child_optional("opensmoke.Table_OmegaD_" + list_fields_[k]))
			{
				tables_omegad_[k].resize(nc_);
				for (unsigned int i = 0; i < nc_; i++)
					tables_omegad_[k][i].resize(nz_);

				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.Table_OmegaD_" + list_fields_[k]));

				for (unsigned int i = 0; i < nc_; i++)
					for (unsigned int j = 0; j < nz_; j++)
						stream >> tables_omegad_[k][i][j];
			}

			// Check fo availability of production/consumption terms
			if (tables_omegap_[k].size() != 0 && tables_omegad_[k].size() != 0)
			{
				available_omega_[k] = true;
				list_available_omega_.push_back(k);
			}
			else if ((tables_omegap_[k].size() != 0 && tables_omegad_[k].size() == 0) ||
				(tables_omegap_[k].size() == 0 && tables_omegad_[k].size() != 0))
			{
				OpenSMOKE::FatalErrorMessage("Production/Consumption terms have to be provided together");
			}


			// Max/min values (mass fractions)
			min_fields_Y_[k] = 1.e16;
			max_fields_Y_[k] = -1.e16;
			mean_fields_Y_[k] = 0.;
			for (unsigned int i = 0; i < nc_; i++)
			{
				auto minmax = std::minmax_element(begin(tables_Y_[k][i]), end(tables_Y_[k][i]));
				if (*minmax.first < min_fields_Y_[k]) min_fields_Y_[k] = *minmax.first;
				if (*minmax.second > max_fields_Y_[k]) max_fields_Y_[k] = *minmax.second;

				mean_fields_Y_[k] += std::accumulate(tables_Y_[k][i].begin(), tables_Y_[k][i].end(), 0.0) / static_cast<double>(tables_Y_[k][i].size());
			}
			mean_fields_Y_[k] /= static_cast<double>(nc_);

			// Max/min values (production terms)
			if (available_omega_[k] == true)
			{
				min_fields_omegap_[k] = 1.e16;
				max_fields_omegap_[k] = -1.e16;
				mean_fields_omegap_[k] = 0.;
				for (unsigned int i = 0; i < nc_; i++)
				{
					auto minmax = std::minmax_element(begin(tables_omegap_[k][i]), end(tables_omegap_[k][i]));
					if (*minmax.first < min_fields_omegap_[k]) min_fields_omegap_[k] = *minmax.first;
					if (*minmax.second > max_fields_omegap_[k]) max_fields_omegap_[k] = *minmax.second;

					mean_fields_omegap_[k] += std::accumulate(tables_omegap_[k][i].begin(), tables_omegap_[k][i].end(), 0.0) / static_cast<double>(tables_omegap_[k][i].size());
				}
				mean_fields_omegap_[k] /= static_cast<double>(nc_);
			}

			// Max/min values (consumption terms)
			if (available_omega_[k] == true)
			{
				min_fields_omegad_[k] = 1.e16;
				max_fields_omegad_[k] = -1.e16;
				mean_fields_omegad_[k] = 0.;
				for (unsigned int i = 0; i < nc_; i++)
				{
					auto minmax = std::minmax_element(begin(tables_omegad_[k][i]), end(tables_omegad_[k][i]));
					if (*minmax.first < min_fields_omegad_[k]) min_fields_omegad_[k] = *minmax.first;
					if (*minmax.second > max_fields_omegad_[k]) max_fields_omegad_[k] = *minmax.second;

					mean_fields_omegad_[k] += std::accumulate(tables_omegad_[k][i].begin(), tables_omegad_[k][i].end(), 0.0) / static_cast<double>(tables_omegad_[k][i].size());
				}
				mean_fields_omegad_[k] /= static_cast<double>(nc_);
			}
		}

		is_active_ = true;
	}

	void LookUpTable_ZY::CheckInput()
	{
		std::cout << "Testing lookup table..." << std::endl;

		// First line
		std::cout << std::left << std::setw(16) << "Z";
		std::cout << std::left << std::setw(16) << "C";
		for (unsigned int k = 0; k < std::min(static_cast<int>(interpolated_Y_.size()), 4); k++)
			std::cout << std::left << std::setw(16) << list_fields_[k];
		std::cout << std::endl;

		// Internal points
		{
			const unsigned int np = 5;
			for (unsigned int i = 0; i < np; i++)
				for (unsigned int j = 0; j < np; j++)
				{
					const double z = min_z_ + (max_z_ - min_z_) / static_cast<double>(np - 1) * static_cast<double>(i);
					const double c = min_overall_c_ + (max_overall_c_ - min_overall_c_) / static_cast<double>(np - 1) * static_cast<double>(j);

					Interpolate(z, c);

					std::cout << std::left << std::setw(16) << z;
					std::cout << std::left << std::setw(16) << c;
					for (unsigned int k = 0; k < std::min(static_cast<int>(interpolated_Y_.size()), 4); k++)
						std::cout << std::left << std::setw(16) << interpolated_Y_[k];

					bool iWarning = false;
					for (unsigned int k = 0; k < interpolated_Y_.size(); k++)
						if ((interpolated_Y_[k] > max_fields_Y_[k]) || (interpolated_Y_[k] < min_fields_Y_[k]))
							iWarning = true;
					if (iWarning == true) std::cout << " WARNING!!!";

					std::cout << std::endl;
				}
		}

		// External points
		{
			const unsigned int np = 5;
			for (unsigned int i = 0; i < np; i++)
				for (unsigned int j = 0; j < np; j++)
				{
					const double z = min_z_ + ((max_z_ + 0.001) - (min_z_ - 0.001)) / static_cast<double>(np - 1) * static_cast<double>(i);
					const double c = min_overall_c_ + ((max_overall_c_ + 0.001) - (min_overall_c_ - 0.001)) / static_cast<double>(np - 1) * static_cast<double>(j);

					Interpolate(z, c);

					std::cout << std::left << std::setw(16) << z;
					std::cout << std::left << std::setw(16) << c;
					for (unsigned int k = 0; k < std::min(static_cast<int>(interpolated_Y_.size()), 4); k++)
						std::cout << std::left << std::setw(16) << interpolated_Y_[k];

					bool iWarning = false;
					for (unsigned int k = 0; k < interpolated_Y_.size(); k++)
						if ((interpolated_Y_[k] > max_fields_Y_[k]) || (interpolated_Y_[k] < min_fields_Y_[k]))
							iWarning = true;
					if (iWarning == true) std::cout << " WARNING!!!";

					std::cout << std::endl;
				}
		}
		std::cout << std::endl;

		// Points directly available in the table
		std::cout << "Checking lookup table on points directly available..." << std::endl;
		{
			std::vector<double> error_mean(interpolated_Y_.size());
			std::fill(error_mean.begin(), error_mean.end(), 0.);
			std::vector<double> error_max(interpolated_Y_.size());
			std::fill(error_max.begin(), error_max.end(), -1.e16);

			std::vector<double> error_mean_omegap(interpolated_omegap_.size());
			std::fill(error_mean_omegap.begin(), error_mean_omegap.end(), 0.);
			std::vector<double> error_max_omegap(interpolated_omegap_.size());
			std::fill(error_max_omegap.begin(), error_max_omegap.end(), -1.e16);

			std::vector<double> error_mean_omegad(interpolated_omegad_.size());
			std::fill(error_mean_omegad.begin(), error_mean_omegad.end(), 0.);
			std::vector<double> error_max_omegad(interpolated_omegad_.size());
			std::fill(error_max_omegad.begin(), error_max_omegad.end(), -1.e16);

			for (unsigned int i = 0; i < nz_; i++)
				for (unsigned int j = 0; j < nc_; j++)
				{
					const double z = z_[i];
					const double ctilde = ctilde_[j];

					InterpolateFromNormalizedProgressVariable(z, ctilde);

					for (unsigned int k = 0; k < interpolated_Y_.size(); k++)
					{
						const double error = std::fabs(interpolated_Y_[k] - tables_Y_[k][j][i]) / std::fabs(mean_fields_Y_[k]);

						if (error > 1.e-2)
							std::cout << "ERROR: Field: " << k << " Point: " << z << " " << ctilde << " Error: " << error << " Table: " << tables_Y_[k][j][i] << " Interpolated: " << interpolated_Y_[k] << std::endl;

						error_mean[k] += error;
						if (error > error_max[k]) error_max[k] = error;
					}

					for (unsigned int k = 0; k < interpolated_omegap_.size(); k++)
						if (tables_omegap_[k].size() != 0)
						{
							const double error = std::fabs(interpolated_omegap_[k] - tables_omegap_[k][j][i]) / std::fabs(mean_fields_omegap_[k]);

							if (error > 1.e-2)
								std::cout << "ERROR: Field: " << k << " Point: " << z << " " << ctilde << " Error: " << error << " Table: " << tables_omegap_[k][j][i] << " Interpolated: " << interpolated_omegap_[k] << std::endl;

							error_mean_omegap[k] += error;
							if (error > error_max_omegap[k]) error_max_omegap[k] = error;
						}

					for (unsigned int k = 0; k < interpolated_omegad_.size(); k++)
						if (tables_omegad_[k].size() != 0)
						{
							const double error = std::fabs(interpolated_omegad_[k] - tables_omegad_[k][j][i]) / std::fabs(mean_fields_omegad_[k]);

							if (error > 1.e-2)
								std::cout << "ERROR: Field: " << k << " Point: " << z << " " << ctilde << " Error: " << error << " Table: " << tables_omegad_[k][j][i] << " Interpolated: " << interpolated_omegad_[k] << std::endl;

							error_mean_omegad[k] += error;
							if (error > error_max_omegad[k]) error_max_omegad[k] = error;
						}
				}

			for (unsigned int k = 0; k < interpolated_Y_.size(); k++)
				error_mean[k] /= static_cast<double>(nz_ * nc_);

			for (unsigned int k = 0; k < interpolated_omegap_.size(); k++)
				error_mean_omegap[k] /= static_cast<double>(nz_ * nc_);

			for (unsigned int k = 0; k < interpolated_omegad_.size(); k++)
				error_mean_omegad[k] /= static_cast<double>(nz_ * nc_);

			std::cout << " * Summary of errors (mass fractions)" << std::endl;
			for (unsigned int k = 0; k < interpolated_Y_.size(); k++)
				std::cout << "   " << list_fields_[k] << " " << error_mean[k] << " " << error_max[k] << std::endl;

			std::cout << " * Summary of errors (production terms)" << std::endl;
			for (unsigned int k = 0; k < interpolated_omegap_.size(); k++)
				if (tables_omegap_[k].size() != 0)
					std::cout << "   " << list_fields_[k] << " " << error_mean_omegap[k] << " " << error_max_omegap[k] << std::endl;

			std::cout << " * Summary of errors (consumption terms)" << std::endl;
			for (unsigned int k = 0; k < interpolated_omegad_.size(); k++)
				if (tables_omegad_[k].size() != 0)
					std::cout << "   " << list_fields_[k] << " " << error_mean_omegad[k] << " " << error_max_omegad[k] << std::endl;
		}
		std::cout << std::endl;
	}

	void LookUpTable_ZY::Summary()
	{
		std::cout << "Progress variable definition..." << std::endl;

		std::cout << "Regular species" << std::endl;
		for (unsigned int i = 0; i < Y_def_indices_.size(); i++)
		{
			std::cout << " * " << thermodynamicsMap_.NamesOfSpecies()[Y_def_indices_[i]] << "\t"
				<< Y_def_weights_[i] * thermodynamicsMap_.MWs()[Y_def_indices_[i]] << "\t"
				<< Y_def_weights_[i] << std::endl;
		}

		std::cout << "PAHs (soot precursors)" << std::endl;
		{
			std::cout << " * PAH12 " << Y_def_alpha_PAH12_ << "\t(" << indices_pah12_.size() << " species)" << std::endl;
			std::cout << " * PAH34 " << Y_def_alpha_PAH34_ << "\t(" << indices_pah34_.size() << " species)" << std::endl;
			std::cout << " * PAHLP " << Y_def_alpha_PAHLP_ << "\t(" << indices_pahlp_.size() << " species)" << std::endl;
		}

		std::cout << "Soot" << std::endl;
		{
			std::cout << " * SP    " << Y_def_alpha_SP_ << "\t(" << indices_sp_.size() << " species)" << std::endl;
			std::cout << " * AGG   " << Y_def_alpha_AGG_ << "\t(" << indices_agg_.size() << " species)" << std::endl;
		}

		std::cout << std::endl;
		std::cout << "Mixture fraction space:  " << min_z_ << " " << max_z_ << " (" << nz_ << ")" << std::endl;
		std::cout << "Progress variable space: " << min_overall_c_ << " " << max_overall_c_ << " (" << nc_ << ")" << std::endl;
		std::cout << "Zf*: " << ZstarFuel_ << " Zox*: " << ZstarOx_ << std::endl;

		std::cout << "Min/Max/Mean Values (mass fractions)" << std::endl;
		for (unsigned int k = 0; k < list_fields_.size(); k++)
			std::cout << list_fields_[k] << ": " << min_fields_Y_[k] << " " << max_fields_Y_[k] << " " << mean_fields_Y_[k] << std::endl;

		std::cout << "Min/Max/Mean Values (production terms)" << std::endl;
		for (unsigned int k = 0; k < list_fields_.size(); k++)
			if (tables_omegap_[k].size() != 0)
				std::cout << list_fields_[k] << ": " << min_fields_omegap_[k] << " " << max_fields_omegap_[k] << " " << mean_fields_omegap_[k] << std::endl;

		std::cout << "Min/Max/Mean Values (consumption terms)" << std::endl;
		for (unsigned int k = 0; k < list_fields_.size(); k++)
			if (tables_omegad_[k].size() != 0)
				std::cout << list_fields_[k] << ": " << min_fields_omegad_[k] << " " << max_fields_omegad_[k] << " " << mean_fields_omegad_[k] << std::endl;

	}

	double LookUpTable_ZY::MinimumLocalProgressVariable(const double z)
	{
		double Cmin = c_min_[nz_ - 1];
		{
			const double query_z = std::min(std::max(z, 0.), 1.);
			const unsigned int iz = std::floor((query_z - min_z_) / dz_);

			if (iz < nz_ - 1)
			{
				const double Q1 = c_min_[iz];
				const double Q2 = c_min_[iz + 1];
				Cmin = Q1 + (Q2 - Q1) / dz_ * (query_z - z_[iz]);
			}
		}

		return Cmin;
	}

	double LookUpTable_ZY::MaximumLocalProgressVariable(const double z)
	{
		double Cmax = c_max_[nz_ - 1];
		{
			const double query_z = std::max(std::min(z, 1.), 0.);
			const unsigned int iz = std::floor((query_z - min_z_) / dz_);

			if (iz < nz_ - 1)
			{
				const double Q1 = c_max_[iz];
				const double Q2 = c_max_[iz + 1];
				Cmax = Q1 + (Q2 - Q1) / dz_ * (query_z - z_[iz]);
			}
		}

		return Cmax;
	}

	void LookUpTable_ZY::InterpolateFromNormalizedProgressVariable(const double z, const double ctilde)
	{
		// Minimum/Maximum local progress variable
		const double Cmin = MinimumLocalProgressVariable(z);
		const double Cmax = MaximumLocalProgressVariable(z);
		const double c = ctilde * (Cmax - Cmin) + Cmin;

		Interpolate(z, c);
	}

	void LookUpTable_ZY::Interpolate(const double z, const double c)
	{
		if (z == query_z_ && c == query_c_)
			return;

		query_z_ = z;
		query_c_ = c;

		// Minimum/Maximum local progress variable
		const double Cmin = MinimumLocalProgressVariable(query_z_);
		const double Cmax = MaximumLocalProgressVariable(query_z_);

		// In case of anomalous points
		if (Cmax - Cmin <= 0.)
		{
			const unsigned int iz = std::max(std::min(static_cast<int>(std::floor((query_z_ - min_z_) / dz_)), static_cast<int>(nz_ - 1)), 0);

			for (unsigned int k = 0; k < list_fields_.size(); k++)
				interpolated_Y_[k] = tables_Y_[k][0][iz];

			for (unsigned int k = 0; k < list_available_omega_.size(); k++)
			{
				const unsigned int j = list_available_omega_[k];
				interpolated_omegap_[j] = tables_omegap_[j][0][iz];
			}

			for (unsigned int k = 0; k < list_available_omega_.size(); k++)
			{
				const unsigned int j = list_available_omega_[k];
				interpolated_omegad_[j] = tables_omegad_[j][0][iz];
			}
		}

		else
		{
			const double query_z = std::max(std::min(query_z_, 1.), 0.);
			const double query_ctilde = std::max(std::min((query_c_ - Cmin) / (Cmax - Cmin), 1.), 0.);

			const unsigned int iz = std::floor((query_z - min_z_) / dz_);
			const unsigned int ic = std::floor((query_ctilde - 0.) / dctilde_);

			// Regular points
			if ((iz < nz_ - 1) && (ic < nc_ - 1))
			{
				// Interpolation coefficients
				const double z_minus_z1 = query_z - z_[iz];
				const double z2_minus_z = z_[iz + 1] - query_z;
				const double c_minus_c1 = query_ctilde - ctilde_[ic];
				const double c2_minus_c = ctilde_[ic + 1] - query_ctilde;
				const double coeff = 1. / (dz_ * dctilde_);

				// Mass fractions
				for (unsigned int k = 0; k < list_fields_.size(); k++)
				{
					// Quadrature points
					const double Q11 = tables_Y_[k][ic][iz];
					const double Q21 = tables_Y_[k][ic][iz + 1];
					const double Q12 = tables_Y_[k][ic + 1][iz];
					const double Q22 = tables_Y_[k][ic + 1][iz + 1];

					// Interpolation
					interpolated_Y_[k] = coeff * (z2_minus_z * (Q11 * c2_minus_c + Q12 * c_minus_c1) +
						z_minus_z1 * (Q21 * c2_minus_c + Q22 * c_minus_c1));
				}

				// Production terms
				for (unsigned int k = 0; k < list_available_omega_.size(); k++)
				{
					const unsigned int j = list_available_omega_[k];

					// Quadrature points
					const double Q11 = tables_omegap_[j][ic][iz];
					const double Q21 = tables_omegap_[j][ic][iz + 1];
					const double Q12 = tables_omegap_[j][ic + 1][iz];
					const double Q22 = tables_omegap_[j][ic + 1][iz + 1];

					// Interpolation
					interpolated_omegap_[j] = coeff * (z2_minus_z * (Q11 * c2_minus_c + Q12 * c_minus_c1) +
						z_minus_z1 * (Q21 * c2_minus_c + Q22 * c_minus_c1));
				}

				// Consumption terms
				for (unsigned int k = 0; k < list_available_omega_.size(); k++)
				{
					const unsigned int j = list_available_omega_[k];

					// Quadrature points
					const double Q11 = tables_omegad_[j][ic][iz];
					const double Q21 = tables_omegad_[j][ic][iz + 1];
					const double Q12 = tables_omegad_[j][ic + 1][iz];
					const double Q22 = tables_omegad_[j][ic + 1][iz + 1];

					// Interpolation
					interpolated_omegad_[j] = coeff * (z2_minus_z * (Q11 * c2_minus_c + Q12 * c_minus_c1) +
						z_minus_z1 * (Q21 * c2_minus_c + Q22 * c_minus_c1));
				}
			}

			else if ((iz < nz_ - 1) && (ic >= nc_ - 1))
			{
				const double z1 = z_[iz];

				// Mass fractions
				for (unsigned int k = 0; k < list_fields_.size(); k++)
				{
					const double Q1 = tables_Y_[k][nc_ - 1][iz];
					const double Q2 = tables_Y_[k][nc_ - 1][iz + 1];

					interpolated_Y_[k] = Q1 + (Q2 - Q1) / dz_ * (query_z - z1);
				}

				// Production terms
				for (unsigned int k = 0; k < list_available_omega_.size(); k++)
				{
					const unsigned int j = list_available_omega_[k];

					const double Q1 = tables_omegap_[j][nc_ - 1][iz];
					const double Q2 = tables_omegap_[j][nc_ - 1][iz + 1];

					interpolated_omegap_[j] = Q1 + (Q2 - Q1) / dz_ * (query_z - z1);
				}

				// Consumption terms
				for (unsigned int k = 0; k < list_available_omega_.size(); k++)
				{
					const unsigned int j = list_available_omega_[k];

					const double Q1 = tables_omegad_[j][nc_ - 1][iz];
					const double Q2 = tables_omegad_[j][nc_ - 1][iz + 1];

					interpolated_omegad_[j] = Q1 + (Q2 - Q1) / dz_ * (query_z - z1);
				}

			}

			else if ((iz >= nz_ - 1) && (ic < nc_ - 1))
			{
				const double c1 = ctilde_[ic];

				// Mass fractions
				for (unsigned int k = 0; k < list_fields_.size(); k++)
				{
					const double Q1 = tables_Y_[k][ic][nz_ - 1];
					const double Q2 = tables_Y_[k][ic + 1][nz_ - 1];

					interpolated_Y_[k] = Q1 + (Q2 - Q1) / dctilde_ * (query_ctilde - c1);
				}

				// Production terms
				for (unsigned int k = 0; k < list_available_omega_.size(); k++)
				{
					const unsigned int j = list_available_omega_[k];

					const double Q1 = tables_omegap_[j][ic][nz_ - 1];
					const double Q2 = tables_omegap_[j][ic + 1][nz_ - 1];

					interpolated_omegap_[j] = Q1 + (Q2 - Q1) / dctilde_ * (query_ctilde - c1);
				}

				// Consumption terms
				for (unsigned int k = 0; k < list_available_omega_.size(); k++)
				{
					const unsigned int j = list_available_omega_[k];

					const double Q1 = tables_omegad_[j][ic][nz_ - 1];
					const double Q2 = tables_omegad_[j][ic + 1][nz_ - 1];

					interpolated_omegad_[j] = Q1 + (Q2 - Q1) / dctilde_ * (query_ctilde - c1);
				}
			}

			else if ((iz >= nz_ - 1) && (ic >= nc_ - 1))
			{
				// Mass fractions
				for (unsigned int k = 0; k < list_fields_.size(); k++)
					interpolated_Y_[k] = tables_Y_[k][nc_ - 1][nz_ - 1];

				// Production terms
				for (unsigned int k = 0; k < list_available_omega_.size(); k++)
				{
					const unsigned int j = list_available_omega_[k];
					interpolated_omegap_[j] = tables_omegap_[j][nc_ - 1][nz_ - 1];
				}

				// Consumption terms
				for (unsigned int k = 0; k < list_available_omega_.size(); k++)
				{
					const unsigned int j = list_available_omega_[k];
					interpolated_omegad_[j] = tables_omegad_[j][nc_ - 1][nz_ - 1];
				}
			}

			else
			{
				std::cout << "WARNING: Unmapped point: " << query_z_ << " " << query_c_ << std::endl;
			}
		}
	}

	void LookUpTable_ZY::Interpolate(const unsigned int k, const double z, const double c)
	{
		query_z_ = z;
		query_c_ = c;

		// Minimum/Maximum local progress variable
		const double Cmin = MinimumLocalProgressVariable(query_z_);
		const double Cmax = MaximumLocalProgressVariable(query_z_);

		// In case of anomalous points
		if (Cmax - Cmin <= 0.)
		{
			const unsigned int iz = std::max(std::min(static_cast<int>(std::floor((query_z_ - min_z_) / dz_)), static_cast<int>(nz_ - 1)), 0);

			// Mass fractions
			interpolated_Y_[k] = tables_Y_[k][0][iz];

			// Production terms
			interpolated_omegap_[k] = tables_omegap_[k][0][iz];

			// Consumption terms
			interpolated_omegad_[k] = tables_omegad_[k][0][iz];
		}

		else
		{
			const double query_z = std::max(std::min(query_z_, 1.), 0.);
			const double query_ctilde = std::max(std::min((query_c_ - Cmin) / (Cmax - Cmin), 1.), 0.);

			const unsigned int iz = std::floor((query_z - min_z_) / dz_);
			const unsigned int ic = std::floor((query_ctilde - 0.) / dctilde_);

			// Regular points
			if ((iz < nz_ - 1) && (ic < nc_ - 1))
			{
				// Interpolation coefficients
				const double z_minus_z1 = query_z - z_[iz];
				const double z2_minus_z = z_[iz + 1] - query_z;
				const double c_minus_c1 = query_ctilde - ctilde_[ic];
				const double c2_minus_c = ctilde_[ic + 1] - query_ctilde;
				const double coeff = 1. / (dz_ * dctilde_);

				// Mass fractions
				{
					const double Q11 = tables_Y_[k][ic][iz];
					const double Q21 = tables_Y_[k][ic][iz + 1];
					const double Q12 = tables_Y_[k][ic + 1][iz];
					const double Q22 = tables_Y_[k][ic + 1][iz + 1];

					interpolated_Y_[k] = coeff * (z2_minus_z * (Q11 * c2_minus_c + Q12 * c_minus_c1) +
						z_minus_z1 * (Q21 * c2_minus_c + Q22 * c_minus_c1));
				}

				// Production terms
				{
					const double Q11 = tables_omegap_[k][ic][iz];
					const double Q21 = tables_omegap_[k][ic][iz + 1];
					const double Q12 = tables_omegap_[k][ic + 1][iz];
					const double Q22 = tables_omegap_[k][ic + 1][iz + 1];

					interpolated_omegap_[k] = coeff * (z2_minus_z * (Q11 * c2_minus_c + Q12 * c_minus_c1) +
						z_minus_z1 * (Q21 * c2_minus_c + Q22 * c_minus_c1));
				}

				// Consumption terms
				{
					const double Q11 = tables_omegad_[k][ic][iz];
					const double Q21 = tables_omegad_[k][ic][iz + 1];
					const double Q12 = tables_omegad_[k][ic + 1][iz];
					const double Q22 = tables_omegad_[k][ic + 1][iz + 1];

					interpolated_omegad_[k] = coeff * (z2_minus_z * (Q11 * c2_minus_c + Q12 * c_minus_c1) +
						z_minus_z1 * (Q21 * c2_minus_c + Q22 * c_minus_c1));
				}
			}

			else if ((iz < nz_ - 1) && (ic >= nc_ - 1))
			{
				const double z1 = z_[iz];

				// Mass fractions
				{
					const double Q1 = tables_Y_[k][nc_ - 1][iz];
					const double Q2 = tables_Y_[k][nc_ - 1][iz + 1];

					interpolated_Y_[k] = Q1 + (Q2 - Q1) / dz_ * (query_z - z1);
				}

				// Production terms
				{
					const double Q1 = tables_omegap_[k][nc_ - 1][iz];
					const double Q2 = tables_omegap_[k][nc_ - 1][iz + 1];

					interpolated_omegap_[k] = Q1 + (Q2 - Q1) / dz_ * (query_z - z1);
				}

				// Consumption terms
				{
					const double Q1 = tables_omegad_[k][nc_ - 1][iz];
					const double Q2 = tables_omegad_[k][nc_ - 1][iz + 1];

					interpolated_omegad_[k] = Q1 + (Q2 - Q1) / dz_ * (query_z - z1);
				}

			}

			else if ((iz >= nz_ - 1) && (ic < nc_ - 1))
			{
				const double c1 = ctilde_[ic];

				// Mass fractions
				{
					const double Q1 = tables_Y_[k][ic][nz_ - 1];
					const double Q2 = tables_Y_[k][ic + 1][nz_ - 1];

					interpolated_Y_[k] = Q1 + (Q2 - Q1) / dctilde_ * (query_ctilde - c1);
				}

				// Production terms
				{
					const double Q1 = tables_omegap_[k][ic][nz_ - 1];
					const double Q2 = tables_omegap_[k][ic + 1][nz_ - 1];

					interpolated_omegap_[k] = Q1 + (Q2 - Q1) / dctilde_ * (query_ctilde - c1);
				}

				// Consumption terms
				{
					const double Q1 = tables_omegad_[k][ic][nz_ - 1];
					const double Q2 = tables_omegad_[k][ic + 1][nz_ - 1];

					interpolated_omegad_[k] = Q1 + (Q2 - Q1) / dctilde_ * (query_ctilde - c1);
				}
			}

			else if ((iz >= nz_ - 1) && (ic >= nc_ - 1))
			{
				// Mass fractions
				interpolated_Y_[k] = tables_Y_[k][nc_ - 1][nz_ - 1];

				// Production terms
				interpolated_omegap_[k] = tables_omegap_[k][nc_ - 1][nz_ - 1];

				// Consumption terms
				interpolated_omegad_[k] = tables_omegad_[k][nc_ - 1][nz_ - 1];
			}

			else
			{
				std::cout << "WARNING: Unmapped point: " << query_z_ << " " << query_c_ << std::endl;
			}
		}
	}



	int LookUpTable_ZY::Index(const std::string name) const
	{
		const auto itr = std::find(list_fields_.begin(), list_fields_.end(), name);
		if (itr != list_fields_.cend())
			return std::distance(list_fields_.begin(), itr);

		OpenSMOKE::FatalErrorMessage("Species/Class not available: " + name);

		return -1;
	}

	bool LookUpTable_ZY::available_omega(const std::string name) const
	{
		const auto itr = std::find(list_fields_.begin(), list_fields_.end(), name);
		if (itr != list_fields_.cend())
		{
			const int index = std::distance(list_fields_.begin(), itr);
			return available_omega_[index];
		}
		else
		{
			return false;
		}
	}

	double LookUpTable_ZY::Y(const std::string name) const
	{
		auto itr = std::find(list_fields_.begin(), list_fields_.end(), name);
		if (itr != list_fields_.cend())
			return interpolated_Y_[std::distance(list_fields_.begin(), itr)];
		else return 0.;
	}

	double LookUpTable_ZY::OmegaP(const std::string name) const
	{
		auto itr = std::find(list_fields_.begin(), list_fields_.end(), name);
		if (itr != list_fields_.cend())
			return interpolated_omegap_[std::distance(list_fields_.begin(), itr)];
		else return 0.;
	}

	double LookUpTable_ZY::OmegaD(const std::string name) const
	{
		auto itr = std::find(list_fields_.begin(), list_fields_.end(), name);
		if (itr != list_fields_.cend())
			return interpolated_omegad_[std::distance(list_fields_.begin(), itr)];
		else return 0.;
	}

	double LookUpTable_ZY::ReconstructMixtureFraction(const double* Y)
	{
		// Mole fractions
		double mw;
		Eigen::VectorXd X(thermodynamicsMap_.NumberOfSpecies());
		thermodynamicsMap_.MoleFractions_From_MassFractions(X.data(), mw, Y);

		// Reconstruct mixture fractions
		const double ZC = X.dot(thermodynamicsMap_.atomic_composition().col(jC_)) * WC_ / mw;
		const double ZH = X.dot(thermodynamicsMap_.atomic_composition().col(jH_)) * WH_ / mw;
		const double ZO = X.dot(thermodynamicsMap_.atomic_composition().col(jO_)) * WO_ / mw;

		const double Zstar = 2. * ZC / WC_ + 0.50 * ZH / WH_ - ZO / WO_;

		const double Z = (Zstar - ZstarOx_) / (ZstarFuel_ - ZstarOx_);

		return Z;
	}

	double LookUpTable_ZY::ReconstructProgressVariable(const double* Y)
	{
		double C = 0.;

		// Regular species
		for (unsigned int i = 0; i < Y_def_indices_.size(); i++)
			C += Y_def_weights_[i] * Y[Y_def_indices_[i]];

		// PAH12 (if available)
		for (unsigned int i = 0; i < indices_pah12_.size(); i++)
			C += Y_def_alpha_PAH12_ * Y[indices_pah12_[i]] / thermodynamicsMap_.MWs()[indices_pah12_[i]];

		// PAH34 (if available)
		for (unsigned int i = 0; i < indices_pah34_.size(); i++)
			C += Y_def_alpha_PAH34_ * Y[indices_pah34_[i]] / thermodynamicsMap_.MWs()[indices_pah34_[i]];

		// PAHLP (if available)
		for (unsigned int i = 0; i < indices_pahlp_.size(); i++)
			C += Y_def_alpha_PAHLP_ * Y[indices_pahlp_[i]] / thermodynamicsMap_.MWs()[indices_pahlp_[i]];

		// SP (if available)
		for (unsigned int i = 0; i < indices_sp_.size(); i++)
			C += Y_def_alpha_SP_ * Y[indices_sp_[i]] / thermodynamicsMap_.MWs()[indices_sp_[i]];

		// AGG (if available)
		for (unsigned int i = 0; i < indices_agg_.size(); i++)
			C += Y_def_alpha_AGG_ * Y[indices_agg_[i]] / thermodynamicsMap_.MWs()[indices_agg_[i]];

		return C;
	}

}
