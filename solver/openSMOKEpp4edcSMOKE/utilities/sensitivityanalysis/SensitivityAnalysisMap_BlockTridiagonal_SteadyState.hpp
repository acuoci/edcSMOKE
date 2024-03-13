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

namespace OpenSMOKE
{
	SensitivityMap_BlockTridiagonal_SteadyState::SensitivityMap_BlockTridiagonal_SteadyState(KineticsMap_CHEMKIN& kinetics, const unsigned int number_of_equations, const unsigned int block_dimension) :
	kinetics_(kinetics),
	number_of_equations_(number_of_equations),
	block_dimension_(block_dimension)
	{
		// Number of blocks
		number_of_blocks_ = number_of_equations_ / block_dimension_;

		// Default values
		energy_type_ = CONSTANT_PRESSURE_SYSTEM;
		sensitivity_type_ = PhysicalConstants::SENSITIVITY_FREQUENCY_FACTOR;


		number_of_species_ = boost::lexical_cast<unsigned int>(kinetics.NamesOfSpecies().size());
		number_of_reactions_ = boost::lexical_cast<unsigned int>(kinetics.NumberOfReactions());

		index_of_species_ = 0;				// index of first species
		index_of_temperature_ = -1;			// negative value means that there is no temperature equation
		index_of_mass_flow_rate_ = -1;		// negative value means that there is no density equation

		if (sensitivity_type_ == PhysicalConstants::SENSITIVITY_KINETIC_CONSTANT)
			number_of_parameters_ = kinetics.NumberOfReactions();
		else if (sensitivity_type_ == PhysicalConstants::SENSITIVITY_FREQUENCY_FACTOR)
			number_of_parameters_ = kinetics.NumberOfReactions() + kinetics.NumberOfFallOffReactions() + kinetics.NumberOfCABRReactions();

		// Kinetic parameters
		ChangeDimensions(number_of_parameters_, &parameters_, true);
	
		// CPU time
		cpuTimeSingleFactorization_ = 0.;
		cpuTimeFactorization_ = 0.;
		cpuTimeSingleAssembling_ = 0.;
		cpuTimeAssembling_ = 0.;
		cpuTimeSingleSolution_ = 0.;
		cpuTimeSolution_ = 0.;
	}

	void SensitivityMap_BlockTridiagonal_SteadyState::SetSensitivityType(const PhysicalConstants::sensitivity_type type)
	{
		sensitivity_type_ = type;
		if (sensitivity_type_ == PhysicalConstants::SENSITIVITY_KINETIC_CONSTANT)
			number_of_parameters_ = kinetics_.NumberOfReactions();
		else if (sensitivity_type_ == PhysicalConstants::SENSITIVITY_FREQUENCY_FACTOR)
			number_of_parameters_ = kinetics_.NumberOfReactions() + kinetics_.NumberOfFallOffReactions() + kinetics_.NumberOfCABRReactions();

		ChangeDimensions(number_of_parameters_, &parameters_, true);
	}

	void SensitivityMap_BlockTridiagonal_SteadyState::SetIndexOfTemperature(const unsigned int index)
	{
		index_of_temperature_ = index;
	}

	void SensitivityMap_BlockTridiagonal_SteadyState::SetIndexOfMassFlowRate(const unsigned int index)
	{
		index_of_mass_flow_rate_ = index;
	}

	void SensitivityMap_BlockTridiagonal_SteadyState::SetIndicesOfSpecies(const Eigen::VectorXi& indices_species)
	{
		indices_species_ = indices_species;
	}

	void SensitivityMap_BlockTridiagonal_SteadyState::SetEnergyEquationType(const EnergyEquationType type)
	{
		energy_type_ = type;
	}

	void SensitivityMap_BlockTridiagonal_SteadyState::Calculate(const Eigen::VectorXd& T, const double P_Pa, const std::vector<Eigen::VectorXd>& X, OpenSMOKE::OpenSMOKEBandMatrixDouble& J, const std::vector<Eigen::VectorXd>& scaling_Jp)
	{
		Eigen::MatrixXd sensitivity_coeffs_;

		// Allocate memory: sensitivity_temperature
		sensitivity_temperature_.resize(number_of_blocks_, number_of_parameters_);
		sensitivity_temperature_.setZero();

		// Allocate memory: sensitivity_mass_flow_rate
		sensitivity_mass_flow_rate_.resize(number_of_blocks_, number_of_parameters_);	
		sensitivity_mass_flow_rate_.setZero();

		// Allocate memory: sensitivity_species
		if (indices_species_.size() != 0)
		{
			sensitivity_species_ = new Eigen::MatrixXd[indices_species_.size()];
			for (int i = 0; i < indices_species_.size(); i++)
			{
				sensitivity_species_[i].resize(number_of_blocks_, number_of_parameters_);
				sensitivity_species_[i].setZero();
			}
		}

		// Factorization of matrix J
		J.Factorize();

		// Calculation of block size
		double MB_MAX = 256.;															// max RAM available
		unsigned int NP_BLOCK = int(MB_MAX*1.e6 / (8.*number_of_blocks_*number_of_species_));	// maximum number of parameters per block
		unsigned int N_BLOCKS = number_of_parameters_ / NP_BLOCK;													// number of blocks

		for (unsigned int n = 1; n <= N_BLOCKS + 1; n++)
		{
			const int iStart = NP_BLOCK*(n - 1) + 1;
			const int iEnd = std::min(NP_BLOCK*n, number_of_parameters_);
			const int block_size = iEnd - iStart + 1;

			if (block_size > 0)
			{
				std::cout << "Sensitivity Block #" << n << "/" << N_BLOCKS + 1 << "(size: " << block_size << ")" << std::endl;

				sensitivity_coeffs_.resize(number_of_equations_, block_size);
				sensitivity_coeffs_.setZero();

				const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

				for (unsigned k = 2; k <= number_of_blocks_ - 1; k++)
				{
					OpenSMOKE::OpenSMOKEVectorDouble mole_fractions_local(number_of_species_);
					OpenSMOKE::OpenSMOKEVectorDouble c_local(number_of_species_);
					OpenSMOKE::OpenSMOKEVectorDouble R_local(number_of_species_);
					OpenSMOKE::OpenSMOKEVectorDouble Jp_species_(number_of_species_);
					OpenSMOKE::OpenSMOKEVectorDouble b_local(block_dimension_);

					const double cTot = P_Pa / T(k - 1) / PhysicalConstants::R_J_kmol;
					for (unsigned int i = 1; i <= number_of_species_; i++)
					{
						mole_fractions_local[i] = X[k - 1](i - 1);
						c_local[i] = cTot*mole_fractions_local[i];
					}

					// Set maps
					kinetics_.thermodynamics().SetTemperature(T(k - 1));
					kinetics_.thermodynamics().SetPressure(P_Pa);
					kinetics_.SetTemperature(T(k - 1));
					kinetics_.SetPressure(P_Pa);

					// Update reaction and formation rates
					kinetics_.ReactionRates(c_local.GetHandle());
					kinetics_.FormationRates(R_local.GetHandle());

					// Assembling rhs
					for (int j = iStart; j <= iEnd; j++)
					{
						double Jp_T_;
						if (index_of_temperature_ == -1)
							kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, j, c_local.GetHandle(), Jp_species_.GetHandle(), parameters_[j]);
						else
							kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, energy_type_, j, c_local.GetHandle(), mole_fractions_local.GetHandle(), Jp_species_.GetHandle(), Jp_T_, parameters_[j]);

						for (unsigned int i = 1; i <= number_of_species_; i++)
							b_local[i] = -Jp_species_[i] * scaling_Jp[k - 1](i - 1);

						if (index_of_temperature_ != -1)
							b_local[index_of_temperature_ + 1] = -Jp_T_*scaling_Jp[k - 1](index_of_temperature_);

						if (index_of_mass_flow_rate_ != -1)
							b_local[index_of_mass_flow_rate_ + 1] = 0.;

						const unsigned j_local = j - iStart + 1;
						for (unsigned int i = 1; i <= block_dimension_; i++)
							sensitivity_coeffs_( (k - 1)*block_dimension_ + i-1, j_local-1 ) = b_local[i];
					}
				}
				const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

				cpuTimeSingleAssembling_ = tend - tstart;
				cpuTimeAssembling_ += cpuTimeSingleAssembling_;


				// Solution
				{
					const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

					J.Solve(block_size, sensitivity_coeffs_.data());

					const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

					cpuTimeSingleSolution_ = tend - tstart;
					cpuTimeSolution_ += cpuTimeSingleSolution_;
				}

				// From global to local
				{
					// Temperature
					if (index_of_temperature_ != -1)
					{
						for (unsigned k = 1; k <= number_of_blocks_; k++)
						{
							int k_global = (k - 1)*block_dimension_ + index_of_temperature_ + 1;
							for (int j = iStart; j <= iEnd; j++)
								sensitivity_temperature_(k-1,j-1) = sensitivity_coeffs_(k_global-1, j - iStart + 1-1);
						}
					}

					// Mass flow rate
					if (index_of_mass_flow_rate_ != -1)
					{
						for (unsigned k = 1; k <= number_of_blocks_; k++)
						{
							int k_global = (k - 1)*block_dimension_ + index_of_mass_flow_rate_ + 1;
							for (int j = iStart; j <= iEnd; j++)
								sensitivity_mass_flow_rate_(k-1,j-1) = sensitivity_coeffs_(k_global-1, j - iStart + 1-1);
						}
					}

					for (int i = 0; i < indices_species_.size(); i++)
					{
						for (unsigned k = 1; k <= number_of_blocks_; k++)
						{
							int k_global = (k - 1)*block_dimension_ + indices_species_[i] + 1;
							for (int j = iStart; j <= iEnd; j++)
								sensitivity_species_[i](k-1,j-1) = sensitivity_coeffs_(k_global-1, j - iStart + 1-1);
						}
					}
					
				}
			} // if block size is greater than 0
		} // end cycle on blocks
	}

	void SensitivityMap_BlockTridiagonal_SteadyState::SaveOnXMLFile(const std::string folder_name)
	{
		unsigned int number_of_variables = indices_species_.size();

		if (index_of_temperature_ != -1)	number_of_variables++;
		if (index_of_mass_flow_rate_ != -1)	number_of_variables++;

		// Main file
		{
			boost::filesystem::path file_name = folder_name;
			file_name /= "Sensitivities.xml";
			std::ofstream fXML;
			fXML.open(file_name.c_str(), std::ios::out);
			fXML.setf(std::ios::scientific);

			fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			fXML << "<opensmoke version=\"0.1a\">" << std::endl;

			fXML << "<variables>" << std::endl;
			fXML << number_of_variables << std::endl;
		
			for (int i = 0; i < indices_species_.size(); i++)
				fXML << kinetics_.NamesOfSpecies()[indices_species_[i]] << " " << i << " " << indices_species_[i]+1 << std::endl;
			if (index_of_temperature_ != -1)
				fXML << "temperature" << " " << indices_species_.size() << " " << index_of_temperature_+1 << std::endl;
			if (index_of_mass_flow_rate_ != -1)
				fXML << "mass-flow-rate" << " " << indices_species_.size() + 1 << " " << index_of_mass_flow_rate_+1 << std::endl;
			fXML << "</variables>" << std::endl;

			fXML << "<n-parameters>" << std::endl;
			fXML << parameters_.Size() << std::endl;
			fXML << "</n-parameters>" << std::endl;
			fXML << "<points>" << std::endl;
			fXML << number_of_blocks_ << std::endl;
			fXML << "</points>" << std::endl;
			fXML << "<constant-parameters>" << std::endl;
			for (int i = 1; i <= parameters_.Size(); i++)
				fXML << parameters_[i] << std::endl;
			fXML << "</constant-parameters>" << std::endl;
			fXML << "</opensmoke>" << std::endl;

			fXML.close();
		}
	
		// Temperature
		if (index_of_temperature_ != -1)
		{
			boost::filesystem::path file_name = folder_name;
			file_name /= "Sensitivities.temperature.xml";
			std::ofstream fXML;
			fXML.open(file_name.c_str(), std::ios::out);
			fXML.setf(std::ios::scientific);

			fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			fXML << "<opensmoke version=\"0.1a\">" << std::endl;
			fXML << "<coefficients>" << std::endl;
			for (unsigned int j = 1; j <= number_of_blocks_; j++)
			{
				for (int i = 1; i <= parameters_.Size(); i++)
					fXML << std::left << std::setw(16) << sensitivity_temperature_(j-1, i-1);

				fXML << std::endl;
			}
			fXML << "</coefficients>" << std::endl;
			fXML << "</opensmoke>" << std::endl;

			fXML.close();
		}

		// Species
		for (int z = 0; z < indices_species_.size(); z++)
		{
			boost::filesystem::path file_name = folder_name;
			file_name /= ("Sensitivities." + kinetics_.NamesOfSpecies()[indices_species_[z]] + ".xml");
			std::ofstream fXML;
			fXML.open(file_name.c_str(), std::ios::out);
			fXML.setf(std::ios::scientific);

			fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			fXML << "<opensmoke version=\"0.1a\">" << std::endl;
			fXML << "<coefficients>" << std::endl;
			for (unsigned int j = 1; j <= number_of_blocks_; j++)
			{
				for (int i = 1; i <= parameters_.Size(); i++)
					fXML << std::left << std::setw(16) << sensitivity_species_[z](j-1,i-1);

				fXML << std::endl;
			}
			fXML << "</coefficients>" << std::endl;
			fXML << "</opensmoke>" << std::endl;

			fXML.close();
		}

		// Mass flow rate
		if (index_of_mass_flow_rate_ != -1)
		{
			boost::filesystem::path file_name = folder_name;
			file_name /= "Sensitivities.mass-flow-rate.xml";
			std::ofstream fXML;
			fXML.open(file_name.c_str(), std::ios::out);
			fXML.setf(std::ios::scientific);

			fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			fXML << "<opensmoke version=\"0.1a\">" << std::endl;
			fXML << "<coefficients>" << std::endl;
			for (unsigned int j = 1; j <= number_of_blocks_; j++)
			{
				for (int i = 1; i <= parameters_.Size(); i++)
					fXML << std::left << std::setw(16) << sensitivity_mass_flow_rate_(j-1,i-1);

				fXML << std::endl;
			}
			fXML << "</coefficients>" << std::endl;
			fXML << "</opensmoke>" << std::endl;

			fXML.close();
		}
	}
}





