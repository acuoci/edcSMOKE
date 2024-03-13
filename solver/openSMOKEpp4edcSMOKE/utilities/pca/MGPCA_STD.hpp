/*----------------------------------------------------------------------*\
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

#include "math/OpenSMOKEUtilities.h"

namespace OpenSMOKE
{	
	template< typename ThermoMap >
	MGPCA_STD<ThermoMap>::MGPCA_STD(const ThermoMap& thermo, const boost::filesystem::path& folder_name, const MGPCA_Type type)
	{
		type_ = type;

		// Removed species
		{
			std::cout << "MGPCA (STD): Reading names of removed species" << std::endl;
			boost::filesystem::path file_name = folder_name / "removed_species";
			std::ifstream fRemoved;
			fRemoved.open(file_name.string().c_str(), std::ios::in);
			
			for(;;)
			{
				std::string dummy;
				fRemoved >> dummy;
				if (dummy == "")
					break;

				unsigned int index = thermo.IndexOfSpecies(dummy);
				indices_of_removed_species_.push_back(index-1);

				if (fRemoved.eof())
					break;
			}
			fRemoved.close();

			number_of_removed_species_ = indices_of_removed_species_.size();
		}

		// Retained and unretained species
		{
			std::cout << "MGPCA (STD): Reading names of retained species" << std::endl;
			boost::filesystem::path file_name = folder_name / "ret_names_Flamelets.out";
			std::ifstream fRetainedNames;
			fRetainedNames.open(file_name.string().c_str(), std::ios::in);
			
			for(;;)
			{
				std::string dummy;
				fRetainedNames >> dummy;
				if (dummy == "")
					break;

				unsigned int index = thermo.IndexOfSpecies(dummy);
				indices_of_retained_species_.push_back(index-1);

				if (fRetainedNames.eof())
					break;
			}
			fRetainedNames.close();

			number_of_species_ = thermo.NumberOfSpecies();
			number_of_retained_species_ = indices_of_retained_species_.size();
			number_of_unretained_species_ = (number_of_species_-number_of_removed_species_) - number_of_retained_species_;

			for(unsigned int i=0;i<thermo.NumberOfSpecies();i++)
			{
				bool unretained_species = true;
				for(unsigned int j=0;j<number_of_retained_species_;j++)
					if (i==indices_of_retained_species_[j])
					{
						unretained_species = false;
						break;
					}

				for(unsigned int j=0;j<number_of_removed_species_;j++)
					if (i==indices_of_removed_species_[j])
					{
						unretained_species = false;
						break;
					}

				if (unretained_species == true)
					indices_of_unretained_species_.push_back(i);
			}

			// Tracked species
			{
				number_of_tracked_species_ = number_of_retained_species_+number_of_unretained_species_;
				indices_of_tracked_species_.resize(number_of_tracked_species_);
				unsigned int count = 0;
				for(unsigned int i=0;i<number_of_retained_species_;i++)
					indices_of_tracked_species_[count++] = indices_of_retained_species_[i];
				for(unsigned int i=0;i<number_of_unretained_species_;i++)
					indices_of_tracked_species_[count++] = indices_of_unretained_species_[i];
				std::sort(indices_of_tracked_species_.begin(), indices_of_tracked_species_.end());

				if (type_ == MGPCA_Standard)
				{
					number_of_mixture_fraction_points_ = 1;
				}
				else if (type_ == MGPCA_Standard_LocalScaling)
				{
					boost::filesystem::path		file_name_averages = folder_name / "X_cent";
					std::ifstream fInput;
					fInput.open(file_name_averages.string().c_str(), std::ios::in);
			
					unsigned count = 0;
					for(;;)
					{
						std::string dummy;
						fInput >> dummy;
						if (dummy == "")
							break;

						count++;

						if (fInput.eof())
							break;
					}
					fInput.close();

					number_of_mixture_fraction_points_ = count / number_of_species_;

					unsigned int check = count % number_of_species_;
					if (check != 0)
						OpenSMOKE::FatalErrorMessage("The X_cent file is not consistent with the other PCA files!");
					
					std::cout << "Number of entry points:            " << count << std::endl;
					std::cout << "Number of mixture fraction points: " << number_of_mixture_fraction_points_ << std::endl;
					
				}

				average_retained_.resize(number_of_mixture_fraction_points_);
				average_unretained_.resize(number_of_mixture_fraction_points_);
				gamma_retained_.resize(number_of_retained_species_);
				gamma_unretained_.resize(number_of_retained_species_);

				current_average_retained_.resize(number_of_retained_species_);
				current_average_unretained_.resize(number_of_unretained_species_);

				for (unsigned int i=0; i<number_of_mixture_fraction_points_; i++)
				{
					average_retained_[i].resize(number_of_retained_species_);
					average_unretained_[i].resize(number_of_unretained_species_);
				}
			}

			x_retained_.resize(number_of_retained_species_);
			x_unretained_.resize(number_of_unretained_species_);
		}

		{
			std::cout << "MGPCA: Scaling Factors" << std::endl;
			
			boost::filesystem::path file_name_gamma = folder_name / "gamma.out";
			std::ifstream fGamma;
			fGamma.open(file_name_gamma.string().c_str(), std::ios::in);

			for(unsigned int i=0;i<number_of_tracked_species_;i++)
			{
				double dummy_gamma;
				fGamma >> dummy_gamma;
				
				unsigned int index = indices_of_tracked_species_[i];

				for(unsigned int j=0;j<number_of_retained_species_;j++)
					if ( index == indices_of_retained_species_[j])
						gamma_retained_[j] = dummy_gamma;
				
				for(unsigned int j=0;j<number_of_unretained_species_;j++)
					if (index == indices_of_unretained_species_[j])
						gamma_unretained_[j] = dummy_gamma;
			}
		}

		if (type_ == MGPCA_Standard)
		{
			std::cout << "MGPCA: Centering Factors" << std::endl;
			
			boost::filesystem::path file_name_averages = folder_name / "X_ave";
			std::ifstream fAverages;
			fAverages.open(file_name_averages.string().c_str(), std::ios::in);

			for(unsigned int i=0;i<number_of_tracked_species_;i++)
			{
				double dummy_averages;
				fAverages >> dummy_averages;
				
				unsigned int index = indices_of_tracked_species_[i];

				for(unsigned int j=0;j<number_of_retained_species_;j++)
					if ( index == indices_of_retained_species_[j])
						average_retained_[0][j] = dummy_averages;
				
				for(unsigned int j=0;j<number_of_unretained_species_;j++)
					if (index == indices_of_unretained_species_[j])
						average_unretained_[0][j] = dummy_averages;
			}
		}
		else if (type_ == MGPCA_Standard_LocalScaling)
		{
			std::cout << "MGPCA: Centering Factors" << std::endl;

			std::ifstream fAverages;
			boost::filesystem::path file_name_averages = folder_name / "X_cent";

			fAverages.open(file_name_averages.string().c_str(), std::ios::in);
		
			for(unsigned int k=0;k<number_of_mixture_fraction_points_;k++)
			{
				for(unsigned int i=0;i<number_of_species_;i++)
				{
					double dummy_average;
					fAverages >> dummy_average;
				
					unsigned int index = i;

					for(unsigned int j=0;j<number_of_retained_species_;j++)
						if ( index == indices_of_retained_species_[j])
							average_retained_[k][j] = dummy_average;
				
					for(unsigned int j=0;j<number_of_unretained_species_;j++)
						if (index == indices_of_unretained_species_[j])
							average_unretained_[k][j] = dummy_average;
				}
			}
		}

		// Projection matrix
		{
			std::cout << "MGPCA (STD): Reading projection matrix" << std::endl;
			boost::filesystem::path file_name = folder_name / "B_proj";
			std::ifstream fB;
			fB.open(file_name.string().c_str(), std::ios::in);

			B_.resize(number_of_retained_species_, number_of_unretained_species_);
			for(unsigned int i=0;i<number_of_retained_species_;i++)
				for(unsigned int j=0;j<number_of_unretained_species_;j++)
					fB >> B_(i,j);
		}
		
		// Summary
		{
			std::cout << "Removed species (with 0-based index)" << std::endl;
			for(unsigned int i=0;i<number_of_removed_species_;i++)
				std::cout << thermo.NamesOfSpecies()[indices_of_removed_species_[i]] << "\t" << indices_of_removed_species_[i] << std::endl;

			std::cout << "Retained species (with 0-based index)" << std::endl;
			for(unsigned int i=0;i<number_of_retained_species_;i++)
				std::cout	<< thermo.NamesOfSpecies()[indices_of_retained_species_[i]] << "\t" 
							<< indices_of_retained_species_[i]							<< "\t"
							<< average_retained_[0][i]										<< "\t"
							<< gamma_retained_[i]										<< "\t"
							<< std::endl;

			std::cout << "Un-Retained species (with 0-based index)" << std::endl;
			for(unsigned int i=0;i<number_of_unretained_species_;i++)
				std::cout	<< thermo.NamesOfSpecies()[indices_of_unretained_species_[i]]	<< "\t" 
							<< indices_of_unretained_species_[i] 							<< "\t" 
							<< average_unretained_[0][i]										<< "\t"
							<< gamma_unretained_[i]											<< "\t"
							<< std::endl;
		}
	}

	template< typename ThermoMap >
	void MGPCA_STD<ThermoMap>::SetNormalization(const bool normalization)
	{
		normalization_ = normalization;
	}

	template< typename ThermoMap >
	void MGPCA_STD<ThermoMap>::Reconstruct(Eigen::VectorXd& x)
	{
		if (normalization_ == true)
		{
			// Centering and scaling
			double sum_retained = 0.;
			for(unsigned int i=0;i<number_of_retained_species_;i++)
			{
				sum_retained += x(indices_of_retained_species_[i]);
				x_retained_(i) = ( x(indices_of_retained_species_[i]) - average_retained_[0](i) ) / gamma_retained_(i);
			}

			if (sum_retained < 1.)
			{
				// Reconstruction
				for(unsigned int j=0;j<number_of_unretained_species_;j++)
				{
					x_unretained_(j) = 0.;
					for(unsigned int i=0;i<number_of_retained_species_;i++)
						x_unretained_(j) += x_retained_(i)*B_(i,j);
				}
		
				// Centering and scaling
				double sum_unretained_uncorrected = 0.;
				for(unsigned int i=0;i<number_of_unretained_species_;i++)
				{
					x(indices_of_unretained_species_[i]) = std::max(0.,x_unretained_(i)*gamma_unretained_(i) + average_unretained_[0](i));
					sum_unretained_uncorrected+=x(indices_of_unretained_species_[i]);
				}

				// Normalization
				if (sum_unretained_uncorrected>0.)
				{
					const double ratio = (1.-sum_retained)/sum_unretained_uncorrected;
					for(unsigned int i=0;i<number_of_unretained_species_;i++)
						x(indices_of_unretained_species_[i]) *= ratio;
				}
				
			//	double ratio = 1./(sum_retained+sum_unretained_uncorrected);
			//	for(unsigned int i=0;i<number_of_species_;i++)
			//		x(i) *= ratio;
			}
			else
			{
				for(unsigned int i=0;i<number_of_unretained_species_;i++)
					x(indices_of_unretained_species_[i]) = 0.;
			}

		}
		else
		{
			// Centering and scaling
			for(unsigned int i=0;i<number_of_retained_species_;i++)
				x_retained_(i) = ( x(indices_of_retained_species_[i]) - average_retained_[0](i) ) / gamma_retained_(i);

			// Reconstruction
			for(unsigned int j=0;j<number_of_unretained_species_;j++)
			{
				x_unretained_(j) = 0.;
				for(unsigned int i=0;i<number_of_retained_species_;i++)
					x_unretained_(j) += x_retained_(i)*B_(i,j);
			}
		
			// Centering and scaling
			for(unsigned int i=0;i<number_of_unretained_species_;i++)
				x(indices_of_unretained_species_[i]) = std::max(0.,x_unretained_(i)*gamma_unretained_(i) + average_unretained_[0](i));
		}
	}

	template< typename ThermoMap >
	void MGPCA_STD<ThermoMap>::Reconstruct(const double z, Eigen::VectorXd& x)
	{
		if (z<=0.)
		{
			current_average_retained_ = average_retained_[0];
			current_average_unretained_ = average_unretained_[0];
		}
		else if (z>=1.)
		{
			current_average_retained_ = average_retained_[number_of_mixture_fraction_points_-1];
			current_average_unretained_ = average_unretained_[number_of_mixture_fraction_points_-1];
		}
		else
		{
			const double delta = 1./double(number_of_mixture_fraction_points_-1);
			unsigned int i = (unsigned int)(z/delta);
			const double ratio = (z-i*delta)/delta;

			for(unsigned int j=0;j<number_of_retained_species_;j++)
				current_average_retained_(j) = average_retained_[i](j) + ratio*(average_retained_[i+1](j)-average_retained_[i](j));
			for(unsigned int j=0;j<number_of_unretained_species_;j++)
				current_average_unretained_(j) = average_unretained_[i](j) + ratio*(average_unretained_[i+1](j)-average_unretained_[i](j));

			
		}

		if (normalization_ == true)
		{
			// Centering and scaling
			double sum_retained = 0.;
			for(unsigned int i=0;i<number_of_retained_species_;i++)
			{
				sum_retained += x(indices_of_retained_species_[i]);
				x_retained_(i) = ( x(indices_of_retained_species_[i]) - current_average_retained_(i) ) / gamma_retained_(i);
			}

			if (sum_retained < 1.)
			{
				// Reconstruction
				for(unsigned int j=0;j<number_of_unretained_species_;j++)
				{
					x_unretained_(j) = 0.;
					for(unsigned int i=0;i<number_of_retained_species_;i++)
						x_unretained_(j) += x_retained_(i)*B_(i,j);
				}
		
				// Centering and scaling
				double sum_unretained_uncorrected = 0.;
				for(unsigned int i=0;i<number_of_unretained_species_;i++)
				{
					x(indices_of_unretained_species_[i]) = std::max(0.,x_unretained_(i)*gamma_unretained_(i) + current_average_unretained_(i));
					sum_unretained_uncorrected+=x(indices_of_unretained_species_[i]);
				}

				// Normalization
				if (sum_unretained_uncorrected>0.)
				{
					const double ratio = (1.-sum_retained)/sum_unretained_uncorrected;
					for(unsigned int i=0;i<number_of_unretained_species_;i++)
						x(indices_of_unretained_species_[i]) *= ratio;
				}
			}
			else
			{
				for(unsigned int i=0;i<number_of_unretained_species_;i++)
					x(indices_of_unretained_species_[i]) = 0.;
			}

		}
		else
		{
			// Centering and scaling
			for(unsigned int i=0;i<number_of_retained_species_;i++)
				x_retained_(i) = ( x(indices_of_retained_species_[i]) - current_average_retained_(i) ) / gamma_retained_(i);

			// Reconstruction
			for(unsigned int j=0;j<number_of_unretained_species_;j++)
			{
				x_unretained_(j) = 0.;
				for(unsigned int i=0;i<number_of_retained_species_;i++)
					x_unretained_(j) += x_retained_(i)*B_(i,j);
			}
		
			// Centering and scaling
			for(unsigned int i=0;i<number_of_unretained_species_;i++)
				x(indices_of_unretained_species_[i]) = std::max(0.,x_unretained_(i)*gamma_unretained_(i) + current_average_unretained_(i));
		}
	}
}
