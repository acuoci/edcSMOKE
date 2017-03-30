/*-----------------------------------------------------------------------*\
|                  _       _____ __  __  ____  _  ________                |
|                 | |     / ____|  \/  |/ __ \| |/ /  ____|               |
|          ___  __| | ___| (___ | \  / | |  | | ' /| |__                  |
|         / _ \/ _` |/ __|\___ \| |\/| | |  | |  < |  __|                 |
|        |  __/ (_| | (__ ____) | |  | | |__| | . \| |____                |
|         \___|\__,_|\___|_____/|_|  |_|\____/|_|\_\______|               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, M.R. Malik, Z. Li, A. Parente                      |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|   Contacts: Mohammad Rafi Malik, Zhiyi Li, Alessandro Parente           |
|   Aero-Thermo-Mechanical Department                                     |
|   Université Libre de Bruxelles                                         |
|   Avenue F. D. Roosevelt 50, 1050 Bruxelles (Belgium)                   |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of edcSMOKE solver.                                 |
|                                                                         |
|	License                                                           |
|                                                                         |
|   Copyright(C) 2017-2014 A. Cuoci, A. Parente                           |
|   edcSMOKE is free software: you can redistribute it and/or modify      |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   edcSMOKE is distributed in the hope that it will be useful,           |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with edcSMOKE. If not, see <http://www.gnu.org/licenses/>.      |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include <queue>

namespace OpenSMOKE
{
	DRG::DRG(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
			 OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML) :

			 thermodynamicsMapXML_(*thermodynamicsMapXML),
			 kineticsMapXML_(*kineticsMapXML)
	{
		epsilon_ = 1.e-2;

		NS_ = thermodynamicsMapXML_.NumberOfSpecies();
		NR_ = kineticsMapXML_.NumberOfReactions();

		ChangeDimensions(NR_, &rNet_, true);
		r_.resize(NS_, NS_);
		

		important_species_.resize(NS_);
                important_reactions_.resize(NR_);

		number_important_species_ = NS_;
		number_unimportant_species_ = 0;
		number_unimportant_reactions_ = 0;

		ResetCounters();

		// Build a full matrix of net stoichiometric coefficients nu = nuB - nuF
		Eigen::MatrixXd nu_(NR_, NS_);						
		{
			// Be careful: eigen vectors and matrices are 0-index based
			// Be careful: if the kinetic scheme is large, this matrix, since it is full, can be very memory expensive
			//             Example: 10^3 species, 10^4 reactions = size of the matrix 10^7 elements!
			//             This is the reason why we store stoichiometric matrices in sparse format.
			//             Of course te advantage of having a full matrix, is that you access the elements directly, without
			//             using iterators and pointers, as reported above
			nu_.setZero();

			// Loop over all the reactions (product side)
			for (int k = 0; k < kineticsMapXML_.stoichiometry().stoichiometric_matrix_products().outerSize(); ++k)
			{
				// Loop over all the non-zero stoichiometric coefficients (product side) of reaction k
				for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML_.stoichiometry().stoichiometric_matrix_products(), k); it; ++it)
				{
					nu_(it.row(), it.col()) += it.value();
				}
			}

			// Loop over all the reactions (product side)
			for (int k = 0; k < kineticsMapXML_.stoichiometry().stoichiometric_matrix_reactants().outerSize(); ++k)
			{
				// Loop over all the non-zero stoichiometric coefficients (product side) of reaction k
				for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML_.stoichiometry().stoichiometric_matrix_reactants(), k); it; ++it)
				{
					nu_(it.row(), it.col()) -= it.value();
				}
			}
		}

		// Build the delta matrix (dense matrix) used by the DRG method
		Eigen::MatrixXd delta_(NR_, NS_);					
		{
			for (unsigned int i = 0; i < NR_; i++)
			{
				for (unsigned int j = 0; j < NS_; j++)
					delta_(i, j) = (nu_(i, j) == 0) ? 0 : 1;
			}
		}

		// Sparse Algebra
		{			
			// Sparse net stoichiometric matrix
			{
				nu_sparse_.resize(NS_,NR_);

				typedef Eigen::Triplet<double> T;
				std::vector<T> tripletList;
				tripletList.reserve(NR_*4);
				for (unsigned int i = 0; i < NR_; i++)
					for (unsigned int j = 0; j < NS_; j++)
					{
  						if ( nu_(i,j) != 0.)
  						tripletList.push_back(T(j,i,nu_(i,j)));
					}
		
				nu_sparse_.setFromTriplets(tripletList.begin(), tripletList.end());
			}

			// Sparse delta matrix, 1 means the species is involved in the reaction, (NR x NS)
			delta_sparse_.resize(NS_,NR_);
			{
				typedef Eigen::Triplet<double> T;
				std::vector<T> tripletList;
				tripletList.reserve(NR_*4);
				for (unsigned int i = 0; i < NR_; i++)
					for (unsigned int j = 0; j < NS_; j++)
					{
  						if ( delta_(i,j) != 0.)
  						tripletList.push_back(T(j,i,1.));
					}
		
				delta_sparse_.setFromTriplets(tripletList.begin(), tripletList.end());
			}

			// Nu times delta
			{
				nu_times_delta_.resize(NS_);

				for (unsigned int k = 0; k < NS_; k++)
				{
					nu_times_delta_[k].resize(NS_, NR_);

					typedef Eigen::Triplet<double> T;
					std::vector<T> tripletList;
					tripletList.reserve(NR_*4);
					for (unsigned int i = 0; i < NR_; i++)
						for (unsigned int j = 0; j < NS_; j++)
						{
							const double prod = nu_(i,k) * delta_(i,j);
  							if ( prod != 0.)
  								tripletList.push_back(T(j,i, prod));
						}
		
					nu_times_delta_[k].setFromTriplets(tripletList.begin(), tripletList.end());
				}
			}
		}
	}

	void DRG::ResetCounters()
	{
		counter_analysis_ = 0;

		cpuTimeCumulative_PairWiseErrorMatrix_ = 0.;
		//cpuTimeCumulative_PairWiseErrorMatrix_Step1_ = 0.;
		cpuTimeCumulative_PairWiseErrorMatrix_Step2_ = 0.;
		cpuTimeCumulative_ParsePairWiseErrorMatrix_ = 0.;

		number_important_species_cumulative_ = 0;
		number_important_reactions_cumulative_ = 0;
		number_important_species_min_ = 1e6;
		number_important_species_max_ = 0;
		number_important_reactions_min_ = 1e6;
		number_important_reactions_max_ = 0;
	}

	void DRG::SetKeySpecies(const std::vector<std::string> names_key_species)
	{
		index_key_species_.resize(names_key_species.size());
		for (unsigned int i = 0; i < names_key_species.size(); i++)
			index_key_species_[i] = thermodynamicsMapXML_.IndexOfSpecies(names_key_species[i]) - 1;
	}

	void DRG::SetKeySpecies(const std::vector<unsigned int> key_species)
	{
		index_key_species_ = key_species;
	}

	void DRG::SetEpsilon(const double epsilon)
	{
		epsilon_ = epsilon;
	}

	void DRG::Analysis(const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c)
	{
		// Updating counter
		counter_analysis_++;

		// Construction pair wise matrix
		{
			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

			PairWiseErrorMatrix(T, P_Pa, c);

			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			const double cpuTimeLocal_PairWiseErrorMatrix = tEnd-tStart;
			cpuTimeCumulative_PairWiseErrorMatrix_ += cpuTimeLocal_PairWiseErrorMatrix;
			cpuTimeAverage_PairWiseErrorMatrix_ = cpuTimeCumulative_PairWiseErrorMatrix_/static_cast<double>(counter_analysis_);
		}

		// Parsing pair wise matrix
		{
			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

			ParsePairWiseErrorMatrix();

			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			const double cpuTimeLocal_ParsePairWiseErrorMatrix = tEnd-tStart;
			cpuTimeCumulative_ParsePairWiseErrorMatrix_ += cpuTimeLocal_ParsePairWiseErrorMatrix;
			cpuTimeAverage_ParsePairWiseErrorMatrix_ = cpuTimeCumulative_ParsePairWiseErrorMatrix_/static_cast<double>(counter_analysis_);
		}
	}

	void DRG::PairWiseErrorMatrix(const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c)
	{
		// Now we know T, P and composition. 
		// We have to pass those data to the thermodynamic and kinetic maps
		kineticsMapXML_.SetTemperature(T);
		kineticsMapXML_.SetPressure(P_Pa);
		thermodynamicsMapXML_.SetTemperature(T);
		thermodynamicsMapXML_.SetPressure(P_Pa);

		// Now we can calculate (internally) the reaction rates concentrations are needed
		kineticsMapXML_.ReactionRates(c.GetHandle());
		kineticsMapXML_.GiveMeReactionRates(rNet_.GetHandle());	// [kmol/m3/s]

		// Calculate the pair-wise error matrix
		{
			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

			{
				r_.setZero();

				Eigen::VectorXd numerator_(NS_);
				Eigen::VectorXd denominator_(NS_);

				// Denominator			
				denominator_.setConstant(0.);		
				for (int k=0; k<nu_sparse_.outerSize(); ++k)
	  			{
					const double rnet = rNet_[k+1];
					for (Eigen::SparseMatrix<double>::InnerIterator it(nu_sparse_,k); it; ++it)
	  				{
	    					denominator_[it.row()] += std::fabs(it.value() * rnet);
	  				}
				}

				// Numerator
				for (int i = 0; i < NS_; i++)
				{
					numerator_.setConstant(0.);
					for (int k=0; k<nu_times_delta_[i].outerSize(); ++k)
	  				{
						const double rnet = rNet_[k+1];
						for (Eigen::SparseMatrix<double>::InnerIterator it(nu_times_delta_[i],k); it; ++it)
	  					{
	    						numerator_[it.row()] += std::fabs(it.value() * rnet);
	  					}
					}

					for (int j = 0; j < NS_; j++)
						r_(i,j) = numerator_[j]/(1.e-64+denominator_[i]);
				}
			}

			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			const double cpuTimeLocal_PairWiseErrorMatrix_Step2_ = tEnd-tStart;
			cpuTimeCumulative_PairWiseErrorMatrix_Step2_ += cpuTimeLocal_PairWiseErrorMatrix_Step2_;
			cpuTimeAverage_PairWiseErrorMatrix_Step2_ = cpuTimeCumulative_PairWiseErrorMatrix_Step2_/static_cast<double>(counter_analysis_);
		}
        }
        
	void DRG::ParsePairWiseErrorMatrix()
	{
		// Reset important species and important reactions
		important_species_.assign(NS_,false);
		important_reactions_.assign(NR_,true);
     
		// Initialize the queue with key-species
		std::queue <int> Q;
		for ( int i=0; i<index_key_species_.size(); i++)
		{
			Q.push(index_key_species_[i]);
			important_species_[index_key_species_[i]]=true;
		}
           
		// DFS with rAB
		while (!Q.empty())
		{
			for ( int k=0; k<NS_; k++)
			{
				if (important_species_[k] == false)                    
				{
					if (r_(Q.front(), k) > epsilon_)
					{                            
						important_species_[k] = true;
						Q.push(k);
					}
				}
			}
			Q.pop();
		}	
                
		// Important reactions
		for (int k=0; k<delta_sparse_.outerSize(); ++k)
  		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(delta_sparse_,k); it; ++it)
  			{
				if (important_species_[it.row()] == false)
				{
					important_reactions_[k] = false;
				}
  			}
		}

		// Count important species and reactions
		number_important_species_ = std::count (important_species_.begin(), important_species_.end(), true);
		number_unimportant_species_ = NS_ - number_important_species_;
		number_unimportant_reactions_ = std::count (important_reactions_.begin(), important_reactions_.end(), false);

		// Vector containing the indices of important species (zero-based)
		{
			indices_important_species_.resize(number_important_species_);
			indices_unimportant_species_.resize(number_unimportant_species_);

			unsigned int count_important = 0;
			unsigned int count_unimportant = 0;
			for(unsigned int k = 0; k < NS_; k++)
			{
				if (important_species_[k] == true)
				{
					indices_important_species_[count_important] = k;
					count_important++;
				}
				else
				{
					indices_unimportant_species_[count_unimportant] = k;
					count_unimportant++;
				}
			}
                }

		// Vector containing the indices of unimportant species (zero-based)
		{
			indices_unimportant_reactions_.resize(number_unimportant_reactions_);
			unsigned int count = 0;
			for(unsigned int k = 0; k < NR_; k++)
				if (important_reactions_[k] == false)
				{
					indices_unimportant_reactions_[count] = k;
					count++;
				}
                }
		
		number_important_species_cumulative_ += number_important_species_;
		number_important_reactions_cumulative_ += (NR_-number_unimportant_reactions_);
		number_important_species_average_ = number_important_species_cumulative_/static_cast<double>(counter_analysis_);
		number_important_reactions_average_ = number_important_reactions_cumulative_/static_cast<double>(counter_analysis_);

		number_important_species_max_ = (number_important_species_ > number_important_species_max_) ? number_important_species_ : number_important_species_max_;
		number_important_species_min_ = (number_important_species_ < number_important_species_min_) ? number_important_species_ : number_important_species_min_;

		number_important_reactions_max_ = ((NR_-number_unimportant_reactions_) > number_important_reactions_max_) ? (NR_-number_unimportant_reactions_) : number_important_reactions_max_;
		number_important_reactions_min_ = ((NR_-number_unimportant_reactions_) < number_important_reactions_min_) ? (NR_-number_unimportant_reactions_) : number_important_reactions_min_;
	}

	void DRG::WriteCpuTimesOnTheScreen()
	{
		std::cout << "DRG Analysis" << std::endl;
		std::cout << " * Pair-Wise Error Matrix (overall):         " << cpuTimeAverage_PairWiseErrorMatrix_ << std::endl;
		std::cout << " * Pair-Wise Error Matrix (construction):    " << cpuTimeAverage_PairWiseErrorMatrix_Step2_ << std::endl;
		std::cout << " * Pair-Wise Error Matrix (parsing):         " << cpuTimeAverage_ParsePairWiseErrorMatrix_ << std::endl;
		std::cout << " * Number important species (avg : min : max : tot):   " << number_important_species_average_ << " : " << number_important_species_min_ << " : " << number_important_species_max_ << " : " << NS_ << std::endl;
		std::cout << " * Number important reactions (avg : min : max : tot): " << number_important_reactions_average_ << " : " << number_important_reactions_min_ << " : " << number_important_reactions_max_ << " : " << NR_ << std::endl;			
	}
}
