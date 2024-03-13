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

#include "Grammar_OnTheFlyCEMA.h"

#if OPENSMOKE_USE_MKL == 1
#include "mkl.h"
#include "mkl_lapacke.h"
#elif OPENSMOKE_USE_OPENBLAS == 1
#include "cblas.h"
#include "lapacke.h"
#endif

namespace OpenSMOKE
{
	template <typename T> int sgn(T val) 
	{
		return (T(0) < val) - (val < T(0));
	}

	double cem_for_plots(const double lambda_exp_real)
	{
		return sgn(lambda_exp_real)*std::log10(1. + std::fabs(lambda_exp_real));
	}

	OnTheFlyCEMA::OnTheFlyCEMA(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								const boost::filesystem::path path_output) :

		thermodynamicsMap_(thermodynamicsMap),
		kineticsMap_(kineticsMap)
	{
		is_active_ = false;
		is_print_on_file_active_ = false;
		print_ei_ = false;
		print_pi_ = false;

		path_output_ = path_output;
		N_add_conservative_modes_ = 0;
		lapack_mode_ = true;
	}

	void OnTheFlyCEMA::Setup(const bool lapack_mode, const int additional_conservative_modes, const std::vector<std::string>& output_species, const std::vector<std::string> output_reactions)
	{
		is_active_ = true;
		lapack_mode_ = lapack_mode;
		N_add_conservative_modes_ = additional_conservative_modes;

		// Set required species
		{
			if (output_species[0] == "ALL" || output_species[0] == "all")
			{
				print_ei_ = true;
			}
			else if (output_species[0] == "NONE" || output_species[0] == "none")
			{
				print_ei_ = false;
			}
			else
			{
				print_ei_ = true;
				indices_of_output_species_.resize(output_species.size());
				for (unsigned int i = 0; i < output_species.size(); i++)
					indices_of_output_species_[i] = thermodynamicsMap_.IndexOfSpecies(output_species[i]);
			}
		}

		// Set required reactions
		{
			if (output_reactions[0] == "ALL" || output_reactions[0] == "all")
			{
				print_pi_ = true;
			}
			else if (output_reactions[0] == "NONE" || output_reactions[0] == "none")
			{
				print_pi_ = false;
			}
			else
			{
				print_pi_ = true;
				indices_of_output_reactions_.resize(output_reactions.size());
				for (unsigned int i = 0; i < output_reactions.size(); i++)
					indices_of_output_reactions_[i] = boost::lexical_cast<unsigned int>(output_reactions[i]);
			}
		}

		// Print on file
		if (print_ei_ == true || print_pi_ == true)
			is_print_on_file_active_ = true;

		MemoryAllocation();
	}

	void OnTheFlyCEMA::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		is_active_ = true;
		is_print_on_file_active_ = true;
		print_ei_ = true;
		print_pi_ = true;

		Grammar_OnTheFlyCEMA grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@Print") == true)
			dictionary.ReadBool("@Print", is_print_on_file_active_);

		if (dictionary.CheckOption("@LapackMode") == true)
			dictionary.ReadBool("@LapackMode", lapack_mode_);

		if (dictionary.CheckOption("@Species") == true)
		{
			std::vector<std::string> output_species;
			dictionary.ReadOption("@Species", output_species);

			if (output_species[0] == "ALL" || output_species[0] == "all")
			{
				print_ei_ = true;
			}
			else if (output_species[0] == "NONE" || output_species[0] == "none")
			{
				print_ei_ = false;
			}
			else
			{
				print_ei_ = true;
				indices_of_output_species_.resize(output_species.size());
				for(unsigned int i=0;i<output_species.size();i++)
					indices_of_output_species_ [i] = thermodynamicsMap_.IndexOfSpecies(output_species[i]);
			}
		}

		if (dictionary.CheckOption("@Reactions") == true)
		{
			std::vector<std::string> output_reactions;
			dictionary.ReadOption("@Reactions", output_reactions);

			if (output_reactions[0] == "ALL" || output_reactions[0] == "all")
			{
				print_pi_ = true;
			}
			else if (output_reactions[0] == "NONE" || output_reactions[0] == "none")
			{
				print_pi_ = false;
			}
			else
			{
				print_pi_ = true;
				indices_of_output_reactions_.resize(output_reactions.size());
				for (unsigned int i = 0; i<output_reactions.size(); i++)
					indices_of_output_reactions_[i] = boost::lexical_cast<unsigned int>(output_reactions[i]);
			}
		}

		if (dictionary.CheckOption("@AdditionalConservativeModes") == true)
			dictionary.ReadInt("@AdditionalConservativeModes", N_add_conservative_modes_);

		MemoryAllocation();
	}

	void OnTheFlyCEMA::MemoryAllocation()
	{
		// Basic variables
		NE_ = static_cast<unsigned int>(thermodynamicsMap_.elements().size());
		NS_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();
		N_ = (NS_ + 1) - (NE_ + 1 + N_add_conservative_modes_);
		Nstar_ = NS_ + 1;

		// Jacobian matrices
		Jomega_.resize(Nstar_, Nstar_);
		if (lapack_mode_ == false)
			JomegaTranspose_.resize(Nstar_, Nstar_);

		// Right and left eigenvectors associated with CEM
		ae_.resize(Nstar_);
		be_.resize(Nstar_);

		// Explosive and participation indices
		EI_.resize(Nstar_);
		PI_.resize(NR_);

		// Summary
		std::cout << std::endl;
		std::cout << "-------------------------------------------------------------" << std::endl;
		std::cout << "              CHEMICAL EXPLOSIVE MODE ANALYSIS               " << std::endl;
		std::cout << "-------------------------------------------------------------" << std::endl;
		std::cout << " * Number of elements:					  " << NE_ << std::endl;
		std::cout << " * Number of additional conservative modes: " << N_add_conservative_modes_ << std::endl;
		std::cout << " * Number of species:                       " << NS_ << std::endl;
		std::cout << " * Number of reactions:                     " << NR_ << std::endl;
		std::cout << " * Number of effective modes:               " << N_ << std::endl;
		std::cout << "-------------------------------------------------------------" << std::endl;
		std::cout << std::endl;
	}

	void OnTheFlyCEMA::CloseOutputFiles()
	{
		fCEMA_.close();
	}

	void OnTheFlyCEMA::PrepareOutputFiles()
	{
		if (is_print_on_file_active_ == true)
		{
			if (indices_of_output_species_.size() != 0)
			{
				widths_of_output_species_.resize(indices_of_output_species_.size());
				for (unsigned int i = 0; i < indices_of_output_species_.size(); i++)
					widths_of_output_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap_.NamesOfSpecies()[indices_of_output_species_[i] - 1], NS_);
			}
			else
			{
				widths_of_output_species_.resize(NS_);
				for (unsigned int i = 0; i < NS_; i++)
					widths_of_output_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap_.NamesOfSpecies()[i], NS_);
			}

			const boost::filesystem::path file_name = path_output_ / "CEMA.out";
			fCEMA_.open(file_name.c_str(), std::ios::out);
			fCEMA_.setf(std::ios::scientific);

			unsigned int counter = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fCEMA_, "t[s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fCEMA_, "T[K]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fCEMA_, "P[Pa]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fCEMA_, "lambda[1/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fCEMA_, "lambda_plt[1/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fCEMA_, "CEM[1-" + boost::lexical_cast<std::string>(cem_index_) + "]", counter);

			if (print_ei_ == true)
			{
				OpenSMOKE::PrintTagOnASCIILabel(20, fCEMA_, "EI_T", counter);

				if (indices_of_output_species_.size() != 0)
				{
					for (unsigned int i = 0; i < indices_of_output_species_.size(); i++)
						OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fCEMA_, "EI_" + thermodynamicsMap_.NamesOfSpecies()[indices_of_output_species_[i] - 1], counter);
				}
				else
				{
					for (unsigned int i = 0; i < NS_; i++)
						OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fCEMA_, "EI_" + thermodynamicsMap_.NamesOfSpecies()[i], counter);
				}
			}

			if (print_pi_ == true)
			{
				if (indices_of_output_reactions_.size() != 0)
				{
					for (unsigned int i = 0; i < indices_of_output_reactions_.size(); i++)
						OpenSMOKE::PrintTagOnASCIILabel(20, fCEMA_, "PI_" + boost::lexical_cast<std::string>(indices_of_output_reactions_[i]), counter);
				}
				else
				{
					for (unsigned int i = 0; i < NR_; i++)
						OpenSMOKE::PrintTagOnASCIILabel(20, fCEMA_, "PI_" + boost::lexical_cast<std::string>(i + 1), counter);
				}
			}

			fCEMA_ << std::endl;
		}
	}

	void reconstruct_eigenvectors(const int n, double* wi, double* vr, Eigen::MatrixXcd& eig)
	{
		for (int i = 0; i < n; i++)
		{
			int j = 0;
			while (j < n)
			{
				if (wi[j] == (double)0.0)
				{
					eig.real()(i,j) = vr[i + j*n];
					eig.imag()(i,j) = 0.;
					j++;
				}
				else
				{
					eig.real()(i,j) = vr[i + j*n];
					eig.imag()(i,j) = vr[i + (j + 1)*n];
					
					eig.real()(i,j+1) = vr[i + j*n];
					eig.imag()(i,j+1) = -vr[i + (j + 1)*n];
					
					j += 2;
				}
			}
		}
	}

	void OnTheFlyCEMA::CalculateChemicalExplosiveMode(const std::vector<size_t>& lambda_indices_sorted, const std::vector<double>& lambda_real)
	{
		cem_index_  = static_cast<unsigned int>(lambda_indices_sorted[0]);
		cem_lambda_ = lambda_real[cem_index_];
	}

	void OnTheFlyCEMA::CalculateExplosionIndices()
	{
		Eigen::VectorXcd ae_times_be = ae_.cwiseProduct(be_);
		Eigen::VectorXd  ae_times_be_modulus(Nstar_);
		for (unsigned int i = 0; i < Nstar_; i++)
			ae_times_be_modulus(i) = std::sqrt(ae_times_be[i].real()*ae_times_be[i].real() + ae_times_be[i].imag()*ae_times_be[i].imag());
		const double sum = ae_times_be_modulus.sum();

		EI_ = ae_times_be_modulus / sum;
	}


	void OnTheFlyCEMA::CalculateParticipationIndices(const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c)
	{
		OpenSMOKE::OpenSMOKEVectorDouble r(NR_);
		Eigen::VectorXcd r_matrix(NR_);

		// Formation rates
		{
			thermodynamicsMap_.SetTemperature(T);
			thermodynamicsMap_.SetPressure(P_Pa);
			kineticsMap_.SetTemperature(T);
			kineticsMap_.SetPressure(P_Pa);

			kineticsMap_.ReactionRates(c.GetHandle());
			kineticsMap_.ReactionRates(r.GetHandle());

			for (unsigned int i = 0; i < NR_; i++)
				r_matrix(i) = r[i + 1];
		}

		// Participation indices
		{
			Eigen::MatrixXcd be_matrix(1, NS_);
			for (unsigned int i = 0; i < NS_; i++)
				be_matrix(0, i) = be_(i);

			// Product
			Eigen::MatrixXcd P_matrix = be_matrix*kineticsMap_.stoichiometry().stoichiometric_matrix_reactants().transpose();

			Eigen::VectorXcd P(NR_);
			for (unsigned int i = 0; i < NR_; i++)
				P[i] = P_matrix(0, i);

			Eigen::VectorXcd P_times_R = P.cwiseProduct(r_matrix);
			Eigen::VectorXd  P_times_R_modulus(NR_);
			for (unsigned int i = 0; i < NR_; i++)
				P_times_R_modulus(i) = std::sqrt(P_times_R[i].real()*P_times_R[i].real() + P_times_R[i].imag()*P_times_R[i].imag());
			const double sum = P_times_R_modulus.sum();

			for (unsigned int i = 0; i < NR_; i++)
				PI_(i) = P_times_R_modulus(i) / sum;
		}
	}

	void OnTheFlyCEMA::SearchConservativeModes(const std::vector<double>& lambda_real, const std::vector<double>& lambda_mod, std::vector<double>& lambda_real_cleaned, std::vector<size_t>& lambda_original_associations)
	{
		std::vector<size_t> sorted_indices(Nstar_);
		sorted_indices = SortAndTrackIndicesIncreasing(lambda_mod);
		std::vector<size_t> elements_to_remove = std::vector<size_t>(sorted_indices.begin(), sorted_indices.begin() + NE_);

		lambda_real_cleaned = lambda_real;

		std::sort(elements_to_remove.begin(), elements_to_remove.end());
		std::reverse(elements_to_remove.begin(), elements_to_remove.end());

		for (unsigned int i = 0; i<elements_to_remove.size(); i++)
		{
			lambda_real_cleaned.erase(lambda_real_cleaned.begin() + elements_to_remove[i]);
			lambda_original_associations.erase(lambda_original_associations.begin() + elements_to_remove[i]);
		}
	}

	void OnTheFlyCEMA::SortVectors(const std::vector<double>& lambda_real_cleaned, const std::vector<size_t>& lambda_original_associations, std::vector<size_t>& lambda_original_associations_sorted)
	{
		std::vector<size_t> provisional_index(N_);
		provisional_index = SortAndTrackIndicesDecreasing(lambda_real_cleaned);

		for (unsigned int i = 0; i < N_; i++)
			lambda_original_associations_sorted[i] = lambda_original_associations[provisional_index[i]];
	}

	void OnTheFlyCEMA::Calculate(const double t, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEMatrixDouble& J)
	{
		if (lapack_mode_ == true)
			CalculateLapack(t, T, P_Pa, c, J);
		else
			CalculateEigen(t, T, P_Pa, c, J);
	}

	void OnTheFlyCEMA::CalculateEigen(const double t, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEMatrixDouble& J)
	{
		// Populate Jacobian matrices
		for (unsigned int i = 0; i<Nstar_; i++)
			for (unsigned int j = 0; j < Nstar_; j++)
			{
				Jomega_(i, j) = J[i + 1][j + 1];
				JomegaTranspose_(i, j) = J[j + 1][i + 1];
			}

		{
			// Local variables
			std::vector<double> right_lambda_real(Nstar_);
			std::vector<double> right_lambda_imag(Nstar_);
			std::vector<double> right_lambda_mod(Nstar_);
			std::vector<size_t> right_lambda_original_associations_sorted(Nstar_);

			// Calculate right eigenvectors/eigenvalues
			Eigen::EigenSolver<Eigen::MatrixXd> eigenJomega(Jomega_);
			for (unsigned int i = 0; i < Nstar_; i++)
			{
				right_lambda_real[i] = eigenJomega.eigenvalues()[i].real();
				right_lambda_imag[i] = eigenJomega.eigenvalues()[i].imag();
				right_lambda_mod[i] = std::sqrt(right_lambda_real[i] * right_lambda_real[i] + right_lambda_imag[i] * right_lambda_imag[i]);
			}

			{
				std::vector<double> right_lambda_real_cleaned(Nstar_);
				std::vector<size_t> right_lambda_original_associations(Nstar_);
				for (unsigned int i = 0; i < Nstar_; i++)
					right_lambda_original_associations[i] = i;

				// Search conservative modes
				SearchConservativeModes(right_lambda_real, right_lambda_mod, right_lambda_real_cleaned, right_lambda_original_associations);

				// Sort vectors (right eigenvectors)
				SortVectors(right_lambda_real_cleaned, right_lambda_original_associations, right_lambda_original_associations_sorted);

			}

			// Right eigenvector ae (associated to the explosive mode)
			{
				cem_index_ = static_cast<unsigned int>(right_lambda_original_associations_sorted[0]);
				cem_lambda_ = right_lambda_real[cem_index_];
				ae_ = eigenJomega.eigenvectors().col(cem_index_);
			}
		}

		{
			// Local variables
			std::vector<double> left_lambda_real(Nstar_);
			std::vector<double> left_lambda_imag(Nstar_);
			std::vector<double> left_lambda_mod(Nstar_);
			std::vector<size_t> left_lambda_original_associations_sorted(Nstar_);

			// Calculate left eigenvectors/eigenvalues
			Eigen::EigenSolver<Eigen::MatrixXd> eigenJomegaTranspose(JomegaTranspose_);
			for (unsigned int i = 0; i < Nstar_; i++)
			{
				left_lambda_real[i] = eigenJomegaTranspose.eigenvalues()[i].real();
				left_lambda_imag[i] = eigenJomegaTranspose.eigenvalues()[i].imag();
				left_lambda_mod[i] = std::sqrt(left_lambda_real[i] * left_lambda_real[i] + left_lambda_imag[i] * left_lambda_imag[i]);
			}

			{
				std::vector<double> left_lambda_real_cleaned(Nstar_);
				std::vector<size_t> left_lambda_original_associations(Nstar_);
				for (unsigned int i = 0; i < Nstar_; i++)
					left_lambda_original_associations[i] = i;

				// Search conservative modes
				SearchConservativeModes(left_lambda_real, left_lambda_mod, left_lambda_real_cleaned, left_lambda_original_associations);

				// Sort vectors (left eigenvectors)
				SortVectors(left_lambda_real_cleaned, left_lambda_original_associations, left_lambda_original_associations_sorted);
			}

			// Left eigenvector b
			{
				const unsigned int left_cem_index = static_cast<unsigned int>(left_lambda_original_associations_sorted[0]);
				const double left_cem_lambda = left_lambda_real[left_cem_index];
				be_ = eigenJomegaTranspose.eigenvectors().row(left_cem_index);

				// Difference with right calculations
				const double relative_error = std::fabs(left_cem_lambda - cem_lambda_) / (1.e-12 + 0.50*(left_cem_lambda + cem_lambda_));
				if (relative_error > 1.e-3)
				{
					std::cout << " * Explosive Eigenvalue (right side): " << cem_lambda_ << std::endl;
					std::cout << " * Explosive Eigenvalue (left side):  " << left_cem_lambda << std::endl;
					std::cout << " * Relative error (%):                " << relative_error*100. << std::endl;
					OpenSMOKE::FatalErrorMessage("Something was wrong with CEMA: left and right explosive eigenvalues do not match each other.");
				}
			}
		}

		// Explosion Indices (EI)
		CalculateExplosionIndices();

		// Participation Indices (PI)
		CalculateParticipationIndices(T, P_Pa, c);
	}

	void OnTheFlyCEMA::CalculateLapack(const double t, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEMatrixDouble& J)
	{
		// Local variables
		double* vl = new double[Nstar_*Nstar_];
		double* vr = new double[Nstar_*Nstar_];
		std::vector<double> lambda_real(Nstar_);
		std::vector<double> lambda_imag(Nstar_);
		std::vector<double> lambda_mod(Nstar_);
		std::vector<size_t> lambda_original_associations_sorted(Nstar_);

		// Populate Jacobian matrix
		for (unsigned int i = 0; i<Nstar_; i++)
			for (unsigned int j = 0; j < Nstar_; j++)
				Jomega_(i, j) = J[i + 1][j + 1];

		// Calculate eigenvectors/eigenvalues
		#if OPENSMOKE_USE_MKL == 1
		{
			int info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'V', 'V', Nstar_, Jomega_.data(), Nstar_, lambda_real.data(), lambda_imag.data(), vl, Nstar_, vr, Nstar_);

			if (info != 0)
				OpenSMOKE::FatalErrorMessage("OnTheFlyCEMA::Calculate: LAPACKE_dgeev failure");
		}
		#else
			OpenSMOKE::FatalErrorMessage("OnTheFlyCEMA::Calculate: available only if compilation was carried out with Intel MKL");
		#endif

		for (unsigned int i = 0; i < Nstar_; i++)
			lambda_mod[i] = std::sqrt(lambda_real[i]*lambda_real[i] + lambda_imag[i]*lambda_imag[i]);

		// Search for conservative modes and sort eigenvalues
		{
			std::vector<double> lambda_real_cleaned(Nstar_);
			std::vector<size_t> lambda_original_associations(Nstar_);
			for (unsigned int i = 0; i < Nstar_; i++)
				lambda_original_associations[i] = i;

			// Search conservative modes
			SearchConservativeModes(lambda_real, lambda_mod, lambda_real_cleaned, lambda_original_associations);

			// Sort vectors (right eigenvectors)
			SortVectors(lambda_real_cleaned, lambda_original_associations, lambda_original_associations_sorted);
		}

		// Chemical Explosive Mode (CEM)
		CalculateChemicalExplosiveMode(lambda_original_associations_sorted, lambda_real);

		// Eigenvectors associated to the CEM
		{
			// EigenVectors
			Eigen::MatrixXcd right_eigenvectors(Nstar_, Nstar_);
			Eigen::MatrixXcd left_eigenvectors(Nstar_, Nstar_);
			reconstruct_eigenvectors(Nstar_, lambda_imag.data(), vr, right_eigenvectors);
			reconstruct_eigenvectors(Nstar_, lambda_imag.data(), vl, left_eigenvectors);

			ae_ = right_eigenvectors.col(cem_index_);
			be_ = left_eigenvectors.col(cem_index_);
		}

		// Explosion Indices (EI)
		CalculateExplosionIndices();

		// Participation Indices (PI)
		CalculateParticipationIndices(T, P_Pa, c);
	}
	
	void OnTheFlyCEMA::WriteOnFile(const double t, const double T, const double P_Pa)
	{
		// Print on file
		if (is_print_on_file_active_ == true)
		{
			// Basic variables
			fCEMA_ << std::setw(20) << std::left << t;
			fCEMA_ << std::setw(20) << std::left << T;
			fCEMA_ << std::setw(20) << std::left << P_Pa;

			// CEM: Chemical explosive mode
			fCEMA_ << std::setw(20) << std::left << cem_lambda_;
			fCEMA_ << std::setw(20) << std::left << cem_for_plots(cem_lambda_);
			fCEMA_ << std::setw(20) << std::left << cem_index_+1;

			// Explosive indices
			if (print_ei_ == true)
			{
				fCEMA_ << std::setw(20) << std::left << EI_(Nstar_ - 1);

				if (indices_of_output_species_.size() != 0)
				{
					for (unsigned int i = 0; i < indices_of_output_species_.size(); i++)
						fCEMA_ << std::setw(widths_of_output_species_[i]) << std::left << EI_(indices_of_output_species_[i] - 1);
				}
				else
				{
					for (unsigned int i = 0; i < NS_; i++)
						fCEMA_ << std::setw(widths_of_output_species_[i]) << std::left << EI_(i);
				}
			}

			// Participation indices
			if (print_pi_ == true)
			{
				if (indices_of_output_reactions_.size() != 0)
				{
					for (unsigned int i = 0; i < indices_of_output_reactions_.size(); i++)
						fCEMA_ << std::setw(20) << std::left << PI_(indices_of_output_reactions_[i] - 1);
				}
				else
				{
					for (unsigned int i = 0; i < NR_; i++)
						fCEMA_ << std::setw(20) << std::left << PI_(i);
				}
			}

			fCEMA_ << std::endl;
		}
	}
}

