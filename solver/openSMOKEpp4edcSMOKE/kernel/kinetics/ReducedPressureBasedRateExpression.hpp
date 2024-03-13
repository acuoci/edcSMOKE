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

#include "KineticsUtilityFunctions.h"

namespace OpenSMOKE
{

	ReducedPressureBasedRateExpression::ReducedPressureBasedRateExpression()
	{
	}
	
	ReducedPressureBasedRateExpression::ReducedPressureBasedRateExpression(const ReducedPressureBasedRateExpression& orig)
	{
	}
	
	//ReducedPressureBasedRateExpression::~ReducedPressureBasedRateExpression()
	//{
	//}

	void ReducedPressureBasedRateExpression::ErrorMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  ReducedPressureBasedRateExpression"	<< std::endl;
		std::cout << "Error:  " << message							<< std::endl;
		std::cout << "Press a key to continue... "					<< std::endl;
		getchar();
		exit(-1);
	}

	void ReducedPressureBasedRateExpression::WarningMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  ReducedPressureBasedRateExpression"	<< std::endl;
		std::cout << "Warning:  "	<< message						<< std::endl;
		std::cout << "Press a key to continue... "					<< std::endl;
		getchar();
	}

	void ReducedPressureBasedRateExpression::SetKineticParameters(const double A, const double beta, const double E)
	{
		A_ = A;
		beta_ = beta;
		E_ = E;
	}

	void ReducedPressureBasedRateExpression::SetStoichiometry(const std::vector<unsigned int>& reactant_nu_indices, const std::vector<double>& reactant_nu, const std::vector<unsigned int>& product_nu_indices, const std::vector<double>& product_nu)
	{
		reactant_nu_indices_ = reactant_nu_indices;
		reactant_nu_ = reactant_nu;
		product_nu_indices_ = product_nu_indices;
		product_nu_ = product_nu;
	}

	void ReducedPressureBasedRateExpression::SetReactionOrders(const std::vector<unsigned int>& reactant_lambda_indices, const std::vector<double>& reactant_lambda, const std::vector<unsigned int>& product_lambda_indices, const std::vector<double>& product_lambda)
	{
		reactant_lambda_indices_ = reactant_lambda_indices;
		reactant_lambda_ = reactant_lambda;
		product_lambda_indices_ = product_lambda_indices;
		product_lambda_ = product_lambda;
	}

	void ReducedPressureBasedRateExpression::Setup(const std::vector<std::string>& lines)
	{
		// Check if the mixture rule is linear or non-linear
		{
			std::string line = lines[1];
			if (OpenSMOKE_Utilities::ReadReactionKeyWordPlusWord("RPBMR", line) == "LINEAR")			non_linear_mixture_rule_ = false;
			else if (OpenSMOKE_Utilities::ReadReactionKeyWordPlusWord("RPBMR", line) == "NONLINEAR")	non_linear_mixture_rule_ = true;
			else ErrorMessage("Unrecognized RPBMR type. Available options: LINEAR | NONLINEAR");
		}

		// Analyze collider lines
		{
			std::vector<unsigned int> lines_with_collider;

			// Look for lines with colliders
			for (unsigned int k = 2; k < lines.size(); k++)
			{
				std::string line = lines[k];
				std::string keyword;
				OpenSMOKE_Utilities::LookForKeyWord(line, keyword);
				if (keyword == "COLL")
					lines_with_collider.push_back(k);
			}

			// Analyze line with collider
			std::vector< std::vector<std::string> > lines_kmus;

			for (unsigned int k = 0; k < lines_with_collider.size(); k++)
			{
				// Identify parameters for Lambda0
				std::string line = lines[lines_with_collider[k]];
				std::vector<std::string> words;
				bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusWords("RPBMR. ", line, 4, words);
				if (tag == false) 
					ErrorMessage("Error in line: " + line);

				// Identify groups of lines for kMu
				std::vector<std::string> lines_kmu;
				if (k == lines_with_collider.size() - 1)
				{
					for (unsigned int i = lines_with_collider[k] + 1; i < lines.size(); i++)
						lines_kmu.push_back(lines[i]);
				}
				else
				{
					for (unsigned int i = lines_with_collider[k] + 1; i < lines_with_collider[k + 1]; i++)
						lines_kmu.push_back(lines[i]);
				}

				if (words[0] == "DEFAULT")
				{
					lambda0_A_default_ = boost::lexical_cast<double>(words[1]);		//!< frequency factor default collider [a.u.]
					lambda0_Beta_default_ = boost::lexical_cast<double>(words[2]);	//!< frequency factor default collider [-]
					lambda0_E_default_ = boost::lexical_cast<double>(words[3]);		//!< frequency factor default collider [J/kmol]

					// Third-body efficiencies
					{
						std::string line = lines[lines_with_collider[k]+1];
						boost::algorithm::trim(line);
						typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_slash;
						boost::char_separator<char> sep_slash("/");
						tokenizer_slash tokens(line, sep_slash);
						const std::size_t n = std::distance(tokens.begin(), tokens.end());

						if (n % 2 != 0 || line.at(line.size() - 1) != '/')
							ErrorMessage("Error in third-body efficiencies: " + line);

						tokenizer_slash::iterator tok_slash = tokens.begin();
						for (std::size_t i = 1; i <= n / 2; i++)
						{
							std::string name = *tok_slash;
							boost::algorithm::trim(name);
							++tok_slash;
							std::string number = *(tok_slash);
							boost::algorithm::trim(number);
							++tok_slash;

							tb_names_.push_back(name);
							tb_efficiencies_.push_back(boost::lexical_cast<double>(number));
						}
					}

					// Analyze data for kMu
					lines_kmu.erase(lines_kmu.begin());

					std::cout << "Collider " << words[0] << std::endl;
					for (unsigned int i=0;i< lines_kmu.size(); i++)
						std::cout << lines_kmu[i] << std::endl;

					// Analyze
					AnalyzeKineticConstantLines(-1, lines_kmu);
				}
				else
				{
					names_colliders_.push_back(words[0]);							//!< names of explicit colliders	
					lambda0_A_.push_back(boost::lexical_cast<double>(words[1]));	//!< frequency factor default collider [a.u.]
					lambda0_Beta_.push_back(boost::lexical_cast<double>(words[2]));	//!< frequency factor default collider [-]
					lambda0_E_.push_back(boost::lexical_cast<double>(words[3]));	//!< frequency factor default collider [J/kmol]

					// Analyze data for kMu
					std::cout << "Collider " << words[0] << std::endl;
					for (unsigned int i = 0; i < lines_kmu.size(); i++)
						std::cout << lines_kmu[i] << std::endl;

					// Analyze
					std::cout << "TORMV Transfer" << std::endl;
					lines_kmus.push_back(lines_kmu);	
				}
			}

			// Final checks
			// TODO

			// Final operations
			N_ = names_colliders_.size();
			NTB_ = tb_names_.size();

			// Summary on the screen
				std::cout << "DEFAULT" << " " << lambda0_A_default_ << " " << lambda0_Beta_default_ << " " << lambda0_E_default_ << std::endl;
			for (unsigned int i = 0; i < N_; i++)
				std::cout << names_colliders_[i] << " " << lambda0_A_[i] << " " << lambda0_Beta_[i] << " " << lambda0_E_[i] << std::endl;

			for (unsigned int i = 0; i < NTB_; i++)
				std::cout << tb_names_[i] << " " << tb_efficiencies_[i] << std::endl;

			std::cout << "Allocation" << std::endl;
			kappaMu_type_.resize(N_);
			kappaMu_PLOG_.resize(N_);
			kappaMu_CHEB_.resize(N_);
			kappaMu_TROE_.resize(N_);
			kappaMu_TAB_.resize(N_);

			for (unsigned int i = 0; i < N_; i++)
			{
				std::cout << "Analyze " << i << std::endl;
				AnalyzeKineticConstantLines(i, lines_kmus[i]);
			}
			
		}
	}

	void ReducedPressureBasedRateExpression::AnalyzeKineticConstantLines(const int index, const std::vector<std::string>& lines)
	{
		// Recognize reaction type
		TypeRateCoefficients type_constant = TypeRateCoefficients::UNDEFINED;

		for (unsigned int j = 0; j < lines.size(); j++)
		{
			std::string line = lines[j];
			std::string keyword;
			OpenSMOKE_Utilities::LookForKeyWord(line, keyword);

			if (keyword == "PLOG")
			{
				type_constant = TypeRateCoefficients::PLOG;
				break;
			}
			else if (keyword == "CHEB")
			{
				type_constant = TypeRateCoefficients::CHEBISHEV;
				break;
			}
			else if (keyword == "TROE")
			{
				type_constant = TypeRateCoefficients::TROE;
				break;
			}
			else if (keyword == "TABULATION")
			{
				type_constant = TypeRateCoefficients::TABULATION;
				break;
			}
		}


		if (type_constant == TypeRateCoefficients::PLOG)
		{
			std::vector<double> plog_coefficients;

			for (unsigned int j = 0; j < lines.size(); j++)
			{
				std::string line = lines[j];

				std::string keyword;
				OpenSMOKE_Utilities::LookForKeyWord(line, keyword);

				if (keyword == "PLOG")
				{
					std::vector<double> coefficients;
					bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Pressure Dependence through Logarithmic Interpolation Rate Expressions (PLOG). ", line, 4, coefficients);

					if (tag == false)
						ErrorMessage("Error in PLOG reaction at line: " + line);

					for (unsigned int i = 0; i < coefficients.size(); i++)
						plog_coefficients.push_back(coefficients[i]);
				}
				else
					ErrorMessage("Unexpected option for a PLOG reaction at line: " + line);
			}

			// TODO (conversion factors)
			const double conversion_factor_A = 1.;
			const double conversion_factor_E = 1.;
			plog_coefficients.push_back(conversion_factor_A);
			plog_coefficients.push_back(conversion_factor_E);

			// Creating PLOG reaction
			if (index == -1)
			{
				kappaMu_type_default_ = TypeRateCoefficients::PLOG;
				kappaMu_PLOG_default_.Setup(plog_coefficients);
			}
			else
			{
				kappaMu_type_[index] = TypeRateCoefficients::PLOG;
				kappaMu_PLOG_[index].Setup(plog_coefficients);
			}
		}
		else if (type_constant == TypeRateCoefficients::TROE)
		{
			std::vector<double> troe_coefficients;
			std::vector<double> low_coefficients;
			std::vector<double> inf_coefficients(3);
			inf_coefficients[0] = A_;
			inf_coefficients[1] = beta_;
			inf_coefficients[2] = E_;

			for (unsigned int j = 0; j < lines.size(); j++)
			{
				std::string line = lines[j];

				std::string keyword;
				OpenSMOKE_Utilities::LookForKeyWord(line, keyword);

				if (keyword == "TROE")
				{
					bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Pressure-dependent reaction (TROE). ", line, 3, 4, troe_coefficients);
					if (tag == false)
						ErrorMessage("Error in TROE reaction at line: " + line);
				}
				else if (keyword == "LOW")
				{
					std::vector<double> coefficients;
					bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Low-pressure limit (LOW). ", line, 3, low_coefficients);

					if (tag == false)
						ErrorMessage("Error in LOW reaction at line: " + line);
				}
			}

			// TODO: conversions


			// Creating TROE reaction
			if (index == -1)
			{
				kappaMu_type_default_ = TypeRateCoefficients::TROE;
				kappaMu_TROE_default_.Setup(troe_coefficients, low_coefficients, inf_coefficients);
			}
			else
			{
				kappaMu_type_[index] = TypeRateCoefficients::TROE;
				kappaMu_TROE_[index].Setup(troe_coefficients, low_coefficients, inf_coefficients);
			}
		}
		else if (type_constant == TypeRateCoefficients::CHEBISHEV)
		{
			std::vector<double> chebyshev_coefficients;
			std::vector<double> chebyshev_pressure_limits;			//!< pressure limits
			std::vector<double> chebyshev_temperature_limits;		//!< temperature limits

			for (unsigned int j = 0; j < lines.size(); j++)
			{
				std::string line = lines[j];

				std::string keyword;
				OpenSMOKE_Utilities::LookForKeyWord(line, keyword);

				// Check for PCHEB
				if (keyword == "PCHEB")
				{
					const bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Pressure limit for Chebyshev Polynomial Rate Expression (PCHEB). ", line, 2, chebyshev_pressure_limits);
					if (tag == false)
						ErrorMessage("Error in PCHEB option at line: " + line);
				}
				else if (keyword == "TCHEB")
				{
					const bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Temperature limit for Chebyshev Polynomial Rate Expression (TCHEB). ", line, 2, chebyshev_temperature_limits);
					if (tag == false)
						ErrorMessage("Error in TCHEB option at line: " + line);
				}
				else if (keyword == "CHEB")
				{
					std::vector<double> coefficients;
					bool tag = OpenSMOKE_Utilities::ReadCoefficients("Chebishev Polynomial Rate Expressions (CHEB). ", line, coefficients);

					if (tag == false)
						ErrorMessage("Error in CHEB reaction at line: " + line);

					for (unsigned int i = 0; i < coefficients.size(); i++)
						chebyshev_coefficients.push_back(coefficients[i]);
				}
				else
					ErrorMessage("Unexpected option for a CHEB reaction at line: " + line);
			}

			// Checking coefficients
			unsigned int NxM = int(chebyshev_coefficients[0]) * int(chebyshev_coefficients[1]);
			if (NxM != chebyshev_coefficients.size() - 2)
			{
				std::cout << "Error in the definition of the Chebyshev coefficients." << std::endl;
				std::cout << "Expected coefficients: " << NxM << " - Given: " << chebyshev_coefficients.size() - 2 << std::endl;
				ErrorMessage("Error in CHEB reaction");
			}

			if (chebyshev_pressure_limits.size() == 0)
			{
				chebyshev_pressure_limits.resize(2);
				chebyshev_pressure_limits[0] = 0.001;
				chebyshev_pressure_limits[1] = 100.;
			}
			if (chebyshev_temperature_limits.size() == 0)
			{
				chebyshev_temperature_limits.resize(2);
				chebyshev_temperature_limits[0] = 300.;
				chebyshev_temperature_limits[1] = 2500.;
			}

			// TODO (conversion factors)
			const double conversion_factor_A = 1.;
			chebyshev_coefficients.push_back(conversion_factor_A);

			// Creating CHEB reaction
			if (index == -1)
			{
				kappaMu_type_default_ = TypeRateCoefficients::CHEBISHEV;
				kappaMu_CHEB_default_.Setup(chebyshev_coefficients, chebyshev_pressure_limits, chebyshev_temperature_limits);
			}
			else
			{
				kappaMu_type_[index] = TypeRateCoefficients::CHEBISHEV;
				kappaMu_CHEB_[index].Setup(chebyshev_coefficients, chebyshev_pressure_limits, chebyshev_temperature_limits);
			}
		}
		else if (type_constant == TypeRateCoefficients::TABULATION)
		{
		}
	}

	void ReducedPressureBasedRateExpression::Setup(std::vector<double> coefficients)
	{

	}

	void ReducedPressureBasedRateExpression::CalculateActivityCoefficients(const double T, const double* X, double& f_default, std::vector<double>& f)
	{
		// TODO
	}


	double ReducedPressureBasedRateExpression::KineticConstant(const double T, const double P, const double* X)
	{
		// Activity coefficients
		double f_default = 1.;
		std::vector<double> f(N_);
		std::fill(f.begin(), f.end(), 1.0);

		// Calculate activity coefficients in case of non-linear mixture rule
		if (non_linear_mixture_rule_ == true)
			CalculateActivityCoefficients(T, X, f_default, f);

		// Absolute values of least negative chemically significant eigenvalue in the low-pressure limit (a.u.)
		const double lambda0_default = lambda0_A_default_ * std::pow(T, lambda0_Beta_default_) * std::exp(-lambda0_E_default_ / PhysicalConstants::R_J_kmol / T);
		std::vector<double> lambda0(N_);
		for (unsigned int i=0;i<N_;i++)
			lambda0[i] = lambda0_A_[i] * std::pow(T, lambda0_Beta_[i]) * std::exp(-lambda0_E_[i] / PhysicalConstants::R_J_kmol / T);

		// Effective "default" mole fraction
		double X_default = 1.;
		for (unsigned int i = 0; i < N_; i++)
			X_default -= X[i_colliders_[i]];
		for (unsigned int i = 0; i < NTB_; i++)
		{
			const unsigned int j = tb_indices_[i];
			X_default -= X[j];
			X_default += X[j] * tb_efficiencies_[i];
		}

		// Absolute pressures and fractional contributions to the reduced pressure
		double numerator = f_default * lambda0_default * X_default;
		for (unsigned int i = 0; i < N_; i++)
			numerator += f[i] * lambda0[i] * X[i_colliders_[i]];	
		
		const double Pi_default = P * numerator / lambda0_default;
		std::vector<double> Pi(N_);
		for (unsigned int i = 0; i < N_; i++)
			Pi[i] = P * numerator / lambda0[i];

		const double Xtilde_default = f_default * lambda0_default * X_default / numerator ;
		std::vector<double> Xtilde(N_);
		for (unsigned int i = 0; i < N_; i++)
			Xtilde[i] =  f[i] * lambda0[i] * X[i_colliders_[i]] / numerator;

		// Linear mixture rule
		double kappaMu_default = 0.;
		if (kappaMu_type_default_ == TypeRateCoefficients::TROE)
			kappaMu_default = kappaMu_TROE_default_.KineticConstant(T, Pi_default);
		else if (kappaMu_type_default_ == TypeRateCoefficients::PLOG)
			kappaMu_default = kappaMu_PLOG_default_.KineticConstant(T, Pi_default);
		else if (kappaMu_type_default_ == TypeRateCoefficients::CHEBISHEV)
			kappaMu_default = kappaMu_CHEB_default_.KineticConstant(T, Pi_default);
		else if (kappaMu_type_default_ == TypeRateCoefficients::TABULATION)
			kappaMu_default = kappaMu_TAB_default_.KineticConstant(T, Pi_default);
		
		std::vector<double> kappaMu(N_);
		for (unsigned int i = 0; i < N_; i++)
		{
			if (kappaMu_type_[i] == TypeRateCoefficients::TROE)
				kappaMu[i] = kappaMu_TROE_[i].KineticConstant(T, Pi[i]);
			else if (kappaMu_type_[i] == TypeRateCoefficients::PLOG)
				kappaMu[i] = kappaMu_PLOG_[i].KineticConstant(T, Pi[i]);
			else if (kappaMu_type_[i] == TypeRateCoefficients::CHEBISHEV)
				kappaMu[i] = kappaMu_CHEB_[i].KineticConstant(T, Pi[i]);
			else if (kappaMu_type_[i] == TypeRateCoefficients::TABULATION)
				kappaMu[i] = kappaMu_TAB_[i].KineticConstant(T, Pi[i]);
		}
		
		double kappa = kappaMu_default * Xtilde_default;
		for (unsigned int i = 0; i < N_; i++)
			kappa += kappaMu[i] * Xtilde[i];

		return kappa;
	}

	void ReducedPressureBasedRateExpression::ReadFromASCIIFile(std::istream& fInput)
	{
		// TODO
	}

	void ReducedPressureBasedRateExpression::WriteOnASCIIFileOldStyle(std::ostream& fOutput) const
	{
		// TODO
	}

	void ReducedPressureBasedRateExpression::WriteShortSummaryOnASCIIFile(std::ostream& fOutput, const double conversion_factor_A) const
	{
		// TODO
	}

	void ReducedPressureBasedRateExpression::WriteCHEMKINReactionData(std::stringstream& reaction_data, const double conversion_factor_A) const
	{
		// TODO
	}
}
