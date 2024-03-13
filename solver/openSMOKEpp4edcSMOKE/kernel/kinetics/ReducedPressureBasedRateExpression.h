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


#ifndef OpenSMOKE_ReducedPressureBasedRateExpression_H
#define	OpenSMOKE_ReducedPressureBasedRateExpression_H

#include "math/PhysicalConstants.h"
#include <vector>
#include <string.h>

#include "TroeRateExpression.h"
#include "ChebyshevPolynomialRateExpression.h"
#include "PressureLogarithmicRateExpression.h"

class TabulatedRateExpression
{
public:

	/**
	*@brief Evaluates the kinetic constant
	*/
	double KineticConstant(const double T, const double P) { return 0.;  }
};

namespace OpenSMOKE
{
	//!  A class to manage a pressure-base mixture rules rate expression
	/*!
		 See: https://github.com/Simple-ape/MixtureRules
	*/

	class ReducedPressureBasedRateExpression
	{

	private:

		enum TypeRateCoefficients { UNDEFINED, TROE, PLOG, CHEBISHEV, TABULATION };

	public:

		/**
		* Default constructor
		*/
		ReducedPressureBasedRateExpression();

		/**
		* Default copy constructor
		*/
		ReducedPressureBasedRateExpression(const ReducedPressureBasedRateExpression& orig);

		/**
		* Default destructor
		*/
		//virtual ~ReducedPressureBasedRateExpression();

		void SetKineticParameters(const double A, const double beta, const double E);

		void SetStoichiometry(const std::vector<unsigned int>& reactant_nu_indices, const std::vector<double>& reactant_nu, const std::vector<unsigned int>& product_nu_indices, const std::vector<double>& product_nu);

		void SetReactionOrders(const std::vector<unsigned int>& reactant_lambda_indices, const std::vector<double>& reactant_lambda, const std::vector<unsigned int>& product_lambda_indices, const std::vector<double>& product_lambda);

		/**
		*@brief Prepares all the data to be used from the CHEMKIN lines
		*/
		void Setup(const std::vector<std::string>& lines);

		/**
		*@brief Prepares all the data to be used for evaluating the reaction rate
		*/
		void Setup(std::vector<double> coefficients);

		/**
		*@brief Analyzes the CHEMKIN lines describing the kinetic constant parameters
		*/
		void AnalyzeKineticConstantLines(const int index, const std::vector<std::string>& lines);

		/**
		*@brief Calculates the activity coefficients (only in case of non-linear mixture rule)
		*/
		void CalculateActivityCoefficients(const double T, const double* X, double& f_default, std::vector<double>& f);

		/**
		*@brief Evaluates the kinetic constant
		*/
		double KineticConstant(const double T, const double P, const double* X);

		/**
		*@brief Writes on a file the data about the pressure table
		*/
		void WriteOnASCIIFileOldStyle(std::ostream& fOutput) const;

		/**
		*@brief Reads from a file the data about the reaction
		*/		
		void ReadFromASCIIFile(std::istream& fInput);

		/**
		*@brief Writes a short summary on file
		*/	
		void WriteShortSummaryOnASCIIFile(std::ostream& fOutput, const double conversion_factor_A) const;

		/**
		*@brief Writes reaction data in CHEMKIN format
		*/
		void WriteCHEMKINReactionData(std::stringstream& reaction_data, const double conversion_factor_A) const;
                
        private:
            
        /**
		*@brief Returns an error message
		*/
		void ErrorMessage(const std::string message);

		/**
		*@brief Returns a warning message
		*/
		void WarningMessage(const std::string message);

	private:

		bool non_linear_mixture_rule_;						//!< non-linear mixture rule

		unsigned int N_;									//!< number of explicit colliders
		std::vector<unsigned int>	i_colliders_;			//!< indices of explicit colliders (0-based)
		std::vector<std::string>	names_colliders_;		//!< names of explicit colliders (0-based)

		unsigned int NTB_;									//!< number of non-unitary third-body efficiencies
		std::vector<double>			tb_efficiencies_;		//!< third-body non-unitary efficiencies
		std::vector<unsigned int>	tb_indices_;			//!< third-body non-unitary efficiencies species indices (0-based)
		std::vector<std::string>	tb_names_;				//!< third-body non-unitary efficiencies species names

		// Absolute values of least negative chemically significant eigenvalue in the low-pressure limit
		std::vector<double> lambda0_A_;						//!< frequency-factor explicit colliders [a.u.]
		std::vector<double> lambda0_Beta_;					//!< temperature exponent for explicit colliders [-]
		std::vector<double> lambda0_E_;						//!< activation energy explicit colliders [J/kmol]
		double lambda0_A_default_;							//!< frequency factor default collider [a.u.]
		double lambda0_Beta_default_;						//!< frequency factor default collider [-]
		double lambda0_E_default_;							//!< frequency factor default collider [J/kmol]

		std::vector<TypeRateCoefficients> kappaMu_type_;	//!< 
		TypeRateCoefficients kappaMu_type_default_;			//!<

		std::vector<PressureLogarithmicRateExpression>	kappaMu_PLOG_;
		std::vector<ChebyshevPolynomialRateExpression>	kappaMu_CHEB_;
		std::vector<TroeRateExpression>					kappaMu_TROE_;
		std::vector<TabulatedRateExpression>			kappaMu_TAB_;

		PressureLogarithmicRateExpression				kappaMu_PLOG_default_;
		ChebyshevPolynomialRateExpression				kappaMu_CHEB_default_;
		TroeRateExpression								kappaMu_TROE_default_;
		TabulatedRateExpression							kappaMu_TAB_default_;

		double A_;
		double beta_;
		double E_;

		std::vector<unsigned int> reactant_nu_indices_;
		std::vector<double> reactant_nu_;
		std::vector<unsigned int> product_nu_indices_;
		std::vector<double> product_nu_;

		std::vector<unsigned int> reactant_lambda_indices_;
		std::vector<double> reactant_lambda_;
		std::vector<unsigned int> product_lambda_indices_;
		std::vector<double> product_lambda_;
	};

}

#include "ReducedPressureBasedRateExpression.hpp"

#endif	/* OpenSMOKE_ReducedPressureBasedRateExpression_H */

