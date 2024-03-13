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


#ifndef OpenSMOKE_TroeRateExpression_H
#define	OpenSMOKE_TroeRateExpression_H

#include "math/PhysicalConstants.h"
#include <vector>
#include <string.h>
#include <Eigen/Dense>

namespace OpenSMOKE
{
	//!  A class to manage a Troe Rate Expression
	/*!
		 This class provides the interface and the tools to manage a Troe Rate Expression
		 (TROE in CHEMKIN standard)
	*/

	class TroeRateExpression
	{
	public:

		/**
		* Default constructor
		*/
		TroeRateExpression();

		/**
		* Default copy constructor
		*/
		TroeRateExpression(const TroeRateExpression& orig);

		/**
		* Default destructor
		*/
		//virtual ~TroeRateExpression();

		/**
		*@brief Prepares all the data to be used for evaluating the reaction rate
		*/
		void Setup(const std::vector<double>& troe_coefficients, const std::vector<double>& low_coefficients, const std::vector<double>& inf_coefficients);

		/**
		*@brief Evaluates the kinetic constant
		*/
		double KineticConstant(const double T, const double P);

		/**
		*@brief Evaluates the kinetic constant
		*/
		void FitttingArrheniusLaw(std::ostream& fOutput);

		/**
		*@brief Writes on a file the data about the Troe rate expression
		*/
        void WriteStatus(std::ostream& fOutput) const;

		/**
		*@brief Reads from a file the data about the Troe rate expression
		*/		
		void ReadFromASCIIFile(std::istream& fInput);

		/**
		*@brief Writes a short summary on file
		*/	
		void WriteShortSummaryOnASCIIFile(std::ostream& fOutput) const;
                
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

		std::vector<double> troe_coefficients_;
		std::vector<double> low_coefficients_;
		std::vector<double> inf_coefficients_;
	};

}

#include "TroeRateExpression.hpp"

#endif	/* OpenSMOKE_TroeRateExpression_H */

