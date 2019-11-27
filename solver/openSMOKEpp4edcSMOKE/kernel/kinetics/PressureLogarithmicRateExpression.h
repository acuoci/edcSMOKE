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
|	License                                                               |
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


#ifndef OpenSMOKE_PressureLogarithmicRateExpression_H
#define	OpenSMOKE_PressureLogarithmicRateExpression_H

#include "math/PhysicalConstants.h"
#include <vector>
#include <string.h>

namespace OpenSMOKE
{
	//!  A class to manage a Pressure Dependence through logarithmic interpolation rate expression
	/*!
		 This class provides the interface and the tools to manage a Pressure Dependence through 
		 logarithmic interpolation rate expression (PLOG in CHEMKIN standard)
		 (see Section 3.6.3 of CHEMKIN manual)
	*/

	class PressureLogarithmicRateExpression
	{
	public:

		/**
		* Default constructor
		*/
		PressureLogarithmicRateExpression();

		/**
		* Default copy constructor
		*/
		PressureLogarithmicRateExpression(const PressureLogarithmicRateExpression& orig);

		/**
		* Default destructor
		*/
		//virtual ~PressureLogarithmicRateExpression();

		/**
		*@brief Prepares all the data to be used for evaluating the reaction rate
		*/
		void Setup(std::vector<double> coefficients_);

		/**
		*@brief Evaluates the kinetic constant
		*/
		double KineticConstant(const double T, const double P);

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

		unsigned int N;										//!< number of points along the pressure axis
		std::vector< std::vector<double> > lnA_;			//!< logarithm (natural) of factor frequencies
		std::vector< std::vector<double> > Beta_;			//!< temperature exponents
		std::vector< std::vector<double> > E_over_R_;		//!< normalized activation energy [K]
		std::vector< std::vector<double> > p_;				//!< list of pressure points [Pa]
		std::vector< std::vector<double> > lnp_;			//!< logarithm of pressure points (for efficiency reasons only)
	};

}

#include "PressureLogarithmicRateExpression.hpp"

#endif	/* OpenSMOKE_PressureLogarithmicRateExpression_H */

