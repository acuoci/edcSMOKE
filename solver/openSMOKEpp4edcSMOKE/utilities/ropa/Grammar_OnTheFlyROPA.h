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

#ifndef OpenSMOKE_Grammar_OnTheFlyROPA_H
#define	OpenSMOKE_Grammar_OnTheFlyROPA_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
	class Grammar_OnTheFlyROPA : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Species", 
																OpenSMOKE::VECTOR_STRING, 
																"Species for which the ROPA is turned on", 
																true) );

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ReferenceSpecies",
																OpenSMOKE::SINGLE_STRING,
																"Reference species used to calculated the conversion",
																true));
                                                                
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ReferenceLiquidMass",
																OpenSMOKE::SINGLE_BOOL,
																"Conversions are calculated using the Liquid Mass as reference",
																false));                                                                

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Conversions", 
																OpenSMOKE::VECTOR_DOUBLE, 
																"List of conversions of main fuel at which the ROPA has to be calculated", 
																false, 
																"none",
																"none",
																"@Times @NumberOfSteps" ) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Times", 
																OpenSMOKE::VECTOR_STRING, 
																"List of times at which the the ROPA has to be calculated", 
																false,
																"none",
																"none",
																"@Conversions @NumberOfSteps"));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfSteps", 
																OpenSMOKE::SINGLE_INT, 
																"Frequency with which the ROPA has to be calculated (default: 25)", 
																false,
																"none",
																"none",
																"@Times @Conversions"));
		
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Threshold", 
																OpenSMOKE::SINGLE_DOUBLE, 
																"Threshold: (default 3%) ", 
																false) );

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MergeForwardAndBackwardReactions",
																OpenSMOKE::SINGLE_BOOL,
																"The forward and backward contributions are merged together before performing the ROPA (default: false)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CompactOutput",
																OpenSMOKE::SINGLE_BOOL,
																"Only species whose amount is larger than a minimum threshold are written on the output file (default: true)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@WriteXML",
																OpenSMOKE::SINGLE_BOOL,
																"Write the ROPA results also in a xml file (default: false)",
																false));
		}
	};
}

#endif	/* OpenSMOKE_Grammar_OnTheFlyROPA_H */

