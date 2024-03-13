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

#ifndef OpenSMOKE_Grammar_OnTheFlyPostProcessing_H
#define	OpenSMOKE_Grammar_OnTheFlyPostProcessing_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
	class Grammar_OnTheFlyPostProcessing : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FormationRates", 
																OpenSMOKE::VECTOR_STRING, 
																"Species for which the formation rates are required (units are speciied trough the @FormationRatesUnits option)", 
																false) );

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FormationRatesUnits",
																OpenSMOKE::SINGLE_STRING,
																"Units for formation rates: mass (kg/m3/s) | moles (kmol/m3/s) (default: moles)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ReactionRates",
																OpenSMOKE::VECTOR_STRING,
																"Reactions for which the reaction rates are required (reactions start from 1, units in kmol/m3/s)",
																false));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ClassGrouping", 
																OpenSMOKE::SINGLE_PATH, 
																"Path to the file containing the list of species and families to carry out the class grouping", 
																false) );
		}
	};
}

#endif	/* OpenSMOKE_Grammar_OnTheFlyPostProcessing_H */

