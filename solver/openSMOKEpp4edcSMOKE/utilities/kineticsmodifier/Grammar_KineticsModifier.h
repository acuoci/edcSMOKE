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
|   Copyright(C) 2019  Alberto Cuoci                                      |
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

#ifndef OpenSMOKE_Grammar_KineticsModifier_H
#define	OpenSMOKE_Grammar_KineticsModifier_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
	class Grammar_KineticsModifier : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@A", 
																OpenSMOKE::VECTOR_STRING, 
																"List of reaction indices and frequency factors (e.g. 3 1e5 12 1e8)", 
																false) );
			
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EoverR", 
																OpenSMOKE::VECTOR_STRING, 
																"List of reaction indices and activation energies (e.g. 3 20000. 9 15000.)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@n", 
																OpenSMOKE::VECTOR_STRING, 
																"List of reaction indices and temperature exponents (e.g. 3 1.8 12 0.7)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ThirdBody", 
																OpenSMOKE::VECTOR_STRING, 
																"List of reaction indices, species and third-body efficiencies (e.g. 3 H2O 2.3  12 CO2 0.9)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FrozenSpecies", 
																OpenSMOKE::VECTOR_STRING, 
																"List of species to be kept frozen during the simulation (e.g. BIN1A BIN2B)", 
																false) );
		}
	};
}

#endif	/* OpenSMOKE_Grammar_KineticsModifier_H */

