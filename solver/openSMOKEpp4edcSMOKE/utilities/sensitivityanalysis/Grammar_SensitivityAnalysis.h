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

#ifndef OpenSMOKE_Grammar_SensitivityAnalysis_H
#define	OpenSMOKE_Grammar_SensitivityAnalysis_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{

	class Grammar_SensitivityAnalysis : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type", 
																OpenSMOKE::SINGLE_STRING, 
																"Type of sensitivity analysis: arrhenius-parameters | kinetic-constants", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FullPivoting", 
																OpenSMOKE::SINGLE_BOOL, 
																"Full pivoting vs Partial pivoting (LU decomposition)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@LinearAlgebra", 
																OpenSMOKE::SINGLE_STRING, 
																"Linear algebra package: Eigen | BzzMath", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SubSteps", 
																OpenSMOKE::SINGLE_INT, 
																"Number of substeps when performing sensitivity analysis (default 2)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Species", 
																OpenSMOKE::VECTOR_STRING, 
																"List of species for which the sensitivity coefficients will be written", 
																false) );															
		}
	};
}

#endif	/* OpenSMOKE_Grammar_SensitivityAnalysis_H */

