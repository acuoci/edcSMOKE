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
|   This file is part of OpenSMOKE++ Suite.                               |
|                                                                         |
|   Copyright(C) 2015, 2014, 2013  Alberto Cuoci                          |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_Grammar_Grid1D_H
#define OpenSMOKE_Grammar_Grid1D_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"
#include <Eigen/Dense>

namespace OpenSMOKE
{
	class Grammar_Grid1D : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type",
				OpenSMOKE::SINGLE_STRING,
				"Type of grid: centered | database | equispaced | liquid-pool",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Length",
				OpenSMOKE::SINGLE_MEASURE,
				"Length of computational domain",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InitialPoints",
				OpenSMOKE::SINGLE_INT,
				"Initial number of points",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Center",
				OpenSMOKE::SINGLE_MEASURE,
				"Center position (default: middle point)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Width",
				OpenSMOKE::SINGLE_MEASURE,
				"Width of flame (default: 0.2*length)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FixedPoint",
				OpenSMOKE::SINGLE_MEASURE,
				"Coordinate of fixed point (default: center of the grid)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaxPoints",
				OpenSMOKE::SINGLE_INT,
				"Maximum number of allowed points (default 300)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaxAdaptivePoints",
				OpenSMOKE::SINGLE_INT,
				"Maximum number of points that can be added per grid adaptation (default 10)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@GradientCoefficient",
				OpenSMOKE::SINGLE_DOUBLE,
				"Controls the maximum gradient allowed between two points (default 0.1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CurvatureCoefficient",
				OpenSMOKE::SINGLE_DOUBLE,
				"Controls the maximum curvature allowed between two points (default 0.5)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Threshold",
				OpenSMOKE::SINGLE_DOUBLE,
				"Controls the minimum amount of each species to be considered active in the adaptive process (default: 1e-7)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RegridPoints",
				OpenSMOKE::SINGLE_INT,
				"Points used for regrid operation (default: 20)",
				false));
		}
	};
}

#endif /* OpenSMOKE_Grammar_Grid1D_H */
