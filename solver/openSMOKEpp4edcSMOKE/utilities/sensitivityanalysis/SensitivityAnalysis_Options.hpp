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

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
	class Grammar_SensitivityAnalysis_Options : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type", 
																OpenSMOKE::SINGLE_STRING, 
																"Type of sensitivity analysis: arrhenius-parameters | kinetic-constants", 
																false) );			

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SubSteps", 
																OpenSMOKE::SINGLE_INT, 
																"Number of substeps when performing sensitivity analysis (default 2)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Species", 
																OpenSMOKE::VECTOR_STRING, 
																"List of species for which the sensitivity coefficients will be written", 
																true) );	

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DenseSolver",
																OpenSMOKE::SINGLE_STRING,
																"Dense solver: Eigen",
																false,
																"none",
																"none", 
																"@SparseSolver"));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DenseFullPivoting", 
																OpenSMOKE::SINGLE_BOOL, 
																"Full pivoting vs Partial pivoting (LU decomposition)", 
																false) );

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SparseSolver",
																OpenSMOKE::SINGLE_STRING,
																"Sparse solver: none | EigenBiCGSTAB | EigenGMRES | EigenDGMRES |EigenSparseLU | Pardiso | SuperLUSerial | UMFPack | LIS",
																false,
																"none",
																"none",
																"@DenseSolver"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SparsePreconditioner",
																OpenSMOKE::SINGLE_STRING,
																"Preconditioner: diagonal | ILUT",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SparsePreconditionerDropTolerance",
																OpenSMOKE::SINGLE_DOUBLE,
																"Preconditioner drop tolerance (default 1e-6)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SparsePreconditionerFillFactor",
																OpenSMOKE::SINGLE_INT,
																"Preconditioner fill factor (default 10)",
																false));
		}
	};

	SensitivityAnalysis_Options::SensitivityAnalysis_Options()
	{
		dense_solver_type_ = OpenSMOKE::SOLVER_DENSE_EIGEN;
		dense_decomposition_type_ = OpenSMOKE::DENSE_DECOMPOSITION_FULL_PIVOTING_LU;

		sparse_solver_type_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU;
		sparse_preconditioner_type_ = OpenSMOKE::PRECONDITIONER_SPARSE_ILUT;
		sparse_preconditioner_drop_tolerance_ = 1.e-6;
		sparse_preconditioner_fill_factor_ = 10;
		

		sensitivity_type_ = PhysicalConstants::SENSITIVITY_FREQUENCY_FACTOR;
		number_of_substeps_ = 2;
		list_of_species_.resize(0);			
	}

	void SensitivityAnalysis_Options::TransferOptions(OpenSMOKE::SensitivityMap& map)
	{
		map.SetNumberOfSubSteps(number_of_substeps_);
		map.SetSensitivityType(sensitivity_type_);
		map.SetDenseSolverType(dense_solver_type_);
		map.SetDenseDecompositionType(dense_decomposition_type_);
		map.SetSparseSolverType(sparse_solver_type_);
		map.SetSparsePreconditionerType(sparse_preconditioner_type_);
		map.SetSparsePreconditionerDropTolerance(sparse_preconditioner_drop_tolerance_);
		map.SetSparsePreconditionerFillFactor(sparse_preconditioner_fill_factor_);
	}

	void SensitivityAnalysis_Options::SetSensitivityType(const PhysicalConstants::sensitivity_type type)
	{
		sensitivity_type_ = type;
	}

	void SensitivityAnalysis_Options::SetDenseSolverType(const DenseSolverType type)
	{
		dense_solver_type_ = type;
		sparse_solver_type_ = OpenSMOKE::SOLVER_SPARSE_NONE;
	}

	void SensitivityAnalysis_Options::SetSparseSolverType(const SparseSolverType type)
	{
		sparse_solver_type_ = type;
		dense_solver_type_ = OpenSMOKE::SOLVER_DENSE_NONE;
	}

	void SensitivityAnalysis_Options::SetDenseDecompositionType(const DenseDecompositionType type)
	{
		dense_decomposition_type_ = type;
	}

	void SensitivityAnalysis_Options::SetSparsePreconditionerType(const SparsePreconditionerType type)
	{
		sparse_preconditioner_type_ = type;
	}

	void SensitivityAnalysis_Options::SetNumberOfSubSteps(const unsigned int number_of_substeps)
	{
		number_of_substeps_ = number_of_substeps;
	}

	void SensitivityAnalysis_Options::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_SensitivityAnalysis_Options grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@Type") == true)
		{
			std::string value;
			dictionary.ReadString("@Type", value);
			if (value == "arrhenius-parameters")
				sensitivity_type_ = PhysicalConstants::SENSITIVITY_FREQUENCY_FACTOR;
			else if (value == "kinetic-constants")
				sensitivity_type_ = PhysicalConstants::SENSITIVITY_KINETIC_CONSTANT;
			else
				OpenSMOKE::FatalErrorMessage("Unknown sensitivity analysis type: " + value );
		}

		if (dictionary.CheckOption("@DenseSolver") == true)
		{
			std::string value;
			dictionary.ReadString("@DenseSolver", value);
			if (value == "Eigen")
				dense_solver_type_ = SOLVER_DENSE_EIGEN;
			else
				OpenSMOKE::FatalErrorMessage("Unknown dense solver type: " + value );

			sparse_solver_type_ = OpenSMOKE::SOLVER_SPARSE_NONE;
		}

		if (dictionary.CheckOption("@DenseFullPivoting") == true)
		{
			bool value;
			dictionary.ReadBool("@DenseFullPivoting", value);
			if (value == true)
				dense_decomposition_type_ = DENSE_DECOMPOSITION_FULL_PIVOTING_LU;
			else
				dense_decomposition_type_ = DENSE_DECOMPOSITION_PARTIAL_PIVOTING_LU;
		}

		if (dictionary.CheckOption("@SparseSolver") == true)
		{
			std::string name;
			dictionary.ReadString("@SparseSolver", name);

			if (name == "EigenBiCGSTAB")		sparse_solver_type_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB;
			else if (name == "EigenGMRES")		sparse_solver_type_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES;
			else if (name == "EigenDGMRES")		sparse_solver_type_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES;
			else if (name == "EigenSparseLU")	sparse_solver_type_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU;
			else if (name == "Pardiso")
			{
				#if OPENSMOKE_USE_MKL == 0
					FatalErrorMessage("OpenSMOKE++ was built without the MKL support. Please select a different @SparseSolver option");
				#endif
					sparse_solver_type_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO;
			}
			else if (name == "SuperLUSerial")
			{
				#if OPENSMOKE_USE_SUPERLU_SERIAL == 0
					FatalErrorMessage("OpenSMOKE++ was built without the SuperLU support. Please select a different @SparseSolver option");
				#endif
					sparse_solver_type_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL;
			}
			else if (name == "UMFPack")
			{
				#if OPENSMOKE_USE_UMFPACK == 0
					FatalErrorMessage("OpenSMOKE++ was built without the UMFPACK support. Please select a different @SparseSolver option");
				#endif
					sparse_solver_type_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK;
			}
			else if (name == "LIS")
			{
				#if OPENSMOKE_USE_LIS == 0
					FatalErrorMessage("OpenSMOKE++ was built without the LIS support. Please select a different @SparseSolver option");
				#endif
					sparse_solver_type_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS;
			}
			else
			{
				FatalErrorMessage("Unknown Sparse Solver type: " + name);
			}

			dense_solver_type_ = OpenSMOKE::SOLVER_DENSE_NONE;
		}

		if (dictionary.CheckOption("@SparsePreconditioner") == true)
		{
			std::string name;
			dictionary.ReadString("@SparsePreconditioner", name);

			if (name == "diagonal")		sparse_preconditioner_type_ = OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL;
			else if (name == "ILUT")	sparse_preconditioner_type_ = OpenSMOKE::PRECONDITIONER_SPARSE_ILUT;

			else FatalErrorMessage("Unknown Sparse Preconditioner type: " + name);
		}

		if (dictionary.CheckOption("@SparsePreconditionerDropTolerance") == true)
			dictionary.ReadDouble("@SparsePreconditionerDropTolerance", sparse_preconditioner_drop_tolerance_);

		if (dictionary.CheckOption("@SparsePreconditionerFillFactor") == true)
			dictionary.ReadInt("@SparsePreconditionerFillFactor", sparse_preconditioner_fill_factor_);

		if (dictionary.CheckOption("@SubSteps") == true)
			dictionary.ReadInt("@SubSteps", number_of_substeps_);

		if (dictionary.CheckOption("@Species") == true)
			dictionary.ReadOption("@Species", list_of_species_);
	}

}
