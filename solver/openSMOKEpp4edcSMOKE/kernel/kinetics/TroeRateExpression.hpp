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

namespace OpenSMOKE
{

	TroeRateExpression::TroeRateExpression()
	{
	}
	
	TroeRateExpression::TroeRateExpression(const TroeRateExpression& orig)
	{
	}
	
	//TroeRateExpression::~TroeRateExpression()
	//{
	//}

	void TroeRateExpression::ErrorMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  TroeRateExpression"	<< std::endl;
		std::cout << "Error:  " << message							<< std::endl;
		std::cout << "Press a key to continue... "					<< std::endl;
		getchar();
		exit(-1);
	}

	void TroeRateExpression::WarningMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  TroeRateExpression"	<< std::endl;
		std::cout << "Warning:  "	<< message						<< std::endl;
		std::cout << "Press a key to continue... "					<< std::endl;
		getchar();
	}

	void TroeRateExpression::Setup(const std::vector<double>& troe_coefficients, const std::vector<double>& low_coefficients, const std::vector<double>& inf_coefficients)
	{
		troe_coefficients_ = troe_coefficients;
		low_coefficients_ = low_coefficients;
		inf_coefficients_ = inf_coefficients;
	}

	double TroeRateExpression::KineticConstant(const double T, const double P)
	{
		return 0.;
	}

	void TroeRateExpression::WriteStatus(std::ostream& fOutput) const
	{
		fOutput << std::setw(9) << " "; fOutput << "Troe Parameters" << std::endl;
		fOutput << std::setw(9) << " "; fOutput << std::scientific << "a    " << troe_coefficients_[0] << std::endl;
		fOutput << std::setw(9) << " "; fOutput << std::scientific << "T*** " << troe_coefficients_[1] << std::endl;
		fOutput << std::setw(9) << " "; fOutput << std::scientific << "T*   " << troe_coefficients_[2] << std::endl;

		if (troe_coefficients_.size() == 3)
			fOutput << std::setw(9) << " "; fOutput << std::scientific << "T**  " << 0. << std::endl;

		if (troe_coefficients_.size() == 4)
			fOutput << std::setw(9) << " "; fOutput << std::scientific << "T**  " << troe_coefficients_[3] << std::endl;

		fOutput << std::setw(9) << " "; fOutput << "Inf coefficients" << std::endl;
		fOutput << std::setw(9) << " "; fOutput << std::scientific << "A    " << inf_coefficients_[0] << std::endl;
		fOutput << std::setw(9) << " "; fOutput << std::scientific << "n    " << inf_coefficients_[1] << std::endl;
		fOutput << std::setw(9) << " "; fOutput << std::scientific << "E    " << inf_coefficients_[2] << std::endl;

		fOutput << std::setw(9) << " "; fOutput << "Low coefficients" << std::endl;
		fOutput << std::setw(9) << " "; fOutput << std::scientific << "A    " << low_coefficients_[0] << std::endl;
		fOutput << std::setw(9) << " "; fOutput << std::scientific << "n    " << low_coefficients_[1] << std::endl;
		fOutput << std::setw(9) << " "; fOutput << std::scientific << "E    " << low_coefficients_[2] << std::endl;
	}

	void TroeRateExpression::ReadFromASCIIFile(std::istream& fInput)
	{
		// TODO
	}

	void TroeRateExpression::WriteShortSummaryOnASCIIFile(std::ostream& fOutput) const
	{
		// TODO
	}

	void TroeRateExpression::FitttingArrheniusLaw(std::ostream& fOutput)
	{
		// TODO
	}
}
