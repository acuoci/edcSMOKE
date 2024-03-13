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
|   Copyright(C) 2022  Alberto Cuoci                                      |
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

// STL 
#include <vector>
#include <iostream>
#include <string>

namespace OpenSMOKE
{

	void MulticomponentTransportLibrary::Initialize(const std::string file_transport, const std::string file_thermo, const std::vector<std::string> names, const double* M)
	{
		// Default values
		exact_ = true;
		Fickr_ = 1;
		FickM_ = 0;
		FickThr_ = 0.1;
		Soretr_ = 1;
		SwapPolicy_ = -2;
		use_thermal_conductivity_ = true;

		// Read transport file
		std::cout << " * Multicomponent diffusion (MuTLib): Reading transport file..." << std::endl;
		if (!kineticdata_.ReadFile(file_transport.c_str()))
		{
			std::cout << "Error reading transport file" << std::endl;
			exit(0);
		}

		// Read thermodynamic file
		std::cout << " * Multicomponent diffusion (MuTLib): Reading thermodynamic file..." << std::endl;
		if (!thermodata_.ReadFile(file_thermo.c_str()))
		{
			std::cout << "Error reading thermodynamic file" << std::endl;
			exit(0);
		}

		// Set species names and properties
		std::cout << " * Multicomponent diffusion (MuTLib): Adding species to the mixture..." << std::endl;
		N_ = names.size();
		for (unsigned int i = 0; i < N_; i++)
		{
			const std::string name = names[i];
			const int indexK = kineticdata_.CheckName(name);
			const int indexT = thermodata_.CheckName(name);
			if (indexK == -1 || indexT == -1)
			{
				std::cout << "   Error, species " << name << " not found" << std::endl;
				exit(0);
			}

			species_.AddSpecie(name, M[i], 1. / static_cast<double>(N_), indexK, indexT);
		}

		//Sets up the transport class
		std::cout << " * Multicomponent diffusion (MuTLib): Setup transport..." << std::endl;
		species_.temp = 298.15;		// temperature in K
		species_.press = 101325.;	// pressure in Pa
		transport_.SetUp(&species_, &kineticdata_, &thermodata_);
	}

	void MulticomponentTransportLibrary::Initialize(const OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMapXML, const OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXML)
	{
		// Default values
		exact_ = true;
		Fickr_ = 1;
		FickM_ = 0;
		FickThr_ = 0.1;
		Soretr_ = 1;
		SwapPolicy_ = -2;
		use_thermal_conductivity_ = true;

		// Check if raw transport parameters are available
		if (transportMapXML.is_raw_transport_parameters_available() == false)
			OpenSMOKE::FatalErrorMessage("The kinetic mechanism was preprocessd without exporting the raw transport parameters of species (Raw-Transport-Parameters leaf in kinetics.xml)");

		// Read transport data
		std::cout << " * Multicomponent diffusion (MuTLib): Importing transport file..." << std::endl;
		{
			std::vector<int> shape_factor(thermodynamicsMapXML.NumberOfSpecies());
			std::vector<double> epsilon_over_kb(thermodynamicsMapXML.NumberOfSpecies());
			std::vector<double> sigma(thermodynamicsMapXML.NumberOfSpecies());
			std::vector<double> mu(thermodynamicsMapXML.NumberOfSpecies());
			std::vector<double> alfa(thermodynamicsMapXML.NumberOfSpecies());
			std::vector<double> zRot298(thermodynamicsMapXML.NumberOfSpecies());
			for (unsigned int i = 0; i < thermodynamicsMapXML.NumberOfSpecies(); i++)
			{
				shape_factor[i] = static_cast<int>(transportMapXML.raw_transport_parameters()[i][0]);
				epsilon_over_kb[i] = transportMapXML.raw_transport_parameters()[i][1];
				sigma[i] = transportMapXML.raw_transport_parameters()[i][2];
				mu[i] = transportMapXML.raw_transport_parameters()[i][3];
				alfa[i] = transportMapXML.raw_transport_parameters()[i][4];
				zRot298[i] = transportMapXML.raw_transport_parameters()[i][5];
			}

			kineticdata_.ImportCoefficients(	thermodynamicsMapXML.NamesOfSpecies(),
												shape_factor, epsilon_over_kb, sigma, mu, alfa, zRot298);
		}

		// Read thermodynamic data
		std::cout << " * Multicomponent diffusion (MuTLib): Importing thermodynamic file..." << std::endl;
		{
			std::vector<std::vector<double>> CoeffLowT;
			std::vector<std::vector<double>> CoeffHighT;
			std::vector<std::vector<double>> TRange;
			for (unsigned int i = 0; i < thermodynamicsMapXML.NumberOfSpecies(); i++)
			{
				std::vector<double> low(7);
				thermodynamicsMapXML.NASA_LowT(i, low.data());
				CoeffLowT.push_back(low);

				std::vector<double> high(7);
				thermodynamicsMapXML.NASA_HighT(i, high.data());
				CoeffHighT.push_back(high);

				std::vector<double> t(3);
				thermodynamicsMapXML.NASA_TRange(i, t.data());
				TRange.push_back(t);
			}

			thermodata_.ImportCoefficients(thermodynamicsMapXML.NamesOfSpecies(), CoeffLowT, CoeffHighT, TRange);
		}

		// Set species names and properties
		std::cout << " * Multicomponent diffusion (MuTLib): Adding species to the mixture..." << std::endl;
		N_ = thermodynamicsMapXML.NamesOfSpecies().size();
		for (unsigned int i = 0; i < N_; i++)
		{
			const std::string name = thermodynamicsMapXML.NamesOfSpecies()[i];
			const int indexK = kineticdata_.CheckName(name);
			const int indexT = thermodata_.CheckName(name);
			if (indexK == -1 || indexT == -1)
			{
				std::cout << "   Error, species " << name << " not found" << std::endl;
				exit(0);
			}

			species_.AddSpecie(name, thermodynamicsMapXML.MWs()[i], 1. / static_cast<double>(N_), indexK, indexT);
		}

		//Sets up the transport class
		std::cout << " * Multicomponent diffusion (MuTLib): Setup transport..." << std::endl;
		species_.temp = 298.15;		// temperature in K
		species_.press = 101325.;	// pressure in Pa
		transport_.SetUp(&species_, &kineticdata_, &thermodata_);
	}

	void MulticomponentTransportLibrary::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_MulticomponentTransportLibrary grammar;
		dictionary.SetGrammar(grammar);

		// Read type
		if (dictionary.CheckOption("@SwapPolicy") == true)
		{
			std::string policy;
			dictionary.ReadString("@SwapPolicy", policy);

			if (policy == "most-abundant")
				SetSwapPolicy(-2);
			else if (policy == "last")
				SetSwapPolicy(-1);
			else
			{
				bool found = false;
				for (unsigned int i = 0; i < N_; i++)
					if (species_.GetName(i) == policy)
					{
						SetSwapPolicy(i);
						found = true;
						break;
					}
				if (found == false)
					OpenSMOKE::FatalErrorMessage("The species name indicated in @SwapPolicy does not exist. Please check.");
			}
		}

		// Read Neumann steps for Fick
		{
			int value = 1;
			if (dictionary.CheckOption("@FickNeumannSteps") == true)
				dictionary.ReadInt("@FickNeumannSteps", value);
			Fickr_ = static_cast<unsigned int>(value);
		}

		// Read Neumann steps for Soret
		{
			int value = 1;
			if (dictionary.CheckOption("@SoretNeumannSteps") == true)
				dictionary.ReadInt("@SoretNeumannSteps", value);
			Soretr_ = static_cast<unsigned int>(value);
		}

		// Read exact flag
		{
			bool flag = true;
			if (dictionary.CheckOption("@Exact") == true)
				dictionary.ReadBool("@Exact", flag);
			exact_ = flag;
		}

		// Use thermal conductivity
		{
			bool flag = false;
			if (dictionary.CheckOption("@UseMulticomponentThermalConductivity") == true)
				dictionary.ReadBool("@UseMulticomponentThermalConductivity", flag);
			use_thermal_conductivity_ = flag;
		}

	}

	void MulticomponentTransportLibrary::SetupMixture(const double T, const double P, const double* X)
	{
		species_.temp = T;		// temperature in K
		species_.press = P;		// pressure in Pa
		for (unsigned int i = 0; i < N_; i++)
			species_.SetMolFrac(i, X[i]);

		transport_.Calculate(exact_, Fickr_, FickM_, FickThr_, Soretr_, SwapPolicy_);
	}

	void MulticomponentTransportLibrary::GetFickFluxes(const Eigen::VectorXd& nablaX, Eigen::VectorXd& jfick)
	{
		jfick = transport_.GetFickFluxes(nablaX);
	}


	void MulticomponentTransportLibrary::GetSoretFluxes(const double nablaLogT, Eigen::VectorXd& jsoret)
	{
		jsoret = transport_.GetSoretFluxes(nablaLogT);
	}


	void MulticomponentTransportLibrary::SetExact(const bool flag)
	{
		exact_ = flag;
	}


	void MulticomponentTransportLibrary::SetFickNeumannSteps(const unsigned int n)
	{
		Fickr_ = n;
	}


	void MulticomponentTransportLibrary::SetSoretNeumannSteps(const unsigned int n)
	{
		Soretr_ = n;
	}


	void MulticomponentTransportLibrary::SetSwapPolicy(const int policy)
	{
		SwapPolicy_ = policy;
	}


	void MulticomponentTransportLibrary::SetUserThermalConductivity(const bool flag)
	{
		use_thermal_conductivity_ = flag;
	}

}
