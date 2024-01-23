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

#include "math/OpenSMOKEUtilities.h"
#include "TransportPropertiesMap_CHEMKIN.h"
#include "preprocessing/CollisionIntegralMatrices.hpp"

#if OPENSMOKE_USE_MKL == 1
#include "mkl.h"
#endif

namespace OpenSMOKE
{
	const double TransportPropertiesMap_CHEMKIN::threshold_ = 1.e-14;

	TransportPropertiesMap_CHEMKIN::TransportPropertiesMap_CHEMKIN(const unsigned int nSpecies)
	{
		this->nspecies_ = nSpecies;
		temperature_lambda_must_be_recalculated_ = true;
		temperature_eta_must_be_recalculated_ = true;
		temperature_gamma_must_be_recalculated_ = true;
		temperature_teta_must_be_recalculated_ = true;
		pressure_gamma_must_be_recalculated_ = true;

		this->T_ = this->P_ = 0.;
		this->T_old_ = this->P_old_ = 0.;

		this->species_bundling_ = false;
		is_lennard_jones_available_ = false;
		is_raw_transport_parameters_available_ = false;
		stefanmaxwell_skip_policy_ = STEFANMAXWELL_SKIP_POLICY_MOST_ABUNDANT;
		stefanmaxwell_skip_ = this->nspecies_ - 1;

		MemoryAllocation();
	}

	TransportPropertiesMap_CHEMKIN::TransportPropertiesMap_CHEMKIN(boost::property_tree::ptree& ptree)
	{
		temperature_lambda_must_be_recalculated_ = true;
		temperature_eta_must_be_recalculated_ = true;
		temperature_gamma_must_be_recalculated_ = true;
		temperature_teta_must_be_recalculated_ = true;
		pressure_gamma_must_be_recalculated_ = true;
		this->T_ = this->P_ = 0.;
		this->T_old_ = this->P_old_ = 0.;

		this->species_bundling_ = false;

		ImportSpeciesFromXMLFile(ptree);
		ImportViscosityModelFromXMLFile(ptree);
		MemoryAllocation();
		ImportCoefficientsFromXMLFile(ptree);

		stefanmaxwell_skip_policy_ = STEFANMAXWELL_SKIP_POLICY_MOST_ABUNDANT;
		stefanmaxwell_skip_ = this->nspecies_ - 1;
	}

	TransportPropertiesMap_CHEMKIN::TransportPropertiesMap_CHEMKIN(const TransportPropertiesMap_CHEMKIN& rhs)
	{
		CopyFromMap(rhs);
	}

	TransportPropertiesMap_CHEMKIN::~TransportPropertiesMap_CHEMKIN()
	{
		delete[] this->M;
		delete[] this->fittingLambda;
		delete[] this->fittingEta;
		delete[] this->fittingTeta;
		delete[] this->fittingGamma;
		delete[] this->fittingGammaSelfDiffusion;

		if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE)
		{
			delete[] MWRatio1over4;
			delete[] phi_eta_sup;
			delete[] phi_eta_inf;
			delete[] sqrtEta;
			delete[] usqrtEta;
		}
		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING)
		{
			delete[] sqrtMWRatio_inf;
			delete[] sqrtMWRatio_sup;
			delete[] sqrtMW;
		}

		delete[] sum_diffusion_coefficients;
		delete[] x_corrected;

		if (this->species_bundling_ == true)
		{
			delete[] this->bundling_fittingGammaSelfDiffusion_;
			delete[] this->bundling_fittingGamma_;
			delete[] this->bundling_gammaSpecies_;
			delete[] this->bundling_gammaSpeciesSelfDiffusion_;
		}

		iThermalDiffusionRatios_.clear();
	}

	void TransportPropertiesMap_CHEMKIN::SetTemperature(const double& T)
	{
		this->T_old_ = this->T_;
		this->T_ = T;

		const double epsilon = 0.;
		if (std::fabs(this->T_ - this->T_old_) / this->T_ > epsilon)
		{
			temperature_lambda_must_be_recalculated_ = true;
			temperature_eta_must_be_recalculated_ = true;
			temperature_gamma_must_be_recalculated_ = true;
			temperature_teta_must_be_recalculated_ = true;
		}
	}

	void TransportPropertiesMap_CHEMKIN::SetPressure(const double& P)
	{
		this->P_old_ = this->P_;
		this->P_ = P;

		if (std::fabs(this->P_ - this->P_old_) / this->P_ > 1.e-14)
		{
			pressure_gamma_must_be_recalculated_ = true;
		}
	}

	void TransportPropertiesMap_CHEMKIN::CopyFromMap(const TransportPropertiesMap_CHEMKIN& rhs)
	{
		this->nspecies_ = rhs.nspecies_;
		this->species_bundling_ = rhs.species_bundling_;

		MemoryAllocation();

		this->count_species_thermal_diffusion_ratios_ = rhs.count_species_thermal_diffusion_ratios_;
		this->sumK = rhs.sumK;
		this->iThermalDiffusionRatios_ = rhs.iThermalDiffusionRatios_;

		for (unsigned int i = 0; i < this->nspecies_; i++)
			this->M[i] = rhs.M[i];

		for (unsigned int i = 0; i < this->nspecies_ * 4; i++)
			this->fittingLambda[i] = rhs.fittingLambda[i];

		for (unsigned int i = 0; i < this->nspecies_ * 4; i++)
			this->fittingEta[i] = rhs.fittingEta[i];

		for (unsigned int i = 0; i < this->nspecies_ * 4; i++)
			this->fittingGammaSelfDiffusion[i] = rhs.fittingGammaSelfDiffusion[i];

		for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2 * 4; i++)
			this->fittingGamma[i] = rhs.fittingGamma[i];

		this->fittingTeta = new double[4 * count_species_thermal_diffusion_ratios_*this->nspecies_];
		for (unsigned int i = 0; i < this->count_species_thermal_diffusion_ratios_*(this->nspecies_) * 4; i++)
			this->fittingTeta[i] = rhs.fittingTeta[i];

		if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE)
		{
			for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2; i++)
				this->MWRatio1over4[i] = rhs.MWRatio1over4[i];

			for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2; i++)
				this->phi_eta_sup[i] = rhs.phi_eta_sup[i];

			for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2; i++)
				this->phi_eta_inf[i] = rhs.phi_eta_inf[i];

			for (unsigned int i = 0; i < this->nspecies_; i++)
				this->sqrtEta[i] = rhs.sqrtEta[i];

			for (unsigned int i = 0; i < this->nspecies_; i++)
				this->usqrtEta[i] = rhs.usqrtEta[i];
		}
		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING)
		{
			for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2; i++)
				this->sqrtMWRatio_inf[i] = rhs.sqrtMWRatio_inf[i];

			for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2; i++)
				this->sqrtMWRatio_sup[i] = rhs.sqrtMWRatio_sup[i];

			for (unsigned int i = 0; i < this->nspecies_; i++)
				this->sqrtMW[i] = rhs.sqrtMW[i];
		}

		if (this->species_bundling_ == true)
		{
			this->bundling_number_groups_ = rhs.bundling_number_groups_;
			this->bundling_reference_species_ = rhs.bundling_reference_species_;
			this->bundling_groups_ = rhs.bundling_groups_;
			this->bundling_species_group_ = rhs.bundling_species_group_;

			this->bundling_fittingGammaSelfDiffusion_ = new double[4 * this->nspecies_];
			this->bundling_fittingGamma_ = new double[this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2 * 4];
			this->bundling_gammaSpecies_ = new double[this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2];
			this->bundling_gammaSpeciesSelfDiffusion_ = new double[this->bundling_number_groups_];

			for (unsigned int i = 0; i < 4 * this->nspecies_; i++)
				this->bundling_fittingGammaSelfDiffusion_[i] = rhs.bundling_fittingGammaSelfDiffusion_[i];

			for (unsigned int i = 0; i < this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2 * 4; i++)
				this->bundling_fittingGamma_[i] = rhs.bundling_fittingGamma_[i];

			for (unsigned int i = 0; i < this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2; i++)
				this->bundling_gammaSpecies_[i] = rhs.bundling_gammaSpecies_[i];

			for (unsigned int i = 0; i < this->bundling_number_groups_; i++)
				this->bundling_gammaSpeciesSelfDiffusion_[i] = rhs.bundling_gammaSpeciesSelfDiffusion_[i];
		}

		this->tetaSpecies_.resize(this->nspecies_*iThermalDiffusionRatios_.size());
		this->tetaSpecies_.setZero();
	}

	void TransportPropertiesMap_CHEMKIN::MemoryAllocation()
	{
		// viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING;
		// viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_MATHUR_SAXENA;
		// viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE; (default)

		sum_threshold_ = 1. + threshold_*double(this->nspecies_);

		M = new double[this->nspecies_];
		fittingLambda = new double[this->nspecies_ * 4];
		fittingEta = new double[this->nspecies_ * 4];
		fittingGammaSelfDiffusion = new double[this->nspecies_ * 4];
		fittingGamma = new double[this->nspecies_*(this->nspecies_ - 1) / 2 * 4];

		// fittingTeta is not inizialized here, because the number of species for which the thermal diffusion ratios is 
		// important is not yet known

		sumK.resize(this->nspecies_);
		sumK.setZero();

		if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE)
		{
			MWRatio1over4 = new double[this->nspecies_*(this->nspecies_ - 1) / 2];
			phi_eta_sup = new double[this->nspecies_*(this->nspecies_ - 1) / 2];
			phi_eta_inf = new double[this->nspecies_*(this->nspecies_ - 1) / 2];

			sqrtEta = new double[this->nspecies_];
			usqrtEta = new double[this->nspecies_];


			// Provisional implementation for 1+M model
			// It can be improved
			{
				mMWRatio1over4.resize(this->nspecies_);
				mphi_eta_sup.resize(this->nspecies_);
				for (unsigned int i=0;i<this->nspecies_;i++)
				{
					mMWRatio1over4[i].resize(this->nspecies_);
					mphi_eta_sup[i].resize(this->nspecies_);
				}
			}
		}
		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING)
		{
			sqrtMWRatio_inf = new double[this->nspecies_*(this->nspecies_ - 1) / 2];
			sqrtMWRatio_sup = new double[this->nspecies_*(this->nspecies_ - 1) / 2];
			sqrtMW = new double[this->nspecies_];
		}

		sum_diffusion_coefficients = new double[this->nspecies_];
		x_corrected = new double[this->nspecies_];

		this->lambdaSpecies_.resize(this->nspecies_);
		this->lambdaSpecies_.setZero();

		this->etaSpecies_.resize(this->nspecies_);
		this->etaSpecies_.setZero();

		this->gammaSpecies_.resize(this->nspecies_*(this->nspecies_ - 1) / 2);
		this->gammaSpecies_.setZero();

		// [Bertrand Naud] 1+M model
		{
			this->gammaSpeciesSelfDiffusion_.resize(this->nspecies_);
			this->gammaSpeciesSelfDiffusion_.setZero();

			epsij_over_kb_.resize(this->nspecies_);
			for (unsigned int i=0;i<this->nspecies_;i++)
			{
				epsij_over_kb_[i].resize(this->nspecies_);
			}
			fParker_298_.resize(this->nspecies_);
		}

		M_0511_ = new double[this->nspecies_];
		M_0489_ = new double[this->nspecies_];


		// Thermal diffusion ratios are allocated after reading the properties
	}

	void TransportPropertiesMap_CHEMKIN::SetCoefficients(const unsigned int j, const double* coefficients)
	{
		// Molecular weight
		M[j] = coefficients[0];

		// Thermal conductivity
		{
			unsigned int i = j * 4;
			fittingLambda[i++] = coefficients[1];
			fittingLambda[i++] = coefficients[2];
			fittingLambda[i++] = coefficients[3];
			fittingLambda[i++] = coefficients[4];
		}

		// Dynamic viscosity
		{
			unsigned int i = j * 4;
			fittingEta[i++] = coefficients[5];
			fittingEta[i++] = coefficients[6];
			fittingEta[i++] = coefficients[7];
			fittingEta[i++] = coefficients[8];
		}

		// Mass diffusion coefficients
		{
			unsigned int countFitting = 0;
			for (unsigned int i = 1; i <= j; i++)
				countFitting += this->nspecies_ - i;
			countFitting *= 4;
			unsigned int countCoefficients = 9 + 4 * (j + 1);
			for (unsigned int k = j + 1; k < this->nspecies_; k++)
			for (unsigned int i = 1; i <= 4; i++)
				fittingGamma[countFitting++] = coefficients[countCoefficients++];
		}

		// Thermal diffusion coefficients
		{
			if (M[j] <= PhysicalConstants::MAX_MW_THERMALDIFFUSION_RATIOS)
				iThermalDiffusionRatios_.push_back(j + 1);

			unsigned int i = j * 4 * this->nspecies_;
			unsigned int count = 8 + 4 * this->nspecies_ + 1;

			for (unsigned int k = 0; k < this->nspecies_; k++)
			for (unsigned int m = 1; m <= 4; m++)
				fittingTeta[i++] = coefficients[count++];
		}

		if (j == this->nspecies_ - 1)
			CompleteInitialization();
	}

	void TransportPropertiesMap_CHEMKIN::ImportCoefficientsFromASCIIFile(std::istream& fInput)
	{
		{
			fInput >> count_species_thermal_diffusion_ratios_;
			fittingTeta = new double[4 * count_species_thermal_diffusion_ratios_*this->nspecies_];
		}

		std::string tag;
		for (unsigned int j = 0; j < this->nspecies_; j++)
		{
			// Molecular weight
			fInput >> M[j];

			// Thermal conductivity
			{
				unsigned int i = j * 4;
				fInput >> fittingLambda[i++];
				fInput >> fittingLambda[i++];
				fInput >> fittingLambda[i++];
				fInput >> fittingLambda[i++];
			}

			// Dynamic viscosity
			{
				unsigned int i = j * 4;
				fInput >> fittingEta[i++];
				fInput >> fittingEta[i++];
				fInput >> fittingEta[i++];
				fInput >> fittingEta[i++];
			}

			// [Bertrand Naud] Mass Self diffusion coefficients
			{
				unsigned int i = j * 4;
				fInput >> fittingGammaSelfDiffusion[i++];
				fInput >> fittingGammaSelfDiffusion[i++];
				fInput >> fittingGammaSelfDiffusion[i++];
				fInput >> fittingGammaSelfDiffusion[i++];
			}

			// Mass diffusion coefficients
			{
				unsigned int countFitting = 0;
				for (unsigned int i = 1; i <= j; i++)
					countFitting += this->nspecies_ - i;
				countFitting *= 4;
				for (unsigned int k = j + 1; k < this->nspecies_; k++)
				for (unsigned int i = 1; i <= 4; i++)
					fInput >> fittingGamma[countFitting++];
			}

			// Thermal diffusion coefficients
			{
				if (M[j] <= PhysicalConstants::MAX_MW_THERMALDIFFUSION_RATIOS)
				{
					iThermalDiffusionRatios_.push_back(j + 1);

					std::size_t i = 4 * this->nspecies_*(iThermalDiffusionRatios_.size() - 1);

					for (unsigned int k = 0; k < this->nspecies_; k++)
					for (unsigned int m = 1; m <= 4; m++)
						fInput >> fittingTeta[i++];
				}
			}
		}

		// Check that everything is read properly
		fInput >> tag;
		CheckForFatalError(tag == "E");

		// Complete the initialization
		CompleteInitialization();
	}

	void TransportPropertiesMap_CHEMKIN::ImportLennardJonesCoefficientsFromASCIIFile(std::istream& fInput)
	{
		is_lennard_jones_available_ = true;
		mu_.resize(this->nspecies_);
		sigma_.resize(this->nspecies_);
		epsilon_over_kb_.resize(this->nspecies_);

		for (unsigned int j = 0; j < this->nspecies_; j++)
		{
			fInput >> mu_[j];					// [kg]
			fInput >> sigma_[j];				// [m]
			fInput >> epsilon_over_kb_[j];		// [K]
		}
	}

	void TransportPropertiesMap_CHEMKIN::ImportRawTransportParametersFromASCIIFile(std::istream& fInput)
	{
		is_raw_transport_parameters_available_ = true;
		raw_transport_parameters_.resize(this->nspecies_);

		for (unsigned int j = 0; j < this->nspecies_; j++)
		{
			std::vector<double> transport(6);
			fInput >> transport[0];	transport[0] = static_cast<int>(transport[0]);
			fInput >> transport[1];
			fInput >> transport[2];
			fInput >> transport[3];
			fInput >> transport[4];
			fInput >> transport[5];
			raw_transport_parameters_[j] = transport;	
		}

		// [Bertrand Naud] 1+M diffusion model
		for(unsigned int j=0;j<this->nspecies_;j++)
		{
			const double eps_over_kB_j = raw_transport_parameters_[j][1];
			const double sigma_j = raw_transport_parameters_[j][2];
			const double dipmu_j = raw_transport_parameters_[j][3];
			const double alpha_j = raw_transport_parameters_[j][4];
			for(unsigned int k=0;k<=j;k++)
			{
				const double eps_over_kB_k = raw_transport_parameters_[k][1];
				const double sigma_k = raw_transport_parameters_[k][2];
				const double dipmu_k = raw_transport_parameters_[k][3];
				const double alpha_k = raw_transport_parameters_[k][4];
				if (dipmu_j < 1.e-20 && dipmu_k > 1.e-20)
				{
					//--- k is polar, j is nonpolar
					const double alphaStar = alpha_j / boost::math::pow<3>(sigma_j);
					const double muStar_k = dipmu_k / std::sqrt(eps_over_kB_k*boost::math::pow<3>(sigma_k));
					const double xi = 1. + 0.25*alphaStar * boost::math::pow<2>(muStar_k) * std::sqrt(eps_over_kB_k/eps_over_kB_j);
					//sigmaij[k][j] = 0.5 * (sigma_j+sigma_k) * std::pow(xsi, -1./6.);
					//sigmaij[j][k] = r_sgm[k][j];
					epsij_over_kb_[k][j] = std::sqrt(eps_over_kB_j*eps_over_kB_k) * xi*xi;
					epsij_over_kb_[j][k] = epsij_over_kb_[k][j];
				}
				else if (dipmu_j > 1.e-20 && dipmu_k < 1.e-20)
				{
					//--- j is polar, k is nonpolar
					const double alphaStar = alpha_k / boost::math::pow<3>(sigma_k);
					const double muStar_j = dipmu_j / std::sqrt(eps_over_kB_j*boost::math::pow<3>(sigma_j));
					const double xi = 1. + 0.25*alphaStar * boost::math::pow<2>(muStar_j) * std::sqrt(eps_over_kB_j/eps_over_kB_k);
					//sigmaij[k][j] = 0.5 * (sigma_j+sigma_k) * std::pow(xsi, -1./6.);
					//sigmaij[j][k] = r_sgm[k][j];
					epsij_over_kb_[k][j] = std::sqrt(eps_over_kB_j*eps_over_kB_k) * xi*xi;
					epsij_over_kb_[j][k] = epsij_over_kb_[k][j];
				}
				else
				{
					//--- normal case, either both polar or both nonpolar
					//sigmaij[k][j] = 0.5 * (sigma_j + sigma_k);
					//sigmaij[j][k] = r_sgm[k][j];
					epsij_over_kb_[k][j] = std::sqrt(eps_over_kB_j*eps_over_kB_k);
					epsij_over_kb_[j][k] = epsij_over_kb_[k][j];
				}
			}

			// Parker-Brau-Jonkman correction
			const double aux298 = std::sqrt(eps_over_kB_j / 298.);
			fParker_298_[j] = 1. + aux298 * ( PhysicalConstants::ZROTA + aux298 * (PhysicalConstants::ZROTB + aux298 * PhysicalConstants::ZROTC));
		}
	}

	void TransportPropertiesMap_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)
	{
		std::cout << " * Reading transport properties from XML file..." << std::endl;

		// Main transport coefficients
		{
			std::stringstream stream;
			stream.str( ptree.get< std::string >("opensmoke.Transport") );
			ImportCoefficientsFromASCIIFile(stream);
		}
		
		// Lennard-Jones parametes (not stricly needed)
		{
			std::stringstream stream;
			stream.str( ptree.get< std::string >("opensmoke.Lennard-Jones") );
			ImportLennardJonesCoefficientsFromASCIIFile(stream);
		}

		// Raw transport parameters (needed for using the multicomponent transport library mutlib)
		{
			//if (ptree.get_optional<bool>("opensmoke.Raw-Transport-Parameters").is_initialized() == true)
			// Now it is autmatically written by the OpenSMOKE++ PreProcessor, no need of any check
			{
				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.Raw-Transport-Parameters"));
				ImportRawTransportParametersFromASCIIFile(stream);
			}
		}
	}

	unsigned int index_from_mtrix_to_vector(const unsigned int n, const unsigned int i, const unsigned int j)
	{
		unsigned int index = 0;
		for (unsigned int k = 1; k <= i; k++)
			index += (n - k);

		index += j - i - 1;

		return index;
	}

	void TransportPropertiesMap_CHEMKIN::ImportSpeciesBundlingFromXMLFile(boost::property_tree::ptree& ptree, const double epsilon)
	{
		std::cout << " * Reading species bundling from XML file..." << std::endl;

		// Read self diffusion coefficients
		{
			std::stringstream stream;
			stream.str( ptree.get< std::string >("opensmoke.SpeciesBundling.SelfDiffusion") );

			this->bundling_fittingGammaSelfDiffusion_ = new double[4 * this->nspecies_];

			for (unsigned int j = 0; j < this->nspecies_; j++)
			{
				{
					unsigned int i = j * 4;
					stream >> this->bundling_fittingGammaSelfDiffusion_[i++];
					stream >> this->bundling_fittingGammaSelfDiffusion_[i++];
					stream >> this->bundling_fittingGammaSelfDiffusion_[i++];
					stream >> this->bundling_fittingGammaSelfDiffusion_[i++];
				}
			}
		}
		
		// Read bundling
		bool iFound = false;
		BOOST_FOREACH( boost::property_tree::ptree::value_type const& node, ptree.get_child( "opensmoke.SpeciesBundling" ) ) 
		{

			boost::property_tree::ptree subtree = node.second;  

			if( node.first == "Bundling" ) 
			{
				const double epsilon_double = subtree.get<double>("<xmlattr>.epsilon");
				
				std::cout << "Looking for eps=" << epsilon << " - Found eps=" << epsilon_double << std::endl;

				if (iFound == false && std::fabs(epsilon_double - epsilon)/epsilon < 1.e-4)
				{
					this->species_bundling_ = true;

					iFound = true;
					std::cout << "   Maximum error (epsilon): " << epsilon_double << std::endl;

					// Number of groups
					{
						this->bundling_number_groups_ = static_cast<unsigned int>(subtree.get<double>("NumberOfGroups"));
					}

					// Reference species
					{
						std::stringstream stream;
						stream.str( subtree.get< std::string >("ReferenceSpecies") );

						unsigned int n;
						stream >> n;
						this->bundling_reference_species_.resize(n);
						for (unsigned int i = 0; i<n; i++)
							stream >> this->bundling_reference_species_[i];
					}

					// Species in groups
					{
						std::stringstream stream;
						stream.str( subtree.get< std::string >("SpeciesInGroups") );

						this->bundling_groups_.resize(this->bundling_number_groups_);
						for (unsigned k = 0; k < this->bundling_number_groups_; k++)
						{
							unsigned int n;
							stream >> n;
							this->bundling_groups_[k].resize(n);
							for (unsigned int i = 0; i<n; i++)
								stream >> this->bundling_groups_[k][i];
						}
					}

					// Group of species
					{
						std::stringstream stream;
						stream.str( subtree.get< std::string >("GroupOfSpecies") );

						unsigned int n;
						stream >> n;
						this->bundling_species_group_.resize(n);
						for (unsigned int i = 0; i<n; i++)
							stream >> this->bundling_species_group_[i];
					}

					std::cout << "   Number of groups: " << this->bundling_number_groups_ << "/" << this->nspecies_ << std::endl;
					std::cout << "   List of reference species: " << std::endl;
					for (unsigned int i = 0; i < this->bundling_number_groups_; i++)
						std::cout << "    - " << this->bundling_reference_species_[i] << " (" << this->bundling_groups_[i].size() << ")" << std::endl;

					// Diffusion coefficients
					bundling_sum_diffusion_coefficients_.resize(this->bundling_number_groups_);
					bundling_sum_x_groups_.resize(this->bundling_number_groups_);

					this->bundling_fittingGamma_ = new double[this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2 * 4];
					this->bundling_gammaSpecies_ = new double[this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2];
					this->bundling_gammaSpeciesSelfDiffusion_ = new double[this->bundling_number_groups_];

					unsigned int countFitting = 0;
					for (unsigned int i = 0; i < this->bundling_number_groups_; i++)
					for (unsigned int j = i + 1; j < this->bundling_number_groups_; j++)
					{
						unsigned int index_i = this->bundling_reference_species_[i];
						unsigned int index_j = this->bundling_reference_species_[j];
						if (this->bundling_reference_species_[i] > this->bundling_reference_species_[j])
						{
							index_j = this->bundling_reference_species_[i];
							index_i = this->bundling_reference_species_[j];
						}

						unsigned int index = index_from_mtrix_to_vector(this->nspecies_, index_i, index_j) * 4;

						for (unsigned int i = 0; i < 4; i++)
							this->bundling_fittingGamma_[countFitting + i] = fittingGamma[index + i];
						countFitting += 4;

					}
				}
			}
		}

		if (iFound == false)
			ErrorMessage("TransportPropertiesMap_CHEMKIN::ImportSpeciesBundlingFromXMLFile", "The requested epsilon for species bundling was not found");
	}


	void TransportPropertiesMap_CHEMKIN::ImportSpeciesFromXMLFile(boost::property_tree::ptree& ptree)
	{
		this->nspecies_ = ptree.get<unsigned int>("opensmoke.NumberOfSpecies");  

		std::stringstream stream;
		stream.str( ptree.get< std::string >("opensmoke.NamesOfSpecies") );

		std::vector<std::string> names(this->nspecies_);  
		for(unsigned int i=0;i<this->nspecies_;i++)
			stream >> names[i]; 

		index_H2O_ = 0;
		index_CO2_ = 0;
		index_CO_  = 0;
		index_CH4_ = 0;
		index_NH3_ = 0;
		index_NO_  = 0;
		index_N2O_ = 0;


		for (unsigned int i = 0; i < this->nspecies_; i++)
		{
			if (names[i] == "H2O" || names[i] == "h2o")	index_H2O_ = i + 1;
			if (names[i] == "CO2" || names[i] == "co2")	index_CO2_ = i + 1;
			if (names[i] == "CO"  || names[i] == "co")	index_CO_  = i + 1;
			if (names[i] == "CH4" || names[i] == "ch4")	index_CH4_ = i + 1;
			if (names[i] == "NH3" || names[i] == "nh3")	index_NH3_ = i + 1;
			if (names[i] == "NO"  || names[i] == "no")	index_NO_  = i + 1;
			if (names[i] == "N2O" || names[i] == "n2o")	index_N2O_ = i + 1;
		}
	}


	void TransportPropertiesMap_CHEMKIN::ImportViscosityModelFromXMLFile(boost::property_tree::ptree& ptree)
	{
		viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE;

		boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.ViscosityModel");
		if (child)
		{
			const std::string model = ptree.get<std::string>("opensmoke.ViscosityModel"); 

			std::cout << " * User defined viscosity model: " << model << std::endl;

			if (model == "Wilke")			viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE;
			else if (model == "Herning")		viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING;
			else if (model == "MathurSaxena")	viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_MATHUR_SAXENA;
			else ErrorMessage("TransportPropertiesMap_CHEMKIN::ImportViscosityModelFromXMLFile", "Error in reading the viscosity model.");
	}
	}

	void TransportPropertiesMap_CHEMKIN::CompleteInitialization()
	{
		if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE)
		{
			unsigned int i=0;
			for(unsigned int k=0;k<this->nspecies_;k++)
				for(unsigned int j=k+1;j<this->nspecies_;j++)
				{
					phi_eta_sup[i] = 1./std::sqrt( 8.*(1.+M[k]/M[j]) );
					phi_eta_inf[i] = 1./std::sqrt( 8.*(1.+M[j]/M[k]) );
					i++;
				}

			i=0;
			for(unsigned int k=0;k<this->nspecies_;k++)
				for(unsigned int j=k+1;j<this->nspecies_;j++)
					MWRatio1over4[i++] = std::sqrt(std::sqrt(M[j]/M[k]));


			// Provisional implementation for 1+M model
			// It can be improved
			{
				for(unsigned int k=0;k<this->nspecies_;k++)
					for(unsigned int j=0;j<this->nspecies_;j++)
					{
						mMWRatio1over4[j][k] = std::sqrt(std::sqrt(M[j]/M[k]));
						mphi_eta_sup[j][k] = 1./std::sqrt( 8.*(1.+M[k]/M[j]) );
					}
			}
		}

		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING)
		{
			for(unsigned int k=0;k<this->nspecies_;k++)
				sqrtMW[k] = std::sqrt(M[k]);

			unsigned int i=0;
			for(unsigned int k=0;k<this->nspecies_;k++)
				for(unsigned int j=k+1;j<this->nspecies_;j++)
				{
					sqrtMWRatio_sup[i] = std::sqrt(M[j]/M[k]);
					sqrtMWRatio_inf[i++] = std::sqrt(M[k]/M[j]);
				}
		}
		
		this->tetaSpecies_.resize(this->nspecies_*boost::lexical_cast<unsigned int>(iThermalDiffusionRatios_.size()));
		this->tetaSpecies_.setZero();
		
		// Molecular weight functions to be used in the Kuo's model for thermal diffusion coefficients
		for (unsigned int k = 0; k < this->nspecies_; k++)
		{
			M_0511_[k] = std::pow(M[k], 0.511);
			M_0489_[k] = std::pow(M[k], 0.489);
		}
	}

	inline void TransportPropertiesMap_CHEMKIN::lambda()
	{
        if (temperature_lambda_must_be_recalculated_ == true)
        {
			const double logT=std::log(T_);
			const double logT2=logT*logT;
			const double logT3=logT*logT2;

			#if OPENSMOKE_USE_MKL == 0

				const double* k = fittingLambda;
				for (unsigned int j=0;j<this->nspecies_;j++)
				{
					double sum = (*k++);
							sum+= (*k++)*logT;
							sum+= (*k++)*logT2;
							sum+= (*k++)*logT3;
					this->lambdaSpecies_(j)=std::exp(sum);
				}

			#elif OPENSMOKE_USE_MKL == 1
	
				const double* k = fittingLambda;
				for (unsigned int j=0;j<this->nspecies_;j++)
				{
					double sum = (*k++);
							sum+= (*k++)*logT;
							sum+= (*k++)*logT2;
							sum+= (*k++)*logT3;
					this->lambdaSpecies_(j)=sum;
				}
				vdExp( this->nspecies_, this->lambdaSpecies_.data(), this->lambdaSpecies_.data() );
	
			#endif

            temperature_lambda_must_be_recalculated_ = false;
        }
	}

	inline void TransportPropertiesMap_CHEMKIN::eta()
	{
        if (temperature_eta_must_be_recalculated_ == true)
        {
			const double logT=std::log(T_);
			const double logT2=logT*logT;
			const double logT3=logT*logT2;

			#if OPENSMOKE_USE_MKL == 0

				const double* k = fittingEta;
				for (unsigned int j=0;j<this->nspecies_;j++)
				{
					double sum = (*k++);
						   sum+= (*k++)*logT;
						   sum+= (*k++)*logT2;
						   sum+= (*k++)*logT3;
					this->etaSpecies_(j)=std::exp(sum);
				}

			#elif OPENSMOKE_USE_MKL == 1
	
				const double* k = fittingEta;
				for (unsigned int j=0;j<this->nspecies_;j++)
				{
					double sum = (*k++);
						   sum+= (*k++)*logT;
						   sum+= (*k++)*logT2;
						   sum+= (*k++)*logT3;
					this->etaSpecies_(j)=sum;
				}
				vdExp( this->nspecies_, this->etaSpecies_.data(), this->etaSpecies_.data() );
	
			#endif

			temperature_eta_must_be_recalculated_ = false;
		}
	}
	
	inline void TransportPropertiesMap_CHEMKIN::eta(const unsigned int Mmain, const std::vector<unsigned int>& kreorder)
	{
		// In the 1+M formulation we cannot exploit:
		// 1. Intel MKL/oneAPI optimized exponential function
		// 2. Reuse of previously calculated individual viscosities if the temperature does not change
		// We should try to see if there is a way to optimize the code

		const double logT=std::log(T_);
		const double logT2=logT*logT;
		const double logT3=logT*logT2;

		for (unsigned int i=0;i<1+Mmain;i++)
		{
			const unsigned int k = kreorder[i];
			const unsigned int kk = k*4;
			const double sum = fittingEta[kk] + fittingEta[kk+1]*logT + fittingEta[kk+2]*logT2 + fittingEta[kk+3]*logT3;
			this->etaSpecies_(k)=std::exp(sum);
		}
	}
	
	inline void TransportPropertiesMap_CHEMKIN::teta()
	{
		if (temperature_teta_must_be_recalculated_ == true)
		{
			const double T2=T_*T_;
			const double T3=T_*T2;
		
			// Thermal diffusion ratios evaluation
			for (unsigned int i=1;i<=iThermalDiffusionRatios_.size();i++)
			{
				unsigned int j = 4*(i-1)*this->nspecies_;
				unsigned int m = (i-1)*this->nspecies_;

				double* ptFitting = &fittingTeta[j];

				for (unsigned int k=1;k<=this->nspecies_;k++)
				{
					double sum  = (*ptFitting++);
						   sum += (*ptFitting++)*T_;
						   sum += (*ptFitting++)*T2;
						   sum += (*ptFitting++)*T3;
					this->tetaSpecies_[m++] = sum; 

				}
			}
                
			temperature_teta_must_be_recalculated_ = false;
		}
	}

	inline void TransportPropertiesMap_CHEMKIN::gamma()
	{
		Mijkl_must_be_recalculated_ = false;

		if ( temperature_gamma_must_be_recalculated_ == false && pressure_gamma_must_be_recalculated_  == true )
		{
			Mijkl_must_be_recalculated_ = true;

            		const double multiplier = P_/P_old_;
			double* D = this->gammaSpecies_.data();
			for (unsigned int j=0;j<this->nspecies_;j++)
			for (unsigned int k=j+1;k<this->nspecies_;k++)
			{
				*D++ *= multiplier;
			}

			// [Bertrand Naud] Self diffusion coefficients
			double* Dself = this->gammaSpeciesSelfDiffusion_.data();
			for (unsigned int j = 0; j < this->nspecies_; j++)
			{
				*Dself++ *= multiplier;
			}
			
			pressure_gamma_must_be_recalculated_ = false;
		}
            
        if ( temperature_gamma_must_be_recalculated_ == true || pressure_gamma_must_be_recalculated_ == true )
        {

			Mijkl_must_be_recalculated_ = true;

			// Only the upper hals of this matrix is evaluated (the main diagonalis not evaluated)
			// Indeed the matrix is symmetric and the mixture rule is able to exploit this kind of symmetry

			const double P_bar = P_/100000.;
			const double logT=std::log(T_);
			const double logT2=logT*logT;
			const double logT3=logT*logT2;

			#if OPENSMOKE_USE_MKL == 0

			const double* d = fittingGamma;
			double* D = this->gammaSpecies_.data();
			for (unsigned int j=1;j<=this->nspecies_;j++)
				for (unsigned int k=j+1;k<=this->nspecies_;k++)
				{
					double sum = (*d++);
							sum+= (*d++)*logT;
							sum+= (*d++)*logT2;
							sum+= (*d++)*logT3;
					*D++ = P_bar/std::exp(sum);
				}

			#elif OPENSMOKE_USE_MKL == 1
	
			const double lnP_bar = std::log(P_bar);
			const double* d = fittingGamma;
			double* D = this->gammaSpecies_.data();
			for (unsigned int j=0;j<this->nspecies_;j++)
				for (unsigned int k=j+1;k<this->nspecies_;k++)
				{
					double sum = (*d++);
							sum+= (*d++)*logT;
							sum+= (*d++)*logT2;
							sum+= (*d++)*logT3;
					*D++ = lnP_bar-sum;
				}

			vdExp( this->nspecies_*(this->nspecies_-1)/2, this->gammaSpecies_.data(), this->gammaSpecies_.data() );

			#endif

			// [Bertrand Naud] Self diffusion coefficients
			const double* dself = fittingGammaSelfDiffusion;
			double* Dself = this->gammaSpeciesSelfDiffusion_.data();
			for (unsigned int j = 0; j <this->nspecies_; j++)
			{
				double sum = (*dself++);
				sum += (*dself++)*logT;
				sum += (*dself++)*logT2;
				sum += (*dself++)*logT3;
				*Dself++ = P_bar / std::exp(sum);
			}

            temperature_gamma_must_be_recalculated_ = false;
            pressure_gamma_must_be_recalculated_ = false;
        }
	}

	inline void TransportPropertiesMap_CHEMKIN::bundling_gamma()
	{
		if (temperature_gamma_must_be_recalculated_ == false && pressure_gamma_must_be_recalculated_ == true )
		{
			const double multiplier = P_ / P_old_;

			double* D = this->bundling_gammaSpecies_;
			for (unsigned int j = 0; j<this->bundling_number_groups_; j++)
			for (unsigned int k = j + 1; k<this->bundling_number_groups_; k++)
			{
				*D++ *= multiplier;
			}

			double* Dself = this->bundling_gammaSpeciesSelfDiffusion_;
			for (unsigned int j = 0; j < this->bundling_number_groups_; j++)
			{
				*Dself++ *= multiplier;
			}

			pressure_gamma_must_be_recalculated_ = false;
		}

		if ( temperature_gamma_must_be_recalculated_ == true || pressure_gamma_must_be_recalculated_ == true)
		{
			// Only the upper hals of this matrix is evaluated (the main diagonalis not evaluated)
			// Indeed the matrix is symmetric and the mixture rule is able to exploit this kind of symmetry

			const double P_bar = P_ / 100000.;
			const double logT = std::log(T_);
			const double logT2 = logT*logT;
			const double logT3 = logT*logT2;

			// Self diffusion coefficients
			double* Dself = this->bundling_gammaSpeciesSelfDiffusion_;
			for (unsigned int j = 0; j <this->bundling_number_groups_; j++)
			{
				unsigned int index = this->bundling_reference_species_[j];
				const double* d = &this->bundling_fittingGammaSelfDiffusion_[index * 4];

				double sum = (*d++);
				sum += (*d++)*logT;
				sum += (*d++)*logT2;
				sum += (*d++)*logT3;

				*Dself++ = P_bar / std::exp(sum);
			}


			#if OPENSMOKE_USE_MKL == 0

			const double* d = this->bundling_fittingGamma_;
			double* D = this->bundling_gammaSpecies_;
			for (unsigned int j = 1; j <= this->bundling_number_groups_; j++)
			for (unsigned int k = j + 1; k <= this->bundling_number_groups_; k++)
			{
				double sum = (*d++);
				sum += (*d++)*logT;
				sum += (*d++)*logT2;
				sum += (*d++)*logT3;
				*D++ = P_bar / std::exp(sum);
			}

			#elif OPENSMOKE_USE_MKL == 1

			const double lnP_bar = std::log(P_bar);
			const double* d = this->bundling_fittingGamma_;
			double* D = this->bundling_gammaSpecies_;
			for (unsigned int j = 0; j<this->bundling_number_groups_; j++)
			for (unsigned int k = j + 1; k<this->bundling_number_groups_; k++)
			{
				double sum = (*d++);
				sum += (*d++)*logT;
				sum += (*d++)*logT2;
				sum += (*d++)*logT3;
				*D++ = lnP_bar - sum;
			}

			vdExp(this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2, this->bundling_gammaSpecies_, this->bundling_gammaSpecies_);

			#endif

			temperature_gamma_must_be_recalculated_ = false;
			pressure_gamma_must_be_recalculated_ = false;
		}
	}

	double TransportPropertiesMap_CHEMKIN::lambdaMix(const double* moleFractions)
	{
		// Calcolo della conducibilita della miscela
		// Formula di Mathur, Todor, Saxena - Molecular Physics 52:569 (1967)

		const double sum1 = Dot(this->nspecies_, moleFractions, this->lambdaSpecies_.data());
		const double sum2 = UDot(this->nspecies_, moleFractions, this->lambdaSpecies_.data());
		const double lambdamix = 0.50 * ( sum1 + 1./sum2 );
		
		return lambdamix;
	}


	// [Bertrand Naud] 1+M model
	double TransportPropertiesMap_CHEMKIN::lambdaMult(const double* moleFractions, const unsigned int Mmain, const std::vector<unsigned int>& kreorder)
	{
		double* a10_k = this->a10_.data();
		double condtr = 0.;
		for (unsigned int kk=0;kk<1+Mmain;kk++)
		{
			const unsigned int k=kreorder[kk];
			condtr -= moleFractions[k] * (*a10_k);
			a10_k++;
		}

		double* M1001_jk = this->M1001_.data();
		double* inv_M0101_k = this->inv_M0101_.data();
		double condin = 0.;
		for (unsigned int kk=0;kk<1+Mmain;kk++)
		{
			const unsigned int k=kreorder[kk];
			double sum1 = 0.;
			double* a10_j = this->a10_.data();
			for (unsigned int jj=0;jj<1+Mmain;jj++)
			{
				sum1 += (*M1001_jk)*(*a10_j);
				M1001_jk++;
				a10_j++;
			}
			condin -= moleFractions[k] * (*inv_M0101_k) * (1. - sum1);
			inv_M0101_k++;
		}
		const double lambdamult = 6.25 * P_ / T_ * (condtr+condin);

		return lambdamult;
	}

	double TransportPropertiesMap_CHEMKIN::etaMix(const double* moleFractions)
	{
		for (unsigned int k = 0; k < this->nspecies_; k++)
			sumK(k) = moleFractions[k];

		// Wilke - Journal of Chemical Physics 18:517 (1950)
		// Modified by Bird, Stewart, Lightfoot - Transport phenomena (1960)
		// Available in: Reid, Prausnitz, Poling - The properties of gases and liquids, p. 407
		if(viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE)
		{
			#if OPENSMOKE_USE_MKL == 0
			
				for(unsigned int k=0;k<this->nspecies_;k++)
				{
					sqrtEta[k] = std::sqrt(this->etaSpecies_(k));
					usqrtEta[k] = 1./sqrtEta[k];
				}

			#elif OPENSMOKE_USE_MKL == 1
			
				vdSqrt(this->nspecies_, this->etaSpecies_.data(), sqrtEta);
				vdInv(this->nspecies_, sqrtEta, usqrtEta);
			
			#endif

			const double* ptMWRatio1over4=MWRatio1over4;
			const double* ptphi_eta_sup=phi_eta_sup;
			const double* ptphi_eta_inf=phi_eta_inf;

			for (unsigned int k = 0; k < this->nspecies_; k++)
				for (unsigned int j = k + 1; j < this->nspecies_; j++)
				{
					double delta_phi = sqrtEta[k] * usqrtEta[j] * (*ptMWRatio1over4++);	// F.(49)
					sumK(k) += moleFractions[j] * (*ptphi_eta_sup++)*(1. + delta_phi)*(1. + delta_phi);
					sumK(j) += moleFractions[k] * (*ptphi_eta_inf++)*(1. + 1. / delta_phi)*(1. + 1. / delta_phi);
				}

			double etamix = 0.;
			for(unsigned int k=0;k<this->nspecies_;k++)
				etamix += moleFractions[k]*this->etaSpecies_(k)/sumK(k);				// F.(48)

			return etamix;
		}

		// Herning and Zipperer
		// Available in: Reid, Prausnitz, Poling - The properties of gases and liquids, p. 410
		// This is a simplified formula, less accurate, but much faster
		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING)
		{
			double sum = 0;
			for (unsigned int k = 0; k < this->nspecies_; k++)
				sum += moleFractions[k] * sqrtMW[k];

			double etamix = 0.;
			for (unsigned int k = 0; k < this->nspecies_; k++)
				etamix += moleFractions[k] * this->etaSpecies_(k) * sqrtMW[k];				// F.(48)
			etamix /= sum;

			return etamix;
		}
		// Mathur and Saxena
		// Molecular Physics 52:569 (1967)
		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_MATHUR_SAXENA)
		{
			const double sum1 = Dot(this->nspecies_, moleFractions, this->etaSpecies_.data());
			const double sum2 = UDot(this->nspecies_, moleFractions, this->etaSpecies_.data());
			const double etamix = 0.50 * (sum1 + 1. / sum2);
			return etamix;
		}
	
		return 0.;
	}

	// New function added for 1+M model
	double TransportPropertiesMap_CHEMKIN::etaMix(const double* moleFractions, const unsigned int Mmain, const std::vector<unsigned int>& kreorder)
	{
		// Wilke - Journal of Chemical Physics 18:517 (1950)
		// Modified by Bird, Stewart, Lightfoot - Transport phenomena (1960)
		// Available in: Reid, Prausnitz, Poling - The properties of gases and liquids, p. 407
		if(viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE)
		{
			// Precalculation of square roots of viscosities of species (to save computational time)
			for(unsigned int i=0;i<1+Mmain;i++)
			{
				const unsigned int k=kreorder[i];

				sqrtEta[k] = std::sqrt(this->etaSpecies_(k));
				usqrtEta[k] = 1./sqrtEta[k];
			}
			
			// Calculation of mixture viscosity using 1+M species only
			double etamix = 0.;
			for (unsigned int kk=0; kk<1+Mmain; kk++)
			{
				const unsigned int k=kreorder[kk];

				double sumXPhi = 0.;
				for (unsigned int ll=0; ll<1+Mmain; ll++)
				{
					const unsigned int l=kreorder[ll];

					const double one_plus_delta_phi = 1. + sqrtEta[k] * usqrtEta[l] * (mMWRatio1over4[l][k]);

					// BUG
					//sumXPhi += moleFractions[l] * (one_plus_delta_phi*one_plus_delta_phi) *  mphi_eta_sup[k][l];
					sumXPhi += moleFractions[l] * (one_plus_delta_phi*one_plus_delta_phi) *  mphi_eta_sup[l][k];
				}

				etamix += moleFractions[k]*this->etaSpecies_(k)/sumXPhi;
			}

			return etamix;
		}
		else
		{
			ErrorMessage("TransportPropertiesMap_CHEMKIN", "the 1+M model is compatible with the Wilke model only");
		}
	
		return 0.;
	}

	void TransportPropertiesMap_CHEMKIN::gammaMix(double* gammamix, const double* moleFractions)
	{	
		// Reset
		for(unsigned int k=0;k<this->nspecies_;k++)
			sum_diffusion_coefficients[k] = 0.;

		// Adjust mole fractions
		for(unsigned int i=0;i<this->nspecies_;i++)
			x_corrected[i] = (moleFractions[i]+threshold_)/sum_threshold_;

		double MWmix = 0.;
		for(unsigned int i=0;i<this->nspecies_;i++)
			MWmix += x_corrected[i]*M[i];
		
		// a. Evaluating Mass Diffusion coefficients (mixture averaged)
		const double *d = this->gammaSpecies_.data();
		for(unsigned int k=0;k<this->nspecies_;k++)
			for(unsigned int j=k+1;j<this->nspecies_;j++)
			{
				sum_diffusion_coefficients[j] += x_corrected[k] * (*d);
				sum_diffusion_coefficients[k] += x_corrected[j] * (*d);
				d++;
			}
			
		// b. Evaluating Mass Diffusion coefficients (mixture averaged)
		for(unsigned int k=0;k<this->nspecies_;k++)
			gammamix[k] = (MWmix - x_corrected[k]*M[k]) / (MWmix*sum_diffusion_coefficients[k]);


		// TODO
		/*
		Eigen::MatrixXd GammaBar;
		StefanMaxwellGammaBar(GammaBar, moleFractions);

		for (unsigned int k = 0; k < this->nspecies_; k++)
			std::cout << k << " " << gammamix[k] << std::endl;

		Eigen::VectorXd jsm(this->nspecies_);	// Stefan-Maxwell fluxes (kg/m2/s)
		std::cout << "Sum: " << jsm.sum() << std::endl;

		getchar();
		*/
	}

	void TransportPropertiesMap_CHEMKIN::bundling_gammaMix(double* gammamix, const double* moleFractions)
	{
		// Reset
		for (unsigned int k = 0; k<this->bundling_number_groups_; k++)
			this->bundling_sum_diffusion_coefficients_[k] = 0.;

		// Adjust mole fractions
		for (unsigned int i = 0; i<this->nspecies_; i++)
			x_corrected[i] = (moleFractions[i] + threshold_) / sum_threshold_;

		double MWmix = 0.;
		for (unsigned int i = 0; i<this->nspecies_; i++)
			MWmix += x_corrected[i] * M[i];

		for (unsigned int k = 0; k < this->bundling_number_groups_; k++)
		{
			this->bundling_sum_x_groups_[k] = 0.;
			for (unsigned int j = 0; j < this->bundling_groups_[k].size(); j++)
				this->bundling_sum_x_groups_[k] += x_corrected[this->bundling_groups_[k][j]];
		}
			
		// a. Evaluating Mass Diffusion coefficients (mixture averaged)
		const double *d = this->bundling_gammaSpecies_;
		for (unsigned int k = 0; k<this->bundling_number_groups_; k++)
		for (unsigned int j = k + 1; j<this->bundling_number_groups_; j++)
		{
			this->bundling_sum_diffusion_coefficients_[j] += this->bundling_sum_x_groups_[k] * (*d);
			this->bundling_sum_diffusion_coefficients_[k] += this->bundling_sum_x_groups_[j] * (*d);
			d++;
		}

		// b. Evaluating Mass Diffusion coefficients (mixture averaged)
		for (unsigned int k = 0; k < this->nspecies_; k++)
		{
			const unsigned int group = this->bundling_species_group_[k];
			const double correction = (this->bundling_sum_x_groups_[group] - x_corrected[k])*this->bundling_gammaSpeciesSelfDiffusion_[group];

			gammamix[k] =	(MWmix - x_corrected[k] * M[k]) / 
							(MWmix* (this->bundling_sum_diffusion_coefficients_[group] + correction));
		}
	}


	// [Bertrand Naud] 1+M model
	void TransportPropertiesMap_CHEMKIN::gammaMult(double* gammamult, const double* moleFractions, const unsigned int Mmain, const std::vector<unsigned int>& kreorder)
	{
		//--- set binary diffusion coefficients
		Eigen::MatrixXd Dij(1+Mmain,this->nspecies_);
		Dij.setZero();
		const double* d = this->gammaSpecies_.data();
		double* dself = this->gammaSpeciesSelfDiffusion_.data();
		// i<j:   indice D(i,j) = indice D(0,1) + j + i*N - (i+1)*(i+2)/2
		// i>j:   indice D(j,i) = indice D(0,1) + i + j*N - (j+1)*(j+2)/2
		for (unsigned int ii=0; ii<1+Mmain; ii++)
		{
			const unsigned int i=kreorder[ii];
			const int indice_Dij_minus_j=i*(this->nspecies_) - (i+1)*(i+2)/2;
			for (unsigned int jj=0;jj<this->nspecies_;jj++)
			{
				const unsigned int j = kreorder[jj];
				if(i == j)
				{
					Dij(ii,jj) = 1. / *(dself+i);
				}
				else
				{
					// if(i < j) { Dij = 1. / *( d + j + i*(this->nspecies_) - (i+1)*(i+2)/2 ); }
					if(i < j) { Dij(ii,jj) = 1. / *( d + j + indice_Dij_minus_j ); }
					else      { Dij(ii,jj) = 1. / *( d + i + j*(this->nspecies_) - (j+1)*(j+2)/2 ); }
				}
			}
		}

		const unsigned int K=kreorder[0];

		if (Mmain == 0)
		{
			// Dilute species
			for (unsigned int i1=1;i1<this->nspecies_;i1++)
			{
				const unsigned int i=kreorder[i1];
				for (unsigned int j=0;j<i;j++)
				{
					unsigned int l = i*nspecies_ + j;
					gammamult[l] = 0.;
				}

				// j = i
				unsigned int l = i*nspecies_ + i;
				gammamult[l] = Dij(0,i1);

				for (unsigned int j=i+1;j<this->nspecies_;j++)
				{
					unsigned int l = i*nspecies_ + j;
					gammamult[l] = 0.;
				}
			}
		}
		else if (Mmain == 1)
		{

			// Store D_iK.c_ij (main species only)
			// used in thermodiffusion problem, in "solve_thermodiffusion_onePlusM" to approximate vector "a10"
			this->DiKcij_.resize(Mmain*Mmain);
			this->DiKcij_.setZero();
			double* DiKcij = this->DiKcij_.data();

			const unsigned int i=kreorder[1];
			unsigned int l = i*nspecies_ + i;
			gammamult[l] = Dij(0,1) / (1+moleFractions[i]*(M[i]/M[K]-1));
			(*DiKcij) = gammamult[l];

			// dilute species j
			for (unsigned int j1=2; j1<this->nspecies_; j1++)
			{
				const unsigned int j=kreorder[j1];
				l = i*nspecies_ + j;
				gammamult[l] = - moleFractions[i]*(M[j]/M[K]-Dij(0,1)/Dij(1,j1))
									/ ( (1+moleFractions[i]*(M[i]/M[K]-1)) * ((1-moleFractions[i])/Dij(0,j1) + moleFractions[i]/Dij(1,j1)) );
			}
			// main species j=K
			l = i*nspecies_ + K;
			gammamult[l] = 0.;

			// dilute species
			for (unsigned int i1=2;i1<this->nspecies_;i1++)
			{
				const unsigned int i=kreorder[i1];
				for (unsigned int j=0;j<i;j++)
				{
					l = i*nspecies_ + j;
					gammamult[l] = 0.;
				}

				// j = i
				l = i*nspecies_ + i;
				gammamult[l] = 1 / ((1-moleFractions[kreorder[1]])/Dij(0,i1) + moleFractions[kreorder[1]]/Dij(1,i1));
				for (unsigned int j=i+1;j<this->nspecies_;j++)
				{
					l = i*nspecies_ + j;
					gammamult[l] = 0.;
				}
			}
		}
		else
		{
			Eigen::MatrixXd inv_one_plus_A11(Mmain,Mmain);
			inv_one_plus_A11.setZero();
			Eigen::MatrixXd A12(Mmain,this->nspecies_ -Mmain-1);
			A12.setZero();
			Eigen::VectorXd inv_one_plus_A22(this->nspecies_ -Mmain-1);
			inv_one_plus_A22.setZero();

			// set [1+A_11] and [A_12]
			for (unsigned int i1=1;i1<1+Mmain;i1++)
			{
				const unsigned int i = kreorder[i1];
				const unsigned int ii = i1-1;
				//double Xi = moleFractions[i];
				double Xi = std::max(1e-30,moleFractions[i]);

				// [1+A_11]
				inv_one_plus_A11(ii,ii) = 1.;
				for (unsigned int j1=1;j1<1+Mmain;j1++)
				{
					const unsigned int j = kreorder[j1];
					const unsigned int jj = j1-1;
					double ratio_Dij = Dij(0,i1)/Dij(i1,j1);
					//double Xj = moleFractions[j];
					double Xj = std::max(1e-30,moleFractions[j]);
					inv_one_plus_A11(ii,ii) += Xj * (ratio_Dij - 1.);
					inv_one_plus_A11(ii,jj) += Xi * (M[j]/M[K] - ratio_Dij) * Dij(0,j1)/Dij(0,i1);
				}

				// [A_12] and [A_22]
				for (unsigned int j1=1+Mmain;j1<this->nspecies_;j1++)
				{
					const unsigned int j = kreorder[j1];
					const unsigned int jj = j1-Mmain-1;
					inv_one_plus_A22(jj) += Xi * (Dij(0,j1)/Dij(i1,j1) - 1.);
					A12(ii,jj) = Xi * (M[j]/M[K] - Dij(0,i1)/Dij(i1,j1)) * Dij(0,j1)/Dij(0,i1);
				}
			}
			// set 1/[1+A_22]
			for (unsigned int ii=0;ii<this->nspecies_ -Mmain-1;ii++)
			{
				inv_one_plus_A22(ii) = 1. / (1. + inv_one_plus_A22(ii));
			}


			// Invert MxM matrix [1+A_11]
			{ 
				// THIS COULD BE A SUBROUTINE CALLED "INVERT_MATRIX" 
				if (Mmain == 0)
				{
					// nothing to do
				}
				else if (Mmain == 1)
				{
					inv_one_plus_A11(0,0) = 1./inv_one_plus_A11(0,0);
				}
				else if (Mmain == 2)
				{
					double inv_det = 1./(inv_one_plus_A11(0,0)*inv_one_plus_A11(1,1) - inv_one_plus_A11(0,1)*inv_one_plus_A11(1,0));
					double aa =   inv_one_plus_A11(1,1) * inv_det;
					double bb = - inv_one_plus_A11(0,1) * inv_det;
					double cc = - inv_one_plus_A11(1,0) * inv_det;
					double dd =   inv_one_plus_A11(0,0) * inv_det;
					inv_one_plus_A11(0,0) = aa;
					inv_one_plus_A11(0,1) = bb;
					inv_one_plus_A11(1,0) = cc;
					inv_one_plus_A11(1,1) = dd;
				}
				else // Lapack matrix inversion
				{
					//--- Linpack: dgefa ------------------------------------------
					std::vector<int>ipvt(Mmain); 
					// find pivot
					int info = -1;
					int idamax = 0;
					for (unsigned int k=0;k<Mmain-1;k++)
					{
						idamax = 0;
						double dmax = std::abs(inv_one_plus_A11(k,k));
						for (unsigned int i=1;i<Mmain-k;i++)
						{
							double xmag = std::abs(inv_one_plus_A11(k+i,k));
							if (xmag > dmax)
							{
								idamax = i;
								dmax = xmag;
							}
						}
						int l = idamax + k;

						ipvt[k] = l;
						if (inv_one_plus_A11(l,k) != 0.)
						{
							// interchange if necessary
							if (l!=k)
							{
								double tt = inv_one_plus_A11(l,k);
								inv_one_plus_A11(l,k) = inv_one_plus_A11(k,k);
								inv_one_plus_A11(k,k) = tt;
							}
							// compute multipliers
							double tt = -1. / inv_one_plus_A11(k,k);
							for (unsigned int i=1;i<Mmain-k;i++)
							{
								inv_one_plus_A11(k+i,k) = tt*inv_one_plus_A11(k+i,k);
							}
							// row elimination with column indexing
							for (unsigned int j=k+1;j<Mmain;j++)
							{
								tt = inv_one_plus_A11(l,j);
								if (l!=k)
								{
									inv_one_plus_A11(l,j) = inv_one_plus_A11(k,j);
									inv_one_plus_A11(k,j) = tt;
								}
								for (unsigned int i=1;i<Mmain-k;i++)
								{
									inv_one_plus_A11(k+i,j) = inv_one_plus_A11(k+i,j) + tt*inv_one_plus_A11(k+i,k);
								}
							}
						}
						else
						{
							info=k;
						}
					}
					ipvt[Mmain-1] = Mmain-1;
					if (inv_one_plus_A11(Mmain-1,Mmain-1) == 0.) info = Mmain-1;
					if (info!=-1) std::cout << "Error inverting matrix (dgefa), info = " << info << std::endl;
					//-------------------------------------------------------------
					//--- Linpack: dgedi ------------------------------------------
					// compute inverse(u)
					for (unsigned int k=0;k<Mmain;k++)
					{
						inv_one_plus_A11(k,k) = 1. / inv_one_plus_A11(k,k);
						double tt = - inv_one_plus_A11(k,k);
						for (unsigned int i=0;i<k;i++)
						{
							inv_one_plus_A11(i,k) = tt*inv_one_plus_A11(i,k);
						}
						for (unsigned int j=k+1;j<Mmain;j++)
						{
							tt = inv_one_plus_A11(k,j);
							inv_one_plus_A11(k,j) = 0.;
							for (unsigned int i=0;i<=k;i++)
							{
								inv_one_plus_A11(i,j) = inv_one_plus_A11(i,j) + tt*inv_one_plus_A11(i,k);
							}
						}
					}

					// form inverse(u)*inverse(l)
					std::vector<double>work(Mmain); 
					for (unsigned int kb=1;kb<=Mmain-1;kb++)
					{
						int k = Mmain - kb - 1;
						for (unsigned int i=k+1;i<Mmain;i++)
						{
							work[i] = inv_one_plus_A11(i,k);
							inv_one_plus_A11(i,k) = 0.;
						}
						for (unsigned int j=k+1;j<Mmain;j++)
						{
							double tt = work[j];
							for (unsigned int i=0;i<Mmain;i++)
							{
								inv_one_plus_A11(i,k) = inv_one_plus_A11(i,k) + tt*inv_one_plus_A11(i,j);
							}
						}
						int l = ipvt[k];
						if (l!=k)
						{
							for (unsigned int i=0;i<Mmain;i++)
							{
								double dtemp = inv_one_plus_A11(i,k);
								inv_one_plus_A11(i,k) = inv_one_plus_A11(i,l);
								inv_one_plus_A11(i,l) = dtemp;
							}
						}
					}
				}
			} // Closing Invert matrix

			// Store D_iK.c_ij (main species only)
			// used in thermodiffusion problem, in "solve_thermodiffusion_onePlusM" to approximate vector "a10"
			this->DiKcij_.resize(Mmain*Mmain);
			this->DiKcij_.setZero();
			double* DiKcij = this->DiKcij_.data();

			// main species i
			for (unsigned int ii=0;ii<Mmain;ii++)
			{
				const unsigned int i1=ii+1;
				const unsigned int i=kreorder[i1];
				// main species j
				for (unsigned int jj=0;jj<Mmain;jj++)
				{
					const unsigned int j1=jj+1;
					const unsigned int j=kreorder[j1];
					unsigned int l = i*nspecies_ + j;
					gammamult[l] = Dij(0,i1) * inv_one_plus_A11(ii,jj);
					(*DiKcij) = gammamult[l];
					DiKcij++;
				}
				// dilute species j
				for (unsigned int jj=0; jj<this->nspecies_ -Mmain-1; jj++)
				{
					const unsigned int j1=jj+Mmain+1;
					const unsigned int j=kreorder[j1];
					unsigned int l = i*nspecies_ + j;
					gammamult[l] = 0.;
					for (unsigned int kk=0;kk<Mmain;kk++)
					{
						gammamult[l] -= Dij(0,i1) * inv_one_plus_A11(ii,kk) * A12(kk,jj);
					}
					gammamult[l] = gammamult[l] * inv_one_plus_A22(jj);
				}

		 		// main species j=K
				unsigned int l = i*nspecies_ + K;
				gammamult[l] = 0.;
			}

			// dilute species
			for (unsigned int ii=0;ii<this->nspecies_ -Mmain-1;ii++)
			{
				const unsigned int i1=ii+Mmain+1;
				const unsigned int i=kreorder[i1];
				for (unsigned int j=0;j<i;j++)
				{
					unsigned int l = i*nspecies_ + j;
					gammamult[l] = 0.;
				}

				// j = i
				unsigned int l = i*nspecies_ + i;
				gammamult[l] = Dij(0,i1) * inv_one_plus_A22(ii);

				for (unsigned int j=i+1;j<this->nspecies_;j++)
				{
					unsigned int l = i*nspecies_ + j;
					gammamult[l] = 0.;
				}
			}
		}

		// Main species K
		// BUG
		/*
		for (unsigned int ii=1;ii<1+Mmain;ii++)
		{
			const unsigned int i=kreorder[ii];
			// main species j
			double sum1 = 0.;
			for (unsigned int jj=1;jj<this->nspecies_;jj++)
			{
				const unsigned int j=kreorder[jj];
				unsigned int l = j*nspecies_ + i;
				sum1 += M[j]*gammamult[l];
			}
			unsigned int l = K*nspecies_ + i;
			gammamult[l] = - sum1 / M[K];
		}
		for (unsigned int ii=1+Mmain;ii<this->nspecies_;ii++)
		{
			const unsigned int i=kreorder[ii];
			unsigned int l = i*nspecies_ + i;
			gammamult[K*nspecies_+i] = - M[i]*gammamult[l] / M[K];
		}
		*/

		// D_Ki, main species i
		for (unsigned int ii=1;ii<1+Mmain;ii++)
		{
			const unsigned int i=kreorder[ii];
			double sum1 = 0.;
			for (unsigned int jj=1;jj<1+Mmain;jj++)
			{
				const unsigned int j=kreorder[jj];
				unsigned int l = j*nspecies_ + i;
				sum1 += M[j]*gammamult[l];
			}
			unsigned int l = K*nspecies_ + i;
			gammamult[l] = - sum1 / M[K];
		}
		// D_Ki, dilute species i
		for (unsigned int ii=1+Mmain;ii<this->nspecies_;ii++)
		{
			const unsigned int i=kreorder[ii];
			double sum1 = 0.;
			for (unsigned int jj=1;jj<1+Mmain;jj++)
			{
				const unsigned int j=kreorder[jj];
				unsigned int l = j*nspecies_ + i;
				sum1 += M[j]*gammamult[l];
			}
			unsigned int l = i*nspecies_ + i;
			sum1 += M[i]*gammamult[l];
			l = K*nspecies_ + i;
			gammamult[l] = - sum1 / M[K];
		}


		// D_KK
		unsigned int l = K*nspecies_ + K;
		gammamult[l] = 0.;
	}

	void TransportPropertiesMap_CHEMKIN::solve_thermodiffusion_onePlusM(const double* moleFractions, const double* cpspecies, const unsigned int Mmain, const std::vector<unsigned int>& kreorder)
	{
	//if (Mijkl_must_be_recalculated_ == true)
	//{
		//--- set binary diffusion coefficients
		Eigen::MatrixXd Dij(1+Mmain,1+Mmain);
		Dij.setZero();

		const double* d = this->gammaSpecies_.data();
		double* dself = this->gammaSpeciesSelfDiffusion_.data();
		// i<j:   indice D(i,j) = indice D(0,1) + j + i*N - (i+1)*(i+2)/2
		// i>j:   indice D(j,i) = indice D(0,1) + i + j*N - (j+1)*(j+2)/2
		for (unsigned int ii=0; ii<1+Mmain; ii++)
		{
			const unsigned int i=kreorder[ii];
			const int indice_Dij_minus_j=i*(this->nspecies_) - (i+1)*(i+2)/2;
			for (unsigned int jj=0;jj<1+Mmain;jj++)
			{
				const unsigned int j = kreorder[jj];
				if(i == j)
				{
					Dij(ii,jj) = 1. / *(dself+i);
				}
				else
				{
					// if(i < j) { Dij = 1. / *( d + j + i*(this->nspecies_) - (i+1)*(i+2)/2 ); }
					if(i < j) { Dij(ii,jj) = 1. / *( d + j + indice_Dij_minus_j ); }
					else      { Dij(ii,jj) = 1. / *( d + i + j*(this->nspecies_) - (j+1)*(j+2)/2 ); }
				}
			}
		}

		// temperature dependent (crot/kB)/zrot, cint/kB and Dij/(Dint,ij)
		Eigen::VectorXd crot_over_kB_div_zrot(1+Mmain);
		crot_over_kB_div_zrot.setZero();
		Eigen::VectorXd cint_over_kB(1+Mmain);
		cint_over_kB.setZero();
		Eigen::MatrixXd Dij_div_Dijint(1+Mmain,1+Mmain);
		Dij_div_Dijint.setZero();

		double deltaprime = 2985. / std::sqrt(T_*T_*T_);
		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			unsigned int i=kreorder[ii];
			int nlin = raw_transport_parameters_[i][0];
			if (nlin == 0)
			{
				crot_over_kB_div_zrot(ii) = 0.;
				cint_over_kB(ii) = 0.;
			}
			else
			{
				// Parker-Brau-Jonkman correction
				double eps_over_kB = raw_transport_parameters_[i][1];
				double aux = std::sqrt(eps_over_kB / T_);
				double fParker_T = 1. + aux * ( PhysicalConstants::ZROTA + aux * (PhysicalConstants::ZROTB + aux * PhysicalConstants::ZROTC));
				double zrot = fParker_298_[i] / fParker_T;
				double zrot298 = raw_transport_parameters_[i][5];
				// problem when zrot298=0 (e.g. for OH)
				if (zrot298 > 1.)
					{zrot = zrot298 * zrot;}
				if (nlin == 1)
					{crot_over_kB_div_zrot(ii) = 1.0/zrot;}
				else if (nlin == 2)
					{crot_over_kB_div_zrot(ii) = 1.5/zrot;}

				double cp_div_R = cpspecies[i] / 8314.34; // R=8314.34 J/(kmol*K)
				cint_over_kB(ii) = cp_div_R - 2.5;
			}

			for (unsigned int jj=0;jj<1+Mmain;jj++)
			{
				Dij_div_Dijint(ii,jj) = 1.;
			}
			// polar species
			double dipmu = raw_transport_parameters_[i][3];
			if (dipmu > 1.e-20)
				{Dij_div_Dijint(ii,ii) = 1. + deltaprime;}
		}

		Eigen::MatrixXd f1(1+Mmain,1+Mmain);
		f1.setZero();
		Eigen::MatrixXd f2(1+Mmain,1+Mmain);
		f2.setZero();
		Eigen::VectorXd f3(1+Mmain);
		f3.setZero();
		Eigen::MatrixXd f0010(1+Mmain,1+Mmain);
		f0010.setZero();
		Eigen::MatrixXd f1010_2(1+Mmain,1+Mmain);
		f1010_2.setZero();
		Eigen::MatrixXd f1010_3(1+Mmain,1+Mmain);
		f1010_3.setZero();
		Eigen::MatrixXd f1001(1+Mmain,1+Mmain);
		f1001.setZero();
		Eigen::MatrixXd f0101(1+Mmain,1+Mmain);
		f0101.setZero();

		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			const unsigned int i=kreorder[ii];

			// symmetric factors
			for (unsigned int jj=0;jj<ii+1;jj++)
			{
				const unsigned int j=kreorder[jj];

				// f1 and f2
				f1(ii,jj) = 1./(Dij(ii,jj) * (M[i]+M[j]));
				f2(ii,jj) = f1(ii,jj) * M[i]*M[j]/(M[i]+M[j]);

				// temperature dependent Astar, Bstar, Cstar
				double TSLOG = std::log(T_/epsij_over_kb_[i][j]);
				double T1 = TSLOG;
				double T2 = TSLOG*T1;
				double T3 = TSLOG*T2;
				double T4 = TSLOG*T3;
				double T5 = TSLOG*T4;
				double T6 = TSLOG*T5;
				double ASTAR = CollisionIntegralMatrices::FITASTAR_1
						+ CollisionIntegralMatrices::FITASTAR_2*T1
						+ CollisionIntegralMatrices::FITASTAR_3*T2
						+ CollisionIntegralMatrices::FITASTAR_4*T3
						+ CollisionIntegralMatrices::FITASTAR_5*T4
						+ CollisionIntegralMatrices::FITASTAR_6*T5
						+ CollisionIntegralMatrices::FITASTAR_7*T6;
				double BSTAR = CollisionIntegralMatrices::FITBSTAR_1
						+ CollisionIntegralMatrices::FITBSTAR_2*T1
						+ CollisionIntegralMatrices::FITBSTAR_3*T2
						+ CollisionIntegralMatrices::FITBSTAR_4*T3
						+ CollisionIntegralMatrices::FITBSTAR_5*T4
						+ CollisionIntegralMatrices::FITBSTAR_6*T5
						+ CollisionIntegralMatrices::FITBSTAR_7*T6;
				double CSTAR = CollisionIntegralMatrices::FITCSTAR_1
						+ CollisionIntegralMatrices::FITCSTAR_2*T1
						+ CollisionIntegralMatrices::FITCSTAR_3*T2
						+ CollisionIntegralMatrices::FITCSTAR_4*T3
						+ CollisionIntegralMatrices::FITCSTAR_5*T4
						+ CollisionIntegralMatrices::FITCSTAR_6*T5
						+ CollisionIntegralMatrices::FITCSTAR_7*T6;

				f0010(ii,jj) = 2.5 - 3.*CSTAR;
				f1010_2(ii,jj) = 6.25 - 3.*BSTAR;
				f1010_3(ii,jj) = 4. * ASTAR * (1. + 5./(3.*boost::math::constants::pi<double>()) * (crot_over_kB_div_zrot(ii)+crot_over_kB_div_zrot(jj)));
				f1001(ii,jj) = 10./boost::math::constants::pi<double>() * ASTAR;
			}
		}

		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			for (unsigned int jj=ii+1;jj<1+Mmain;jj++)
			{
				f1(ii,jj) = f1(jj,ii);
				f2(ii,jj) = f2(jj,ii);
				f0010(ii,jj)   = f0010(jj,ii);
				f1010_2(ii,jj) = f1010_2(jj,ii);
				f1010_3(ii,jj) = f1010_3(jj,ii);
				f1001(ii,jj)   = f1001(jj,ii);
			}
		}

		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			const unsigned int i=kreorder[ii];
			// f3 and f0101
			int nlin = raw_transport_parameters_[i][0];
			if (nlin == 0)
			{
				f3(ii) = 0.;
				for (unsigned int jj=0;jj<1+Mmain;jj++) {f0101(ii,jj) = 0.;}
			}
			else
			{
				f3(ii) = M[i] * crot_over_kB_div_zrot(ii) / cint_over_kB(ii);
				for (unsigned int jj=0;jj<1+Mmain;jj++)
				{
					const unsigned int j=kreorder[jj];
					f0101(ii,jj) = -6.25/(cint_over_kB(ii) * Dij(ii,jj)) * (Dij_div_Dijint(ii,jj) + (0.24/M[j]) * f1001(ii,jj)*f3(ii));
				}
			}
		}
		//Mijkl_must_be_recalculated_ = false;
		//}


		Eigen::VectorXd a10(1+Mmain);
		a10.setZero();
		Eigen::MatrixXd M0010(1+Mmain,1+Mmain);
		M0010.setZero();
		Eigen::MatrixXd M1000(1+Mmain,1+Mmain);
		M1000.setZero();
		Eigen::MatrixXd M1010(1+Mmain,1+Mmain);
		M1010.setZero();
		Eigen::MatrixXd M1001(1+Mmain,1+Mmain);
		M1001.setZero();
		Eigen::MatrixXd M0110(1+Mmain,1+Mmain);
		M0110.setZero();
		Eigen::VectorXd inv_M0101(1+Mmain);
		inv_M0101.setZero();

		// M0010 and M1000
		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			const unsigned int i=kreorder[ii];
			for (unsigned int jj=0;jj<1+Mmain;jj++)
			{
				const unsigned int j=kreorder[jj];

				const double m0010_ij = f1(ii,jj)*f0010(ii,jj);
				M0010(ii,jj) += moleFractions[i]*M[i]*m0010_ij;
				M0010(ii,ii) -= moleFractions[j]*M[j]*m0010_ij; // delta_ij contribution
				M1000(jj,ii) += moleFractions[j]*M[i]*m0010_ij;
				M1000(ii,ii) -= moleFractions[j]*M[j]*m0010_ij; // delta_ij contribution
			}
		}

		const unsigned int Kmain=kreorder[0];
		for (unsigned int ii=1;ii<1+Mmain;ii++)
		{
			const unsigned int i=kreorder[ii];
			for (unsigned int jj=0;jj<1+Mmain;jj++)
			{
				M0010(ii,jj) -= (moleFractions[i]*M[i])/(moleFractions[Kmain]*M[Kmain]) * M0010(0,jj);
			}
		}

		// M1010
		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			const unsigned int i=kreorder[ii];
			for (unsigned int jj=0;jj<1+Mmain;jj++)
			{
				const unsigned int j=kreorder[jj];
				const double m1010_ij_2 = f1010_2(ii,jj);
				const double m1010_ij_3 = f1010_3(ii,jj);

				M1010(ii,jj) += moleFractions[i]*f2(ii,jj)*(7.5+m1010_ij_2-m1010_ij_3);
				M1010(ii,ii) -= moleFractions[j]*f2(ii,jj)*((M[i]/M[j])*7.5+(M[j]/M[i])*m1010_ij_2+m1010_ij_3); // delta_ij contribution
			}
		}

		// M1001 and M0110
		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			const unsigned int i=kreorder[ii];
			for (unsigned int jj=0;jj<1+Mmain;jj++)
			{
				const unsigned int j=kreorder[jj];
				const double m1001_ij = f1(ii,jj) * f1001(ii,jj);

				M1001(ii,jj) += moleFractions[i]*m1001_ij;
				M1001(ii,ii) += moleFractions[j]*m1001_ij; // delta_ij contribution
			}
		}
		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			for (unsigned int jj=0;jj<1+Mmain;jj++)
			{
				M0110(ii,jj) = f3(ii) * M1001(ii,jj);
				M1001(ii,jj) *= f3(jj);
			}
		}

		// [M0101]⁻¹
		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			const unsigned int i=kreorder[ii];
			int nlin = raw_transport_parameters_[i][0];
			if (nlin == 0)
				{inv_M0101(ii) = 0.;}
			else
			{
				double sum1 = 0.;
				for (unsigned int kk=0;kk<1+Mmain;kk++)
				{
					const unsigned int k=kreorder[kk];
					sum1 += moleFractions[k] * f0101(ii,kk);
				}
				inv_M0101(ii) = 1. / sum1;
			}
		}

		//==================================
		//=== Solve thermodiffusion problem
		//==================================
		// Neumann series approximation at order 0
		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			double sum1 = 0.;
			for (unsigned int kk=0;kk<1+Mmain;kk++)
			{
				sum1 += M0110(kk,ii) * inv_M0101(kk);
			}
			a10(ii) = (1. - sum1) / M1010(ii,ii);
		}

		// Neumann series approximation at order 'norder=3'
		Eigen::VectorXd approx(1+Mmain);
		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			approx(ii) = a10(ii);
		}
		Eigen::VectorXd sum1(1+Mmain);
		sum1.setZero();
		Eigen::VectorXd u1(1+Mmain);
		u1.setZero();
		Eigen::VectorXd v1(1+Mmain);
		v1.setZero();
		Eigen::VectorXd v2(1+Mmain);
		v2.setZero();
		const unsigned int norder=3;
		for (unsigned int n=1;n<=norder;n++)
		{
			for (unsigned int ii=0;ii<1+Mmain;ii++)
			{
				sum1(ii) = 0.;
				for (unsigned int jj=0;jj<1+Mmain;jj++)
				{
					sum1(ii) -= M1010(jj,ii) * approx(jj);
				}
			}

			for (unsigned int ll=1;ll<1+Mmain;ll++)
			{
				u1(ll) = 0.;
				for (unsigned int jj=0;jj<1+Mmain;jj++)
				{
					u1(ll) += M1000(jj,ll) * approx(jj);
				}
			}
			double* DkKckl = this->DiKcij_.data();
			for (unsigned int kk=1;kk<1+Mmain;kk++)
			{
				v1(kk) = 0.;
				const unsigned int k=kreorder[kk];
				if (moleFractions[k] > 1.e-20)
				{
					for (unsigned int ll=1;ll<1+Mmain;ll++)
					{
						const unsigned int l=kreorder[ll];
						v1(kk) += moleFractions[l] * (*DkKckl) * u1(ll);
						DkKckl++;
					}
					v1(kk) += v1(kk)/moleFractions[k];
				}
				else
				{
					for (unsigned int ll=1;ll<1+Mmain;ll++)
					{
						DkKckl++;
					}
				}
			}

			for (unsigned int kk=0;kk<1+Mmain;kk++)
			{
				v2(kk) = 0.;
				for (unsigned int jj=0;jj<1+Mmain;jj++)
				{
					v2(kk) += M1001(jj,kk) * approx(jj);
				}
			}

			for (unsigned int ii=0;ii<1+Mmain;ii++)
			{
				double sum2 = 0.;
				for (unsigned int kk=1;kk<1+Mmain;kk++)
				{
					sum2 += M0010(kk,ii) * v1(kk);
				}

				double sum3 = 0.;
				for (unsigned int kk=0;kk<1+Mmain;kk++)
				{
					sum3 += M0110(kk,ii) * v2(kk)*inv_M0101(kk);
				}

				approx(ii) += (sum1(ii)+sum2+sum3)/M1010(ii,ii);
			}

			for (unsigned int ii=0;ii<1+Mmain;ii++)
			{
				a10(ii) += approx(ii);
			}
		}

		// Store the thermodiffusion arrays to be used in "lambdaMult" and "tetaMult"
		this->a10_.resize(1+Mmain);
		this->a10_.setZero();
		this->M1000_.resize(Mmain,1+Mmain);
		this->M1000_.setZero();
		this->M1001_.resize(1+Mmain,1+Mmain);
		this->M1001_.setZero();
		this->inv_M0101_.resize(1+Mmain);
		this->inv_M0101_.setZero();
		double* a10save = this->a10_.data();
		double* M1000save = this->M1000_.data();
		double* M1001save = this->M1001_.data();
		double* inv_M0101save = this->inv_M0101_.data();
		for (unsigned int ii=0;ii<1+Mmain;ii++)
		{
			*a10save = a10(ii);
			a10save++;
			*inv_M0101save = inv_M0101(ii);
			inv_M0101save++;
			for (unsigned int jj=0;jj<1+Mmain;jj++)
			{
				*M1001save = M1001(jj,ii);
				M1001save++;
			}
		}
		for (unsigned int ii=1;ii<1+Mmain;ii++)
		{
			for (unsigned int jj=0;jj<1+Mmain;jj++)
			{
				*M1000save = M1000(jj,ii);
				M1000save++;
			}
		}
	}

	void TransportPropertiesMap_CHEMKIN::ThermalDiffusionCoefficients(double* GammaSoret, const double T, const double* X, const double* Y)
	{
		double sum1 = 0.;
		for (unsigned int j = 0; j < this->nspecies_; j++)
			sum1 += M_0511_[j] * X[j];

		double sum2 = 0.;
		for (unsigned int j = 0; j < this->nspecies_; j++)
			sum2 += M_0489_[j] * X[j];

		const double C1 = -2.59e-7 * std::pow(T, 0.659);
		const double C2 = sum1 / sum2;

		for (unsigned int j = 0; j < this->nspecies_; j++)
			GammaSoret[j] = C1 * (M_0511_[j] * X[j] / sum1 - Y[j]) * C2;

	}

	void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
	{
		unsigned int numRows = matrix.rows() - 1;
		unsigned int numCols = matrix.cols();

		if (rowToRemove < numRows)
			matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

		matrix.conservativeResize(numRows, numCols);
	}

	void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
	{
		unsigned int numRows = matrix.rows();
		unsigned int numCols = matrix.cols() - 1;

		if (colToRemove < numCols)
			matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.rightCols(numCols - colToRemove);

		matrix.conservativeResize(numRows, numCols);
	}

	unsigned int maxElement(const double* v, const unsigned int n)
	{
		double maxv = v[0];
		unsigned int maxi = 0;

		for (unsigned int i = 1; i < n; i++)
		{
			if (v[i] > maxv)
			{
				maxv = v[i];
				maxi = i;
			}
		}
	
		return maxi;
	}

	void TransportPropertiesMap_CHEMKIN::SetupStefanMaxwellFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, const std::vector<std::string>& names)
	{
		Grammar_StefanMaxwell grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@SwapPolicy") == true)
		{
			std::string policy;
			dictionary.ReadString("@SwapPolicy", policy);
			StefanMaxwellSkipPolicy(policy, names);
		}
	}

	void TransportPropertiesMap_CHEMKIN::StefanMaxwellSkipPolicy(const std::string policy, const std::vector<std::string>& names)
	{
		if (policy == "most-abundant")
		{
			stefanmaxwell_skip_policy_ = STEFANMAXWELL_SKIP_POLICY_MOST_ABUNDANT;
		}
		else if (policy == "last")
		{
			stefanmaxwell_skip_policy_ = STEFANMAXWELL_SKIP_POLICY_POLICY_FIXED;
			stefanmaxwell_skip_ = this->nspecies_ - 1;
		}
		else
		{
			std::vector<std::string>::const_iterator it = std::find(names.begin(), names.end(), policy);
			if (it != names.end())
			{
				stefanmaxwell_skip_policy_ = STEFANMAXWELL_SKIP_POLICY_POLICY_FIXED;
				stefanmaxwell_skip_ = std::distance(names.begin(), it);
			}
			else
			{
				ErrorMessage(	"TransportPropertiesMap_CHEMKIN::StefanMaxwellSkipPolicy",
								"Species " + policy + " is not available in the kinetic mechanism!");
			}
		}
	}

	void TransportPropertiesMap_CHEMKIN::StefanMaxwellGammaBar(Eigen::MatrixXd& GammaBar, const double* X)
	{
		// Assembling binary matrix explicitly (the terms on the main diagonal are not necessary)
		Eigen::MatrixXd uD(this->nspecies_, this->nspecies_);
		uD.setZero();

		// The gammaSpecies_ field is actually the inverse of binary diffusion coefficients
		const double* d = this->gammaSpecies_.data();
		for (unsigned int i = 0; i < this->nspecies_; i++)
			for (unsigned int j = i + 1; j < this->nspecies_; j++)
			{
				uD(i, j) = *d;
				uD(j, i) = *d;
				d++;
			}

		// Adjust mole fractions in order to have a minimum amount of each species in the kinetic mechanism
		for (unsigned int i = 0; i < this->nspecies_; i++)
			x_corrected[i] = (X[i] + threshold_) / sum_threshold_;

		// Corrected molecular weight
		double MWmix = 0.;
		for (unsigned int i = 0; i < this->nspecies_; i++)
			MWmix += x_corrected[i] * M[i];

		// Species to skip
		if (stefanmaxwell_skip_policy_ == STEFANMAXWELL_SKIP_POLICY_MOST_ABUNDANT)
			stefanmaxwell_skip_ = maxElement(x_corrected, this->nspecies_);

		// Assembling matrix A
		Eigen::MatrixXd A(this->nspecies_, this->nspecies_);
		A.setZero();
		{
			// 1. Diagonal
			for (unsigned int i = 0; i < this->nspecies_; i++)
			{
				if (i != stefanmaxwell_skip_)
				{
					double sum = 0.;
					for (unsigned int j = 0; j < this->nspecies_; j++)
						if (j != i)
							sum += x_corrected[j] * uD(i, j) * MWmix / M[i];

					A(i, i) = -(x_corrected[i] * uD(i, stefanmaxwell_skip_) * MWmix / M[stefanmaxwell_skip_] + sum);
				}
			}

			// 2. Extra-diagonal
			for (unsigned int i = 0; i < this->nspecies_; i++)
				for (unsigned int j = 0; j < this->nspecies_; j++)
				{
					if (i != j && i != stefanmaxwell_skip_ && j != stefanmaxwell_skip_)
						A(i, j) = x_corrected[i] * (uD(i, j) * MWmix / M[j] - uD(i, stefanmaxwell_skip_) * MWmix / M[stefanmaxwell_skip_]);

				}
		}

		// Assembling matrix B
		Eigen::MatrixXd B(this->nspecies_, this->nspecies_);
		B.setZero();
		{
			// 1. Diagonal
			for (unsigned int i = 0; i < this->nspecies_; i++)
			{
				if (i != stefanmaxwell_skip_)
					B(i, i) = -(x_corrected[i] * MWmix / M[stefanmaxwell_skip_] + (1. - x_corrected[i]) * MWmix / M[i]);
			}

			// 2. Extra-diagonal
			for (unsigned int i = 0; i < this->nspecies_; i++)
				for (unsigned int j = 0; j < this->nspecies_; j++)
				{
					if (i != j && i != stefanmaxwell_skip_ && j != stefanmaxwell_skip_)
					{
						B(i, j) = x_corrected[i] * (MWmix / M[j] - MWmix / M[stefanmaxwell_skip_]);
					}
				}
		}

		bool debug_message = false;

		debug_message = false;
		if (debug_message == true)
		{
			std::cout << stefanmaxwell_skip_ << std::endl;

			std::cout << MWmix << std::endl;
			std::cout << std::endl;

			std::cout << "uD" << std::endl;
			std::cout << uD << std::endl;
			std::cout << std::endl;

			std::cout << "A" << std::endl;
			std::cout << A << std::endl;
			std::cout << std::endl;
			std::cout << "B" << std::endl;
			std::cout << B << std::endl;
			std::cout << std::endl;
		}

		// Remove row and column corresponding to the selected species for forcing the closure condition
		removeRow(A, stefanmaxwell_skip_);
		removeColumn(A, stefanmaxwell_skip_);
		removeRow(B, stefanmaxwell_skip_);
		removeColumn(B, stefanmaxwell_skip_);


		// Ordinary diffusion coefficients (without zero rows and columns)
		Eigen::MatrixXd G = A.partialPivLu().solve(B);

		debug_message = false;
		if (debug_message == true)
		{
			std::cout << "A" << std::endl;
			std::cout << A << std::endl;
			std::cout << std::endl;

			std::cout << "B" << std::endl;
			std::cout << B << std::endl;
			std::cout << std::endl;

			std::cout << "GammaBar(B)" << std::endl;
			std::cout << G << std::endl;
			std::cout << std::endl;
		}

		// Add row and column of zeros corresponding to the species to skip
		GammaBar.resize(this->nspecies_, this->nspecies_);
		GammaBar.setZero();
		for (unsigned int k = 0; k < this->stefanmaxwell_skip_; k++)
			for (unsigned int j = 0; j < this->stefanmaxwell_skip_; j++)
				GammaBar(k, j) = G(k, j);
		for (unsigned int k = 0; k < this->stefanmaxwell_skip_; k++)
			for (unsigned int j = this->stefanmaxwell_skip_+1; j < this->nspecies_; j++)
				GammaBar(k, j) = G(k, j-1);
		for (unsigned int k = this->stefanmaxwell_skip_ + 1; k < this->nspecies_; k++)
			for (unsigned int j = 0; j < this->stefanmaxwell_skip_; j++)
				GammaBar(k, j) = G(k-1, j);
		for (unsigned int k = this->stefanmaxwell_skip_ + 1; k < this->nspecies_; k++)
			for (unsigned int j = this->stefanmaxwell_skip_ + 1; j < this->nspecies_; j++)
				GammaBar(k, j) = G(k-1, j-1);

		debug_message = false;
		if (debug_message == true)
		{
			std::cout << "GammaBar(A)" << std::endl;
			std::cout << GammaBar << std::endl;
			std::cout << std::endl;

			std::cout << "D(k,k)" << std::endl;
			for (unsigned int k = 0; k < this->nspecies_; k++)
				std::cout << k << " " << GammaBar(k, k) << std::endl;
			std::cout << std::endl;
		}
	}

	void TransportPropertiesMap_CHEMKIN::StefanMaxwellFickFluxes(double* jsm,	const Eigen::MatrixXd& GammaBarL, const double rhoL, 
																				const Eigen::MatrixXd& GammaBarR, const double rhoR, const double* nablaY)
	{
		double sum_jsm = 0.;
		for (unsigned int i = 0; i < this->nspecies_; i++)
		{
			// The species to skip has a diffusion coefficient automatically set equal to 0
			// Thus there is no need to explicitly skip it
			jsm[i] = 0.;
			for (unsigned int j = 0; j < this->nspecies_; j++)
				jsm[i] += -0.50*(rhoL*GammaBarL(i,j)+ rhoR*GammaBarR(i,j)) * nablaY[j];

			// Sum the fluxes
			sum_jsm += jsm[i];
		}

		// Closure on the species to skip
		jsm[stefanmaxwell_skip_] = -sum_jsm;
	}

	void TransportPropertiesMap_CHEMKIN::StefanMaxwellSoretFluxes(double* jsm,
		const double rhoL, const double rhoR,
		const Eigen::VectorXd& GammaTL, const Eigen::VectorXd& GammaTR, 
		const double TL, const double TR, const double nablaT)
	{
		double sum_jsm = 0.;	
		for (unsigned int i = 0; i < this->nspecies_; i++)
		{
			jsm[i] = -0.50 * (GammaTL(i) / TL + GammaTR(i) / TR) * nablaT;
			sum_jsm += jsm[i];
		}

		// Closure on the species to skip
		jsm[stefanmaxwell_skip_] = -sum_jsm;
	}

	void TransportPropertiesMap_CHEMKIN::StefanMaxwellFluxes(	double* jsm, 
																	const Eigen::MatrixXd& GammaBarL, const double rhoL, 
																	const Eigen::MatrixXd& GammaBarR, const double rhoR, const double* nablaY, 
																	const Eigen::VectorXd& GammaTL, const Eigen::VectorXd& GammaTR, const double TL, const double TR, const double nablaT)
	{
		double sum_jsm = 0.;
		for (unsigned int i = 0; i < this->nspecies_; i++)
		{
			// The species to skip has a diffusion coefficient automatically set equal to 0
			// Thus there is no need to explicitly skip it
			jsm[i] = 0.;
			for (unsigned int j = 0; j < this->nspecies_; j++)
				jsm[i] += -0.50 * (rhoL * GammaBarL(i, j) + rhoR * GammaBarR(i, j)) * nablaY[j];

			// In case the Soret effect is accounted for
			jsm[i] += -0.50*(GammaTL(i)/TL+GammaTR(i)/TR) * nablaT;

			// Sum the fluxes
			sum_jsm += jsm[i];
		}

		// Closure on the species to skip
		jsm[stefanmaxwell_skip_] = -sum_jsm;
	}

	void TransportPropertiesMap_CHEMKIN::tetaMix(double* tetamix, const double* moleFractions)
	{
		for (unsigned int j=0;j<this->nspecies_;j++)
			tetamix[j] = 0.;

		for (unsigned int k=0;k<iThermalDiffusionRatios_.size();k++)
		{
			unsigned int iSpecies = iThermalDiffusionRatios_[k]-1;
			unsigned int i = k*this->nspecies_;
			for (unsigned int j=0;j<this->nspecies_;j++)
			{
				tetamix[iSpecies] += this->tetaSpecies_[i++]*moleFractions[iSpecies]*moleFractions[j];
			}
		}
	}

	// [Bertrand Naud] 1+M model
	inline void TransportPropertiesMap_CHEMKIN::reset_iThermalDiffusionRatios()
	{
		iThermalDiffusionRatios_.clear();
		iThermalDiffusionRatios_.resize(this->nspecies_);
		for (unsigned int i=0;i<this->nspecies_;i++)
		{
			iThermalDiffusionRatios_[i] = i+1;
		}
	}

	// [Bertrand Naud] 1+M model
	void TransportPropertiesMap_CHEMKIN::tetaMult(double* tetamix, const double* moleFractions, const unsigned int Mmain, const std::vector<unsigned int>& kreorder)
	{
		for (unsigned int j=0;j<this->nspecies_;j++)
			tetamix[j] = 0.;

		// main species j (other than Kmain)
		double* M1000_kj = this->M1000_.data();
		double sumcoef = 0.;
		for (unsigned int jj=1;jj<1+Mmain;jj++)
		{
			// main species j
			double sum1 = 0.;
			double* a10_k = this->a10_.data();
			for (unsigned int kk=0;kk<1+Mmain;kk++)
			{
				sum1 += (*M1000_kj)*(*a10_k);
				M1000_kj++;
				a10_k++;
			}
			const unsigned int j=kreorder[jj];
			tetamix[j] = 2.5 * moleFractions[j] * sum1;
			sumcoef += M[j] * tetamix[j];
		}

		unsigned int Kmain = kreorder[0];
		tetamix[Kmain] = - sumcoef / M[Kmain];
	}

	double TransportPropertiesMap_CHEMKIN::kPlanckMix(const double* moleFractions)
	{
		// References
		// * Nakamura and Shindo, Effects of radiation heat loss on laminar premixed ammonia/air flames
		//   Proceedings of the Combustion Institute, 37(2), p. 1741-1748 (2019), DOI: 10.1016/j.proci.2018.06.138
		// * Riviere and Soufiani, Updated band model parameters for H2O, CO2, CH4 and CO radiation at high temperature
		//   International Journal of Heat and Mass Transfer, 55(13-14), p. 3349-3358 (2012), DOI: 10.1016/j.ijheatmasstransfer.2012.03.019

		// Constants
		const double one_over_T = 1. / this->T_;
		const double T = this->T_;

		// 1. Water [1/m/bar]
		double K_H2O = 0.;
		if (T<2500.)
			K_H2O = -1.3876E+00 +one_over_T*( 5.0731E+03 +one_over_T*(-2.4595E+06 +one_over_T*( 5.8019E+09 +one_over_T*(-2.3041E+12 +one_over_T*( 3.2463E+14)))));	// Nakamura
		else
			K_H2O =  5.5460E-01 +one_over_T*(-9.7985E+03 +one_over_T*( 6.4090E+07 +one_over_T*(-1.9257E+11 +one_over_T*( 3.1871E+14 +one_over_T*(-2.0677E+17)))));	// Frassoldati
			 

		// 2. Carbon dioxide [1/m/bar]
		double K_CO2 = 0.;
		if (T<860.)		
			K_CO2 = -3.2951E+02 +T*(4.7368E+00 +T*(-2.5009E-02  +T*(6.6558E-05 +T*(-9.4178E-08 +T*(6.7643E-11 +T*(-1.9425E-14))))));	// Nakamura
		else if (T>2500.)		
			K_CO2 = 6.9379E+01 +T*(-8.5726E-02 +T*( 4.6526E-05  +T*(-1.3925E-08 +T*(2.3958E-12 +T*(-2.2295E-16 +T*( 8.7228E-21))))));	// Frassoldati
		else		
			K_CO2 = -8.0566E+00 +T*( 2.6921E-01 +T*(-5.3035E-04  +T*( 4.5006E-07 +T*(-1.9863E-10 +T*( 4.4901E-14 +T*(-4.1215E-18))))));	// Nakamura


		// 3. Carbon monoxide [1/m/bar]
		double K_CO = 0.;
		if (T<865.)		
			K_CO = 3.2799E+01 +T*(-3.7687E-01 +T*( 1.6586E-03  +T*(-3.5452E-06 +T*( 4.0556E-09 +T*(-2.4106E-12 +T*(5.8836E-16))))));	// Nakamura
		else if (T>2500.)		
			K_CO = 3.4400E+00 +T*(-3.0640E-03 +T*( 1.0630E-06  +T*(-1.6720E-10 +T*( 9.9570E-15 +T*( 0.0000E+00 +T*( 0.0000E+00))))));	// Riviere (extrapolated by Frassoldati)
		else		
			K_CO = 1.1310E+00 +T*( 2.2935E-02 +T*(-5.1068E-05  +T*( 4.5899E-08 +T*(-2.1060E-11 +T*( 4.9059E-15 +T*(-4.6155E-19))))));	// Nakamura


		// 4. Methane [1/m/bar]
		double K_CH4 = 0.;
		if (T<2300.)
			K_CH4 = -1.7121E+00 +one_over_T*( 4.0886E+03 +one_over_T*(2.6917E+06 +one_over_T*(-2.9348E+09 +one_over_T*(9.0450E+11 +one_over_T*(-9.7900E+13)))));	// Riviere (regression done by Frassoldati)
		else
			K_CH4 =  2.3325E+00 +one_over_T*(-3.9324E+04 +one_over_T*(2.5824E+08 +one_over_T*(-8.0965E+11 +one_over_T*(1.2509E+15 +one_over_T*(-7.5719E+17)))));	// Riviere (regression done by Frassoldati)
		

		// 5. Ammonia [1/m/bar]
		double K_NH3 = 0.;
		if (T<1100.)
			K_NH3 =  1.1418E+01 +one_over_T*(-3.5146E+04 +one_over_T*( 3.2908E+07 +one_over_T*(-6.9947E+09 +one_over_T*(2.2340E+11 +one_over_T*( 6.1863E+13)))));	// Nakamura
		else if (T>2500.)
			K_NH3 = -2.9706E-02 +one_over_T*( 3.3703E+02 +one_over_T*(-1.2287E+06 +one_over_T*( 1.5001E+09 +one_over_T*(0.0000E+00 +one_over_T*( 0.0000E+00)))));	// Riviere (extrapolated by Frassoldati)
		else
			K_NH3 = -5.0409E-02 +one_over_T*(-8.4977E+02 +one_over_T*( 9.1183E+06 +one_over_T*(-2.9019E+10 +one_over_T*(3.6044E+13 +one_over_T*(-1.2777E+16)))));	// Nakamura


		// 6. NO [1/m/bar]
		double K_NO = 0.;
		if (T<785.)		
			K_NO = 3.0573E+01 +T*(-3.7121E-01 +T*( 1.7629E-03  +T*(-4.1413E-06 +T*( 5.2412E-09  + T*(-3.4510E-12 +T*( 9.3212E-16))))));	// Nakamura
		else if (T>2500.)		
			K_NO = 1.8298E-01 +T*(-4.7886E-05 +T*( 5.3689E-09  +T*( 0.0000E+00 +T*( 0.0000E+00  + T*( 0.0000E+00 +T*( 0.0000E+00))))));	// Riviere (extrapolated by Frassoldati)
		else		
			K_NO = 3.2858E+00 +T*( 4.8995E-03 +T*(-1.9158E-05  +T*( 1.9830E-08 +T*( -9.7583E-12 + T*( 2.3738E-15 +T*(-2.3010E-19))))));	// Nakamura


		// 7. N2O [1/m/bar]
		double K_N2O = 0.;
		if (T<820.)		
			K_N2O = 9.3493E+01 +T*(-1.0433E+00 +T*( 4.8402E-03  +T*(-1.0351E-05 +T*( 1.1505E-08 +T*(-6.5867E-12 +T*( 1.5560E-15))))));	// Nakamura
		else if (T>2500.)		
			K_N2O = 7.7345E-01 +T*(-3.1345E-04 +T*( 3.5233E-08  +T*( 0.0000E+00 +T*( 0.0000E+00 +T*( 0.0000E+00 +T*( 0.0000E+00))))));	// Riviere (extrapolated by Frassoldati)
		else		
			K_N2O = 2.6905E+01 +T*( 1.2492E-01 +T*(-3.7050E-04  +T*( 3.7380E-07 +T*(-1.8389E-10 +T*( 4.4934E-14 +T*(-4.3785E-18))))));	// Nakamura


		// Old values
		/*
		{
			const double uT = 1000. / this->T_;
			const double T = this->T_;

			// 1. Water [1/m/bar]
			const double K_H2O = -0.23093 + uT*(-1.1239 + uT*(9.4153 + uT*(-2.9988 + uT*(0.51382 + uT*(-1.8684e-5)))));

			// 2. Carbon Dioxide [1/m/bar]
			const double K_CO2 = 18.741  +uT*(-121.31+uT*(273.5  +uT*(-194.05 +uT*( 56.31   + uT*(-5.8169)))));

			// 3. Carbon monoxide [1/m/bar]
			double K_CO;
			if (T < 750.)	K_CO = 4.7869 + T*(-0.06953 + T*(2.95775e-4 + T*(-4.25732e-7 + T*2.02894e-10)));
			else            K_CO = 10.09 + T*(-0.01183 + T*(4.7753e-6 + T*(-5.87209e-10 + T*-2.5334e-14)));

			// 4. Methane [1/m/bar]
			const double K_CH4 = 6.6334 + T*(-0.0035686 + T*(1.6682e-08 + T*(2.5611e-10 - 2.6558e-14*T)));
		}
		*/


		// Mole fractions
		const double x_H2O = (index_H2O_>0) ? moleFractions[index_H2O_-1] : 0.;
		const double x_CO2 = (index_CO2_>0) ? moleFractions[index_CO2_-1] : 0.;
		const double x_CO  = (index_CO_>0)  ? moleFractions[index_CO_-1]  : 0.;
		const double x_CH4 = (index_CH4_>0) ? moleFractions[index_CH4_-1] : 0.;
		const double x_NH3 = (index_NH3_>0) ? moleFractions[index_NH3_-1] : 0.;
		const double x_NO  = (index_NO_>0)  ? moleFractions[index_NO_-1]  : 0.;
		const double x_N2O = (index_N2O_>0) ? moleFractions[index_N2O_-1] : 0.;

		// Absorption coefficients
		const double as_H2O = K_H2O*x_H2O;	// [1/m/bar]
		const double as_CO2 = K_CO2*x_CO2;	// [1/m/bar]
		const double as_CO  = K_CO*x_CO;	// [1/m/bar]
		const double as_CH4 = K_CH4*x_CH4;	// [1/m/bar]
		const double as_NH3 = K_NH3*x_NH3;	// [1/m/bar]
		const double as_NO  = K_NO*x_NO;	// [1/m/bar]
		const double as_N2O = K_N2O*x_N2O;	// [1/m/bar]

		// Global absorption coefficient
		const double kPlanck = (as_H2O + as_CO2 + as_CO + as_CH4 + as_NH3 + as_NO + as_N2O) * (this->P_ / 1.e5);	// [1/m]

		// In case we want to add a degree of freedom with respect to pressure
		// const double kPlanck = (as_H2O + as_CO2 + as_CO + as_CH4 + as_NH3 + as_NO + as_N2O) * std::pow(this->P_ / 1.e5, m);	// [1/m]

		return kPlanck;
	}

	double TransportPropertiesMap_CHEMKIN::kCollision(const unsigned int i, const unsigned int j, const double T)
	{
		const double sigma = 0.50*(sigma_[i] + sigma_[j]);
		const double epsilon_over_kb = std::sqrt(epsilon_over_kb_[i]*epsilon_over_kb_[j]);
		const double mu = mu_[i]*mu_[j]/(mu_[i] + mu_[j]);
		const double TStar = T / epsilon_over_kb;

		return	PhysicalConstants::alphaCollision *
				(sigma*sigma)*Omega11(TStar)*std::sqrt(T/mu);
	}

	double TransportPropertiesMap_CHEMKIN::Omega11(const double TStar)
	{
		return	1.16145*std::pow(TStar, -0.14874) +
				0.52487*std::exp(-0.7732*TStar) +
				2.16178*std::exp(-2.437887*TStar);
	}

	double TransportPropertiesMap_CHEMKIN::Omega11(const double TStar, const double DStar)
	{
		return OpenSMOKE::CollisionIntegral11(TStar, DStar);
	}

	void TransportPropertiesMap_CHEMKIN::Test(const int nLoops, const double& T, int* index)
	{
		double lambdamix, etamix;
		Eigen::VectorXd gammamix(this->nspecies_);
		Eigen::VectorXd tetamix(this->nspecies_);
		
		// Composition (mole fractions)
		Eigen::VectorXd x(this->nspecies_);
		for(unsigned int i=0;i<this->nspecies_;i++)
			x(i) = 1./double(this->nspecies_);

		// Loops
		unsigned int speciesLoopsLambda	= nLoops*125;
		unsigned int speciesLoopsEta		= nLoops*125;
		unsigned int speciesLoopsGamma	= nLoops*1;
		unsigned int speciesLoopsTeta    = nLoops*25;
		unsigned int mixtureLoopsLambda	= nLoops*100;
		unsigned int mixtureLoopsEta		= nLoops*1;
		unsigned int mixtureLoopsGamma	= nLoops*1;
		unsigned int mixtureLoopsTeta	= nLoops*5;


		// Times
		double speciesLambdaTime, speciesEtaTime, speciesGammaTime, speciesTetaTime;
		double mixtureLambdaTime, mixtureEtaTime, mixtureGammaTime, mixtureTetaTime;
                
                SetTemperature(T);
                SetPressure(1e5);

		{
			std::cout << "Species Lambda..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=speciesLoopsLambda;k++)
			{
				lambda();
			}
			double tEnd = OpenSMOKEGetCpuTime();
			speciesLambdaTime = tEnd - tStart;
		}

		{
			std::cout << "Species Eta..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=speciesLoopsEta;k++)
			{
				eta();
			}
			double tEnd = OpenSMOKEGetCpuTime();
			speciesEtaTime = tEnd - tStart;
		}

		{
			std::cout << "Species Gamma..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=speciesLoopsGamma;k++)
			{
				gamma();
			}
			double tEnd = OpenSMOKEGetCpuTime();
			speciesGammaTime = tEnd - tStart;
		}

		{
			std::cout << "Species Teta..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=speciesLoopsTeta;k++)
			{
				teta();
			}
			double tEnd = OpenSMOKEGetCpuTime();
			speciesTetaTime = tEnd - tStart;
		}

		{
			std::cout << "Mixture Lambda..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=mixtureLoopsLambda;k++)
			{
				lambdamix = lambdaMix(x.data());
			}
			double tEnd = OpenSMOKEGetCpuTime();
			mixtureLambdaTime = tEnd - tStart;
		}

		{
			std::cout << "Mixture Eta..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=mixtureLoopsEta;k++)
			{
				etamix = etaMix(x.data());
			}
			double tEnd = OpenSMOKEGetCpuTime();
			mixtureEtaTime = tEnd - tStart;
		}

		{
			std::cout << "Mixture Gamma..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=mixtureLoopsGamma;k++)
			{
				gammaMix(gammamix.data(), x.data());
			}
			double tEnd = OpenSMOKEGetCpuTime();
			mixtureGammaTime = tEnd - tStart;
		}

		{
			std::cout << "Mixture Teta..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=mixtureLoopsTeta;k++)
			{
				tetaMix(tetamix.data(), x.data());
			}
			double tEnd = OpenSMOKEGetCpuTime();
			mixtureTetaTime = tEnd - tStart;
		}

		double lambdaTime	= mixtureLambdaTime/double(mixtureLoopsLambda)	+ speciesLambdaTime/double(speciesLoopsLambda);
		double etaTime		= mixtureEtaTime/double(mixtureLoopsEta)	+ speciesEtaTime/double(speciesLoopsEta);
		double gammaTime	= mixtureGammaTime/double(mixtureLoopsGamma)	+ speciesGammaTime/double(speciesLoopsGamma);
		double tetaTime	= mixtureTetaTime/double(mixtureLoopsTeta)	+ speciesTetaTime/double(speciesLoopsTeta);
		double totTime		= lambdaTime + etaTime + gammaTime + tetaTime;

		std::ofstream fBenchmark;
		fBenchmark.open("Benchmark.plus", std::ios::out);
		fBenchmark.setf(std::ios::scientific);

		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;
		fBenchmark << "                                       PROPERTY VALUES                                                   " << std::endl;
		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;

		fBenchmark << "Thermal conductivities     ";
		fBenchmark << this->lambdaSpecies_(index[1]-1)	<< " " << this->lambdaSpecies_(index[2]-1) << " " << this->lambdaSpecies_(index[3]-1) << std::endl; 

		fBenchmark << "Dynamic viscosities        ";
		fBenchmark << this->etaSpecies_(index[1]-1)		<< " " << this->etaSpecies_(index[2]-1) << " " << this->etaSpecies_(index[3]-1) << std::endl; 

		fBenchmark << "Mass diffusivities         ";
		fBenchmark << this->gammaSpecies_(index[1]-1)	<< " " << this->gammaSpecies_(index[2]-1) << " " << this->gammaSpecies_(index[3]-1) << std::endl; 

		fBenchmark << "Therm. diff ratios         ";
		fBenchmark << this->tetaSpecies_(index[1]-1)	<< " " << this->tetaSpecies_(index[2]-1) << " " << this->tetaSpecies_(index[3]-1) << std::endl; 

		fBenchmark << "Thermal conductivity (mix) ";
		fBenchmark << lambdamix << std::endl;

		fBenchmark << "Dynamic viscosity (mix)    ";
		fBenchmark << etamix << std::endl;

		fBenchmark << "Mass diffusivities (mix)   ";
		fBenchmark << gammamix(index[1]-1) << " " << gammamix(index[2]-1) << " " << gammamix(index[3]-1) << std::endl;

		fBenchmark << "Th. diff. ratios (mix)     ";
		fBenchmark << tetamix(index[1]-1) << " " << tetamix(index[2]-1) << " " << tetamix(index[3]-1) << std::endl;
	
		fBenchmark << std::endl;


		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;
		fBenchmark << "                                      CPU TIME DETAILS                                                   " << std::endl;
		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;

		fBenchmark << "Thermal conductivity (T)       ";
		fBenchmark << speciesLoopsLambda	<< " " << speciesLambdaTime << " " << speciesLambdaTime/double(speciesLoopsLambda)*1000.	<< std::endl;
		
		fBenchmark << "Dynamic Viscosity (T)          ";
		fBenchmark << speciesLoopsEta		<< " " << speciesEtaTime	<< " " << speciesEtaTime/double(speciesLoopsEta)*1000.			<< std::endl;
		
		fBenchmark << "Mass Diffusivities (T)         ";		
		fBenchmark << speciesLoopsGamma		<< " " << speciesGammaTime	<< " " << speciesGammaTime/double(speciesLoopsGamma)*1000.		<< std::endl;
		
		fBenchmark << "Thermal diffusion ratios (T)   ";					
		fBenchmark << speciesLoopsTeta		<< " " << speciesTetaTime	<< " " << speciesTetaTime/double(speciesLoopsTeta)*1000.		<< std::endl;
		
		fBenchmark << "Thermal conductivity (mix)     ";	
		fBenchmark << mixtureLoopsLambda	<< " " << mixtureLambdaTime << " " << mixtureLambdaTime/double(mixtureLoopsLambda)*1000.	<< std::endl;
		
		fBenchmark << "Dynamic Viscosity (mix)        ";
		fBenchmark << mixtureLoopsEta		<< " " << mixtureEtaTime	<< " " << mixtureEtaTime/double(mixtureLoopsEta)*1000.			<< std::endl;
		
		fBenchmark << "Mass Diffusivities (mix)       ";
		fBenchmark << mixtureLoopsGamma		<< " " << mixtureGammaTime	<< " " << mixtureGammaTime/double(mixtureLoopsGamma)*1000.		<< std::endl;
		
		fBenchmark << "Thermal diffusion ratios (mix) ";	
		fBenchmark << mixtureLoopsTeta		<< " " << mixtureTetaTime	<< " " << mixtureTetaTime/double(mixtureLoopsTeta)*1000.		<< std::endl;
		
		fBenchmark << std::endl;


		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;
		fBenchmark << "                                      CPU TIME SUMMARY                                                   " << std::endl;
		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;

		fBenchmark << "Thermal conductivity (tot)      ";		
		fBenchmark << lambdaTime*1000.		<< " ms - " << lambdaTime/totTime*100 << " %" << std::endl;

		fBenchmark << "Dynamic Viscosity (tot)         ";
		fBenchmark << etaTime*1000.			<< " ms - " << etaTime/totTime*100 << " %" << std::endl;

		fBenchmark << "Mass Diffusivities (tot)        ";
		fBenchmark << gammaTime*1000.		<< " ms - " << gammaTime/totTime*100 << " %" << std::endl;

		fBenchmark << "Thermal diffusion ratios (tot)  ";
		fBenchmark << tetaTime*1000.		<< " ms - " << tetaTime/totTime*100 << " %" << std::endl;

		fBenchmark << "All the properties  (tot)       ";
		fBenchmark << totTime*1000.			<< " ms - " << totTime/totTime*100 << " %" << std::endl;

		fBenchmark.close();
	}
}

