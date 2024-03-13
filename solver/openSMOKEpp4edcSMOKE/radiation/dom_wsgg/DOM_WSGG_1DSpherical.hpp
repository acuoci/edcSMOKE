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

#include "radiation/wsgg/WSGG.h"
#include "DOM_WSGG_1DSpherical.h"
#include "math/native-nls-solvers/NonLinearSystemSolver"

namespace OpenSMOKE
{
	double Legendre(const unsigned int n, const double x);

	const double DOM_WSGG_1DSpherical::SIGMA_ = 5.670373e-8;
	const double DOM_WSGG_1DSpherical::PI_ = boost::math::constants::pi<double>();
	const double DOM_WSGG_1DSpherical::PI_TIMES_4_ = 4.*boost::math::constants::pi<double>();
	const double DOM_WSGG_1DSpherical::PI_TIMES_2_ = 2.*boost::math::constants::pi<double>();
	const double DOM_WSGG_1DSpherical::SIGMA_OVER_PI_ = DOM_WSGG_1DSpherical::SIGMA_ / DOM_WSGG_1DSpherical::PI_;
	const double DOM_WSGG_1DSpherical::SIGMA_TIMES_4_ = 4.*DOM_WSGG_1DSpherical::SIGMA_;

	DOM_WSGG_1DSpherical::DOM_WSGG_1DSpherical(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model) :
		grid_(*grid),
		wsgg_model_(model)
	{
		np_ = grid_.N();

		if (wsgg_model_ == OpenSMOKE::EXTINCTION_COEFFICIENT_WSGG_SMITH)
			nbands_ = 4;
		else if (wsgg_model_ == OpenSMOKE::EXTINCTION_COEFFICIENT_WSGG_TRUELOVE)
			nbands_ = 4;
		else if (wsgg_model_ == OpenSMOKE::EXTINCTION_COEFFICIENT_WSGG_YIN)
			nbands_ = 5;
		else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_SMITH_SOOT_SINGLEBAND)
			nbands_ = 1;
		else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_SMITH_SOOT)
			nbands_ = 8;
                else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_CASSOL)
			nbands_ = 100;
	}

	double DOM_WSGG_1DSpherical::Iblack(const double T)
	{
		return SIGMA_OVER_PI_*T*T*T*T;
	}

	void DOM_WSGG_1DSpherical::Finalize()
	{
		for (unsigned int i = 0; i < N_; i++)
		{
			const unsigned int j = N_ + i;
			mu_[j] = -mu_[i];
			wp_[j] = wp_[i];
		}

		mu_times_wp_ = new double[n_];
		for (unsigned int i = 0; i < n_; i++)
			mu_times_wp_[i] = mu_[i] * wp_[i];

		alfa_ = new double[n_ + 1];
		alfa_[0] = 0.;
		for (unsigned int i = 1; i < n_; i++)
			alfa_[i] = alfa_[i - 1] - 2.*wp_[i - 1] * mu_[i - 1];
		alfa_[n_] = 0.;

		MemoryAllocation();

		legendre_ = false;
		radiative_equilibrium_ = false;
		is_scattering_ = false;
		is_anisotropic_scattering_ = false;
		epsilon_wall1_ = 1.;
		epsilon_wall2_ = 1.;

		MemoryAllocation();
	}

	void DOM_WSGG_1DSpherical::SetRadiativeEquilibrium(const bool flag)
	{
		radiative_equilibrium_ = flag;
	}

	void DOM_WSGG_1DSpherical::SetLegendre(const bool flag)
	{
		legendre_ = flag;
		if (legendre_ == true)
		{
			Eigen::MatrixXd A_(n_, n_);
			for (unsigned int k = 0; k < n_; k++)
				A_(0, k) = wp_[k];
			for (unsigned int j = 1; j < n_; j++)
			for (unsigned int k = 0; k < n_; k++)
				A_(j, k) = wp_[k] * Legendre(j, mu_[k]);

			ALU_ = new Eigen::FullPivLU<Eigen::MatrixXd>(A_);
		}
	}

	void DOM_WSGG_1DSpherical::MemoryAllocation()
	{
		I_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
		{
			I_[i].resize(n_);
			for (unsigned int j = 0; j < n_; j++)
				I_[i][j].resize(nbands_);
		}

		dI_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
		{
			dI_[i].resize(n_);
			for (unsigned int j = 0; j < n_; j++)
				dI_[i][j].resize(nbands_);
		}

		fI_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
		{
			fI_[i].resize(n_);
			for (unsigned int j = 0; j < n_; j++)
				fI_[i][j].resize(nbands_);
		}

		T_.resize(np_);
		PCO2_Pa_.resize(np_);
		PH2O_Pa_.resize(np_);
                fv_.resize(np_);

		kappa_.resize(nbands_);
		beta_.resize(nbands_);
		a_.resize(nbands_);

		S_anistropic_scattering_.resize(n_);
	//	sigma_.resize(np_);
	//	A1_.resize(np_);
	//	omega_.resize(np_);

		for (unsigned int j = 0; j < nbands_; j++)
		{
			kappa_[j].resize(np_);
			beta_[j].resize(np_);
			a_[j].resize(np_);
			
			kappa_[j].setConstant(0.);
			beta_[j].setConstant(0.);
			a_[j].setConstant(0.);
		}

		T_.setConstant(0.);
		PCO2_Pa_.setConstant(0.);
		PH2O_Pa_.setConstant(0.);

		S_anistropic_scattering_.setConstant(0.);
	//	sigma_.setConstant(0.);
	//	A1_.setConstant(0.);
	//	omega_.setConstant(0.);
	}

	void DOM_WSGG_1DSpherical::UpdateWSGGCoefficients()
	{
		double* kappa = new double[nbands_];
		double* a = new double[nbands_];

		for (unsigned int i = 0; i < np_; i++)
		{
			if (wsgg_model_ == OpenSMOKE::EXTINCTION_COEFFICIENT_WSGG_TRUELOVE)
				WSGG_Truelove(T_(i), PCO2_Pa_(i), PH2O_Pa_(i), kappa, a);
			else if (wsgg_model_ == OpenSMOKE::EXTINCTION_COEFFICIENT_WSGG_SMITH)
				WSGG_Smith(T_(i), PCO2_Pa_(i), PH2O_Pa_(i), kappa, a);
			else if (wsgg_model_ == OpenSMOKE::EXTINCTION_COEFFICIENT_WSGG_YIN)
				WSGG_Yin(T_(i), PCO2_Pa_(i), PH2O_Pa_(i), kappa, a);
			else if (wsgg_model_ == OpenSMOKE::EXTINCTION_COEFFICIENT_WSGG_SMITH_SOOT_SINGLEBAND)
				WSGG_Smith_Soot_Singleband(T_(i), PCO2_Pa_(i), PH2O_Pa_(i), fv_(i), kappa, a);
			else if (wsgg_model_ == OpenSMOKE::EXTINCTION_COEFFICIENT_WSGG_SMITH_SOOT)
				WSGG_Smith_Soot(T_(i), PCO2_Pa_(i), PH2O_Pa_(i), fv_(i), kappa, a);
                        else if (wsgg_model_ == OpenSMOKE::EXTINCTION_COEFFICIENT_WSGG_CASSOL)
				WSGG_Cassol(T_(i), PCO2_Pa_(i), PH2O_Pa_(i), fv_(i), kappa, a);

			for (unsigned int j = 0; j < nbands_; j++)
			{
				kappa_[j](i) = kappa[j];
				a_[j](i) = a[j];
			}
		}

		delete [] kappa;
		delete [] a;
	}

	void DOM_WSGG_1DSpherical::UpdateProperties()
	{
		UpdateWSGGCoefficients();

		for (unsigned int j = 0; j < nbands_; j++)
		{
			for (unsigned int i = 0; i < np_; i++)
			{
				beta_[j](i) = kappa_[j](i);

			//	beta_[j](i) = kappa_[j](i) + sigma_(i);
			//	if (is_scattering_ == true)
			//	{
			//		if (beta_[j](i)>1e-16)
			//			omega_(i) = sigma_(i) / beta_[j](i);
			//	}
			}
		}
	}

	double DOM_WSGG_1DSpherical::kappa(const unsigned int i)
	{
		double sum = 0.;
		for (unsigned int j = 0; j < nbands_; j++)
			sum += a_[j](i)*kappa_[j](i);
		return sum;
	}

	double DOM_WSGG_1DSpherical::a(const unsigned int i)
	{
		double sum = 0.;
		for (unsigned int j = 0; j < nbands_; j++)
			sum += a_[j](i);
		return sum;
	}

	void DOM_WSGG_1DSpherical::SetTemperature(const double T)
	{ 
		T_.setConstant(T); 
	}

	void DOM_WSGG_1DSpherical::SetTemperature(const Eigen::VectorXd& T)
	{
		if (T.size() != np_)
			OpenSMOKE::FatalErrorMessage("Discrete ordinates: wrong size in temperature field");
		T_ = T;
	}

	void DOM_WSGG_1DSpherical::SetPressureCO2(const double PCO2)
	{
		PCO2_Pa_.setConstant(PCO2);
	}

	void DOM_WSGG_1DSpherical::SetPressureCO2(const Eigen::VectorXd& PCO2)
	{
		if (PCO2.size() != np_)
			OpenSMOKE::FatalErrorMessage("Discrete ordinates: wrong size in CO2 pressure field");
		PCO2_Pa_ = PCO2;
	}

	void DOM_WSGG_1DSpherical::SetPressureH2O(const double PH2O)
	{
		PH2O_Pa_.setConstant(PH2O);
	}

	void DOM_WSGG_1DSpherical::SetPressureH2O(const Eigen::VectorXd& PH2O)
	{
		if (PH2O.size() != np_)
			OpenSMOKE::FatalErrorMessage("Discrete ordinates: wrong size in H2O pressure field");
		PH2O_Pa_ = PH2O;
	}

	void DOM_WSGG_1DSpherical::SetSootVolumeFraction(const double fv)
	{
		fv_.setConstant(fv);
	}

	void DOM_WSGG_1DSpherical::SetSootVolumeFraction(const Eigen::VectorXd& fv)
	{
		if (fv.size() != np_)
			OpenSMOKE::FatalErrorMessage("P1: wrong size in soot volume fraction field");

		fv_ = fv;
	}

	double DOM_WSGG_1DSpherical::J_wall1(const unsigned int k)
	{
		return PI_*Iblack_wall1_*a_[k](0) -(1. - epsilon_wall1_) / epsilon_wall1_ * q(0,k);
	}

	double DOM_WSGG_1DSpherical::J_wall2(const unsigned int k)
	{
		return PI_*Iblack_wall2_*a_[k](np_-1) +(1. - epsilon_wall2_) / epsilon_wall2_ * q(np_ - 1,k);
	}

	double DOM_WSGG_1DSpherical::q_from_wall1()
	{
		// Alternative, equivalent way
		/*
		double sum = 0.;
		for (unsigned int i = 0; i < N_; i++)
			sum += mu_times_wp_[i] * I_[0](i);
		return epsilon_wall1_ * (PI_*Iblack_wall1_ + sum);
		*/

		double sum = 0.;
		for (unsigned int j = 0; j < nbands_; j++)
			for (unsigned int i = 0; i < n_; i++)
				sum += mu_times_wp_[i] * I_[0][i](j);

		return sum;
	}

	double DOM_WSGG_1DSpherical::q_from_wall2()
	{
		// Alternative, equivalent way
		/*
		double sum = 0.;
		for (unsigned int i = 0; i < N_; i++)
		{
			const unsigned int j = N_ + i;
			sum += mu_times_wp_[j] * I_[np_ - 1](j);
		}
		return -epsilon_wall2_ * (PI_*Iblack_wall2_ - sum);
		*/

		double sum = 0.;
		for (unsigned int j = 0; j < nbands_; j++)
			for (unsigned int i = 0; i < n_; i++)
				sum += mu_times_wp_[i] * I_[np_-1][i](j);

		return -sum;
	}

	void DOM_WSGG_1DSpherical::SetWallsTemperatures(const double T_wall1, const double T_wall2)
	{
		T_wall1_ = T_wall1;
		Iblack_wall1_ = Iblack(T_wall1_);
		T_wall2_ = T_wall2;
		Iblack_wall2_ = Iblack(T_wall2_);
	}

	void DOM_WSGG_1DSpherical::SetWallsEmissivities(const double epsilon_wall1, const double epsilon_wall2)
	{
		epsilon_wall1_ = epsilon_wall1;
		epsilon_wall2_ = epsilon_wall2;
	}

	double DOM_WSGG_1DSpherical::q(const unsigned int j)
	{
		double sum = 0.;
		for (unsigned int k = 0; k < nbands_; k++)
			sum += q(j,k);

		return sum;
	}

	double DOM_WSGG_1DSpherical::G(const unsigned int j)
	{
		double sum = 0.;
		for (unsigned int k = 0; k < nbands_; k++)
			sum += G(j,k);

		return sum;
	}

	double DOM_WSGG_1DSpherical::q(const unsigned int j, const unsigned int k)
	{
		double q = 0.;
		for (unsigned int i = 0; i < n_; i++)
			q += mu_times_wp_[i] * I_[j][i](k);

		return q;
	}
	
	double DOM_WSGG_1DSpherical::G(const unsigned int j, const unsigned int k)
	{
		double G = 0.;
		for (unsigned int i = 0; i < n_; i++)
			G += wp_[i] * I_[j][i](k);

		return G;
	}

	double DOM_WSGG_1DSpherical::divq(const unsigned int i)
	{
		double sum = 0.;
		for (unsigned int j = 0; j < nbands_; j++)
			sum += kappa_[j](i)*(a_[j](i)*SIGMA_TIMES_4_*std::pow(T_(i), 4.) - G(i,j));
		return sum;
	}

	void DOM_WSGG_1DSpherical::FirstOrderDerivatives()
	{
		Eigen::VectorXd dummy(np_);
		Eigen::VectorXd phi(np_);
		Eigen::VectorXd dphi(np_);

		for (unsigned int j = 0; j < nbands_; j++)
		{
			for (unsigned int i = 0; i < N_; i++)
			{
				for (unsigned int k = 0; k < np_; k++)
					phi(k) = I_[k][i](j);

				grid_.FirstDerivative(dummy, phi, dphi, Grid1DMultiPurpose::DERIVATIVE_FORWARD);

				for (unsigned int k = 0; k < np_; k++)
					dI_[k][i](j) = dphi(k);
			}
		}

		for (unsigned int jj = 0; jj < nbands_; jj++)
		{
			for (unsigned int i = 0; i < N_; i++)
			{
				const unsigned int j = N_ + i;

				for (unsigned int k = 0; k < np_; k++)
					phi(k) = I_[k][j](jj);

				grid_.FirstDerivative(dummy, phi, dphi, Grid1DMultiPurpose::DERIVATIVE_BACKWARD);

				for (unsigned int k = 0; k < np_; k++)
					dI_[k][j](jj) = dphi(k);
			}
		}
	}

	void DOM_WSGG_1DSpherical::SpatialEquations(const double* y, const double t, double* f)
	{
		// Recover main variables
		{
			unsigned int count = 0;
			for (unsigned int i = 0; i < np_; i++)
				for (unsigned int j = 0; j < n_; j++)
					for (unsigned int k = 0; k < nbands_; k++)
						I_[i][j](k) = y[count++];
		}

		// Radiative equilibrium
		if (radiative_equilibrium_ == true)
		{
			for (unsigned int i = 0; i < np_; i++)
				T_(i) = std::pow(G(i) / SIGMA_TIMES_4_, 0.25);
		}

		// First order derivatives
		FirstOrderDerivatives();

		// Local contributions
		if (legendre_ == false)
		{
			for (unsigned int j = 0; j < np_; j++)
			{
				const double rj = grid_.x[j + 1];

				double S_ = Iblack(T_(j));
				/*
				if (is_scattering_ == true)
				{
					S_ += -omega_(j)*Iblack(T_(j)) + omega_(j) / PI_TIMES_4_*G(I_[j]);
					if (is_anisotropic_scattering_ == true)
					{
						const double q_ = q(I_[j]);
						for (unsigned int i = 0; i < n_; i++)
							S_anistropic_scattering_(i) = omega_(j) / PI_TIMES_4_*A1_(j)*q_*mu_[i];
					}
				}
				*/

				// First ordinate
				for (unsigned int k = 0; k < nbands_; k++)
					fI_[j][0](k) = mu_[0] * dI_[j][0](k) + mu_[0] * I_[j][0](k) / rj + (alfa_[1] * I_[j][1](k)) / (2.*rj*wp_[0]) - beta_[k](j)*(S_*a_[k](0) + S_anistropic_scattering_(0) - I_[j][0](k));
				
				// Central ordinates
				for (unsigned int i = 1; i < n_ - 1; i++)
					for (unsigned int k = 0; k < nbands_; k++)
						fI_[j][i](k) = mu_[i] * dI_[j][i](k) + mu_[i] * I_[j][i](k) / rj + (alfa_[i + 1] * I_[j][i + 1](k) - alfa_[i] * I_[j][i - 1](k)) / (2.*rj*wp_[i]) - beta_[k](j)*(S_*a_[k](i-1) + S_anistropic_scattering_(i) - I_[j][i](k));

				// Last ordinate
				for (unsigned int k = 0; k < nbands_; k++)
					fI_[j][n_ - 1](k) = mu_[n_ - 1] * dI_[j][n_ - 1](k) + mu_[n_ - 1] * I_[j][n_ - 1](k) / rj + (-alfa_[n_ - 1] * I_[j][n_ - 2](k)) / (2.*rj*wp_[n_ - 1]) - beta_[k](j)*(S_*a_[k](np_-1) + S_anistropic_scattering_(n_ - 1) - I_[j][n_ - 1](k));
			}
		}
		
		/*
		if (legendre_ == true)
		{
			for (unsigned int j = 0; j < np_; j++)
			{
				Eigen::VectorXd b(n_);
				b(0) = 0.;
				for (unsigned int i = 1; i < n_; i++)
				{
					unsigned int m = i + 1;

					double sum = 0.;
					for (unsigned int k = 0; k < n_; k++)
						sum += wp_[k] * I_[j](k) * (Legendre(m, mu_[k]) - Legendre(m - 2, mu_[k]));
					b(i) = sum * m*(m - 1) / (2.*m - 1.);
				}

				Eigen::VectorXd D(n_);
				D = ALU_->solve(b);

				double	S_ = Iblack(T_(j));
				if (is_scattering_ == true)
				{
					S_ += -omega_(j)*Iblack(T_(j)) + omega_(j) / PI_TIMES_4_*G(I_[j]);
					if (is_anisotropic_scattering_ == true)
					{
						const double q_ = q(I_[j]);
						for (unsigned int i = 0; i < n_; i++)
							S_anistropic_scattering_(i) = omega_(j) / PI_TIMES_4_*A1_(j)*q_*mu_[i];
					}
				}

				const double rj = grid_.x[j + 1];

				fI_[j](0) = mu_[0] * dI_[j](0) + 2.*mu_[0] * I_[j](0) / rj + D(0) / rj - beta_(j)*(S_ + S_anistropic_scattering_(0) - I_[j](0));
				for (unsigned int i = 1; i < n_ - 1; i++)
					fI_[j](i) = mu_[i] * dI_[j](i) + 2.*mu_[i] * I_[j](i) / rj + D(i) / rj - beta_(j)*(S_ + S_anistropic_scattering_(i) - I_[j](i));
				fI_[j](n_ - 1) = mu_[n_ - 1] * dI_[j](n_ - 1) + 2.*mu_[n_ - 1] * I_[j](n_ - 1) / rj + D(n_ - 1) / rj - beta_(j)*(S_ + S_anistropic_scattering_(n_ - 1) - I_[j](n_ - 1));
			}
		}
		*/

		// Corrections for boundary conditions
		for (unsigned int k = 0; k < nbands_; k++)
		{
			for (unsigned int i = 0; i < N_; i++)
			{
				const unsigned int j = N_ + i;

				fI_[0][j](k) = I_[0][j](k) - J_wall1(k) / PI_;
				fI_[np_ - 1][i](k) = I_[np_ - 1][i](k) - J_wall2(k) / PI_;
			}
		}
	
		// Fill residuals
		{
			unsigned count = 0;
			for (unsigned int i = 0; i < np_; i++)
				for (unsigned int j = 0; j < n_; j++)
					for (unsigned int k = 0; k < nbands_; k++)
						f[count++] = fI_[i][j](k);
		}
	}

	class OpenSMOKE_MyNlsSystem_DiscreteOrdinatesWSGG
	{
	public:
		void assignDO(OpenSMOKE::DOM_WSGG_1DSpherical *DO)
		{
			ptDO = DO;
		}

	protected:

		unsigned int ne_;

		void MemoryAllocation()
		{
		}

		virtual void Equations(const Eigen::VectorXd &x, Eigen::VectorXd &f)
		{
			ptDO->SpatialEquations(x.data(), 0, f.data());
		}

		void Jacobian(const Eigen::VectorXd &y, Eigen::MatrixXd &J)
		{
		};

		void Print(const int call, const double t, const double phiW, const Eigen::VectorXd &x, const Eigen::VectorXd &f)
		{
		}

	private:

		OpenSMOKE::DOM_WSGG_1DSpherical *ptDO;
	};

	void DOM_WSGG_1DSpherical::SolveSpatialEquations()
	{
		// Update properties
		UpdateProperties();

		// Total number of equations
		const unsigned int neq = n_*np_*nbands_;

		// Initial guess
		Eigen::VectorXd yInitial(neq);
		unsigned int count = 0;
		for (unsigned int i = 1; i <= np_; i++)
		{
			for (unsigned int j = 1; j <= N_; j++)
				for (unsigned int k = 0; k < nbands_; k++)
					yInitial(count++) = Iblack_wall1_;

			for (unsigned int j = 1; j <= N_; j++)
				for (unsigned int k = 0; k < nbands_; k++)
					yInitial(count++) = Iblack_wall2_;
		}

		// Create the system
		typedef NlsSMOKE::KernelBand<OpenSMOKE_MyNlsSystem_DiscreteOrdinatesWSGG> kernel;
		NlsSMOKE::NonLinearSolver<kernel> nls_solver;

		// Set initial conditions
		nls_solver.assignDO(this);
		nls_solver.SetFirstGuessSolution(yInitial);
		nls_solver.SetBandSizes(2*n_*nbands_-1, 2*n_*nbands_-1);

		// Solve the non linear system
		double timeStart = OpenSMOKEGetCpuTime();
		const int status = nls_solver();
		double timeEnd = OpenSMOKEGetCpuTime();

		if (status >= 0)
		{
			Eigen::VectorXd yFinal(neq);
			Eigen::VectorXd fResiduals(neq);

			// Update the current solution
			{	
				nls_solver.Solution(yFinal, fResiduals);
				unsigned int count = 0;
				for (unsigned int i = 0; i < np_; i++)
					for (unsigned int j = 0; j < n_; j++)
						for (unsigned int k = 0; k < nbands_; k++)
							I_[i][j](k) = yFinal(count++);
			}

			bool verbose = false;
			if (verbose == true)
			{
				std::string message("Nls System successfully solved: ");
				if (status == 0)		message += "Start conditions";
				else if (status == 1)	message += "The maximum number of functions calls has been performed";
				else if (status == 2)	message += "The Newton correction has reached the required precision";
				else if (status == 3)	message += "The Quasi Newton correction has reached the required precision";
				else if (status == 4)	message += "The Gradient has reached the required precision";
				else if (status == 5)	message += "The Objective function phiNew has reached the required precision";
				else if (status == 6)	message += "The Objective function phiW has reached the required precision";
				else if (status == 7)	message += "The Objective function phiNew has reached the required precision but solution is dubious";
				else if (status == 8)	message += "Reached the assigned max value for Newton calls";

				std::cout << message << std::endl;

				std::cout << std::endl;
				std::cout << " * CPU time (s):                   " << timeEnd - timeStart << std::endl;
			/*	std::cout << " * number of iterations:           " << nls_object.IterationCounter() << std::endl;
				std::cout << " * number of functions (gradient): " << nls_object.NumFunctionsForGradient() << std::endl;
				std::cout << " * number of Newtons:              " << nls_object.NumNewtons() << std::endl;
				std::cout << " * number of Jacobians:            " << nls_object.NumNumericalJacobians() << std::endl;
				std::cout << " * number of functions (Jacobian): " << nls_object.NumFunctionsForNumericalJacobian() << std::endl;
				std::cout << " * residuals (norm2):              " << fResiduals.Norm2() << std::endl;
				std::cout << " * residuals (norm1):              " << fResiduals.Norm1() << std::endl;
			*/	std::cout << std::endl;
			}
		}
		else
		{
			std::string message("Discrete ordinates: Dae Solver Error: ");
			if (status == -1)		message += "It has been impossible to reach the solution";
			else if (status == -2)	message += "The search has been stopped";
			else if (status == -3)	message += "The object is not initialized";
			else if (status == -4)	message += "It has been impossible to reach the solution in Restart";
			std::cout << message << std::endl;
			
			exit(-1);
		}
	}

	void DOM_WSGG_1DSpherical::PrintFile(const boost::filesystem::path file_name)
	{
		std::ofstream fOut(file_name.c_str(), std::ios::out);
		fOut.setf(std::ios::scientific);

		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "r[m](1)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "T[K](2)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "k[1/m](3)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "k0[1/m](4)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "k1[1/m](5)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "k2[1/m](6)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "k3[1/m](7)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "a(8)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "a0[](9)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "a1[](10)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "a2[](11)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "a3[](12)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "G[W/m2](13)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "G0[W/m2](14)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "G1[W/m2](15)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "G2[W/m2](16)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "G3[W/m2](17)";

		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "q[W/m2](18)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "divq[W/m3](19)";

		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "dummy(20)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "phi[-](21)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "psi[-](22)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "psi[-](23)";

		unsigned int counter = 24;
		for (unsigned int k = 0; k < nbands_; k++)
		for (unsigned int j = 0; j < n_; j++)
		{
			std::stringstream label_number_direction; label_number_direction << j + 1;
			std::stringstream label_number_spectrum; label_number_spectrum << k;
			std::stringstream label_counter; label_counter << counter++;
			std::string label = "I" + label_number_spectrum.str() + "." + label_number_direction.str() + "(" + label_counter.str() + ")";
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << label;
		}
		fOut << std::endl;

		for (unsigned int i = 0; i < np_; i++)
		{
			double G_ = G(i);
			double q_ = q(i);
			double divq_ = divq(i);

			fOut << std::scientific << std::setw(16) << std::setprecision(6) << grid_.x[i + 1];
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << T_(i);

			fOut << std::scientific << std::setw(16) << std::setprecision(6) << kappa(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << kappa_[0](i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << kappa_[1](i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << kappa_[2](i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << kappa_[3](i);

			fOut << std::scientific << std::setw(16) << std::setprecision(6) << a(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << a_[0](i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << a_[1](i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << a_[2](i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << a_[3](i);

			fOut << std::scientific << std::setw(16) << std::setprecision(6) << G(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << G(i, 0);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << G(i, 1);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << G(i, 2);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << G(i, 3);

			fOut << std::scientific << std::setw(16) << std::setprecision(6) << q_;
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << divq_;

			fOut << std::scientific << std::setw(16) << std::setprecision(6) << 0.;
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << (std::pow(T_(i), 4.) - std::pow(T_wall2_, 4.)) / (std::pow(T_wall1_, 4.) - std::pow(T_wall2_, 4.) + 1e-32);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << q_ / SIGMA_ / (std::pow(T_wall1_, 4.) - std::pow(T_wall2_, 4.) + 1e-32);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << q_ / SIGMA_ / (std::pow(T_(i), 4.) - std::pow(T_wall2_, 4.));

			for (unsigned int k = 0; k < nbands_; k++)
				for (unsigned int j = 0; j < n_; j++)
					fOut << std::scientific << std::setw(16) << std::setprecision(6) << I_[i][j](k);

			fOut << std::endl;
		}
		fOut.close();
	}
}
