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

#include "P1_1D.h"
#include "math/native-nls-solvers/NonLinearSystemSolver"

namespace OpenSMOKE
{
	const double P1_1D::SIGMA_ = 5.670373e-8;
	const double P1_1D::PI_ = boost::math::constants::pi<double>();
	const double P1_1D::PI_TIMES_4_ = 4.*boost::math::constants::pi<double>();
	const double P1_1D::PI_TIMES_2_ = 2.*boost::math::constants::pi<double>();
	const double P1_1D::SIGMA_OVER_PI_ = P1_1D::SIGMA_ / P1_1D::PI_;
	const double P1_1D::SIGMA_TIMES_4_ = 4.*P1_1D::SIGMA_ ;

	P1_1D::P1_1D(Grid1DMultiPurpose* grid) :
		grid_(*grid)
	{
		radiative_equilibrium_ = false;
		epsilon_wall1_ = 1.;
		epsilon_wall2_ = 1.;
		is_scattering_ = false;
		is_anisotropic_scattering_ = false;

		np_ = grid_.N();
		MemoryAllocation();
	}

	void P1_1D::MemoryAllocation()
	{
		T_.resize(np_);
		kappa_.resize(np_);
		sigma_.resize(np_);
		A1_.resize(np_);
		beta_.resize(np_);
		gamma_.resize(np_);
		dummy_.resize(np_);
		G_.resize(np_);
		dG_.resize(np_);
		diffG_.resize(np_);

		kappa_.setConstant(0.);
		sigma_.setConstant(0.);
		A1_.setConstant(0.);
		beta_.setConstant(0.);
		gamma_.setConstant(0.);
		dummy_.setConstant(0.);
		G_.setConstant(0.);
		dG_.setConstant(0.);
		diffG_.setConstant(0.);
	}

	void P1_1D::SetRadiativeEquilibrium(const bool flag)
	{
		radiative_equilibrium_ = flag;
	}

	void P1_1D::SetTemperature(const double T)
	{ 
		T_.setConstant(T); 
	}

	void P1_1D::SetPlanckMeanAbsorptionCoefficient(const double kappa)
	{
		kappa_.setConstant(kappa);
	}

	void P1_1D::SetScatteringCoefficient(const double sigma)
	{
		if (sigma > 1.e-32)
		{
			sigma_.setConstant(sigma);
			is_scattering_ = true;
		}
		else
		{
			sigma_.setConstant(0.);
			is_scattering_ = false;
			is_anisotropic_scattering_ = false;
		}
	}

	void P1_1D::SetAnisotropicScatteringCoefficient(const double A1)
	{
		if (std::fabs(A1) > 1.e-32)
		{
			A1_.setConstant(A1);
			is_anisotropic_scattering_ = true;
		}
		else
		{
			A1_.setConstant(0.);
			is_anisotropic_scattering_ = false;
		}
	}

	void P1_1D::SetTemperature(const Eigen::VectorXd& T)
	{
		if (T.size() != np_)
			OpenSMOKE::FatalErrorMessage("Discrete ordinates: wrong size in temperature field");
		T_ = T;
	}

	void P1_1D::SetPlanckMeanAbsorptionCoefficient(const Eigen::VectorXd& kappa)
	{
		if (kappa.size() != np_)
			OpenSMOKE::FatalErrorMessage("Discrete ordinates: wrong size in Planck mean absorption coefficient field");
		kappa_ = kappa;
	}

	void P1_1D::SetScatteringCoefficient(const Eigen::VectorXd& sigma)
	{
		if (sigma.size() != np_)
			OpenSMOKE::FatalErrorMessage("Discrete ordinates: wrong size in scattering coefficient field");

		if (sigma.maxCoeff() > 1e-32)
		{
			sigma_ = sigma;
			is_scattering_ = true;
		}
		else
		{
			sigma_.setConstant(0.);
			is_scattering_ = false;
			is_anisotropic_scattering_ = false;
		}
	}

	void P1_1D::SetAnisotropicScatteringCoefficient(const Eigen::VectorXd& A1)
	{
		if (A1.size() != np_)
			OpenSMOKE::FatalErrorMessage("Discrete ordinates: wrong size in anisotropic scattering coefficient field");

		if (std::fabs(A1.maxCoeff()) > 1e-32)
		{
			A1_ = A1;
			is_anisotropic_scattering_ = true;
		}
		else
		{
			A1_.setConstant(0.);
			is_anisotropic_scattering_ = false;
		}
	}

	void P1_1D::UpdateProperties()
	{
		for (unsigned int i = 0; i < np_; i++)
		{
			beta_(i) = kappa_(i) + sigma_(i);
			gamma_(i) = 1. / (3.*beta_(i) - A1_(i)*sigma_(i));
		}
	}

	double P1_1D::q_from_wall1()
	{
		return -gamma_(0)*grid_.FirstDerivativeSecondOrderLeftSide(G_);
	}

	double P1_1D::q_from_wall2()
	{
		return  gamma_(np_ - 1)*grid_.FirstDerivativeSecondOrderRightSide(G_);
	}

	void P1_1D::SetWallsTemperatures(const double T_wall1, const double T_wall2)
	{
		T_wall1_ = T_wall1;
		T_wall2_ = T_wall2;
	}

	void P1_1D::SetWallsEmissivities(const double epsilon_wall1, const double epsilon_wall2)
	{
		epsilon_wall1_ = epsilon_wall1;
		epsilon_wall2_ = epsilon_wall2;
	}

	double P1_1D::q(const unsigned int i)
	{
		if (i == 0)
			return -gamma_(0)*grid_.FirstDerivativeSecondOrderLeftSide(G_);	
		else if (i==np_-1)
			return -gamma_(np_ - 1)*grid_.FirstDerivativeSecondOrderRightSide(G_);
		else
			return -gamma_(i)*grid_.FirstDerivativeSecondOrder(G_,i);
		
	}

	double P1_1D::G(const unsigned int i)
	{
		return G_(i);
	}

	double P1_1D::divq(const unsigned int i)
	{
		return kappa_(i)*(SIGMA_TIMES_4_*std::pow(T_(i), 4.) - G_(i));
	}

	void P1_1D::SpatialEquations(const double* y, const double t, double* f)
	{
		for (unsigned int i = 0; i < np_; i++)
			G_(i) = y[i];

		// Radiative equilibrium
		if (radiative_equilibrium_ == true)
		{
			for (unsigned int i = 0; i < np_; i++)
				T_(i) = std::pow(G_(i) / SIGMA_TIMES_4_, 0.25);
		}

		// First derivative
		grid_.FirstDerivativeAccurateBC(dummy_, G_, dG_, Grid1DMultiPurpose::DERIVATIVE_CENTERED);
		grid_.DiffusionTerm(G_, gamma_, diffG_);

		// Internal wall
		f[0] = -(2. - epsilon_wall1_) / epsilon_wall1_ * 2.*gamma_(0) * dG_(0) + G_(0) - SIGMA_TIMES_4_*pow(T_wall1_, 4.);

		// Internal points
		for (unsigned int i = 1; i < np_ - 1; i++)
			f[i] = kappa_(i) * (G_(i) - SIGMA_TIMES_4_*pow(T_(i), 4.)) - diffG_(i);

		// External wall
		f[np_ - 1] = (2. - epsilon_wall2_) / epsilon_wall2_ * 2.*gamma_(np_ - 1)*dG_(np_ - 1) + G_(np_ - 1) - SIGMA_TIMES_4_*pow(T_wall2_, 4.);
	}

	class OpenSMOKE_MyNlsSystem_P1 
	{
		public:
			void assignP1(OpenSMOKE::P1_1D *P1)
			{
				ptP1 = P1;
			}
		
		protected:

			unsigned int ne_;

			void MemoryAllocation()
			{
			}

			virtual void Equations(const Eigen::VectorXd &x, Eigen::VectorXd &f)
			{
				ptP1->SpatialEquations(x.data(), 0, f.data());
			}

			void Jacobian(const Eigen::VectorXd &y, Eigen::MatrixXd &J)
			{
			};

			void Print(const int call, const double t, const double phiW, const Eigen::VectorXd &x, const Eigen::VectorXd &f)
			{
			}

		private:
			OpenSMOKE::P1_1D *ptP1;
	};

	void P1_1D::SolveSpatialEquations()
	{
		// Update properties
		UpdateProperties();

		// Initial guess
		Eigen::VectorXd yInitial(np_);
		for (unsigned int j = 0; j < np_; j++)
			yInitial(j) = SIGMA_TIMES_4_*pow(T_(j), 4.);

		// Create the system
		typedef NlsSMOKE::KernelBand<OpenSMOKE_MyNlsSystem_P1> kernel;
		NlsSMOKE::NonLinearSolver<kernel> nls_solver;

		// Set initial conditions
		nls_solver.assignP1(this);
		nls_solver.SetFirstGuessSolution(yInitial);
		nls_solver.SetBandSizes(2, 2);

		// Solve the non linear system
		double timeStart = OpenSMOKEGetCpuTime();
		const unsigned int status = nls_solver();
		double timeEnd = OpenSMOKEGetCpuTime();

		if (status >= 0)
		{
			Eigen::VectorXd fResiduals(np_);

			// Update the current solution
			{
				Eigen::VectorXd yFinal(np_);
				
				nls_solver.Solution(yFinal, fResiduals);
				for (unsigned int j = 0; j < np_; j++)
					G_(j) = yFinal(j);
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
		/*		std::cout << " * number of iterations:           " << nls_solver. << std::endl;
				std::cout << " * number of functions (gradient): " << nls_solver.NumFunctionsForGradient() << std::endl;
				std::cout << " * number of Newtons:              " << nls_solver.numberOfJacobianFactorizations() << std::endl;
				std::cout << " * number of Jacobians:            " << nls_solver.NumNumericalJacobians() << std::endl;
				std::cout << " * number of functions (Jacobian): " << nls_solver.NumFunctionsForNumericalJacobian() << std::endl;
				std::cout << " * number of functions (Jacobian): " << nls_solver.NumFunctionsForNumericalJacobian() << std::endl;
				std::cout << " * residuals (norm2):              " << fResiduals.Norm2() << std::endl;
				std::cout << " * residuals (norm1):              " << fResiduals.Norm1() << std::endl;
		*/		std::cout << std::endl;
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

	void P1_1D::PrintFile(const boost::filesystem::path file_name)
	{
		std::ofstream fOut(file_name.c_str(), std::ios::out);
		fOut.setf(std::ios::scientific);

		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "r[m](1)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "T[K](2)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "k[1/m](3)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "sigma[1/m](4)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "beta[1/m](5)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "gamma[m](6)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "A1[-](7)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "G[W/m2](8)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "q[W/m2](9)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "divq[W/m3](10)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "dummy(11)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "phi[-](12)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "psi[-](13)";
		fOut << std::scientific << std::setw(16) << std::setprecision(6) << "psi[-](14)";
		fOut << std::endl;

		for (unsigned int i = 0; i < np_; i++)
		{
			double q_ = q(i); 
			double divq_ = divq(i);

			fOut << std::scientific << std::setw(16) << std::setprecision(6) << grid_.x[i+1];
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << T_(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << kappa_(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << sigma_(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << beta_(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << gamma_(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << A1_(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << G_(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << q_;
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << divq_;
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << 0.;
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << (std::pow(T_(i), 4.) - std::pow(T_wall2_, 4.)) / (std::pow(T_wall1_, 4.) - std::pow(T_wall2_, 4.) + 1e-32);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << q_ / SIGMA_ / (std::pow(T_wall1_, 4.) - std::pow(T_wall2_, 4.)+1e-32);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << q_ / SIGMA_ / (std::pow(T_(i), 4.) - std::pow(T_wall2_, 4.));

			fOut << std::endl;
		}
		fOut.close();
	}
}
