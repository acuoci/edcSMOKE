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

#include "SP3_1D.h"
#include "math/native-nls-solvers/NonLinearSystemSolver"

namespace OpenSMOKE
{
	const double SP3_1D::SIGMA_ = 5.670373e-8;
	const double SP3_1D::PI_ = boost::math::constants::pi<double>();
	const double SP3_1D::PI_TIMES_4_ = 4.*boost::math::constants::pi<double>();
	const double SP3_1D::PI_TIMES_2_ = 2.*boost::math::constants::pi<double>();
	const double SP3_1D::SIGMA_OVER_PI_ = SP3_1D::SIGMA_ / SP3_1D::PI_;
	const double SP3_1D::SIGMA_TIMES_4_ = 4.*SP3_1D::SIGMA_ ;

	SP3_1D::SP3_1D(Grid1DMultiPurpose* grid) :
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

	void SP3_1D::MemoryAllocation()
	{
		dummy_.resize(np_);

		T_.resize(np_);
		kappa_.resize(np_);
		sigma_.resize(np_);
		A1_.resize(np_);
		beta_.resize(np_);

		gamma0_.resize(np_);
		gamma2_.resize(np_);

		alfa0_.resize(np_);
		alfa1_.resize(np_);
		alfa2_.resize(np_);
		alfa3_.resize(np_);

		J0_.resize(np_);
		dJ0_.resize(np_);
		diffJ0_.resize(np_);
		J2_.resize(np_);
		dJ2_.resize(np_);
		diffJ2_.resize(np_);

		dummy_.setConstant(0.);
		kappa_.setConstant(0.);
		sigma_.setConstant(0.);
		A1_.setConstant(0.);
		beta_.setConstant(0.);
		J0_.setConstant(0.);
		dJ0_.setConstant(0.);
		diffJ0_.setConstant(0.);
		J2_.setConstant(0.);
		dJ2_.setConstant(0.);
		diffJ2_.setConstant(0.);

		alfa0_.setConstant(1.);
		alfa1_.setConstant(1.);
		alfa2_.setConstant(1.);
		alfa3_.setConstant(1.);
	}

	void SP3_1D::SetRadiativeEquilibrium(const bool flag)
	{
		radiative_equilibrium_ = flag;
	}

	void SP3_1D::SetTemperature(const double T)
	{ 
		T_.setConstant(T); 
	}

	void SP3_1D::SetPlanckMeanAbsorptionCoefficient(const double kappa)
	{
		kappa_.setConstant(kappa);
	}

	void SP3_1D::SetScatteringCoefficient(const double sigma)
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

	void SP3_1D::SetAnisotropicScatteringCoefficient(const double A1)
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

	void SP3_1D::SetTemperature(const Eigen::VectorXd& T)
	{
		if (T.size() != np_)
			OpenSMOKE::FatalErrorMessage("Discrete ordinates: wrong size in temperature field");
		T_ = T;
	}

	void SP3_1D::SetPlanckMeanAbsorptionCoefficient(const Eigen::VectorXd& kappa)
	{
		if (kappa.size() != np_)
			OpenSMOKE::FatalErrorMessage("Discrete ordinates: wrong size in Planck mean absorption coefficient field");
		kappa_ = kappa;
	}

	void SP3_1D::SetScatteringCoefficient(const Eigen::VectorXd& sigma)
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

	void SP3_1D::SetAnisotropicScatteringCoefficient(const Eigen::VectorXd& A1)
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

	void SP3_1D::UpdateProperties()
	{
		const double A0 = 0.;
		const double A2 = 0.;
		const double A3 = 0.;

		for (unsigned int i = 0; i < np_; i++)
		{
			beta_(i) = kappa_(i) + sigma_(i);

			const double omega = sigma_(i) / beta_(i);
			alfa0_(i) = 1. - omega*A0/1.;
			alfa1_(i) = 1. - omega*A1_(i)/3.;
			alfa2_(i) = 1. - omega*A2/5.;
			alfa3_(i) = 1. - omega*A3/7.;

			gamma0_(i) = 1. / (3.*alfa1_(i)*beta_(i));
			gamma2_(i) = 1. / (7.*alfa3_(i)*beta_(i));
		}
	}

	double SP3_1D::q_from_wall1()
	{
		return -PI_TIMES_4_*gamma0_(0)*grid_.FirstDerivativeLeftSide(J0_);
	}

	double SP3_1D::q_from_wall2()
	{
		return  PI_TIMES_4_*gamma0_(np_ - 1)*grid_.FirstDerivativeRightSide(J0_);
	}

	void SP3_1D::SetWallsTemperatures(const double T_wall1, const double T_wall2)
	{
		T_wall1_ = T_wall1;
		T_wall2_ = T_wall2;
	}

	void SP3_1D::SetWallsEmissivities(const double epsilon_wall1, const double epsilon_wall2)
	{
		epsilon_wall1_ = epsilon_wall1;
		epsilon_wall2_ = epsilon_wall2;
	}

	double SP3_1D::q(const unsigned int i)
	{
		if (i == 0)
			return -PI_TIMES_4_*gamma0_(0)*grid_.FirstDerivativeLeftSide(J0_);
		else if (i==np_-1)
			return -PI_TIMES_4_*gamma0_(np_ - 1)*grid_.FirstDerivativeRightSide(J0_);
		else
			return -PI_TIMES_4_*gamma0_(i)*grid_.FirstDerivativeSecondOrder(J0_, i);
		
	}

	double SP3_1D::G(const unsigned int i)
	{
		return PI_TIMES_4_*( J0_(i) - 2./3*J2_(i));
	}

	double SP3_1D::divq(const unsigned int i)
	{
		return kappa_(i)*(SIGMA_TIMES_4_*std::pow(T_(i), 4.) - G(i));
	}

	void SP3_1D::SpatialEquations(const double* y, const double t, double* f)
	{
		// Recover unknowns
		{
			unsigned int count = 0;
			for (unsigned int i = 0; i < np_; i++)
			{
				J0_(i) = y[count++];
				J2_(i) = y[count++];
			}
		}

		// Radiative equilibrium
		if (radiative_equilibrium_ == true)
		{
			for (unsigned int i = 0; i < np_; i++)
				T_(i) = std::pow(G(i) / SIGMA_TIMES_4_, 0.25);
		}

		// First derivative
		grid_.FirstDerivative(dummy_, J0_, dJ0_, Grid1DMultiPurpose::DERIVATIVE_CENTERED);
		grid_.FirstDerivative(dummy_, J2_, dJ2_, Grid1DMultiPurpose::DERIVATIVE_CENTERED);
		grid_.DiffusionTerm(J0_, gamma0_, diffJ0_);
		grid_.DiffusionTerm(J2_, gamma2_, diffJ2_);

		// Equations
		{
			unsigned int count = 0;

			// Internal wall
			{
				const double Ibw = SIGMA_OVER_PI_*pow(T_wall1_, 4.);
				const double Jw = Ibw*PI_ - (1.-epsilon_wall1_)/epsilon_wall1_*q(0);

				f[count++] = 0.50*(J0_(0) - Jw / PI_) - 1. / 8.*J2_(0) - gamma0_(0)*grid_.FirstDerivativeLeftSide(J0_);
				f[count++] = -1./8.*(J0_(0) - Jw / PI_) + 7. / 24.*J2_(0) - gamma2_(0)*grid_.FirstDerivativeLeftSide(J2_);
			}

			// Internal points
			for (unsigned int i = 1; i < np_ - 1; i++)
			{
				const double Ib = SIGMA_OVER_PI_*std::pow(T_(i), 4.);

				f[count++] = alfa0_(i)*(J0_(i) - 2. / 3.*J2_(i) - Ib) - diffJ0_(i) / beta_(i);
				f[count++] = (5. / 3.*alfa2_(i) + 4. / 3.*alfa0_(i))*J2_(i) - 2.*alfa0_(i)*(J0_(i) - Ib) - 3. / beta_(i)*diffJ2_(i);
			}

			// External wall
			{
				const double Ibw = SIGMA_OVER_PI_*pow(T_wall2_, 4.);
				const double Jw = Ibw*PI_ + (1. - epsilon_wall2_) / epsilon_wall2_*q(np_-1);

				f[count++] = 0.50*(J0_(np_ - 1) - Jw / PI_) - 1. / 8.*J2_(np_ - 1) + gamma0_(np_ - 1)*grid_.FirstDerivativeRightSide(J0_);
				f[count++] = -1. / 8.*(J0_(np_ - 1) - Jw / PI_) + 7. / 24.*J2_(np_ - 1) + gamma2_(np_ - 1)*grid_.FirstDerivativeRightSide(J2_);
			}
		}
	}

	class OpenSMOKE_MyNlsSystem_SP3
	{
	public:
		void assignSP3(OpenSMOKE::SP3_1D *SP3)
		{
			ptSP3 = SP3;
		}

	protected:

		unsigned int ne_;

		void MemoryAllocation()
		{
		}

		virtual void Equations(const Eigen::VectorXd &x, Eigen::VectorXd &f)
		{
			ptSP3->SpatialEquations(x.data(), 0, f.data());
		}

		void Jacobian(const Eigen::VectorXd &y, Eigen::MatrixXd &J)
		{
		};

		void Print(const int call, const double t, const double phiW, const Eigen::VectorXd &x, const Eigen::VectorXd &f)
		{
		}

	private:
		OpenSMOKE::SP3_1D *ptSP3;
	};

	void SP3_1D::SolveSpatialEquations()
	{
		unsigned int blockSize = 2;
		unsigned int bandSize = 2*blockSize-1;
		unsigned int neq = blockSize * np_;

		// Update properties
		UpdateProperties();

		// Initial guess
		Eigen::VectorXd yInitial(neq);
		unsigned int count = 0;
		for (unsigned int j = 0; j < np_; j++)
		{
			yInitial(count++) = SIGMA_*pow(T_(j), 4.);
			yInitial(count++) = SIGMA_*pow(T_(j), 4.);
		}

		// Create the system
		typedef NlsSMOKE::KernelBand<OpenSMOKE_MyNlsSystem_SP3> kernel;
		NlsSMOKE::NonLinearSolver<kernel> nls_solver;

		// Set initial conditions
		nls_solver.assignSP3(this);
		nls_solver.SetFirstGuessSolution(yInitial);
		nls_solver.SetBandSizes(bandSize, bandSize);

		// Solve the non linear system
		double timeStart = OpenSMOKEGetCpuTime();
		const unsigned int status = nls_solver();
		double timeEnd = OpenSMOKEGetCpuTime();

		if (status >= 0)
		{
			Eigen::VectorXd fResiduals(neq);

			// Update the current solution
			{
				Eigen::VectorXd yFinal(neq);
				nls_solver.Solution(yFinal, fResiduals);

				unsigned int count = 0;
				for (unsigned int j = 0; j < np_; j++)
				{
					J0_(j) = yFinal(count++);
					J2_(j) = yFinal(count++);
				}
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
/*				std::cout << " * number of iterations:           " << nls_object.IterationCounter() << std::endl;
				std::cout << " * number of functions (gradient): " << nls_object.NumFunctionsForGradient() << std::endl;
				std::cout << " * number of Newtons:              " << nls_object.NumNewtons() << std::endl;
				std::cout << " * number of Jacobians:            " << nls_object.NumNumericalJacobians() << std::endl;
				std::cout << " * number of functions (Jacobian): " << nls_object.NumFunctionsForNumericalJacobian() << std::endl;
				std::cout << " * number of functions (Jacobian): " << nls_object.NumFunctionsForNumericalJacobian() << std::endl;
				std::cout << " * residuals (norm2):              " << fResiduals.Norm2() << std::endl;
				std::cout << " * residuals (norm1):              " << fResiduals.Norm1() << std::endl;
*/				std::cout << std::endl;
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

	void SP3_1D::PrintFile(const boost::filesystem::path file_name)
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
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << gamma0_(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << A1_(i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << G(i);
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
