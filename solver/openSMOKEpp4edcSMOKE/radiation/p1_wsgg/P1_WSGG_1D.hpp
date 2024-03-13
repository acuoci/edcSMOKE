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
#include "P1_WSGG_1D.h"
#include "math/native-nls-solvers/NonLinearSystemSolver"

namespace OpenSMOKE
{
	const double P1_WSGG_1D::SIGMA_ = 5.670373e-8;
	const double P1_WSGG_1D::PI_ = boost::math::constants::pi<double>();
	const double P1_WSGG_1D::PI_TIMES_4_ = 4.*boost::math::constants::pi<double>();
	const double P1_WSGG_1D::PI_TIMES_2_ = 2.*boost::math::constants::pi<double>();
	const double P1_WSGG_1D::SIGMA_OVER_PI_ = P1_WSGG_1D::SIGMA_ / P1_WSGG_1D::PI_;
	const double P1_WSGG_1D::SIGMA_TIMES_4_ = 4.*P1_WSGG_1D::SIGMA_ ;

	P1_WSGG_1D::P1_WSGG_1D(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model) :
		grid_(*grid),
		wsgg_model_(model)
	{
		radiative_equilibrium_ = false;
		epsilon_wall1_ = 1.;
		epsilon_wall2_ = 1.;
		is_scattering_ = false;
		is_anisotropic_scattering_ = false;

		if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_SMITH)
		{
			nbands_ = 4;
			jMin_ = 1;
		}
		else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_TRUELOVE)
		{
			nbands_ = 4;
			jMin_ = 1;
		}
		else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_YIN)
		{
			nbands_ = 5;
			jMin_ = 1;
		}
		else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_SMITH_SOOT)
		{
			nbands_ = 8;
			jMin_ = 1;
		}
		else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_SMITH_SOOT_SINGLEBAND)
		{
			nbands_ = 1;
			jMin_ = 0;
		}
                
                else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_CASSOL)
                {
			nbands_ = 100;
			jMin_ = 0;
                }

		np_ = grid_.N();
		MemoryAllocation();
	}

	void P1_WSGG_1D::MemoryAllocation()
	{
		T_.resize(np_);
		PCO2_Pa_.resize(np_);
		PH2O_Pa_.resize(np_);
                fv_.resize(np_);

		kappa_.resize(nbands_);
		gamma_.resize(nbands_);
		a_.resize(nbands_);
		G_.resize(nbands_);
		dG_.resize(nbands_);
		diffG_.resize(nbands_);

		for (unsigned int i = 0; i < nbands_; i++)
		{
			kappa_[i].resize(np_);
			gamma_[i].resize(np_);
			a_[i].resize(np_);
			G_[i].resize(np_);
			dG_[i].resize(np_);
			diffG_[i].resize(np_);
		}

		for (unsigned int i = 0; i < nbands_; i++)
		{
			kappa_[i].setConstant(0.);
			gamma_[i].setConstant(0.);
			a_[i].setConstant(0.);
			G_[i].setConstant(0.);
			dG_[i].setConstant(0.);
			diffG_[i].setConstant(0.);
		}
	}

	void P1_WSGG_1D::SetRadiativeEquilibrium(const bool flag)
	{
		radiative_equilibrium_ = flag;
	}

	void P1_WSGG_1D::SetTemperature(const double T)
	{ 
		T_.setConstant(T); 
	}

	void P1_WSGG_1D::SetTemperature(const Eigen::VectorXd& T)
	{
		if (T.size() != np_)
			OpenSMOKE::FatalErrorMessage("P1: wrong size in temperature field");
		T_ = T;
	}

	void P1_WSGG_1D::SetPressureCO2(const double PCO2_Pa)
	{
		PCO2_Pa_.setConstant(PCO2_Pa);
	}

	void P1_WSGG_1D::SetPressureCO2(const Eigen::VectorXd& PCO2_Pa)
	{
		if (PCO2_Pa.size() != np_)
			OpenSMOKE::FatalErrorMessage("P1: wrong size in pressure field");
		PCO2_Pa_ = PCO2_Pa;
	}

	void P1_WSGG_1D::SetPressureH2O(const double PH2O_Pa)
	{
		PH2O_Pa_.setConstant(PH2O_Pa);
	}

	void P1_WSGG_1D::SetSootVolumeFraction(const double fv)
	{
		fv_.setConstant(fv);
	}

	void P1_WSGG_1D::SetSootVolumeFraction(const Eigen::VectorXd& fv)
	{
		if (fv.size() != np_)
			OpenSMOKE::FatalErrorMessage("P1: wrong size in soot volume fraction field");

		fv_ = fv;
	}

	void P1_WSGG_1D::SetPressureH2O(const Eigen::VectorXd& PH2O_Pa)
	{
		if (PH2O_Pa.size() != np_)
			OpenSMOKE::FatalErrorMessage("P1: wrong size in pressure field");
		PH2O_Pa_ = PH2O_Pa;
	}

	void P1_WSGG_1D::UpdateWSGGCoefficients()
	{
		double* kappa = new double[nbands_];
		double* a = new double[nbands_];

		for (unsigned int i = 0; i < np_; i++)
		{
			if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_TRUELOVE)
				WSGG_Truelove(T_(i), PCO2_Pa_(i), PH2O_Pa_(i), kappa, a);
			else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_SMITH)
				WSGG_Smith(T_(i), PCO2_Pa_(i), PH2O_Pa_(i), kappa, a);
			else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_YIN)
				WSGG_Yin(T_(i), PCO2_Pa_(i), PH2O_Pa_(i), kappa, a);
			else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_SMITH_SOOT_SINGLEBAND)
				WSGG_Smith_Soot_Singleband(T_(i), PCO2_Pa_(i), PH2O_Pa_(i), fv_(i), kappa, a);
			else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_SMITH_SOOT)
				WSGG_Smith_Soot(T_(i), PCO2_Pa_(i), PH2O_Pa_(i), fv_(i), kappa, a);
                        else if (wsgg_model_ == EXTINCTION_COEFFICIENT_WSGG_CASSOL)
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

	void P1_WSGG_1D::UpdateProperties()
	{
		UpdateWSGGCoefficients();

		for (unsigned int j = jMin_; j < nbands_; j++)
			for (unsigned int i = 0; i < np_; i++)
				gamma_[j](i) = 1. / (3.*kappa_[j](i));
	}

	double P1_WSGG_1D::q_from_wall1()
	{
		if (radiative_equilibrium_ == false)
		{
			double sum = 0.;
			for (unsigned int j = jMin_; j < nbands_; j++)
				sum += -gamma_[j](0)*grid_.FirstDerivativeSecondOrderLeftSide(G_[j]);
			return  sum;
		}
		else
		{
			double sum = 0.;
			for (unsigned int j = jMin_; j < nbands_; j++)
				sum += -gamma_[j](0)*grid_.FirstDerivativeLeftSide(G_[j]);
			return  sum;
		}
	}

	double P1_WSGG_1D::q_from_wall2()
	{
		if (radiative_equilibrium_ == false)
		{
			double sum = 0.;
			for (unsigned int j = jMin_; j < nbands_; j++)
				sum += gamma_[j](np_ - 1)*grid_.FirstDerivativeSecondOrderRightSide(G_[j]);
			return  sum;
		}
		else
		{
			double sum = 0.;
			for (unsigned int j = jMin_; j < nbands_; j++)
				sum += gamma_[j](np_ - 1)*grid_.FirstDerivativeRightSide(G_[j]);
			return  sum;
		}

	}

	void P1_WSGG_1D::SetWallsTemperatures(const double T_wall1, const double T_wall2)
	{
		T_wall1_ = T_wall1;
		T_wall2_ = T_wall2;
	}

	void P1_WSGG_1D::SetWallsEmissivities(const double epsilon_wall1, const double epsilon_wall2)
	{
		epsilon_wall1_ = epsilon_wall1;
		epsilon_wall2_ = epsilon_wall2;
	}

	double P1_WSGG_1D::kappa(const unsigned int i)
	{
		double sum = 0.;
		for (unsigned int j = 0; j < nbands_; j++)
			sum += a_[j](i)*kappa_[j](i);
		return sum;
	}

	double P1_WSGG_1D::a(const unsigned int i)
	{
		double sum = 0.;
		for (unsigned int j = 0; j < nbands_; j++)
			sum += a_[j](i);
		return sum;
	}

	double P1_WSGG_1D::q(const unsigned int i)
	{
		if (radiative_equilibrium_ == false)
		{
			double sum = 0.;
			if (i == 0)
			{
				for (unsigned int j = jMin_; j < nbands_; j++)
					sum += -gamma_[j](0)*grid_.FirstDerivativeSecondOrderLeftSide(G_[j]);
			}
			else if (i == np_ - 1)
			{
				for (unsigned int j = jMin_; j < nbands_; j++)
					sum += -gamma_[j](np_ - 1)*grid_.FirstDerivativeSecondOrderRightSide(G_[j]);
			}
			else
			{
				for (unsigned int j = jMin_; j < nbands_; j++)
					sum += -gamma_[j](i)*grid_.FirstDerivativeSecondOrder(G_[j], i);
			}
			return sum;
		}
		else
		{
			double sum = 0.;
			if (i == 0)
			{
				for (unsigned int j = jMin_; j < nbands_; j++)
					sum += -gamma_[j](0)*grid_.FirstDerivativeLeftSide(G_[j]);
			}
			else if (i == np_ - 1)
			{
				for (unsigned int j = jMin_; j < nbands_; j++)
					sum += -gamma_[j](np_ - 1)*grid_.FirstDerivativeRightSide(G_[j]);
			}
			else
			{
				for (unsigned int j = jMin_; j < nbands_; j++)
					sum += -gamma_[j](i)*grid_.FirstDerivativeSecondOrder(G_[j], i);
			}
			return sum;
		}
	}

	double P1_WSGG_1D::G(const unsigned int i)
	{
		double sum = 0.;
		for (unsigned int j = jMin_; j < nbands_; j++)
			sum += G_[j](i);
		return sum;
	}

	double P1_WSGG_1D::divq(const unsigned int i)
	{
		double sum = 0.;
		for (unsigned int j = jMin_; j < nbands_; j++)
			sum += kappa_[j](i)*(SIGMA_TIMES_4_*a_[j](i)*std::pow(T_(i), 4.) - G_[j](i));
		return sum;
	}

	void P1_WSGG_1D::SpatialEquationsGlobal(const double* y, const double t, double* f)
	{
		// Recovering variables
		{
			unsigned int count = 0;
			for (unsigned int i = 0; i < np_; i++)
			for (unsigned int j = 0; j < nbands_; j++)
				G_[j](i) = y[count++];
		}

		// Radiative equilibrium
		if (radiative_equilibrium_ == true)
		{
			for (unsigned int i = 0; i < np_; i++)
				T_(i) = std::pow(G(i) / SIGMA_TIMES_4_, 0.25);

			UpdateProperties();
		}

		// First derivative
		for (unsigned int j = 0; j < nbands_; j++)
		{
			grid_.FirstDerivative(dummy_, G_[j], dG_[j], Grid1DMultiPurpose::DERIVATIVE_CENTERED);
			grid_.DiffusionTerm(G_[j], gamma_[j], diffG_[j]);
		}
		
		// Equations
		{
			unsigned int count = 0;

			// Internal wall
			for (unsigned int j = 0; j < nbands_; j++)
				f[count++] = -(2. - epsilon_wall1_) / epsilon_wall1_ * 2.*gamma_[j](0) * dG_[j](0) + G_[j](0) - SIGMA_TIMES_4_*a_[j](0)*pow(T_wall1_, 4.);

			// Internal points
			for (unsigned int i = 1; i < np_ - 1; i++)
				for (unsigned int j = 0; j < nbands_; j++)
					f[count++] = kappa_[j](i) * (G_[j](i) - SIGMA_TIMES_4_*a_[j](i)*pow(T_(i), 4.)) - diffG_[j](i);

			// External wall
			for (unsigned int j = 0; j < nbands_; j++)
				f[count++] = (2. - epsilon_wall2_) / epsilon_wall2_ * 2.*gamma_[j](np_ - 1)*dG_[j](np_ - 1) + G_[j](np_ - 1) - SIGMA_TIMES_4_*a_[j](np_-1)*pow(T_wall2_, 4.);
		}
	}

	void P1_WSGG_1D::SpatialEquationsPerBand(const unsigned int j, const double* y, const double t, double* f)
	{
		// Recovering variables
		for (unsigned int i = 0; i < np_; i++)
			G_[j](i) = y[i];

		// Radiative equilibrium
		if (radiative_equilibrium_ == true)
			OpenSMOKE::FatalErrorMessage("P1_WSGG_1D: The segregated approach cannot be used for solving radiative equilibrium problmes");

		// First derivative
		grid_.FirstDerivativeAccurateBC(dummy_, G_[j], dG_[j], Grid1DMultiPurpose::DERIVATIVE_CENTERED);
		grid_.DiffusionTerm(G_[j], gamma_[j], diffG_[j]);

		// Equations
		{
			// Internal wall
			f[0] = -(2. - epsilon_wall1_) / epsilon_wall1_ * 2.*gamma_[j](0) * dG_[j](0) + G_[j](0) - SIGMA_TIMES_4_*a_[j](0)*pow(T_wall1_, 4.);

			// Internal points
			for (unsigned int i = 1; i < np_ - 1; i++)
				f[i] = kappa_[j](i) * (G_[j](i) - SIGMA_TIMES_4_*a_[j](i)*pow(T_(i), 4.)) - diffG_[j](i);

			// External wall
			f[np_-1] = (2. - epsilon_wall2_) / epsilon_wall2_ * 2.*gamma_[j](np_ - 1)*dG_[j](np_ - 1) + G_[j](np_ - 1) - SIGMA_TIMES_4_*a_[j](np_ - 1)*pow(T_wall2_, 4.);
		}
	}

	class OpenSMOKE_MyNlsSystem_P1WSGGGlobal
	{
	public:
		void assignP1WSGG(OpenSMOKE::P1_WSGG_1D *P1WSGG)
		{
			ptP1WSGG = P1WSGG;
		}

	protected:

		unsigned int ne_;

		void MemoryAllocation()
		{
		}

		virtual void Equations(const Eigen::VectorXd &x, Eigen::VectorXd &f)
		{
			ptP1WSGG->SpatialEquationsGlobal(x.data(), 0, f.data());
		}

		void Jacobian(const Eigen::VectorXd &y, Eigen::MatrixXd &J)
		{
		};

		void Print(const int call, const double t, const double phiW, const Eigen::VectorXd &x, const Eigen::VectorXd &f)
		{
		}

	private:
		OpenSMOKE::P1_WSGG_1D *ptP1WSGG;
	};

	class OpenSMOKE_MyNlsSystem_P1WSGGperBand
	{
	public:
		void assignP1WSGG(OpenSMOKE::P1_WSGG_1D *P1WSGG)
		{
			ptP1WSGG = P1WSGG;
		}

		void SetBand(const unsigned int jBand)
		{
			jBand_ = jBand;
		}

	protected:

		unsigned int ne_;

		void MemoryAllocation()
		{
		}

		virtual void Equations(const Eigen::VectorXd &x, Eigen::VectorXd &f)
		{
			ptP1WSGG->SpatialEquationsPerBand(jBand_, x.data(), 0, f.data());
		}

		void Jacobian(const Eigen::VectorXd &y, Eigen::MatrixXd &J)
		{
		};

		void Print(const int call, const double t, const double phiW, const Eigen::VectorXd &x, const Eigen::VectorXd &f)
		{
		}

	private:
		OpenSMOKE::P1_WSGG_1D *ptP1WSGG;
		unsigned int jBand_;
	};

	void P1_WSGG_1D::SolveSpatialEquations()
	{
		if (radiative_equilibrium_ == true)
			SolveSpatialEquationsGlobal();
		else
			SolveSpatialEquationsPerBand();
	}

	void P1_WSGG_1D::SolveSpatialEquationsGlobal()
	{
		unsigned int blockSize = nbands_;
		unsigned int bandSize = 2 * blockSize - 1;
		unsigned int neq = blockSize * np_;

		// Update
		UpdateProperties();

		// Initial guess
		Eigen::VectorXd yInitial(neq);
		unsigned int count = 0;
		for (unsigned int i = 0; i < np_; i++)
			for (unsigned int j = 0; j < nbands_; j++)
				yInitial(count++) = SIGMA_TIMES_4_*a_[j](i)*pow(T_(i), 4.);

		// Create the system
		typedef NlsSMOKE::KernelBand<OpenSMOKE_MyNlsSystem_P1WSGGGlobal> kernel;
		NlsSMOKE::NonLinearSolver<kernel> nls_solver;

		// Set initial conditions
		nls_solver.assignP1WSGG(this);
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
				for (unsigned int i = 0; i < np_; i++)
					for (unsigned int j = 0; j < nbands_; j++)
						G_[j](i) = yFinal(count++);
			}

			bool verbose = true;
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

	void P1_WSGG_1D::SolveSpatialEquationsPerBand()
	{
		// Update
		UpdateProperties();

		// Memory allocation
		Eigen::VectorXd yInitial(np_);
		Eigen::VectorXd fResiduals(np_);
		Eigen::VectorXd yFinal(np_);

		for (unsigned int j = jMin_; j < nbands_; j++)
		{
			for (unsigned int i = 0; i < np_; i++)
				yInitial(i) = SIGMA_TIMES_4_*a_[j](i)*pow(T_(i), 4.);

			// Create the system
			typedef NlsSMOKE::KernelBand<OpenSMOKE_MyNlsSystem_P1WSGGperBand> kernel;
			NlsSMOKE::NonLinearSolver<kernel> nls_solver;

			// Set initial conditions
			nls_solver.assignP1WSGG(this);
			nls_solver.SetBand(j);
			nls_solver.SetFirstGuessSolution(yInitial);
			nls_solver.SetBandSizes(2, 2);

			// Solve the non linear system
			double timeStart = OpenSMOKEGetCpuTime();
			const unsigned int status = nls_solver();
			double timeEnd = OpenSMOKEGetCpuTime();

			if (status >= 0)
			{
				// Update the current solution
				{
					nls_solver.Solution(yFinal, fResiduals);

					for (unsigned int i = 0; i < np_; i++)
						G_[j](i) = yFinal(i);
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
/*					std::cout << " * number of iterations:           " << nls_object.IterationCounter() << std::endl;
					std::cout << " * number of functions (gradient): " << nls_object.NumFunctionsForGradient() << std::endl;
					std::cout << " * number of Newtons:              " << nls_object.NumNewtons() << std::endl;
					std::cout << " * number of Jacobians:            " << nls_object.NumNumericalJacobians() << std::endl;
					std::cout << " * number of functions (Jacobian): " << nls_object.NumFunctionsForNumericalJacobian() << std::endl;
					std::cout << " * number of functions (Jacobian): " << nls_object.NumFunctionsForNumericalJacobian() << std::endl;
					std::cout << " * residuals (norm2):              " << fResiduals.Norm2() << std::endl;
					std::cout << " * residuals (norm1):              " << fResiduals.Norm1() << std::endl;
*/					std::cout << std::endl;
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
	}

	void P1_WSGG_1D::PrintFile(const boost::filesystem::path file_name)
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
		fOut << std::endl;

		for (unsigned int i = 0; i < np_; i++)
		{
			double q_ = q(i); 
			double divq_ = divq(i);

			fOut << std::scientific << std::setw(16) << std::setprecision(6) << grid_.x[i+1];
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
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << G_[0](i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << G_[1](i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << G_[2](i);
			fOut << std::scientific << std::setw(16) << std::setprecision(6) << G_[3](i);
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
