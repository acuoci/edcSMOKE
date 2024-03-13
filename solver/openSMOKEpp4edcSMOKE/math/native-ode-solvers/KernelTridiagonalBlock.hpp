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
|   Copyright(C) 2023  Alberto Cuoci                                      |
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

namespace OdeSMOKE
{
	template <typename ODESystemObject>
	KernelTridiagonalBlock<ODESystemObject>::KernelTridiagonalBlock()
	{
	}

	template <typename ODESystemObject>
	void KernelTridiagonalBlock<ODESystemObject>::ResetKernel()
	{
		solverType_ = OpenSMOKE::SOLVER_DENSE_EIGEN;
		pre_processing_ = false;

		numberOfFunctionCallsForJacobian_ = 0;

		cpuTimeToAssembleJacobian_ = 0.;
		cpuTimeToFactorize_ = 0.;
		cpuTimeToSolveLinearSystem_ = 0.;

		cpuTimeSingleJacobianAssembling_ = 0.;
		cpuTimeSingleFactorization_ = 0.;
		cpuTimeSingleLinearSystemSolution_ = 0.;
	}

	template <typename ODESystemObject>
	void KernelTridiagonalBlock<ODESystemObject>::MemoryAllocationKernel()
	{
		// Allocate memory (if needed) for the DaeSystem object
		this->MemoryAllocation();

		// Internal variables
		f_plus_.resize(this->ne_);
		y_plus_.resize(this->ne_);
		hJ_.resize(this->ne_);
		df_over_dy_.resize(this->ne_);
	}

	template <typename ODESystemObject>
	KernelTridiagonalBlock<ODESystemObject>::~KernelTridiagonalBlock()
	{
		J_->DestroyMat();
		G_->DestroyMat();
	}

	template <typename ODESystemObject>
	void KernelTridiagonalBlock<ODESystemObject>::SetTridiagonalBlockSize(const unsigned int block_size)
	{
		J_ = new OpenSMOKE::OpenSMOKETridiagonalBlockMatrixDouble(this->ne_, block_size);
		G_ = new OpenSMOKE::OpenSMOKETridiagonalBlockMatrixDouble(this->ne_, block_size);
		if (J_ == NULL || G_ == NULL)
			OpenSMOKE::FatalErrorMessage("Memory allocation for tridiagonal-block matrix");

		J_->SetToZero();
		G_->SetToZero();

		// Sparsity data
		J_->GetSparsityData(group_cols_, group_rows_);
	}

	template <typename ODESystemObject>
	void KernelTridiagonalBlock<ODESystemObject>::PreProcessing()
	{
		pre_processing_ = true;
	}

	template <typename ODESystemObject>
	void KernelTridiagonalBlock<ODESystemObject>::SetLinearAlgebraSolver(const std::string linear_algebra_solver)
	{
		if (linear_algebra_solver == "Eigen")
		{
			solverType_ = OpenSMOKE::SOLVER_DENSE_EIGEN;
		}
		else
			OpenSMOKE::ErrorMessage("KernelTridiagonalBlock<ODESystemObject>", "Requested linear algebra is not supported!");
	}

	template <typename ODESystemObject>
	void KernelTridiagonalBlock<ODESystemObject>::JacobianTimesVector(const Eigen::VectorXd& v_in, Eigen::VectorXd* v_out)
	{
		J_->Product(v_in.data(), v_out->data());
	}

	template <typename ODESystemObject>
	void KernelTridiagonalBlock<ODESystemObject>::UserDefinedJacobian(const Eigen::VectorXd& y, const double t)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		// TODO
		// this->Jacobian(y, t, J_);
		OpenSMOKE::ErrorMessage("KernelTridiagonalBlock<ODESystemObject>", "User defined Jacobian is still not supported!");

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleJacobianAssembling_ = tend - tstart;
		cpuTimeToAssembleJacobian_ += cpuTimeSingleJacobianAssembling_;
	}

	template <typename ODESystemObject>
	void KernelTridiagonalBlock<ODESystemObject>::NumericalJacobian(Eigen::VectorXd& y, const double t, const Eigen::VectorXd& f, const double h, const Eigen::VectorXd& e,
									const bool max_constraints, const Eigen::VectorXd& yMax)
	{
		if (pre_processing_ == false)
			PreProcessing();

		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double BETA = 1.e+3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;

		double hf = BETA * std::fabs(h) * OpenSMOKE::ErrorControl(f, e) * double(this->ne_);
		if (hf < 1.e-10)
			hf = 1.;

		// Save the original vector
		df_over_dy_.setZero();
		y_plus_ = y;

		// Loop
		for (int group = 0; group < group_cols_.size(); group++)
		{
			if (max_constraints == false)
			{
				for (int i = 0; i < group_cols_[group].size(); i++)
				{
					const unsigned int j = group_cols_[group][i];

					const double yh = y(j);
					const double hJf = hf / e(j);

					hJ_(j) = ETA2*std::fabs(std::max(yh, 1. / e(j)));
					hJ_(j) = std::max(hJ_(j), hJf);
					hJ_(j) = std::max(hJ_(j), ZERO_DER);
					hJ_(j) = std::min(hJ_(j), 0.001 + 0.001*std::fabs(yh));

					y_plus_(j) += hJ_(j);
				}
			}
			else
			{
				for (int i = 0; i < group_cols_[group].size(); i++)
				{
					const unsigned int j = group_cols_[group][i];

					const double yh = y(j);
					const double hJf = hf / e(j);

					hJ_(j) = ETA2*std::fabs(std::max(yh, 1. / e(j)));
					hJ_(j) = std::max(hJ_(j), hJf);
					hJ_(j) = std::max(hJ_(j), ZERO_DER);
					hJ_(j) = std::min(hJ_(j), 0.001 + 0.001*std::fabs(yh));

					if (yh + hJ_(j) > yMax(j))
						hJ_(j) = -hJ_(j);

					y_plus_(j) += hJ_(j);
				}
			}

			this->Equations(y_plus_, t, f_plus_);
			numberOfFunctionCallsForJacobian_++;

			for (int i = 0; i < group_cols_[group].size(); i++)
			{
				// Set Jacobian
				const unsigned int j = group_cols_[group][i];
				const unsigned int row = group_rows_[group][i][0];
				const unsigned int n = group_rows_[group][i].size();
				
				for (int k=row;k<row+n;k++)
					df_over_dy_(k) = (f_plus_(k) - f(k)) / hJ_(j);

				J_->SetColumn(j, row, n, df_over_dy_.data());

				// Reset
				y_plus_(j) = y(j);
				for (int k=row;k<row+n;k++)
					df_over_dy_(k) = 0.;
			}
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleJacobianAssembling_ = tend - tstart;
		cpuTimeToAssembleJacobian_ += cpuTimeSingleJacobianAssembling_;
	}

	template <typename ODESystemObject>
	void KernelTridiagonalBlock<ODESystemObject>::BuildAndFactorizeMatrixG(const double hr0)
	{
		// 1. Assembling the G matrix
		// ----------------------------------------------------------------------------------------
		J_->CopyTo(G_);			// G_ = J_
		G_->Scale(-hr0);		// G_ = -hr0*J_
		G_->AddIdentity();		// G_ = I - hr0J_

		// 2. Factorizing the G matrix
		// ----------------------------------------------------------------------------------------
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();
		G_->Factorize();
		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleFactorization_ = tend - tstart;
		cpuTimeToFactorize_ += cpuTimeSingleFactorization_;
	}

	template <typename ODESystemObject>
	void KernelTridiagonalBlock<ODESystemObject>::SolveLinearSystem(Eigen::VectorXd& db)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();
		G_->Solve(db.data());
		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleLinearSystemSolution_ = tend - tstart;
		cpuTimeToSolveLinearSystem_ += cpuTimeSingleLinearSystemSolution_;
	}

	template <typename ODESystemObject>
	void KernelTridiagonalBlock<ODESystemObject>::OdeSolverKernelSummary(std::ostream& out)
	{
		double totalCpu = cpuTimeToAssembleJacobian_ + cpuTimeToFactorize_ + cpuTimeToSolveLinearSystem_;
		double totalSingleCpu = cpuTimeSingleFactorization_ + cpuTimeSingleJacobianAssembling_ + cpuTimeSingleLinearSystemSolution_;

		out << std::endl;
		out << "Data for the TridiagonalBlock ODE solver Kernel" << std::endl;
		out << "---------------------------------------------------------------------------------------------------------" << std::endl;
		out << "Number of function calls (only to assemble Jacobian): " << numberOfFunctionCallsForJacobian_ << " (" << numberOfFunctionCallsForJacobian_ / this->ne_ << ")" << std::endl;
		out << "Cumulative CPU time for assembling Jacobian:          " << cpuTimeToAssembleJacobian_ << " (" << cpuTimeToAssembleJacobian_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU time for LU decomposition:             " << cpuTimeToFactorize_ << " (" << cpuTimeToFactorize_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU time for solving the linear system:    " << cpuTimeToSolveLinearSystem_ << " (" << cpuTimeToSolveLinearSystem_ / totalCpu *100. << "%)" << std::endl;
		out << "CPU time for assembling Jacobian:                     " << cpuTimeSingleJacobianAssembling_ << " (" << cpuTimeSingleFactorization_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU time for LU decomposition:                        " << cpuTimeSingleFactorization_ << " (" << cpuTimeSingleJacobianAssembling_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU time for solving the linear system:               " << cpuTimeSingleLinearSystemSolution_ << " (" << cpuTimeSingleLinearSystemSolution_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "---------------------------------------------------------------------------------------------------------" << std::endl;
		out << std::endl;
	}
}
