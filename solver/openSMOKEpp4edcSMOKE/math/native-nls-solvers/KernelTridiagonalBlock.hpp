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

namespace NlsSMOKE
{
	template <typename NLSSystemObject>
	KernelTridiagonalBlock<NLSSystemObject>::KernelTridiagonalBlock()
	{
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::ResetKernel()
	{
		solverType_ = OpenSMOKE::SOLVER_DENSE_EIGEN;
		pre_processing_ = false;

		numberOfSystemCallsForJacobian_ = 0;
		numberOfJacobianFullAssembling_ = 0;
		numberOfJacobianQuasiAssembling_ = 0;
		numberOfJacobianFactorizations_ = 0;
		numberOfLinearSystemSolutions_ = 0;

		cpuTimeJacobianFullAssembling_ = 0.;
		cpuTimeJacobianQuasiAssembling_ = 0.;
		cpuTimeJacobianFactorization_ = 0.;
		cpuTimeLinearSystemSolution_ = 0.;

		cpuTimeSingleJacobianFullAssembling_ = 0.;
		cpuTimeSingleJacobianQuasiAssembling_ = 0.;
		cpuTimeSingleJacobianFactorization_ = 0.;
		cpuTimeSingleLinearSystemSolution_ = 0.;
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::MemoryAllocationKernel()
	{
		// Allocate memory (if needed) for the DaeSystem object
		this->MemoryAllocation();

		// Internal variables
		hJ_.resize(this->ne_);
		x_plus_.resize(this->ne_);
		f_plus_.resize(this->ne_);
		df_over_dx_.resize(this->ne_);

		// Set zero
		hJ_.setZero();
		x_plus_.setZero();
		f_plus_.setZero();
		df_over_dx_.setZero();
	}

	template <typename NLSSystemObject>
	KernelTridiagonalBlock<NLSSystemObject>::~KernelTridiagonalBlock()
	{
		J_->DestroyMat();
		J_factorized_->DestroyMat();
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::SetTridiagonalBlockSize(const unsigned int blockSize)
	{
		J_ = new OpenSMOKE::OpenSMOKETridiagonalBlockMatrixDouble(this->ne_, blockSize);
		if (J_ == NULL)
			OpenSMOKE::FatalErrorMessage("Memory allocation for TridiagonalBlock matrix");

		J_factorized_ = new OpenSMOKE::OpenSMOKETridiagonalBlockMatrixDouble(this->ne_, blockSize);
		if (J_factorized_ == NULL)
			OpenSMOKE::FatalErrorMessage("Memory allocation for TridiagonalBlock matrix");

		J_->SetToZero();
		J_factorized_->SetToZero();

		// Sparsity data
		J_->GetSparsityData(group_cols_, group_rows_);


		numberOfSystemCallsPerJacobian_ = group_cols_.size();
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::PreProcessing()
	{
		pre_processing_ = true;
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::SetLinearAlgebraSolver(const std::string linear_algebra_solver)
	{
		if (linear_algebra_solver == "Eigen")
		{
			solverType_ = OpenSMOKE::SOLVER_DENSE_EIGEN;
		}
		else
			OpenSMOKE::ErrorMessage("KernelTridiagonalBlock<NLSSystemObject>", "Requested linear algebra is not supported!");
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::UserDefinedJacobian(const Eigen::VectorXd& y, const double t)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		// TODO
		// this->Jacobian(y, t, J_);
		OpenSMOKE::ErrorMessage("KernelTridiagonalBlock<NLSSystemObject>", "User defined Jacobian is still not supported!");

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfSystemCallsForJacobian_ += this->pattern_.number_groups();
		numberOfJacobianFullAssembling_++;
		cpuTimeSingleJacobianFullAssembling_ = tend - tstart;
		cpuTimeJacobianFullAssembling_ += cpuTimeSingleJacobianFullAssembling_;
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::NumericalJacobian(const Eigen::VectorXd& x, const Eigen::VectorXd& f, Eigen::VectorXd& x_dimensions, const bool max_constraints, const Eigen::VectorXd& xMax)
	{
		if (pre_processing_ == false)
			PreProcessing();

		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		const double ZERO_DER = 1.e-8;
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);

		// Save the original vector
		df_over_dx_.setZero();
		x_plus_ = x;

		// Loop
		for (int group = 0; group < group_cols_.size(); group++)
		{
			if (max_constraints == false)
			{
				for (int i = 0; i < group_cols_[group].size(); i++)
				{
					const unsigned int j = group_cols_[group][i];

					const double xh = std::fabs(x(j));
					const double xdh = std::fabs(x_dimensions(j));
					hJ_(j) = ETA2*std::max(xh, xdh);
					hJ_(j) = std::max(hJ_(j), ZERO_DER);
					hJ_(j) = std::min(hJ_(j), 0.001 + 0.001*std::fabs(xh));

					x_plus_(j) += hJ_(j);
				}
			}
			else
			{
				for (int i = 0; i < group_cols_[group].size(); i++)
				{
					const unsigned int j = group_cols_[group][i];

					const double xh = std::fabs(x(j));
					const double xdh = std::fabs(x_dimensions(j));
					hJ_(j) = ETA2*std::max(xh, xdh);
					hJ_(j) = std::max(hJ_(j), ZERO_DER);
					hJ_(j) = std::min(hJ_(j), 0.001 + 0.001*std::fabs(xh));

					if (xh + hJ_(j) > xMax(j))
						hJ_(j) = -hJ_(j);

					x_plus_(j) += hJ_(j);
				}
			}

			this->Equations(x_plus_, f_plus_);

			for (int i = 0; i < group_cols_[group].size(); i++)
			{
				// Set Jacobian
				const unsigned int j = group_cols_[group][i];
				const unsigned int row = group_rows_[group][i][0];
				const unsigned int n = group_rows_[group][i].size();

				for (int k=row;k<row+n;k++)
					df_over_dx_(k) = (f_plus_[k] - f[k]) / hJ_(j);

				J_->SetColumn(j, row, n, df_over_dx_.data());

				// Reset
				x_plus_(j) = x(j);
				for (int k=row;k<row+n;k++)
					df_over_dx_(k) = 0.;
			}
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfSystemCallsForJacobian_ += group_cols_.size();
		numberOfJacobianFullAssembling_++;
		cpuTimeSingleJacobianFullAssembling_ = tend - tstart;
		cpuTimeJacobianFullAssembling_ += cpuTimeSingleJacobianFullAssembling_;
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::DoubleProduct(Eigen::VectorXd& gi, Eigen::VectorXd& aux)
	{
		J_->TProduct(aux.data(), gi.data());
		J_->Product(gi.data(), aux.data());
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::QuasiNewton(const Eigen::VectorXd& dxi, const Eigen::VectorXd& dfi)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		{
			// The auxiliary vector named x_plus is used here
			Eigen::VectorXd* normSquared = &x_plus_;
			
			normSquared->setZero();
			for (int group = 0; group < group_cols_.size(); group++)
			{
				for (int i = 0; i < group_cols_[group].size(); i++)
				{
					// Set Jacobian
					const unsigned int j = group_cols_[group][i];
					const unsigned int row = group_rows_[group][i][0];
					const unsigned int n = group_rows_[group][i].size();

					for (int k=row;k<row+n;k++)
						(*normSquared)(k) += dxi(j)*dxi(j);
				}
			}

			// The auxiliary vector named x_plus is used here
			Eigen::VectorXd* sum_vector = &f_plus_;

			(*sum_vector) = dfi;

			for (int group = 0; group < group_cols_.size(); group++)
			{
				for (int i = 0; i < group_cols_[group].size(); i++)
				{
					// Set Jacobian
					const unsigned int j = group_cols_[group][i];
					const unsigned int row = group_rows_[group][i][0];
					const unsigned int n = group_rows_[group][i].size();

					for (int k=row;k<row+n;k++)
						(*sum_vector)(k) -= J_->Get(k,j)*dxi(j);
				}
			}

			const double eps = 1.e-10;
			for (int j = 0; j < static_cast<int>(this->ne_); j++)
				(*sum_vector)(j) /= ((*normSquared)(j) + eps);

			for (int group = 0; group < group_cols_.size(); group++)
			{
				for (int i = 0; i < group_cols_[group].size(); i++)
				{
					// Set Jacobian
					const unsigned int j = group_cols_[group][i];
					const unsigned int row = group_rows_[group][i][0];
					const unsigned int n = group_rows_[group][i].size();

					for (int k=row;k<row+n;k++)
					{
						const double value = J_->Get(k,j) + (*sum_vector)(k)*dxi(j);
						J_->Set(k,j,value);
					}
				}
			}
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfJacobianQuasiAssembling_++;
		cpuTimeSingleJacobianQuasiAssembling_ = tend - tstart;
		cpuTimeJacobianQuasiAssembling_ += cpuTimeSingleJacobianQuasiAssembling_;
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::Factorize()
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		J_->CopyTo(J_factorized_);
		J_factorized_->Factorize();

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfJacobianFactorizations_++;
		cpuTimeSingleJacobianFactorization_ = tend - tstart;
		cpuTimeJacobianFactorization_ += cpuTimeSingleJacobianFactorization_;
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::Solve(Eigen::VectorXd& pi)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		J_factorized_->Solve(pi.data());

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfLinearSystemSolutions_++;
		cpuTimeSingleLinearSystemSolution_ = tend - tstart;
		cpuTimeLinearSystemSolution_ += cpuTimeSingleLinearSystemSolution_;
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::CalculatesNormOfJacobianRows(Eigen::VectorXd& row_norms)
	{
		row_norms.setZero();

		for (int group = 0; group < group_cols_.size(); group++)
		{
			for (int i = 0; i < group_cols_[group].size(); i++)
			{
				// Set Jacobian
				const unsigned int j = group_cols_[group][i];
				const unsigned int row = group_rows_[group][i][0];
				const unsigned int n = group_rows_[group][i].size();

				for (int k=row;k<row+n;k++)
				{
					const double J = J_->Get(k,j);
					row_norms(k) += J*J;
				}
			}
		}

		for (unsigned int i = 0; i < this->ne_; i++)
			row_norms(i) = std::sqrt(row_norms(i));
	}

	template <typename NLSSystemObject>
	void KernelTridiagonalBlock<NLSSystemObject>::NlsSolverKernelSummary(std::ostream& out)
	{
		const double totalCpu = cpuTimeJacobianFullAssembling_ + cpuTimeJacobianQuasiAssembling_ + cpuTimeJacobianFactorization_ + cpuTimeLinearSystemSolution_;
		const double totalSingleCpu = cpuTimeSingleJacobianFullAssembling_ + cpuTimeSingleJacobianQuasiAssembling_ + cpuTimeSingleJacobianFactorization_ + cpuTimeSingleLinearSystemSolution_;

		out << std::endl;
		out << "Data for the TridiagonalBlocked NLS solver Kernel" << std::endl;
		out << "---------------------------------------------------------------------------------------------------------" << std::endl;
		out << "Number of system calls (only to assemble Jacobian):     " << numberOfSystemCallsForJacobian_ << " (" << numberOfSystemCallsPerJacobian_ << ")" << std::endl;
		out << "Number of full Jacobian constructions (from scratch):   " << numberOfJacobianFullAssembling_ << std::endl;
		out << "Number of approximated Jacobian constructions:          " << numberOfJacobianQuasiAssembling_ << std::endl;
		out << "Number of Jacobian factorizations:                      " << numberOfJacobianFactorizations_ << std::endl;
		out << "Number of linear system solutions:                      " << numberOfLinearSystemSolutions_ << std::endl;

		out << "Cumulative CPU for constructing the Jacobian (full):    " << cpuTimeJacobianFullAssembling_ << " (" << cpuTimeJacobianFullAssembling_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU for constructing the Jacobian (approx.): " << cpuTimeJacobianQuasiAssembling_ << " (" << cpuTimeJacobianQuasiAssembling_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU for factorizing the Jacobian:            " << cpuTimeJacobianFactorization_ << " (" << cpuTimeJacobianFactorization_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU for solving the linear system:           " << cpuTimeLinearSystemSolution_ << " (" << cpuTimeLinearSystemSolution_ / totalCpu *100. << "%)" << std::endl;

		out << "CPU for constructing the Jacobian (full):               " << cpuTimeSingleJacobianFullAssembling_ << " (" << cpuTimeSingleJacobianFullAssembling_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU for constructing the Jacobian (approx.):            " << cpuTimeSingleJacobianQuasiAssembling_ << " (" << cpuTimeSingleJacobianQuasiAssembling_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU for factorizing the Jacobian:                       " << cpuTimeSingleJacobianFactorization_ << " (" << cpuTimeSingleJacobianFactorization_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU for solving the linear system:                      " << cpuTimeSingleLinearSystemSolution_ << " (" << cpuTimeSingleLinearSystemSolution_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "---------------------------------------------------------------------------------------------------------" << std::endl;
		out << std::endl;
	}
}
