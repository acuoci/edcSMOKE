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

#include <typeinfo>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iomanip>

#if OPENSMOKE_USE_MKL == 1
	#include "mkl.h"
	#include "mkl_lapacke.h"
	#if defined(_WIN32) || defined(_WIN64) 
	#else
		#include "mm_malloc.h"
	#endif
#elif OPENSMOKE_USE_OPENBLAS == 1
	#include "cblas.h"
	#include "lapacke.h"
	#ifndef __APPLE__
		#include "malloc.h"
	#endif
#endif

#if defined(_WIN32) || defined(_WIN64) 

	extern "C" {	void DGEBLTTRF(int* n, int* nb, double* d, double* b, double* c, double* du2, int* ipiv, int* info); }

	extern "C" {	void DGEBLTTRS(int* n, int* nb, int* NRHS, double* d, double* b, double* c, double* du2, int* ipiv, double* f, int* ldf, int* info); }

#else

	extern "C" {	void dgeblttrf_(int* n, int* nb, double* d, double* b, double* c, double* du2, int* ipiv, int* info); }

	extern "C" {	void dgeblttrs_(int* n, int* nb, int* NRHS, double* d, double* b, double* c, double* du2, int* ipiv, double* f, int* ldf, int* info); }

#endif


namespace OpenSMOKE
{
	#if defined(_WIN32) || defined(_WIN64) 
		// Nothing to declare
		template<typename TT>
		OpenSMOKETridiagonalBlockMatrix<TT>* NewTridiagonalBlockMatrix(const int N, const int NB);
	#else
		template<typename TT>
		OpenSMOKETridiagonalBlockMatrix<TT>* NewTridiagonalBlockMatrix(const int N, const int NB);
	#endif

	template<typename T>
	OpenSMOKETridiagonalBlockMatrix<T>::OpenSMOKETridiagonalBlockMatrix(const int nEquations, const int dimBlock)
	{
		NE = nEquations;
		NB = dimBlock;
		N = NE/NB;

		#if defined(_WIN32) || defined(_WIN64) 
		*this = *NewTridiagonalBlockMatrix<T>(N, NB);
		#else
		*this = *NewTridiagonalBlockMatrix<T>(N, NB);
		#endif
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::AddIdentity()
	{
		const T ONE = 1.;
		for (int i=0;i<NE;i++)
			D[diagonal_indices_[i]] += ONE;
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::AddDiagonal(const T *d)
	{
		for (int i=0;i<NE;i++)
			D[diagonal_indices_[i]] += d[i];
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::SetCentralBlock(const int i, const T *a)
	{
		int count = i*(NB*NB);
		for (int k=0;k<NB*NB;k++)
		{
			D[count] = a[k];
			count++;
		}
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::SetUpperBlock(const int i, const T *a)
	{
		int count = i*(NB*NB);
		for (int k=0;k<NB*NB;k++)
		{
			C[count] = a[k];
			count++;
		}
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::SetLowerBlock(const int i, const T *a)
	{
		int count = i*(NB*NB);
		for (int k=0;k<NB*NB;k++)
		{
			B[count] = a[k];
			count++;
		}
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::SetColumn(const int col, const int row, const int n, const T *a)
	{
		const int jBlock=std::floor(col/NB);
		const int jCol = col%NB;

		for (unsigned int k=row;k<row+n;k++)
		{
			const int iBlock = std::floor(k/NB);
			const int iRow = k%NB;
			
			if (iBlock == jBlock)		// Central
			{
				D[iBlock*NB*NB +jCol*NB + iRow] = a[k];
			}
			else if (iBlock == jBlock-1)	// Upper
			{
				C[iBlock*NB*NB +jCol*NB + iRow] = a[k];
			}
			else if (iBlock == jBlock+1)	// Lower
			{
				B[jBlock*NB*NB +jCol*NB + iRow] = a[k];
			}
		}
	}

	template<typename T>
	T OpenSMOKETridiagonalBlockMatrix<T>::Get(const int row, const int col) const
	{
		const int iBlock = std::floor(row/NB);
		const int jBlock = std::floor(col/NB);
		const int iRow = row%NB;
		const int jCol = col%NB;
			
		if (iBlock == jBlock)		// Central
		{
			return D[iBlock*NB*NB +jCol*NB + iRow];
		}
		else if (iBlock == jBlock-1)	// Upper
		{
			return C[iBlock*NB*NB +jCol*NB + iRow];
		}
		else if (iBlock == jBlock+1)	// Lower
		{
			return B[jBlock*NB*NB +jCol*NB + iRow];
		}
		else
		{
			std::cout << "Impossible to access the required element (tridiagonal-block matrix)" << std::endl;
			std::cout << "Row/Col (global):  " << row << "/" << col << std::endl;
			std::cout << "Row/Col (local):   " << iRow << "/" << jCol << std::endl;
			std::cout << "Matrix size/block: " << N*NB << "/" << NB << std::endl;
			std::cout << "Block indices:     " << iBlock << "/" << jBlock << std::endl;
			OpenSMOKE::FatalErrorMessage("T OpenSMOKETridiagonalBlockMatrix<T>::Get(const int row, const int col) const failure");
			
			return 0;
		}
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::Set(const int row, const int col, const T a)
	{
		const int iBlock = std::floor(row/NB);
		const int jBlock = std::floor(col/NB);
		const int iRow = row%NB;
		const int jCol = col%NB;
			
		if (iBlock == jBlock)		// Central
		{
			D[iBlock*NB*NB +jCol*NB + iRow] = a;
		}
		else if (iBlock == jBlock-1)	// Upper
		{
			C[iBlock*NB*NB +jCol*NB + iRow] = a;
		}
		else if (iBlock == jBlock+1)	// Lower
		{
			B[jBlock*NB*NB +jCol*NB + iRow] = a;
		}
		else
		{
			std::cout << "Impossible to access the required element (tridiagonal-block matrix)" << std::endl;
			std::cout << "Row/Col (global):  " << row << "/" << col << std::endl;
			std::cout << "Row/Col (local):   " << iRow << "/" << jCol << std::endl;
			std::cout << "Matrix size/block: " << N*NB << "/" << NB << std::endl;
			std::cout << "Block indices:     " << iBlock << "/" << jBlock << std::endl;
			OpenSMOKE::FatalErrorMessage("void OpenSMOKETridiagonalBlockMatrix<T>::Set(const int row, const int col, const T a) failure");
		}
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::Print(std::ostream& out)
	{
		for (int j = 0; j < NB*(N*NB); j++)
			out << rows_d_[j] << " " << cols_d_[j] << " " << D[j] << std::endl;
		
		for (int j = 0; j < NB*((N-1)*NB); j++)
			out << rows_b_[j] << " " << cols_b_[j] << " " << B[j] << std::endl;

		for (int j = 0; j < NB*((N-1)*NB); j++)
			out << rows_c_[j] << " " << cols_c_[j] << " " << C[j] << std::endl;
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::Scale(const double c_differential, const double c_algebraic, const int *index)
	{
		for (int j = 0; j < NB*(N*NB); j++)
		{
			if (index[rows_d_[j]] == 1)
				D[j] *= c_differential;
			else
				D[j] *= c_algebraic;
		}
		
		for (int j = 0; j < NB*((N-1)*NB); j++)
		{
			if (index[rows_b_[j]] == 1)
				B[j] *= c_differential;
			else
				B[j] *= c_algebraic;
		}

		for (int j = 0; j < NB*((N-1)*NB); j++)
		{
			if (index[rows_c_[j]] == 1)
				C[j] *= c_differential;
			else
				C[j] *= c_algebraic;
		}
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::SetToZero()
	{
		const T ZERO = 0.;

		for (int j = 0; j < NB*(N*NB); j++)
			D[j] = ZERO;
		for (int j = 0; j < NB*((N-1)*NB); j++)
			B[j] = ZERO;
		for (int j = 0; j < NB*((N-1)*NB); j++)
			C[j] = ZERO;

	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::DestroyMat()
	{
		#if (OPENSMOKE_USE_MKL == 1)
			_mm_free(D);
			D = NULL;
			_mm_free(B);
			B = NULL;
			_mm_free(C);
			C = NULL;
			_mm_free(DU2);
			DU2 = NULL;
			_mm_free(ipiv);
			ipiv = NULL;
			_mm_free(diagonal_indices_);
			diagonal_indices_ = NULL;
			_mm_free(rows_d_);
			rows_d_ = NULL;
			_mm_free(rows_b_);
			rows_b_ = NULL;
			_mm_free(rows_c_);
			rows_c_ = NULL;
			_mm_free(cols_d_);
			cols_d_ = NULL;
			_mm_free(cols_b_);
			cols_b_ = NULL;
			_mm_free(cols_c_);
			cols_c_ = NULL;
		#else
			free(D);
			D = NULL;
			free(B);
			B = NULL;
			free(C);
			C = NULL;
			free(DU2);
			DU2 = NULL;
			free(ipiv);
			ipiv = NULL;
			free(diagonal_indices_);
			diagonal_indices_ = NULL;
			free(rows_d_);
			rows_d_ = NULL;
			free(rows_b_);
			rows_b_ = NULL;
			free(rows_c_);
			rows_c_ = NULL;
			free(cols_d_);
			cols_d_ = NULL;
			free(cols_b_);
			cols_b_ = NULL;
			free(cols_c_);
			cols_c_ = NULL;
		#endif
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::CopyTo(OpenSMOKETridiagonalBlockMatrix<T>* M)
	{
			#if (OPENSMOKE_USE_MKL == 1)

				cblas_dcopy(NB*(N*NB), D, 1, M->D, 1);
				cblas_dcopy(NB*((N-1)*NB), B, 1, M->B, 1);
				cblas_dcopy(NB*((N-1)*NB), C, 1, M->C, 1);

			#elif (OPENSMOKE_USE_OPENBLAS == 1)

				cblas_dcopy(NB*(N*NB), D, 1, M->D, 1);
				cblas_dcopy(NB*((N-1)*NB), B, 1, M->B, 1);
				cblas_dcopy(NB*((N-1)*NB), C, 1, M->C, 1);

			#else

				for (int j = 0; j < NB*(N*NB); j++)
					M->D[j] = D[j];
				for (int j = 0; j < NB*((N-1)*NB); j++)
					M->B[j] = B[j];
				for (int j = 0; j < NB*((N-1)*NB); j++)
					M->C[j] = C[j];

			#endif
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::Scale(const double c)
	{
		#if __APPLE__

			cblas_dscal(NB*(N*NB), c, D, 1);
			cblas_dscal(NB*((N-1)*NB), c, B, 1);
			cblas_dscal(NB*((N-1)*NB), c, C, 1);

		#elif (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)

			const int one = 1;

			{
				const int lmat = NB*(N*NB);
				dscal(&lmat, &c, D, &one);
			}
			{
				const int lmat = NB*((N-1)*NB);
				dscal(&lmat, &c, B, &one);
			}
			{
				const int lmat = NB*((N-1)*NB);
				dscal(&lmat, &c, C, &one);
			}

		#else

			for (int j = 0; j < NB*(N*NB); j++)
				D[j] *= c;
			for (int j = 0; j < NB*((N-1)*NB); j++)
				B[j] *= c;
			for (int j = 0; j < NB*((N-1)*NB); j++)
				C[j] *= c;

		#endif
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::Product(const T *v_in, T *v_out)
	{
		const T ZERO = 0.;

		for (int i=0;i<NE;i++)
			v_out[i] = ZERO;

		for (int i=0;i<NB*(N*NB);i++)
			v_out[rows_d_[i]] += D[i]*v_in[cols_d_[i]];
		
		for (int i=0;i<NB*((N-1)*NB);i++)
			v_out[rows_b_[i]] += B[i]*v_in[cols_b_[i]];
		
		for (int i=0;i<NB*((N-1)*NB);i++)
			v_out[rows_c_[i]] += C[i]*v_in[cols_c_[i]];
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::TProduct(const T *v_in, T *v_out)
	{
		const T ZERO = 0.;

		for (int i=0;i<NE;i++)
			v_out[i] = ZERO;

		for (int i=0;i<NB*(N*NB);i++)
			v_out[cols_d_[i]] += D[i]*v_in[rows_d_[i]];
		
		for (int i=0;i<NB*((N-1)*NB);i++)
			v_out[cols_b_[i]] += B[i]*v_in[rows_b_[i]];
		
		for (int i=0;i<NB*((N-1)*NB);i++)
			v_out[cols_c_[i]] += C[i]*v_in[rows_c_[i]];
	}

	template<typename T>
	void OpenSMOKETridiagonalBlockMatrix<T>::GetSparsityData(std::vector<std::vector<int>>& group_cols, std::vector<std::vector<std::vector<int>>>& group_rows)
	{
		const int width = NB*3;
		int ngroups = std::min(width, N*NB);

		// Recognize columns
		group_cols.resize(ngroups);
		for (unsigned int i=0;i<ngroups;i++)
		{
			group_cols[i].resize(0);
			const int n = ((N*NB)%width == 0) ? static_cast<int>(N*NB/width) : static_cast<int>(N*NB/width) + 1;
			for(int j=0;j<n;j++)
			{
				const int index = j*width + i;
				if (index < N*NB)
					group_cols[i].push_back(index);
			}
		}

		// Checking
		{
			std::vector<bool> cols(N*NB);
			std::fill(cols.begin(), cols.end(), false);
			int nCols = 0;
			for (unsigned int i=0;i<ngroups;i++)
				for (unsigned int j=0;j<group_cols[i].size();j++)
				{
					cols[group_cols[i][j]] = true;
					nCols++;
				}
		
			if (nCols != N*NB)
				OpenSMOKE::FatalErrorMessage("void OpenSMOKETridiagonalBlockMatrix<T>::GetSparsityData(...) failure");

			if (std::count(cols.begin(), cols.end(), false))
				OpenSMOKE::FatalErrorMessage("void OpenSMOKETridiagonalBlockMatrix<T>::GetSparsityData(...) failure");
		}

		group_rows.resize(ngroups);
		for (unsigned int i=0;i<ngroups;i++)
		{
			group_rows[i].resize(group_cols[i].size());
			for (unsigned int j=0;j<group_cols[i].size();j++)
			{
				const int k = group_cols[i][j];

				if (k<NB)
				{
					for (unsigned int m=0;m<2*NB;m++)
						group_rows[i][j].push_back(m);
				}
				else if (k>=(N-1)*NB)
				{
					for (unsigned int m=0;m<2*NB;m++)
						group_rows[i][j].push_back((N-2)*NB+m);
				}
				else
				{
					const int block = std::floor(k/NB);
					for (unsigned int m=0;m<3*NB;m++)
						group_rows[i][j].push_back((block-1)*NB+m);
				}
			}
		}
	}

	template<typename T>
	int OpenSMOKETridiagonalBlockMatrix<T>::Factorize()
	{
		int info;

		#if defined(_WIN32) || defined(_WIN64)
			DGEBLTTRF(&N, &NB, D, B, C, DU2, ipiv, &info);
		#else
			dgeblttrf_(&N, &NB, D, B, C, DU2, ipiv, &info);
		#endif

		if (info != 0)
		{
			if (info == 0)
			{
				// Successful factorization
			}
			else if (info == -1000)
			{
			}
			else if (info > 0)
			{
				std::cout << "Pivot (from 1) " << info << " is equal to zero (singular matrix)." << std::endl;
				std::cout << "The factorization can be not completed." << std::endl;
				std::cout << " (1) Matrix layout:        " << CblasColMajor << std::endl;
				std::cout << " (2) Equations:            " << NE << std::endl;
				std::cout << " (3) Blocks:               " << N << std::endl;
				std::cout << " (4) Block size:           " << NB << std::endl;
				std::cout << " (5) Block index (from 1): " << int(std::floor((info-1)/NB))+1 << std::endl;
				std::cout << " (6) Local row (from 1):   " << (info-1)%NB+1 << std::endl;
			}
			else if (info < 0) 
			{
				std::cout << "Parameter " << -info << " has an illegal value" << std::endl;

				std::cout << " (1) Matrix layout: " << CblasColMajor << std::endl;
				std::cout << " (2) Equations:     " << NE << std::endl;
				std::cout << " (3) Blocks:        " << N << std::endl;
				std::cout << " (4) Block size:    " << NB << std::endl;

				OpenSMOKE::FatalErrorMessage("Factorizing the tridiagonal block linear system: DGEBLTTRF failed");
			}
		}

		return info;
	}

	template<typename T>
	int OpenSMOKETridiagonalBlockMatrix<T>::Solve(T *b)
	{
		int one = 1;
		int ldf = N*NB;
		int info;

		#if defined(_WIN32) || defined(_WIN64)
			DGEBLTTRS(&N, &NB, &one, D, B, C, DU2, ipiv, b, &ldf, &info);
		#else
			dgeblttrs_(&N, &NB, &one, D, B, C, DU2, ipiv, b, &ldf, &info);
		#endif

		if (info != 0)
		{
			if (info == 0)
			{
				// Successful factorization
			}
			else if (info < 0) 
			{
				std::cout << "Parameter " << -info << " has an illegal value" << std::endl;

				std::cout << " (1) Matrix layout: " << CblasColMajor << std::endl;
				std::cout << " (2) Equations:     " << NE << std::endl;
				std::cout << " (3) Blocks:        " << N << std::endl;
				std::cout << " (4) Block size:    " << NB << std::endl;

				OpenSMOKE::FatalErrorMessage("Solving the tridiagonal block linear system: DGEBLTTRS failed");
			}
		}

		return info;
	}

	template<typename T>
	int OpenSMOKETridiagonalBlockMatrix<T>::Solve(const int nrhs, T *b)
	{
		int ldf = N*NB;
		int info;

		#if defined(_WIN32) || defined(_WIN64)
			DGEBLTTRS(&N, &NB, &nrhs, D, B, C, DU2, ipiv, b, &ldf, &info);
		#else
			dgeblttrs_(&N, &NB, &nrhs, D, B, C, DU2, ipiv, b, &ldf, &info);
		#endif

		if (info != 0)
		{
			if (info == 0)
			{
				// Successful factorization
			}
			else if (info < 0) 
			{
				std::cout << "Parameter " << -info << " has an illegal value" << std::endl;

				std::cout << " (1) Matrix layout: " << CblasColMajor << std::endl;
				std::cout << " (2) Equations:     " << NE << std::endl;
				std::cout << " (3) Blocks:        " << N << std::endl;
				std::cout << " (4) Block size:    " << NB << std::endl;

				OpenSMOKE::FatalErrorMessage("Solving the tridiagonal block linear system: DGEBLTTRS failed");
			}
		}

		return info;
	}

	template<typename T>
	OpenSMOKETridiagonalBlockMatrix<T>* NewTridiagonalBlockMatrix(const int N, const int NB)
	{
		OpenSMOKETridiagonalBlockMatrix<T>* A;
		
		#if (OPENSMOKE_USE_MKL == 1)
			A = NULL;
			A = (OpenSMOKETridiagonalBlockMatrix<T>*)_mm_malloc(sizeof *A, 64);
			if (A == NULL) return (NULL);
		#else
			A = NULL;
			A = (OpenSMOKETridiagonalBlockMatrix<T>*)malloc(sizeof *A);
			if (A == NULL) return (NULL);
		#endif
		
		#if (OPENSMOKE_USE_MKL == 1)
	
			A->D = NULL;
			A->D = (T *)_mm_malloc(NB*(N*NB) * sizeof(T),64);
			if (A->D == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL); 
			}
	
			A->B = NULL;
			A->B = (T *)_mm_malloc(NB*((N-1)*NB) * sizeof(T),64);
			if (A->B == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL); 
			}
	
			A->C = NULL;
			A->C = (T *)_mm_malloc(NB*((N-1)*NB) * sizeof(T),64);
			if (A->C == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL); 
			}
	
			A->DU2 = NULL;
			A->DU2 = (T *)_mm_malloc(NB*((N-2)*NB) * sizeof(T),64);
			if (A->DU2 == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL); 
			}
	
			A->ipiv = NULL;
			A->ipiv = (int *)_mm_malloc(NB*N * sizeof(int),64);
			if (A->ipiv == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL);
			}

			A->diagonal_indices_ = NULL;
			A->diagonal_indices_ = (int *)_mm_malloc(NB*N * sizeof(int),64);
			if (A->diagonal_indices_ == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL);
			}

			A->rows_d_ = NULL;
			A->rows_d_ = (int *)_mm_malloc(NB*(N*NB) * sizeof(int),64);
			if (A->rows_d_ == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL);
			}

			A->rows_b_ = NULL;
			A->rows_b_ = (int *)_mm_malloc(NB*((N-1)*NB) * sizeof(int),64);
			if (A->rows_b_ == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL);
			}

			A->rows_c_ = NULL;
			A->rows_c_ = (int *)_mm_malloc(NB*((N-1)*NB) * sizeof(int),64);
			if (A->rows_c_ == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL);
			}

			A->cols_d_ = NULL;
			A->cols_d_ = (int *)_mm_malloc(NB*(N*NB) * sizeof(int),64);
			if (A->cols_d_ == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL);
			}

			A->cols_b_ = NULL;
			A->cols_b_ = (int *)_mm_malloc(NB*((N-1)*NB) * sizeof(int),64);
			if (A->cols_b_ == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL);
			}

			A->cols_c_ = NULL;
			A->cols_c_ = (int *)_mm_malloc(NB*((N-1)*NB) * sizeof(int),64);
			if (A->cols_c_ == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL);
			}
		#else

			A->D = NULL;
			A->D = (T *)malloc(NB*(N*NB) * sizeof(T));
			if (A->D == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

			A->B = NULL;
			A->B = (T *)malloc(NB*((N-1)*NB) * sizeof(T));
			if (A->B == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

			A->C = NULL;
			A->C = (T *)malloc(NB*((N-1)*NB) * sizeof(T));
			if (A->C == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

			A->DU2 = NULL;
			A->DU2 = (T *)malloc(NB*((N-2)*NB) * sizeof(T));
			if (A->DU2 == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

			A->ipiv = NULL;
			A->ipiv = (int *)malloc(NB*N * sizeof(int));
			if (A->ipiv == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

			A->diagonal_indices_ = NULL;
			A->diagonal_indices_ = (int *)malloc(NB*N * sizeof(int));
			if (A->diagonal_indices_ == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

			A->rows_d_ = NULL;
			A->rows_d_ = (int *)malloc(NB*(N*NB) * sizeof(int));
			if (A->rows_d_ == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

			A->rows_b_ = NULL;
			A->rows_b_ = (int *)malloc(NB*((N-1)*NB) * sizeof(int));
			if (A->rows_b_ == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

			A->rows_c_ = NULL;
			A->rows_c_ = (int *)malloc(NB*((N-1)*NB) * sizeof(int));
			if (A->rows_c_ == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

			A->cols_d_ = NULL;
			A->cols_d_ = (int *)malloc(NB*(N*NB) * sizeof(int));
			if (A->cols_d_ == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

			A->cols_b_ = NULL;
			A->cols_b_ = (int *)malloc(NB*((N-1)*NB) * sizeof(int));
			if (A->cols_b_ == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

			A->cols_c_ = NULL;
			A->cols_c_ = (int *)malloc(NB*((N-1)*NB) * sizeof(int));
			if (A->cols_c_ == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}

		#endif
	
		// Matrix properties
		{
			// Number of blocks
			A->N = N;
		
			// Block size
			A->NB = NB;

			// Number of equations
			A->NE = N*NB;

			// Indices of diagonal elements (to be applied to D)
			{
				int count = 0;
				for (int j=0;j<N;j++)
				{
					const int i = j*(NB*NB);
					for (int k=0;k<NB;k++)
					{
						const int index = i + k*(NB+1);
						A->diagonal_indices_[count] = index;
						count++;
					}
				}
			}

			// Rows of elements on D
			{
				int count = 0;
				for (int j=0;j<N;j++)
				{
					int col = j*NB;
					for (int k=0;k<NB;k++)
					{
						int row = j*NB;
						for (int i=0;i<NB;i++)
						{	
							A->rows_d_[count] = row;
							A->cols_d_[count] = col;
							row++;
							count++;
						}
						col++;
					}
				}
			}

			// Rows of elements on B
			{
				int count = 0;
				for (int j=1;j<N;j++)
				{
					int col = (j-1)*NB;
					for (int k=0;k<NB;k++)
					{
						int row = j*NB;
						for (int i=0;i<NB;i++)
						{	
							A->rows_b_[count] = row;
							A->cols_b_[count] = col;
							row++;
							count++;
						}
						col++;
					}
				}
			}

			// Rows of elements on C
			{
				int count = 0;
				for (int j=0;j<N-1;j++)
				{
					int col = (j+1)*NB;
					for (int k=0;k<NB;k++)
					{
						int row = j*NB;
						for (int i=0;i<NB;i++)
						{	
							A->rows_c_[count] = row;
							A->cols_c_[count] = col;
							row++;
							count++;
						}
						col++;
					}
				}
			}
		}	

		return(A);
	}
}
