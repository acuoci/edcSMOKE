/***************************************************************************
 *   Copyright (C) 2023 by Alberto Cuoci                                   *
 *   alberto.cuoci@polimi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef OpenSMOKE_OpenSMOKETridiagonalBlockMatrix_Hpp
#define OpenSMOKE_OpenSMOKETridiagonalBlockMatrix_Hpp

#include "OpenSMOKEStdInclude.h"
#include "OpenSMOKEBaseClass.h"
#include "OpenSMOKEUtilities.h"
#include "OpenSMOKEFunctions.h"

namespace OpenSMOKE
{
	//!  A class for banded or tridiagonal-block (as a special case) matrices
	/*!
		 A class for banded or tridiagonal-block (as a special case) matrices
	*/

	template<typename T>
	class OpenSMOKETridiagonalBlockMatrix
	{
	public:

		/**
		*@brief Constructor for a tridiagonal-block matrix
		*@param dimBlock block dimension
		*/
		OpenSMOKETridiagonalBlockMatrix(const int nEquations, const int dimBlock);

		/**
		*@brief Sets all the coefficients equal to zero
		*/
		void SetToZero();

		/**
		*@brief Add the identity matrix
		*/
		void AddIdentity();

		/**
		*@brief Adds the specified diagonal matrix
		*@param d the diagonal vector to be added
		*/
		void AddDiagonal(const T *d);

		/**
		*@brief Destroys the matrix
		*/
		void DestroyMat();

		/**
		*@brief Copies the current matrix in the B matrix (sizes must be consistent)
		*@param B the matrix where to copy
		*/
		void CopyTo(OpenSMOKETridiagonalBlockMatrix<T>* B);

		/**
		*@brief Scales all the elements of the matrix by the same scalar
		*@param c scalar
		*/
		void Scale(const double c);

		/**
		*@brief Scales all the elements of the matrix by the two different scalars, according to the type of equation
		*@param c_differential scalar to be used by the differential equations
		*@param c_algebraic scalar to be used by the algebraic equations
		*@param index type of reaction (0=algebraic, 1=differential)
		*/
		void Scale(const double c_differential, const double c_algebraic, const int *index);	// TODO

		/**
		*@brief Initializes and returns the vectors defining the sparsity structure to be used for assembling the Jacobian matrix
		*/
		void GetSparsityData(std::vector<std::vector<int>>& group_cols, std::vector<std::vector<std::vector<int>>>& group_rows);

		/**
		*@brief Sets the central block (from 0 to N-1)
		*@param i block index
		*@param a vector of values corresonding to the block (column-based ordering)
		*/
		void SetCentralBlock(const int i, const T *a);

		/**
		*@brief Sets the upper block (from 0 to N-2)
		*@param i block index
		*@param a vector of values corresonding to the block (column-based ordering)
		*/
		void SetUpperBlock(const int i, const T *a);

		/**
		*@brief Sets the lower block (from 0 to N-2)
		*@param i block index
		*@param a vector of values corresonding to the block (column-based ordering)
		*/
		void SetLowerBlock(const int i, const T *a);

		/**
		*@brief Sets a specified column (used for whan the Jacobian matrix of a ODE/DAE system is assembled)
		*@param col target column
		*@param row starting row
		*@param n numer of rows to be assigned
		*@param a vector of values to be assigned (size NE)
		*/
		void SetColumn(const int col, const int row, const int n, const T *a);

		/**
		*@brief Returns the required element  
		*@param row row index (global, from 0)
		*@param col column index (global, from 0)
		*/
		T Get(const int row, const int col) const;

		/**
		*@brief Sets the required element  
		*@param row row index (global, from 0)
		*@param col column index (global, from 0)
		*@param a value to be set
		*/
		void Set(const int row, const int col, const T a);

		/**
		*@brief Multiplies the current matrix time a vector: y = A *x 
		*@param x vector to be multiplied
		*@param y vector where to put the result
		*/
		void Product(const T *x, T *y);

		/**
		*@brief Multiplies the transpose of current matrix time a vector: y = A' * x
		*@param x vector to be multiplied
		*@param y vector where to put the result
		*/
		void TProduct(const T *x, T *y);

		/**
		*@brief LU factorization
		*/
		int Factorize();

		/**
		*@brief Solves the linear system (only after LU factorization)
		*@param b rhs of linear system
		*/
		int Solve(T *b);

		/**
		*@brief Solves the linear system (only after LU factorization)
		*@param nrhs number of right hand sides
		*@param b rhs of linear system
		*/
		int Solve(const int nrhs, T *b);

		/**
		*@brief Factorizes and solves the linear system - Please, do not use it: it is still under testing
		*@param b rhs of linear system
		*/
		int FactorizeAndSolve(T *b);

		/**
		*@brief Prints the matrix on the screen (for diagnostic purposes)
		*@param out output stream
		*/
		void Print(std::ostream& out);

		/**
		*@brief Returns the block size
		*/
		int BlockSize() const { return NB; }

		/**
		*@brief Returns the number of blocks
		*/
		int BlocksNumber() const { return N; }

	public:

		T *D;		//!< central blocks
		T *C;		//!< upper blocks
		T *B;		//!< lower blocks
		int *diagonal_indices_;	
		int *rows_d_;	
		int *rows_b_;
		int *rows_c_;		
		int *cols_d_;	
		int *cols_b_;
		int *cols_c_;
		
	private:

		int N;		//!< number of blocks
		int NB;		//!< block size
		int NE;		//!< number of equations (N*NB)
		T *DU2;		//!< auxiliary matrix
		int *ipiv;	//!< pivot local row indices
		
				
		
	public:

		template<typename TT>
		friend OpenSMOKETridiagonalBlockMatrix<TT>* NewTridiagonalBlockMatrix(const int N, const int NB);
	};

	typedef OpenSMOKETridiagonalBlockMatrix<float>		OpenSMOKETridiagonalBlockMatrixFloat;
	typedef OpenSMOKETridiagonalBlockMatrix<double>		OpenSMOKETridiagonalBlockMatrixDouble;
}

#include "OpenSMOKETridiagonalBlockMatrix.hpp"

#endif	// OpenSMOKE_OpenSMOKETridiagonalBlockMatrix_Hpp

