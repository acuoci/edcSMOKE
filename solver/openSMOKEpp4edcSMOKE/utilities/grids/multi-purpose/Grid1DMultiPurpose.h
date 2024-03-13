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

#ifndef OpenSMOKE_Grid1DMultiPurpose_H
#define	OpenSMOKE_Grid1DMultiPurpose_H

#include "math/OpenSMOKEVector.h"
#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp> 
#include <boost/math/special_functions/pow.hpp> 

namespace OpenSMOKE
{
	void DimensionlessGrids(const unsigned int NPL, const unsigned int NPG, OpenSMOKE::OpenSMOKEVectorDouble& c_L, OpenSMOKE::OpenSMOKEVectorDouble& c_G);

	class Grid1DMultiPurpose
	{
	public:

		enum derivative_type { DERIVATIVE_FORWARD, DERIVATIVE_CENTERED, DERIVATIVE_BACKWARD, DERIVATIVE_UPWIND, DERIVATIVE_CENTERED_ACCURATE };
		enum geometry_type { GEOMETRY_1D_SPHERICAL, GEOMETRY_1D_CYLINDRICAL, GEOMETRY_1D_PLANAR };

	public:

		Grid1DMultiPurpose(const geometry_type geometry);

		unsigned int N() const { return N_; }

		OpenSMOKE::OpenSMOKEVectorDouble x;
		OpenSMOKE::OpenSMOKEVectorDouble dx;

	public:

		void Setup(const unsigned int number_of_points);
		void SetDimensionlessCoordinates(const OpenSMOKE::OpenSMOKEVectorDouble& c, const double cA, const double cB);

		void SetCoordinates(const OpenSMOKE::OpenSMOKEVectorDouble& c);
		void SetCoordinates(const Eigen::VectorXd& c);

		void Update(const OpenSMOKE::OpenSMOKEVectorDouble& new_coordinates);
		void Update(const Eigen::VectorXd& new_coordinates);
		void FirstDerivative(const OpenSMOKE::OpenSMOKEVectorDouble& u, const OpenSMOKE::OpenSMOKEVectorDouble& v, OpenSMOKE::OpenSMOKEVectorDouble& dv_over_dr, const derivative_type type);

		void SecondDerivative(const OpenSMOKE::OpenSMOKEVectorDouble& v, OpenSMOKE::OpenSMOKEVectorDouble& d2v_over_dr2);
		
		
		void DiffusionTerm(const OpenSMOKE::OpenSMOKEVectorDouble& v, const OpenSMOKE::OpenSMOKEVectorDouble& gamma_a, const OpenSMOKE::OpenSMOKEVectorDouble& gamma_b, OpenSMOKE::OpenSMOKEVectorDouble& diffusion_term);
		double IntegralValue(const OpenSMOKE::OpenSMOKEVectorDouble& v);
		double IntegralValue(const OpenSMOKE::OpenSMOKEVectorDouble& v1, const OpenSMOKE::OpenSMOKEVectorDouble& v2);
		void Interpolate(const OpenSMOKE::OpenSMOKEVectorDouble& vOriginal, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& vInterpolated);

		void FirstDerivative(const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& dv_over_dr, const derivative_type type);
		void FirstDerivativeAccurateBC(const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& dv_over_dr, const derivative_type type);

		void SecondDerivative(const Eigen::VectorXd& v, Eigen::VectorXd& d2v_over_dr2);

		double FirstDerivativeLeftSide(const Eigen::VectorXd& v);
		double FirstDerivativeRightSide(const Eigen::VectorXd& v);

		double FirstDerivativeSecondOrderLeftSide(const Eigen::VectorXd& v);
		double FirstDerivativeSecondOrderRightSide(const Eigen::VectorXd& v);
		double FirstDerivativeSecondOrder(const Eigen::VectorXd& v, const unsigned int i);

		void DiffusionTerm(const Eigen::VectorXd& v, const Eigen::VectorXd& gamma, Eigen::VectorXd& diffusion_term);
		void DiffusionTerm(const OpenSMOKE::OpenSMOKEVectorDouble& v, const OpenSMOKE::OpenSMOKEVectorDouble& gamma, OpenSMOKE::OpenSMOKEVectorDouble& diffusion_term);



	private:

		unsigned int N_;
		OpenSMOKE::OpenSMOKEVectorDouble a2;
		OpenSMOKE::OpenSMOKEVectorDouble a3;
		OpenSMOKE::OpenSMOKEVectorDouble a4;
		OpenSMOKE::OpenSMOKEVectorDouble a5;

		void FirstDerivativeFirstOrderBC(const Eigen::VectorXd& v, Eigen::VectorXd& dv_over_dr);
		void FirstDerivativeSecondOrderBC(const Eigen::VectorXd& v, Eigen::VectorXd& dv_over_dr);
		void FirstDerivativeInternalPoints(const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& dv_over_dr, const derivative_type type);

		void DiffusionTermSpherical(const Eigen::VectorXd& v, const Eigen::VectorXd& gamma, Eigen::VectorXd& diffusion_term);
		void DiffusionTermCylindrical(const Eigen::VectorXd& v, const Eigen::VectorXd& gamma, Eigen::VectorXd& diffusion_term);
		void DiffusionTermPlanar(const Eigen::VectorXd& v, const Eigen::VectorXd& gamma, Eigen::VectorXd& diffusion_term);
		void DiffusionTermSpherical(const OpenSMOKE::OpenSMOKEVectorDouble& v, const OpenSMOKE::OpenSMOKEVectorDouble& gamma, OpenSMOKE::OpenSMOKEVectorDouble& diffusion_term);
		void DiffusionTermCylindrical(const OpenSMOKE::OpenSMOKEVectorDouble& v, const OpenSMOKE::OpenSMOKEVectorDouble& gamma, OpenSMOKE::OpenSMOKEVectorDouble& diffusion_term);
		void DiffusionTermPlanar(const OpenSMOKE::OpenSMOKEVectorDouble& v, const OpenSMOKE::OpenSMOKEVectorDouble& gamma, OpenSMOKE::OpenSMOKEVectorDouble& diffusion_term);

		void ErrorMessage(const std::string message);

		static const double PI_;
		static const double PI_4_OVER_3_;

		geometry_type geometry_;
	};
}

#include "Grid1DMultiPurpose.hpp"

#endif // OpenSMOKE_Grid1DMultiPurpose_H




