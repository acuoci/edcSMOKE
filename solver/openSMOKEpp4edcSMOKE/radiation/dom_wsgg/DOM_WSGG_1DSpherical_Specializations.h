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
|	License                                                           |
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

#ifndef OpenSMOKE_DOM_WSGG_1DSpherical_Specializations_H
#define	OpenSMOKE_DOM_WSGG_1DSpherical_Specializations_H

#include "DOM_WSGG_1DSpherical.h"

namespace OpenSMOKE
{
	//!  A class to model radiative heat transfer in 1D geometries according to the S2 (symmetric) Discrete Ordinates Method
	/*!
	Cosines and weights are taken from: Modest, Radiative Heat Transfer, 3rd Edition, Chapter 17
	*/

	class DOM_WSGG_1DSpherical_S2sym : public DOM_WSGG_1DSpherical
	{

	public:

		DOM_WSGG_1DSpherical_S2sym(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model) :
			DOM_WSGG_1DSpherical(grid, model)
		{
			n_ = 2;
			N_ = 1;

			mu_ = new double[n_];
			mu_[0] = -1. / std::sqrt(3.);

			wp_ = new double[n_];
			wp_[0] = 2.*boost::math::constants::pi<double>();

			Finalize();
		}
	};

	//!  A class to model radiative heat transfer in 1D geometries according to the S2 (non-symmetric) Discrete Ordinates Method
	/*!
	Cosines and weights are taken from: Modest, Radiative Heat Transfer, 3rd Edition, Chapter 17
	*/

	class DOM_WSGG_1DSpherical_S2nonsym : public DOM_WSGG_1DSpherical
	{

	public:

		DOM_WSGG_1DSpherical_S2nonsym(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model) :
			DOM_WSGG_1DSpherical(grid, model)
		{
			n_ = 2;
			N_ = 1;

			mu_ = new double[n_];
			mu_[0] = -1./2.;

			wp_ = new double[n_];
			wp_[0] = 2.*boost::math::constants::pi<double>();

			Finalize();
		}
	};

	//!  A class to model radiative heat transfer in 1D geometries according to the S4 Discrete Ordinates Method
	/*!
	Cosines and weights are taken from: Modest, Radiative Heat Transfer, 3rd Edition, Chapter 17
	*/

	class DOM_WSGG_1DSpherical_S4 : public DOM_WSGG_1DSpherical
	{

	public:

		DOM_WSGG_1DSpherical_S4(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model) :
			DOM_WSGG_1DSpherical(grid, model)
		{
			n_ = 4;
			N_ = 2;

			mu_ = new double[n_];
			mu_[0] = -0.2958759;
			mu_[1] = -0.9082483;

			wp_ = new double[n_];
			wp_[0] = 4. / 3. * boost::math::constants::pi<double>();
			wp_[1] = 2. / 3. * boost::math::constants::pi<double>();					

			Finalize();
		}
	};

	//!  A class to model radiative heat transfer in 1D geometries according to the S6 Discrete Ordinates Method
	/*!
	Cosines and weights are taken from: Modest, Radiative Heat Transfer, 3rd Edition, Chapter 17
	*/

	class DOM_WSGG_1DSpherical_S6 : public DOM_WSGG_1DSpherical
	{

	public:

		DOM_WSGG_1DSpherical_S6(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model) :
			DOM_WSGG_1DSpherical(grid, model)
		{
			n_ = 6;
			N_ = 3;

			mu_ = new double[n_];
			mu_[0] = -0.1838670;
			mu_[1] = -0.6950514;
			mu_[2] = -0.9656013;

			wp_ = new double[n_];
			wp_[0] = 2.7382012;
			wp_[1] = 2.9011752;
			wp_[2] = 0.6438068;

			Finalize();
		}
	};

	//!  A class to model radiative heat transfer in 1D geometries according to the S8 Discrete Ordinates Method
	/*!
	Cosines and weights are taken from: Modest, Radiative Heat Transfer, 3rd Edition, Chapter 17
	*/

	class DOM_WSGG_1DSpherical_S8 : public DOM_WSGG_1DSpherical
	{

	public:

		DOM_WSGG_1DSpherical_S8(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model) :
			DOM_WSGG_1DSpherical(grid, model)
		{
			n_ = 8;
			N_ = 4;

			mu_ = new double[n_];
			mu_[0] = -0.14225553242346792;
			mu_[1] = -0.57735026918962018;
			mu_[2] = -0.80400872517751543;
			mu_[3] = -0.97955435121784884;

			wp_ = new double[n_];
			wp_[0] = 3.4436590771700393e-1;
			wp_[1] = 4.2028034773474127e-1;
			wp_[2] = 1.2634158137950924e-1;
			wp_[3] = 1.0901216316874843e-1;

			for (unsigned int i = 0; i < N_; i++)
				wp_[i] *= 2.*PI_;

			Finalize();
		}
	};

	//!  A class to model radiative heat transfer in 1D geometries according to the S10 Discrete Ordinates Method
	/*!
	Cosines and weights are taken from: ???
	*/

	class DOM_WSGG_1DSpherical_S10 : public DOM_WSGG_1DSpherical
	{

	public:

		DOM_WSGG_1DSpherical_S10(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model) :
			DOM_WSGG_1DSpherical(grid, model)
		{
			n_ = 10;
			N_ = 5;

			mu_ = new double[n_];
			mu_[0] = -1.3727193312794930e-1;
			mu_[1] = -5.0468891002891736e-1;
			mu_[2] = -7.0041288408170388e-1;
			mu_[3] = -8.5231773445654380e-1;
			mu_[4] = -9.8097544961666649e-1;

			wp_ = new double[n_];
			wp_[0] = 3.2024697371374683e-1;
			wp_[1] = 3.3536146934149275e-1;
			wp_[2] = 9.5325916625663462e-2;
			wp_[3] = 1.8894253053589230e-1;
			wp_[4] = 6.0123109783207715e-2;

			for (unsigned int i = 0; i < N_; i++)
				wp_[i] *= PI_TIMES_2_;

			Finalize();
		}
	};

	//!  A class to model radiative heat transfer in 1D geometries according to the S12 Discrete Ordinates Method
	/*!
	Cosines and weights are taken from: ???
	*/

	class DOM_WSGG_1DSpherical_S12 : public DOM_WSGG_1DSpherical
	{

	public:

		DOM_WSGG_1DSpherical_S12(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model) :
			DOM_WSGG_1DSpherical(grid, model)
		{
			n_ = 12;
			N_ = 6;

			mu_ = new double[n_];
			mu_[0] = -1.2028966644554817e-1;
			mu_[1] = -4.5363844804142206e-1;
			mu_[2] = -6.3016353371905831e-1;
			mu_[3] = -7.6708820673840594e-1;
			mu_[4] = -8.8303032485016963e-1;
			mu_[5] = -9.8542416871763827e-1;

			wp_ = new double[n_];
			wp_[0] = 2.8269030105021986e-1;
			wp_[1] = 3.1219784576420512e-1;
			wp_[2] = 8.7042725966186285e-2;
			wp_[3] = 1.4291214605783314e-1;
			wp_[4] = 1.2413797501117962e-1;
			wp_[5] = 5.1019006150373108e-2;

			for (unsigned int i = 0; i < N_; i++)
				wp_[i] *= PI_TIMES_2_;

			Finalize();
		}
	};

	//!  A class to model radiative heat transfer in 1D geometries according to the S2 (symmetric) Discrete Ordinates Method
	/*!
	Cosines and weights are generated using the provided utility
	Constraints on moments: 0, 1, 2, 3, 6, 10, 12, 14, 16
	*/

	class DOM_WSGG_1DSpherical_S16 : public DOM_WSGG_1DSpherical
	{

	public:

		DOM_WSGG_1DSpherical_S16(Grid1DMultiPurpose* grid, const OpenSMOKE::ExtinctionCoefficientModel model) :
			DOM_WSGG_1DSpherical(grid, model)
		{
			n_ = 16;
			N_ = 8;

			mu_ = new double[n_];
			mu_[0] = -9.275094395832e-002;
			mu_[1] = -3.844125296049e-001;
			mu_[2] = -5.356708394428e-001;
			mu_[3] = -6.527737000993e-001;
			mu_[4] = -7.518535488069e-001;
			mu_[5] = -8.393175883809e-001;
			mu_[6] = -9.184902119661e-001;
			mu_[7] = -9.913599370510e-001;

			wp_ = new double[n_];
			wp_[0] = 2.239624623267e-001;
			wp_[1] = 3.134916476929e-001;
			wp_[2] = 3.396323130917e-003;
			wp_[3] = 1.863118255452e-001;
			wp_[4] = 7.710156584389e-002;
			wp_[5] = 6.035214457415e-002;
			wp_[6] = 1.037416397381e-001;
			wp_[7] = 3.164239114814e-002;

			for (unsigned int i = 0; i < N_; i++)
				wp_[i] *= PI_TIMES_2_;

			Finalize();
		}
	};
	
}

#endif
