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

namespace OpenSMOKE
{
	const double Grid1DMultiPurpose::PI_ = boost::math::constants::pi<double>();
	const double Grid1DMultiPurpose::PI_4_OVER_3_ = 4. / 3.*PI_;

	void Grid1DMultiPurpose::ErrorMessage(const std::string message)
	{
		std::cout << "Grid: " << message << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(0);
	}

	Grid1DMultiPurpose::Grid1DMultiPurpose(const geometry_type geometry)
	{
		N_ = 0;
		geometry_ = geometry;
	}

	void Grid1DMultiPurpose::Setup(const unsigned int number_of_points)
	{
		N_ = number_of_points;
		ChangeDimensions(N_, &x, true);
		ChangeDimensions(N_ - 1, &dx, true);
		ChangeDimensions(N_ - 1, &a2, true);
		ChangeDimensions(N_ - 1, &a3, true);
		ChangeDimensions(N_ - 1, &a4, true);
		ChangeDimensions(N_ - 1, &a5, true);
	}

	void Grid1DMultiPurpose::SetDimensionlessCoordinates(const OpenSMOKE::OpenSMOKEVectorDouble& c, const double xA, const double xB)
	{
		if (c.Size() != N_)
                  {
                    std::cout << c.Size() << "   " << N_ << std::endl;
			ErrorMessage("The number of points of the grid was changed.");
                  }

		OpenSMOKE::OpenSMOKEVectorDouble new_coordinates(N_);
		for (unsigned int i = 1; i <= N_; i++)
			new_coordinates[i] = xA + c[i] * (xB - xA);
		Update(new_coordinates);
	}

	void Grid1DMultiPurpose::SetCoordinates(const OpenSMOKE::OpenSMOKEVectorDouble& c)
	{
		if (c.Size() != N_)
			ErrorMessage("The number of points of the grid was changed.");
		Update(c);
	}

	void Grid1DMultiPurpose::SetCoordinates(const Eigen::VectorXd& c)
	{
		if (c.size() != N_)
			ErrorMessage("The number of points of the grid was changed.");
		Update(c);
	}

	void Grid1DMultiPurpose::Update(const Eigen::VectorXd& new_coordinates)
	{
		OpenSMOKE::OpenSMOKEVectorDouble xx(N_);
		for (unsigned int i = 0; i < N_; i++)
			xx[i+1] = new_coordinates(i);
		Update(xx);
	}

	void Grid1DMultiPurpose::Update(const OpenSMOKE::OpenSMOKEVectorDouble& new_coordinates)
	{
		x = new_coordinates;
		for (unsigned int i = 1; i <= N_ - 1; i++)
			dx[i] = x[i + 1] - x[i];

		for (unsigned int i = 2; i <= N_ - 1; i++)
		{
			double alpha = dx[i] / dx[i - 1];
			double ualpha = 1. / alpha;
			double u1ualpha = 1. / (1. + ualpha);
			a3[i] = ualpha *  u1ualpha;
			a4[i] = 1. *  u1ualpha;
			a5[i] = alpha *  u1ualpha;
			a2[i] = a4[i] + a5[i];
		}
	}

	void Grid1DMultiPurpose::FirstDerivative(const OpenSMOKE::OpenSMOKEVectorDouble& u, const OpenSMOKE::OpenSMOKEVectorDouble& v, OpenSMOKE::OpenSMOKEVectorDouble& dv_over_dr, const derivative_type type)
	{
		dv_over_dr[1] = (v[2] - v[1]) / dx[1];
		dv_over_dr[N_] = (v[N_] - v[N_ - 1]) / dx[N_ - 1];

		if (type == DERIVATIVE_UPWIND)
		{
			for (unsigned int i = 2; i <= N_ - 1; i++)
				dv_over_dr[i] = (u[i] >= 0.) ? (v[i] - v[i - 1]) / dx[i - 1] : (v[i + 1] - v[i]) / dx[i];
		}
		else if (type == DERIVATIVE_BACKWARD)
		{
			for (unsigned int i = 2; i <= N_ - 1; i++)
				dv_over_dr[i] = (v[i] - v[i - 1]) / dx[i - 1];
		}
		else if (type == DERIVATIVE_FORWARD)
		{
			for (unsigned int i = 2; i <= N_ - 1; i++)
				dv_over_dr[i] = (v[i + 1] - v[i]) / dx[i];
		}
		else if (type == DERIVATIVE_CENTERED)
		{
			for (unsigned int i = 2; i <= N_ - 1; i++)
				dv_over_dr[i] = (v[i + 1] - v[i-1]) / (dx[i - 1] + dx[i]);
		}
		else if (type == DERIVATIVE_CENTERED_ACCURATE)
		{
			for (unsigned int i = 2; i <= N_ - 1; i++)
				dv_over_dr[i] = (a3[i] * (v[i + 1] - v[i]) - a5[i] * (v[i - 1] - v[i])) / dx[i];
		}
	}

	void Grid1DMultiPurpose::FirstDerivativeFirstOrderBC(const Eigen::VectorXd& v, Eigen::VectorXd& dv_over_dr)
	{
		dv_over_dr(0) = (v(1) - v(0)) / dx[1];
		dv_over_dr(N_ - 1) = (v(N_ - 1) - v(N_ - 2)) / dx[N_ - 1];
	}

	void Grid1DMultiPurpose::FirstDerivativeSecondOrderBC(const Eigen::VectorXd& v, Eigen::VectorXd& dv_over_dr)
	{
		dv_over_dr(0) = FirstDerivativeSecondOrderLeftSide(v);
		dv_over_dr(N_ - 1) = FirstDerivativeSecondOrderRightSide(v);
	}

	double Grid1DMultiPurpose::FirstDerivativeLeftSide(const Eigen::VectorXd& v)
	{
		return  (v(1) - v(0)) / dx[1];
	}

	double Grid1DMultiPurpose::FirstDerivativeRightSide(const Eigen::VectorXd& v)
	{
		return (v(N_ - 1) - v(N_ - 2)) / dx[N_ - 1];
	}

	double Grid1DMultiPurpose::FirstDerivativeSecondOrderLeftSide(const Eigen::VectorXd& v)
	{
		const double a = dx[1];
		const double b = dx[2];
		const double ab = a + b;
		return (-v(2)*a*a + v(1)*ab*ab - v(0)*(ab*ab - a*a)) / (a*ab*b);
	}

	double Grid1DMultiPurpose::FirstDerivativeSecondOrderRightSide(const Eigen::VectorXd& v)
	{
		const double a = dx[N_-1];
		const double b = dx[N_-2];
		const double ab = a + b;
		return ( v(N_-3)*a*a -v(N_-2)*ab*ab + v(N_-1)*(ab*ab - a*a)) / (a*ab*b);
	}

	double Grid1DMultiPurpose::FirstDerivativeSecondOrder(const Eigen::VectorXd& v, const unsigned int i)
	{
		return (v(i + 1) - v(i - 1)) / (dx[i] + dx[i + 1]);
	}

	void Grid1DMultiPurpose::FirstDerivative(const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& dv_over_dr, const derivative_type type)
	{
		FirstDerivativeFirstOrderBC(v, dv_over_dr);
		FirstDerivativeInternalPoints(u, v, dv_over_dr, type);
	}

	void Grid1DMultiPurpose::FirstDerivativeAccurateBC(const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& dv_over_dr, const derivative_type type)
	{
		FirstDerivativeSecondOrderBC(v, dv_over_dr);
		FirstDerivativeInternalPoints(u, v, dv_over_dr, type);
	}

	void Grid1DMultiPurpose::FirstDerivativeInternalPoints(const Eigen::VectorXd& u, const Eigen::VectorXd& v, Eigen::VectorXd& dv_over_dr, const derivative_type type)
	{
		if (type == DERIVATIVE_UPWIND)
		{
			for (unsigned int i = 1; i <= N_ - 2; i++)
				dv_over_dr(i) = (u(i) >= 0.) ? (v(i) - v(i - 1)) / dx[i] : (v(i + 1) - v(i)) / dx[i+1];
		}
		else if (type == DERIVATIVE_BACKWARD)
		{
			for (unsigned int i = 1; i <= N_ - 2; i++)
				dv_over_dr(i) = (v(i) - v(i - 1)) / dx[i];
		}
		else if (type == DERIVATIVE_FORWARD)
		{
			for (unsigned int i = 1; i <= N_ - 2; i++)
				dv_over_dr(i) = (v(i + 1) - v(i)) / dx[i + 1];
		}
		else if (type == DERIVATIVE_CENTERED)
		{
			for (unsigned int i = 1; i <= N_ - 2; i++)
				dv_over_dr(i) = (v(i + 1) - v(i-1)) / (dx[i] + dx[i + 1]);
		}
		else if (type == DERIVATIVE_CENTERED_ACCURATE)
		{
			for (unsigned int i = 1; i <= N_ - 2; i++)
				dv_over_dr(i) = (a3[i + 1] * (v(i + 1) - v(i)) - a5[i + 1] * (v(i - 1) - v(i))) / dx[i + 1];
		}
	}

	
	void Grid1DMultiPurpose::SecondDerivative(const OpenSMOKE::OpenSMOKEVectorDouble& v, OpenSMOKE::OpenSMOKEVectorDouble& d2v_over_dr2)
	{
		d2v_over_dr2[1] = 0.;
		d2v_over_dr2[N_] = 0.;

		for (unsigned int i = 2; i <= N_ - 1; i++)
			d2v_over_dr2[i] = 2.*(a4[i] * (v[i + 1] - v[i]) + a5[i] * (v[i - 1] - v[i])) / (dx[i] * dx[i]);
	}

	void Grid1DMultiPurpose::SecondDerivative(const Eigen::VectorXd& v, Eigen::VectorXd& d2v_over_dr2)
	{
		d2v_over_dr2(0) = 0.;
		d2v_over_dr2(N_-1) = 0.;

		for (unsigned int i = 1; i <= N_ - 2; i++)
			d2v_over_dr2(i) = 2.*(a4[i+1] * (v(i + 1) - v(i)) + a5[i+1] * (v(i - 1) - v(i))) / (dx[i+1] * dx[i+1]);
	}

	void Grid1DMultiPurpose::DiffusionTerm(const OpenSMOKE::OpenSMOKEVectorDouble& v, const OpenSMOKE::OpenSMOKEVectorDouble& gamma, OpenSMOKE::OpenSMOKEVectorDouble& diffusion_term)
	{
		if (geometry_ == GEOMETRY_1D_SPHERICAL)
			DiffusionTermSpherical(v, gamma, diffusion_term);
		else if (geometry_ == GEOMETRY_1D_CYLINDRICAL)
			DiffusionTermCylindrical(v, gamma, diffusion_term);
		else if (geometry_ == GEOMETRY_1D_PLANAR)
			DiffusionTermPlanar(v, gamma, diffusion_term);
	}
	void Grid1DMultiPurpose::DiffusionTerm(const Eigen::VectorXd& v, const Eigen::VectorXd& gamma, Eigen::VectorXd& diffusion_term)
	{
		if (geometry_ == GEOMETRY_1D_SPHERICAL)
			DiffusionTermSpherical(v, gamma, diffusion_term);
		else if (geometry_ == GEOMETRY_1D_CYLINDRICAL)
			DiffusionTermCylindrical(v, gamma, diffusion_term);
		else if (geometry_ == GEOMETRY_1D_PLANAR)
			DiffusionTermPlanar(v, gamma, diffusion_term);
	}

	void Grid1DMultiPurpose::DiffusionTermSpherical(const OpenSMOKE::OpenSMOKEVectorDouble& v, const OpenSMOKE::OpenSMOKEVectorDouble& gamma, OpenSMOKE::OpenSMOKEVectorDouble& diffusion_term)
	{
		diffusion_term[1] = 0.;
		diffusion_term[N_] = 0.;

		for (unsigned int i = 2; i <= N_ - 1; i++)
		{
			const double cPlus = 0.50*(x[i + 1] * x[i + 1] * gamma[i + 1] + x[i] * x[i] * gamma[i]);
			const double cMinus = 0.50*(x[i - 1] * x[i - 1] * gamma[i - 1] + x[i] * x[i] * gamma[i]);
			diffusion_term[i] = (cPlus / dx[i] * (v[i + 1] - v[i]) - cMinus / dx[i - 1] * (v[i] - v[i - 1])) * 2. / (dx[i - 1] + dx[i]) / (x[i] * x[i]);
		}
	}

	void Grid1DMultiPurpose::DiffusionTermSpherical(const Eigen::VectorXd& v, const Eigen::VectorXd& gamma, Eigen::VectorXd& diffusion_term)
	{
		diffusion_term(0) = 0.;
		diffusion_term(N_-1) = 0.;

		for (unsigned int i = 1; i <= N_ - 2; i++)
		{
			const double cPlus = 0.50*(x[i + 2] * x[i + 2] * gamma(i + 1) + x[i+1] * x[i+1] * gamma(i));
			const double cMinus = 0.50*(x[i] * x[i] * gamma(i - 1) + x[i+1] * x[i+1] * gamma(i));
			diffusion_term(i) = (cPlus / dx[i+1] * (v(i + 1) - v(i)) - cMinus / dx[i] * (v(i) - v(i - 1))) * 2. / (dx[i] + dx[i+1]) / (x[i+1] * x[i+1]);
		}
	}

	void Grid1DMultiPurpose::DiffusionTermCylindrical(const OpenSMOKE::OpenSMOKEVectorDouble& v, const OpenSMOKE::OpenSMOKEVectorDouble& gamma, OpenSMOKE::OpenSMOKEVectorDouble& diffusion_term)
	{
		diffusion_term[1] = 0.;
		diffusion_term[N_] = 0.;

		for (unsigned int i = 2; i <= N_ - 1; i++)
		{
			const double cPlus = 0.50*(x[i + 1] * gamma[i + 1] + x[i] * gamma[i]);
			const double cMinus = 0.50*(x[i - 1] * gamma[i - 1] + x[i] * gamma[i]);
			diffusion_term[i] = (cPlus / dx[i] * (v[i + 1] - v[i]) - cMinus / dx[i - 1] * (v[i] - v[i - 1])) * 2. / (dx[i - 1] + dx[i]) / x[i];
		}
	}

	void Grid1DMultiPurpose::DiffusionTermCylindrical(const Eigen::VectorXd& v, const Eigen::VectorXd& gamma, Eigen::VectorXd& diffusion_term)
	{
		diffusion_term(0) = 0.;
		diffusion_term(N_ - 1) = 0.;

		for (unsigned int i = 1; i <= N_ - 2; i++)
		{
			const double cPlus = 0.50*(x[i + 2] * gamma(i + 1) + x[i + 1] * gamma(i));
			const double cMinus = 0.50*(x[i] * gamma(i - 1) + x[i + 1] * gamma(i));
			diffusion_term(i) = (cPlus / dx[i + 1] * (v(i + 1) - v(i)) - cMinus / dx[i] * (v(i) - v(i - 1))) * 2. / (dx[i] + dx[i + 1]) / x[i + 1];
		}
	}

	void Grid1DMultiPurpose::DiffusionTermPlanar(const OpenSMOKE::OpenSMOKEVectorDouble& v, const OpenSMOKE::OpenSMOKEVectorDouble& gamma, OpenSMOKE::OpenSMOKEVectorDouble& diffusion_term)
	{
		diffusion_term[1] = 0.;
		diffusion_term[N_] = 0.;

		for (unsigned int i = 2; i <= N_ - 1; i++)
		{
			const double cPlus = 0.50*(gamma[i + 1] + gamma[i]);
			const double cMinus = 0.50*(gamma[i - 1] + gamma[i]);
			diffusion_term[i] = (cPlus / dx[i] * (v[i + 1] - v[i]) - cMinus / dx[i - 1] * (v[i] - v[i - 1])) * 2. / (dx[i - 1] + dx[i]);
		}
	}

	void Grid1DMultiPurpose::DiffusionTermPlanar(const Eigen::VectorXd& v, const Eigen::VectorXd& gamma, Eigen::VectorXd& diffusion_term)
	{
		diffusion_term(0) = 0.;
		diffusion_term(N_ - 1) = 0.;

		for (unsigned int i = 1; i <= N_ - 2; i++)
		{
			const double cPlus = 0.50*(gamma(i + 1) + gamma(i));
			const double cMinus = 0.50*(gamma(i - 1) + gamma(i));
			diffusion_term(i) = (cPlus / dx[i + 1] * (v(i + 1) - v(i)) - cMinus / dx[i] * (v(i) - v(i - 1))) * 2. / (dx[i] + dx[i + 1]);
		}
	}

	void Grid1DMultiPurpose::DiffusionTerm(const OpenSMOKE::OpenSMOKEVectorDouble& v, const OpenSMOKE::OpenSMOKEVectorDouble& gamma_a, const OpenSMOKE::OpenSMOKEVectorDouble& gamma_b, OpenSMOKE::OpenSMOKEVectorDouble& diffusion_term)
	{
		diffusion_term[1] = 0.;
		diffusion_term[N_] = 0.;

		for (unsigned int i = 2; i <= N_ - 1; i++)
		{
			const double cPlus = 0.50*(x[i + 1] * x[i + 1] * gamma_a[i + 1] * gamma_b[i + 1] + x[i] * x[i] * gamma_a[i] * gamma_b[i]);
			const double cMinus = 0.50*(x[i - 1] * x[i - 1] * gamma_a[i - 1] * gamma_b[i - 1] + x[i] * x[i] * gamma_a[i] * gamma_b[i]);
			diffusion_term[i] = (cPlus / dx[i] * (v[i + 1] - v[i]) - cMinus / dx[i - 1] * (v[i] - v[i - 1])) * 2. / (dx[i - 1] + dx[i]) / (x[i] * x[i]);
		}
	}

	double Grid1DMultiPurpose::IntegralValue(const OpenSMOKE::OpenSMOKEVectorDouble& v)
	{
		
		double sum = 0;
		
		if (geometry_ == GEOMETRY_1D_SPHERICAL)
		{
			for (unsigned int i = 1; i < N_; i++)
			{
				const double m = (v[i + 1] - v[i]) / (x[i + 1] - x[i]);
				sum += PI_4_OVER_3_ * (pow(x[i + 1], 3.)*(v[i] - m*x[i] + 3. / 4.*m*x[i + 1]) - pow(x[i], 3.)*(v[i] - m*x[i] + 3. / 4.*m*x[i]));
			}
		}
		else if (geometry_ == GEOMETRY_1D_CYLINDRICAL)
		{
			// TODO
		}
		else if (geometry_ == GEOMETRY_1D_PLANAR)
		{
			// TODO
		}

		return sum;
	}

	double Grid1DMultiPurpose::IntegralValue(const OpenSMOKE::OpenSMOKEVectorDouble& v1, const OpenSMOKE::OpenSMOKEVectorDouble& v2)
	{
		double sum = 0;
		
		if (geometry_ == GEOMETRY_1D_SPHERICAL)
		{
			for (unsigned int i = 1; i < N_; i++)
			{
				const double m = (v1[i + 1] * v2[i + 1] - v1[i] * v2[i]) / (x[i + 1] - x[i]);
				sum += PI_4_OVER_3_ * (pow(x[i + 1], 3.)*(v1[i] * v2[i] - m*x[i] + 3. / 4.*m*x[i + 1]) - pow(x[i], 3.)*(v1[i] * v2[i] - m*x[i] + 3. / 4.*m*x[i]));
			}
		}
		else if (geometry_ == GEOMETRY_1D_CYLINDRICAL)
		{
			// TODO
		}
		else if (geometry_ == GEOMETRY_1D_PLANAR)
		{
			// TODO
		}

		return sum;
	}

    void DimensionlessGrids(const unsigned int NPL, const unsigned int NPG, OpenSMOKE::OpenSMOKEVectorDouble& c_L, OpenSMOKE::OpenSMOKEVectorDouble& c_G, double& c_air)
	{
		OpenSMOKE::OpenSMOKEVectorDouble r_L(NPL);
		OpenSMOKE::OpenSMOKEVectorDouble r_G(NPG);
		OpenSMOKE::OpenSMOKEVectorDouble dr_L(NPL);
		OpenSMOKE::OpenSMOKEVectorDouble dr_G(NPG);

		ChangeDimensions(NPL, &c_L, true);
		ChangeDimensions(NPG, &c_G, true);

		// Starting grid step
		double dr_starting = 1. / 1000.;

		// Liquid phase
		{
			int NF;
			double alpha_L1;
			double alpha_L2;

			if (NPL == 61)
			{
				NF = 11;
				alpha_L1 = 0.725;
				alpha_L2 = 0.925;
			}
			else if (NPL == 51)
			{
				NF = 9;
				alpha_L1 = 0.675;
				alpha_L2 = 0.900;
			}
			else if (NPL == 41)
			{
				NF = 7;
				alpha_L1 = 0.630;
				alpha_L2 = 0.850;
			}
			else if (NPL == 31)
			{
				NF = 5;
				alpha_L1 = 0.500;
				alpha_L2 = 0.800;
			}
			else
			{
				OpenSMOKE::ErrorMessage("Liquid grid", "The possible number of points are: 30, 40, 50, 60");
			}

			// Liquid phase
			{
				r_L[1] = 0.;
				r_L[NPL] = 1.;
				for (unsigned int i = NPL; i >= NPL - NF; i--)
				{
					dr_L[i] = dr_starting / 30.;
					r_L[i - 1] = r_L[i] - dr_L[i];
				}
				for (int i = NPL - (NF + 1); i >= 3; i--)
				{
					const double w = double(i + 1) / double(NPL - NF);
					double alpha = alpha_L1*(1. - w) + w*alpha_L2;
					dr_L[i] = dr_L[i + 1] / alpha;
					r_L[i - 1] = r_L[i] - dr_L[i];
				}

				dr_L[2] = r_L[2] - r_L[1];
				dr_L[1] = 0.;
			}
		}
		
		// Gas grid
		{
			int N1, N2, N3;
			double alpha0, alpha1, alpha2, alpha3;
			double dr0;

			if (NPG == 200)
			{
				N1 = 30;
				N2 = 160;
				N3 = 170;
				alpha0 = 1.03125;
				alpha1 = 1.01125;
				alpha2 = 1.05125;
				alpha3 = 1.11;
				dr0 = 1. / 100.;
			}
			else if (NPG == 250)
			{
				N1 = 30;
				N2 = 200;
				N3 = 220;
				alpha0 = 1.02125;
				alpha1 = 1.006125;
				alpha2 = 1.05125;
				alpha3 = 1.12;
				dr0 = 1. / 50.;
			}
			else if (NPG == 300)
			{
				N1 = 40;
				N2 = 230;
				N3 = 270;
				alpha0 = 1.025125;
				alpha1 = 1.004125;
				alpha2 = 1.045125;
				alpha3 = 1.12;
				dr0 = 1. / 100.;
			}
			else if (NPG == 360)
			{
				N1 = 50;
				N2 = 280;
				N3 = 330;
				alpha0 = 1.02125;
				alpha1 = 1.0034125;
				alpha2 = 1.040125;
				alpha3 = 1.11;
				dr0 = 8./1000.;
			}
			else if (NPG == 410)
			{
				N1 = 50;
				N2 = 320;
				N3 = 370;
				alpha0 = 1.02125;
				alpha1 = 1.0025125;
				alpha2 = 1.035125;
				alpha3 = 1.09;
				dr0 = 0.008;
			}
			else if (NPG == 1000)
			{
				N1 = 250;
				N2 = 500;
				N3 = 250;
				alpha0 = 1.0;
				alpha1 = 1.0;
				alpha2 = 1.0;
				alpha3 = 1.0;
				dr0 = 20. / double(NPG);
			}
			else
			{
				OpenSMOKE::ErrorMessage("Gas grid", "The possible number of points are: 200, 250 300, 360, 410, 1000");
			}

			// Build grid
			{
				r_G[1] = 1.;
				dr_G[1] = dr0;

				for (int i = 2; i <= N1; i++)
				{
					r_G[i] = r_G[i - 1] + dr_G[i - 1];
					dr_G[i] = alpha0*dr_G[i - 1];
				}

				for (int i = N1 + 1; i <= N2; i++)
				{
					r_G[i] = r_G[i - 1] + dr_G[i - 1];
					dr_G[i] = alpha1*dr_G[i - 1];
				}

				for (int i = N2 + 1; i <= N3; i++)
				{
					r_G[i] = r_G[i - 1] + dr_G[i - 1];
					dr_G[i] = alpha2*dr_G[i - 1];
				}

				for (int i = N3 + 1; i <= static_cast<int>(NPG) - 1; i++)
				{
					r_G[i] = r_G[i - 1] + dr_G[i - 1];
					dr_G[i] = alpha3*dr_G[i - 1];
				}
				
				// Common lines

				r_G[NPG] = r_G[NPG - 1] + dr_G[NPG - 1];
				dr_G[NPG] = 0.;
			}

			// Print grid on the screen
			bool print_grid_on_the_screen = false;
			if (print_grid_on_the_screen == true)
			{
				std::cout << "Gas(1)" << std::endl;
				for (int i = 1; i <= N1; i++)
					std::cout << i << " " << r_G[i] << std::endl;

				std::cout << "Gas(2)" << std::endl;
				for (int i = N1 + 1; i <= N2; i++)
					std::cout << i << " " << r_G[i] << std::endl;

				std::cout << "Gas(3)" << std::endl;
				for (int i = N2 + 1; i <= N3; i++)
					std::cout << i << " " << r_G[i] << std::endl;

				std::cout << "Gas(4)" << std::endl;
				for (unsigned int i = N3 + 1; i <= NPG; i++)
					std::cout << i << " " << r_G[i] << std::endl;
			}
		}

		// Final operations
		{
			for (unsigned int i = 1; i <= NPL; i++)
				c_L[i] = r_L(i) / 1.;

			for (unsigned int i = 1; i <= NPG; i++)
				c_G[i] = (r_G[i] - 1.) / (r_G[NPG] - 1.);

			c_air = r_G[NPG];
		}
	}


	void DimensionlessEquispacedGrids(const unsigned int NPL, const unsigned int NPG, const double ratio_radii, OpenSMOKE::OpenSMOKEVectorDouble& c_L, OpenSMOKE::OpenSMOKEVectorDouble& c_G, double& c_air)
	{
		OpenSMOKE::OpenSMOKEVectorDouble r_L(NPL);
		OpenSMOKE::OpenSMOKEVectorDouble r_G(NPG);

		ChangeDimensions(NPL, &c_L, true);
		ChangeDimensions(NPG, &c_G, true);

		const double dr_liquid = 1. / double(NPL - 1.);
		r_L[1] = 0.;
		for (unsigned int i = 2; i <= NPL; i++)
			r_L[i] = r_L[i - 1] + dr_liquid;

		const double dr_gas = (ratio_radii - 1.) / double(NPG - 1.);
		r_G[1] = 1.;
		for (unsigned int i = 2; i <= NPG; i++)
			r_G[i] = r_G[i - 1] + dr_gas;

		for (unsigned int i = 1; i <= NPL; i++)
			c_L[i] = r_L(i) / 1.;

		for (unsigned int i = 1; i <= NPG; i++)
			c_G[i] = (r_G[i] - 1.) / (r_G[NPG] - 1.);

		c_air = r_G[NPG];
	}

	void Interpolate(const OpenSMOKE::OpenSMOKEVectorDouble& xOriginal, const OpenSMOKE::OpenSMOKEVectorDouble& vOriginal, const OpenSMOKE::OpenSMOKEVectorDouble& xInterpolated, OpenSMOKE::OpenSMOKEVectorDouble& vInterpolated)
	{
		const unsigned int np_old = xOriginal.Size();
		const unsigned int np_new = xInterpolated.Size();
		OpenSMOKE::ChangeDimensions(np_new, &vInterpolated, true);
		vInterpolated[1] = vOriginal[1];
		vInterpolated[np_new] = vOriginal[np_old];

		unsigned int jLast = 2;
		for (unsigned int i = 2; i <= np_new - 1; i++)
		{
			for (unsigned int j = jLast; j <= np_old; j++)
			if (xOriginal[j] >= xInterpolated[i])
			{
				jLast = j;
				vInterpolated[i] = vOriginal[j - 1] + (vOriginal[j] - vOriginal[j - 1]) / (xOriginal[j] - xOriginal[j - 1]) * (xInterpolated[i] - xOriginal[j - 1]);
				break;
			}
		}
	}

	void Split(const OpenSMOKE::OpenSMOKEVectorDouble& xold, const unsigned int block, OpenSMOKE::OpenSMOKEVectorDouble& xnew)
	{
		const unsigned int np_old = xold.Size();
		const unsigned int np_new = (np_old - 1)*block + 1;
		OpenSMOKE::ChangeDimensions(np_new, &xnew, true);

		unsigned int count = 1;
		for (unsigned int i = 1; i < np_old; i++)
		for (unsigned int j = 1; j <= block; j++)
			xnew[count++] = xold[i] + (j - 1)*(xold[i + 1] - xold[i]) / double(block);
		xnew[np_new] = xold[np_old];
	}
}
