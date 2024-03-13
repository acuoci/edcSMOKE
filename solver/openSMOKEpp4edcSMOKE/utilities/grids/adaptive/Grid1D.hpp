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
|   License                                                               |
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

#include "math/OpenSMOKEUtilities.h"


namespace OpenSMOKE
{
	
	Grid1D::Grid1D(const int Np, const double x0, const double xF) :
		Np_(Np),
		x0_(x0),
		xF_(xF)
	{
		
		Ni_ = Np_ - 1;
		L_ = xF_-x0_;

		Allocate();

		dxw_.setConstant(L_ / static_cast<double>(Ni_));
		dxw_(0) = 0.;

		fixed_point_ = 0;
		x_fixed_point_ = 0.;

		Build();

		grid_adapter_ = new Adapter_Grid1D();
	}

	Grid1D::Grid1D(const int Np, const double x0, const double xF, const double alpha) :
		Np_(Np),
		x0_(x0),
		xF_(xF)
	{
		Ni_ = Np_ - 1;
		L_ = xF_ - x0_;

		Allocate();

		double sum = 1.;
		for (unsigned int i = 1; i <= Np_-2; i++)
			sum += std::pow(alpha, static_cast<double>(i));

		const double dx0 = L_ / sum;

		dxw_(0) = 0.;
		for (unsigned int i = 1; i < Np_; i++)
			dxw_(i) = dx0*std::pow(alpha, static_cast<double>(i-1));

		fixed_point_ = 0;
		x_fixed_point_ = 0.;

		Build();

		grid_adapter_ = new Adapter_Grid1D();
	}

	Grid1D::Grid1D(const Eigen::VectorXd& coordinates) 
	{
		BuildFromCoordinates(coordinates);
	}

	Grid1D::Grid1D(const std::vector<double>& coordinates)
	{
		Eigen::VectorXd coordinates_(coordinates.size());
		for (unsigned int i = 0; i < coordinates_.size(); i++)
			coordinates_(i) = coordinates[i];

		BuildFromCoordinates(coordinates_);
	}

	void Grid1D::BuildFromCoordinates(const Eigen::VectorXd& coordinates)
	{
		Np_ = static_cast<int>(coordinates.size());
		Ni_ = Np_ - 1;

		x0_ = coordinates(0);
		xF_ = coordinates(Np_ - 1);
		L_ = xF_ - x0_;

		Allocate();

		dxw_(0) = 0.;
		for (unsigned int i = 1; i < Np_; i++)
			dxw_(i) = coordinates(i) - coordinates(i - 1);

		fixed_point_ = 0;
		x_fixed_point_ = 0.;

		Build();

		grid_adapter_ = new Adapter_Grid1D();
	}

	void linear_interpolation(const Eigen::VectorXd& x_old, const std::vector<Eigen::VectorXd>& v_old, const Eigen::VectorXd& x_new, std::vector<Eigen::VectorXd>& v_new)
	{
		unsigned int n_old = x_old.size();
		unsigned int n_new = x_new.size();
		unsigned int n = v_old.size();

		v_new.resize(n);
		for (unsigned int i = 0; i < n; i++)
			v_new[i].resize(n_new);

		// First point
		for (unsigned int k = 0; k < n; k++)
			v_new[k](0) = v_old[k](0);
		
		// Internal points
		for (unsigned int i = 1; i < n_new - 1; i++)
		{
			for (unsigned int j = 0; j < n_old; j++)
			if (x_old(j) >= x_new(i))
			{
				for (unsigned int k = 0; k < n; k++)
					v_new[k](i) = v_old[k](j - 1) + (v_old[k](j) - v_old[k](j - 1)) / (x_old(j) - x_old(j - 1))*(x_new(i) - x_old(j - 1));
				break;
			}
		}

		// Last point
		for (unsigned int k = 0; k < n; k++)
			v_new[k](n_new - 1) = v_old[k](n_old - 1);
	}

	void Grid1D::Double(const std::vector<Eigen::VectorXd>& phi, std::vector<Eigen::VectorXd>& phi_new)
	{
		Eigen::VectorXd x_new(2*x_.size()-1);

		for (int i = 0; i < x_.size(); i++)
			x_new(2 * i) = x_(i);
		for (int i = 1; i < x_new.size() - 1; i += 2)
			x_new(i) = 0.50*(x_new(i - 1) + x_new(i + 1));

		linear_interpolation(x_, phi, x_new, phi_new);
		Update(x_new);

		for (int i = 0; i < x_new.size(); i++)
			std::cout << std::scientific << x_new(i) << std::endl;
	}

	void Grid1D::Refine(const double xA, const double xB, const std::vector<Eigen::VectorXd>& phi, std::vector<Eigen::VectorXd>& phi_new)
	{
		unsigned int iA = 0;
		for (int i = 0; i < x_.size(); i++)
			if (x_(i) >= xA)
			{
				iA = i;
				break;
			}

		unsigned int iB = 0;
		for (int i = iA; i < x_.size(); i++)
			if (x_(i) >= xB)
			{
				iB = i;
				break;
			}

		std::vector<double> x_prov(x_.size());
		for (int i = 0; i < x_.size(); i++)
			x_prov[i] = x_(i);
		for (unsigned int i = iA; i < iB; i++)
			x_prov.push_back(0.50*(x_(i) + x_(i + 1)));
		std::sort(x_prov.begin(), x_prov.end());

		Eigen::VectorXd x_new(x_prov.size());
		for (int i = 0; i < x_new.size(); i++)
			x_new(i) = x_prov[i];

		linear_interpolation(x_, phi, x_new, phi_new);
		Update(x_new);

		for (int i = 0; i < x_new.size(); i++)
			std::cout << std::scientific << x_new(i) << std::endl;
	}


	Adapter_Grid1D_Status Grid1D::Refine(const std::vector<Eigen::VectorXd>& phi, std::vector<Eigen::VectorXd>& phi_new)
	{
		Eigen::VectorXd x_new;
		Adapter_Grid1D_Status status = grid_adapter_->Refine(x_, phi, x_new);

		if (status != OpenSMOKE::NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED)
		{
			linear_interpolation(x_, phi, x_new, phi_new);
			Update(x_new);
		}

		return status;
	}
	
	Adapter_Grid1D_Status Grid1D::Regrid(const unsigned int index_psi, const std::vector<Eigen::VectorXd>& phi, std::vector<Eigen::VectorXd>& phi_new)
	{
		const unsigned int n = phi.size();
		const unsigned int np = phi[0].size();
		Eigen::VectorXd psi = phi[index_psi];
		
		// Regrid
		Eigen::VectorXd x_new;
		grid_adapter_->Regrid(x_, psi, x_new);


		// Check for fixed point
		if (fixed_point_ != 0)
		{
			unsigned int j_fixed_point = 0;

			Eigen::VectorXd saved(n);
			for (unsigned int i = 0; i < n; i++)
				saved(i) = phi[i](fixed_point_);

			const double epsilon = 1e-12;
			for (int j = 1; j < x_new.size(); j++)
			{
				if (x_new(j) < x_fixed_point_*(1. + epsilon) && x_new(j) > x_fixed_point_*(1. - epsilon))
				{
					x_new(j) = x_fixed_point_;
					j_fixed_point = j;
					break;
				}
			}

			if (j_fixed_point == 0)
			{
				for (int j = 1; j < x_new.size(); j++)
				if (x_fixed_point_ > x_new(j) && x_fixed_point_ < x_new(j + 1))
				{
					j_fixed_point = j;
					break;
				}

				x_new.conservativeResize(x_new.size() + 1);
				for (unsigned int j = x_new.size() - 1; j>j_fixed_point+1; j--)
					x_new(j) = x_new(j - 1);
				x_new(j_fixed_point + 1) = x_fixed_point_;
				fixed_point_ = j_fixed_point + 1;
			}

			linear_interpolation(x_, phi, x_new, phi_new);
			for (unsigned int i = 0; i < n; i++)
				phi_new[i](fixed_point_) = saved(i);
		}
		else
		{
			linear_interpolation(x_, phi, x_new, phi_new);
		}

		Update(x_new);

		return REGRID_SUCCESS;
	}
	
	void Grid1D::Update(const Eigen::VectorXd& coordinates)
	{
		Np_ = static_cast<int>(coordinates.size());
		Ni_ = Np_ - 1;

		x0_ = coordinates(0);
		xF_ = coordinates(Np_ - 1);
		L_ = xF_ - x0_;

		Allocate();

		dxw_(0) = 0.;
		for (unsigned int i = 1; i < Np_; i++)
			dxw_(i) = coordinates(i) - coordinates(i - 1);

		Build();

		SetFixedPoint(x_fixed_point_);
	}

	void Grid1D::Allocate()
	{
		x_.resize(Np_);
		dxe_.resize(Np_);
		dxw_.resize(Np_);
		ux_.resize(Np_);
		x2_.resize(Np_);
		udxw_.resize(Np_);
		udxe_.resize(Np_);
		dxc_.resize(Np_);
		udxc_.resize(Np_);
		dxc_over_2_.resize(Np_);
		udxc_over_2_.resize(Np_);
		cenp_.resize(Np_);
		cenc_.resize(Np_);
		cenm_.resize(Np_);
	}

	void Grid1D::Build()
	{
		x_(0) = x0_;
		for (unsigned int i = 1; i < Np_; i++)
			x_(i) = x_(i-1) + dxw_(i);

		for (unsigned int i = 0; i < Np_; i++)
			x2_(i) = x_(i) * x_(i);

		if (x_(0) != 0.) ux_(0) = 1. / x_(0);
		else ux_(0) = 0.;
		for (unsigned int i = 1; i < Np_; i++)
			ux_(i) = 1. / x_(i);

		for (unsigned int i = 1; i < Np_; i++)
			udxw_(i) = 1. / dxw_(i);

		for (unsigned int i = 0; i < Np_ - 1; i++)
		{
			dxe_(i) = x_(i+1) - x_(i);
			udxe_(i) = 1. / dxe_(i);
		}

		for (unsigned int i = 1; i < Ni_; i++)
		{
			dxc_(i) = (x_(i + 1) - x_(i-1));
			udxc_(i) = 1. / dxc_(i);
			dxc_over_2_(i) = 0.50*dxc_(i);
			udxc_over_2_(i) = 1. / dxc_over_2_(i);
		}

		for (unsigned int i = 1; i < Ni_; i++)
		{
			cenp_(i) = dxw_(i) * udxe_(i) * udxc_(i);
			cenc_(i) = (dxe_(i) - dxw_(i))*udxe_(i) * udxw_(i);
			cenm_(i) = dxe_(i) * udxw_(i) * udxc_(i);
		}
	}

	void Grid1D::Derivative(const enum derivative_type type, const Eigen::VectorXd &u, const Eigen::VectorXd &phi, Eigen::VectorXd* dphi)
	{

		if (type == DERIVATIVE_1ST_UPWIND)
		{
			(*dphi)(0) = (phi(1) - phi(0))*udxe_(0);
			for (unsigned int i = 1; i < Ni_; i++)
				(*dphi)(i) = (u(i) <= 0.) ? (-phi(i) + phi(i + 1))*udxe_(i) : (-phi(i - 1) + phi(i))*udxw_(i);
			(*dphi)(Ni_) = (phi(Ni_) - phi(Ni_ - 1))*udxw_(Ni_);
		}

		else if (type == DERIVATIVE_1ST_CENTERED)
		{
			(*dphi)(0) = (phi(1) - phi(0))*udxe_(0);
			for (unsigned int i = 1; i < Ni_; i++)
				(*dphi)(i) = -cenm_(i) * phi(i - 1) + cenc_(i) * phi(i) + cenp_(i) * phi(i + 1);
			(*dphi)(Ni_) = (phi(Ni_) - phi(Ni_ - 1))*udxw_(Ni_);
		}

		else if (type == DERIVATIVE_1ST_BACKWARD)
		{
			(*dphi)(0) = (phi(1) - phi(0))*udxe_(0); 
			for (unsigned int i = 1; i < Np_; i++)
				(*dphi)(i) = (-phi(i - 1) + phi(i))*udxw_(i);
		}

		else if (type == DERIVATIVE_1ST_FORWARD)
		{
			for (unsigned int i = 0; i < Ni_; i++)
				(*dphi)(i) = (-phi(i) + phi(i + 1))*udxe_(i);
			(*dphi)(Ni_) = (phi(Ni_) - phi(Ni_ - 1))*udxw_(Ni_);
		}
	}

	void Grid1D::Derivative(const enum derivative_type type, const Eigen::VectorXd &u, const std::vector<Eigen::VectorXd> &phi, std::vector<Eigen::VectorXd>* dphi)
	{
		if (type == DERIVATIVE_1ST_UPWIND)
		{ 
			for (int k = 0; k < phi[0].size(); k++)
			{
				(*dphi)[0](k) = (phi[1](k)-phi[0](k))*udxe_(0);
				for (unsigned int i = 1; i < Ni_; i++)
					(*dphi)[i](k) = (u(i) <= 0.) ? (-phi[i](k)+phi[i + 1](k))*udxe_(i) : (-phi[i - 1](k)+phi[i](k))*udxw_(i);
				(*dphi)[Ni_](k) = (phi[Ni_](k)-phi[Ni_ - 1](k))*udxw_(Ni_);
			}
		}

		else if (type == DERIVATIVE_1ST_CENTERED)
		{
			for (int k = 0; k < phi[0].size(); k++)
			{
				(*dphi)[0](k) = (phi[1](k) - phi[0](k))*udxe_(0);
				for (unsigned int i = 1; i < Ni_; i++)
					(*dphi)[i](k) = -cenm_(i) * phi[i - 1](k) + cenc_(i) * phi[i](k) + cenp_(i) * phi[i + 1](k);
				(*dphi)[Ni_](k) = (phi[Ni_](k) - phi[Ni_ - 1](k))*udxw_(Ni_);
			}
		}

		else if (type == DERIVATIVE_1ST_BACKWARD)
		{
			for (int k = 0; k < phi[0].size(); k++)
			{
				(*dphi)[0](k) = (phi[1](k) - phi[0](k))*udxe_(0);
				for (unsigned int i = 1; i < Np_; i++)
					(*dphi)[i](k) = (-phi[i - 1](k) + phi[i](k))*udxw_(i);

			}
		}

		else if (type == DERIVATIVE_1ST_FORWARD)
		{
			for (int k = 0; k < phi[0].size(); k++)
			{
				for (unsigned int i = 0; i < Ni_; i++)
					(*dphi)[i](k) = (-phi[i](k)+phi[i + 1](k))*udxe_(i);
				(*dphi)[Ni_](k) = (phi[Ni_](k)-phi[Ni_ - 1](k))*udxw_(Ni_);
			}
		}
	}

	void Grid1D::SecondDerivative(const Eigen::VectorXd &phi, Eigen::VectorXd* d2phi)
	{
		for (unsigned int i = 1; i < Ni_; i++)
			(*d2phi)(i) = ( (phi(i + 1) - phi(i) )*udxe_(i) - (phi(i) - phi(i-1))*udxw_(i)) * udxc_over_2_(i);
	}

	void Grid1D::SecondDerivative(const std::vector<Eigen::VectorXd> &phi, std::vector<Eigen::VectorXd>* d2phi)
	{
		for (int k = 0; k < phi[0].size(); k++)
		{
			(*d2phi)[0](k) = 0.;
			for (unsigned int i = 1; i < Ni_; i++)
				(*d2phi)[i](k) = ((phi[i + 1](k) - phi[i](k))*udxe_(i) - (phi[i](k) - phi[i - 1](k))*udxw_(i)) * udxc_over_2_(i);
			(*d2phi)[Ni_](k) = 0.;
		}
	}

	void Grid1D::SecondDerivative(const Eigen::VectorXd &coeff, const Eigen::VectorXd &phi, Eigen::VectorXd* d2phi)
	{
		for (unsigned int i = 1; i < Ni_; i++)
			(*d2phi)(i) = (	0.50*(coeff(i + 1) + coeff(i))*(phi(i + 1) - phi(i))*udxe_(i) - 
							0.50*(coeff(i - 1) + coeff(i))*(phi(i) - phi(i - 1))*udxw_(i) ) * udxc_over_2_(i);
	}

	void Grid1D::SetFixedPoint(const double x_fixed_point)
	{
		if (x_fixed_point == 0. && fixed_point_ == 0)
			return;
			
		const double epsilon = 1e-15;
		for (unsigned int i = 1; i < Ni_; i++)
		{
			if (x_(i) < x_fixed_point*(1. + epsilon) && x_(i) > x_fixed_point*(1. - epsilon))
			{
				fixed_point_ = i;
				x_fixed_point_ = x_fixed_point;
				return;
			}
		}

		OpenSMOKE::FatalErrorMessage("The grid was not able to find the requested fixed point");
	}

	void Grid1D::ResetFixedPoint()
	{
		fixed_point_ = 0;
		x_fixed_point_ = 0.;
	}

	Grid1D::Grid1D(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, Eigen::VectorXd& w)
	{
		Grammar_Grid1D grammar;
		dictionary.SetGrammar(grammar);

		double length = 0.;
		double x_fixed_point = 0.;
		int points = 0;
		std::string type;

		// @Type (To add)
		{
			if (dictionary.CheckOption("@Type") == true)
			{
				dictionary.ReadString("@Type", type);
				if (type != "centered" && type != "database" && type != "equispaced" && type != "liquid-pool")
					OpenSMOKE::FatalErrorMessage("Unknown grid type. Known types are: centered | database | equispaced | liquid-pool");
			}
		}

		// @Length
		{
			if (dictionary.CheckOption("@Length") == true)
			{
				double value;
				std::string units;
				dictionary.ReadMeasure("@Length", value, units);

				if (units == "m")			length = value;
				else if (units == "cm")		length = value / 100.;
				else if (units == "mm")		length = value / 1000.;
				else OpenSMOKE::FatalErrorMessage("Unknown length units");

				x_fixed_point = length / 2.;
			}
		}

		// @InitialPoints
		{
			if (dictionary.CheckOption("@InitialPoints") == true)
			{
				dictionary.ReadInt("@InitialPoints", points);
				if (points < 7)
					OpenSMOKE::FatalErrorMessage("The minimum number of points is 7.");
			}
		}

		// @FixedPoint
		{
			if (dictionary.CheckOption("@FixedPoint") == true)
			{
				std::string units;
				dictionary.ReadMeasure("@FixedPoint", x_fixed_point, units);
				if (units == "m")		x_fixed_point *= 1.;
				else if (units == "cm")	x_fixed_point *= 1.e-2;
				else if (units == "mm")	x_fixed_point *= 1.e-3;
				else OpenSMOKE::FatalErrorMessage("Wrong units in @FixedPoint option");
			}
		}

		// Setup of grid
		{
			Eigen::VectorXd x(points);
			w.resize(points);

			if (type == "equispaced")
			{
				const double dx = length / double(points - 1);
				const double dw = 1. / double(points - 1);

				x(0) = 0.;
				w(0) = 0.;
				for (int i = 1; i < points; i++)
				{
					x(i) = x(i-1) + dx;
					w(i) = w(i-1) + dw;
				}

				// Fixed point is in the middle of computational domain
				x_fixed_point = x(points / 2);
			}

			else if (type == "centered")
			{
				double center = 0.5*length;
				double width = 0.075*length;

				// @InitialPoints
				{
					if (dictionary.CheckOption("@Center") == true)
					{
						double value;
						std::string units;
						dictionary.ReadMeasure("@Center", value, units);

						if (units == "m")			center = value;
						else if (units == "cm")		center = value / 100.;
						else if (units == "mm")		center = value / 1000.;
						else OpenSMOKE::FatalErrorMessage("Unknown length units");
					}
				}

				// @Width
				{
					if (dictionary.CheckOption("@Width") == true)
					{
						double value;
						std::string units;
						dictionary.ReadMeasure("@Width", value, units);

						if (units == "m")			width = value;
						else if (units == "cm")		width = value / 100.;
						else if (units == "mm")		width = value / 1000.;
						else OpenSMOKE::FatalErrorMessage("Unknown length units");
					}
				}

				const double x_left = center - 0.50*width;
				const double x_right = center + 0.50*width;

				x(0) = 0.;
				x(1) = x_left / 2.;
				x(2) = x_left;

				x(points - 1) = length;
				x(points - 2) = (x_right + length) / 2.;
				x(points - 3) = x_right;

				int remaining_points = points - 6;

				if (remaining_points % 2)
				{
					int i_center = (points - 1) / 2;
					x(i_center) = center;

					int n_intervals = (remaining_points - 1) / 2 + 1;

					double alpha = 1.5;
					double sum = 1.;
					for (int i = 1; i < n_intervals; i++)
						sum += std::pow(alpha, i);
					double dx = (center - x_left) / sum;


					for (int i = 1; i < n_intervals; i++)
					{
						x(i_center + i) = x(i_center + i - 1) + dx*std::pow(alpha, i - 1);
						x(i_center - i) = x(i_center - i + 1) - dx*std::pow(alpha, i - 1);
					}
				}
				else
				{
					int i_center = points / 2;
					x(i_center) = center;

					// Left
					{
						int n_intervals = remaining_points / 2 + 1;

						double alpha = 2.5;
						double sum = 1.;
						for (int i = 1; i < n_intervals; i++)
							sum += pow(alpha, i);
						double dx = (center - x_left) / sum;

						for (int i = 1; i < n_intervals; i++)
							x(i_center - i) = x(i_center - i + 1) - dx*pow(alpha, i - 1);
					}

					// Right
					{
						int n_intervals = remaining_points / 2;

						double alpha = 2.5;
						double sum = 1.;
						for (int i = 1; i < n_intervals; i++)
							sum += pow(alpha, i);
						double dx = (center - x_left) / sum;

						for (int i = 1; i < n_intervals; i++)
							x(i_center + i) = x(i_center + i - 1) + dx*pow(alpha, i - 1);
					}
				}

				w(0) = 0.;
				w(1) = 0.;
				w(2) = 0.;
				w(points - 1) = 1.;
				w(points - 2) = 1.;
				w(points - 3) = 1.;

				for (int i = 3; i < points - 3; i++)
					w(i) = w(2) + 1. / (x(points - 3) - x(2))*(x(i) - x(2));
			}

			else if (type == "database")
			{
				if (points == 12)
				{
					x << 0.000000E+00, 2.450200E-01, 4.802400E-01, 4.900400E-01, 4.933600E-01, 4.966800E-01, 5.000000E-01, 5.033200E-01, 5.066400E-01, 5.099600E-01, 7.549800E-01, 1.000000E+00;
					w << 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 1.666667E-01, 3.333333E-01, 5.000000E-01, 6.666667E-01, 8.333333E-01, 1.000000E+00, 1.000000E+00, 1.000000E+00;
				}
				else if (points == 16)
				{
					x << 0.000000E+00, 1.225100E-01, 2.450200E-01, 3.626300E-01, 4.802400E-01, 4.900400E-01, 4.933600E-01, 4.966800E-01, 5.000000E-01, 5.033200E-01, 5.066400E-01, 5.099600E-01, 6.324700E-01, 7.549800E-01, 8.774900E-01, 1.000000E+00;
					w << 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 1.666667E-01, 3.333333E-01, 5.000000E-01, 6.666667E-01, 8.333333E-01, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00;
				}
				else if (points == 20)
				{
					x << 0.000000E+00, 1.225100E-01, 2.450200E-01, 3.626300E-01, 4.214350E-01, 4.508375E-01, 4.802400E-01, 4.900400E-01, 4.933600E-01, 4.966800E-01, 5.000000E-01, 5.033200E-01, 5.066400E-01, 5.099600E-01, 6.324700E-01, 7.549800E-01, 8.774900E-01, 9.387450E-01, 9.693725E-01, 1.000000E+00;
					w << 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 1.666667E-01, 3.333333E-01, 5.000000E-01, 6.666667E-01, 8.333333E-01, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00;
				}
				else
					OpenSMOKE::FatalErrorMessage("Database grids are available only for a number of points of: 12, 16, 20");

				x *= length;
			}

			else if (type == "liquid-pool")
			{

				const double grid_alpha_fuel = 1.10;
				const double grid_alpha_ox = 1.10;
				const double grid_points_fraction = 0.80;
				const double grid_distance_fraction = 0.60;
				
				const int NpA = static_cast<int>(grid_points_fraction*static_cast<double>(points));
				const int NiA = NpA-1;
				const int NpB = points-NiA;
				const int NiB = NpB-1;

				Eigen::VectorXd dxw(points);
				dxw.setZero();

				double sum=1.;
				for(unsigned int i=2;i<=NiA;i++)
					sum += std::pow(grid_alpha_fuel, static_cast<double>(i-1));

				double dx1=length* grid_distance_fraction /sum;
				dxw(0)=0.;
				for(unsigned int i=2;i<=NpA;i++)
					dxw(i-1)=dx1*std::pow(grid_alpha_fuel, static_cast<double>(i-2));

				sum=1.;
				for(unsigned int i=2;i<=NiB;i++)
					sum += std::pow(grid_alpha_ox, static_cast<double>(i-1));

				dx1=length*(1.- grid_distance_fraction)/sum;
				dxw(points-1)=dx1;
				for(unsigned int i=points-1;i>=NpA+1;i--)
				{
					const unsigned int j=points-i;
					dxw(i-1)=dx1*std::pow(grid_alpha_ox, j);
				}

				x(0) = 0.;
				w(0) = 0.;
				for (int i = 1; i < points; i++)
				{
					x(i) = x(i-1) + dxw(i);
					w(i) = w(i-1) + dxw(i)/length;
				}

				// Fixed point is in the middle of computational domain
				x_fixed_point = x(points / 2);
			}

			BuildFromCoordinates(x);
			SetFixedPoint(x_fixed_point);
		}

		// @MaxPoints
		{
			if (dictionary.CheckOption("@MaxPoints") == true)
			{
				int max_points;
				dictionary.ReadInt("@MaxPoints", max_points);
				grid_adapter_->SetMaxPoints(max_points);
			}
		}

		// @MaxAdaptivePoints
		{
			if (dictionary.CheckOption("@MaxAdaptivePoints") == true)
			{
				int max_adaptive_points;
				dictionary.ReadInt("@MaxAdaptivePoints", max_adaptive_points);
				grid_adapter_->SetMaxPointsToBeAdded(max_adaptive_points);
			}
		}

		// @GradientCoefficient
		{
			if (dictionary.CheckOption("@GradientCoefficient") == true)
			{
				double gradient_coefficient;
				dictionary.ReadDouble("@GradientCoefficient", gradient_coefficient);
				grid_adapter_->SetCoefficientGradient(gradient_coefficient);
			}
		}

		// @CurvatureCoefficient
		{
			if (dictionary.CheckOption("@CurvatureCoefficient") == true)
			{
				double curvature_coefficient;
				dictionary.ReadDouble("@CurvatureCoefficient", curvature_coefficient);
				grid_adapter_->SetCoefficientCurvature(curvature_coefficient);
			}
		}

		// @Threshold
		{
			if (dictionary.CheckOption("@Threshold") == true)
			{
				double threshold;
				dictionary.ReadDouble("@Threshold", threshold);
				grid_adapter_->SetThreshold(threshold);
			}
		}

		// @RegridPoints
		{
			if (dictionary.CheckOption("@RegridPoints") == true)
			{
				int regrid_points_;
				dictionary.ReadInt("@RegridPoints", regrid_points_);
				grid_adapter_->SetNumberOfRegridPoints(regrid_points_);
			}
		}

		
	}
}
