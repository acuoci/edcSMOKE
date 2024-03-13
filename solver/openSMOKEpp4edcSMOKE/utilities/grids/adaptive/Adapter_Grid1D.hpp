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

namespace OpenSMOKE
{
	Adapter_Grid1D::Adapter_Grid1D()
	{
		coeff_grad_ = 0.1;
		coeff_curv_ = 0.5;
		max_points_ = 300;
		max_points_to_be_added_ = 10;
		threshold_ = 1.e-7;
		fraction_active_points_ = 0.60;
		ratio_gradient_curvature_ = 1.50;
		regrid_points_ = 20;
	}

	Adapter_Grid1D_Status Adapter_Grid1D::Regrid(const Eigen::VectorXd& x, const Eigen::VectorXd& phi, Eigen::VectorXd& x_new)
	{
		x_new.resize(regrid_points_);
		const unsigned int n_old = x.size();

		const double r0 = 1. - fraction_active_points_;
		const double r1 = fraction_active_points_*ratio_gradient_curvature_ / (ratio_gradient_curvature_ + 1.);
		const double r2 = fraction_active_points_ - r1;

		double tv1 = 0.;
		for (unsigned int i = 1; i < n_old; i++)
			tv1 += std::fabs(phi(i) - phi(i - 1));

		double tv2 = 0.;
		for (unsigned int i = 1; i < n_old - 1; i++)
			tv2 += std::fabs((phi(i + 1) - phi(i)) / (x(i + 1) - x(i)) - (phi(i) - phi(i - 1)) / (x(i) - x(i - 1)));

		const double L = x(n_old - 1) - x(0);
		const double b1 = r1*L / (tv1*(1. - fraction_active_points_));
		const double b2 = r2*L / (tv2*(1. - fraction_active_points_));

		Eigen::VectorXd weights(n_old);
		weights(0) = 0.;
		for (unsigned int i = 1; i < n_old; i++)
		{
			const double dx = x(i) - x(i - 1);
			weights(i) = dx + +weights(i - 1) +
				b1*std::fabs(phi(i) - phi(i - 1)) +
				b2*std::fabs((phi(i + 1) - phi(i)) / (x(i + 1) - x(i)) - (phi(i) - phi(i - 1)) / (x(i) - x(i - 1)));
		}

		for (unsigned int i = 0; i < n_old; i++)
			weights(i) /= weights(n_old - 1);

		x_new(0) = x(0);
		x_new(regrid_points_ - 1) = x(n_old - 1);
		double deta = 1. / double(regrid_points_ - 1);

		int istart = 1;
		for (unsigned int j = 1; j < n_old - 1; j++)
		{
			double etaj = (j - 1 + 1)*deta;
			for (unsigned int i = istart; i < n_old; i++)
			{
				if (etaj < weights(i))
				{
					const double delta = (etaj - weights(i - 1)) / (weights(i) - weights(i - 1));
					x_new(j) = x(i - 1) + (x(i) - x(i - 1))*delta;
					break;
				}
				else
				{
					istart = i;
				}
			}
		}

		return REGRID_SUCCESS;
	}



	Adapter_Grid1D_Status Adapter_Grid1D::Refine(const Eigen::VectorXd& x, const std::vector<Eigen::VectorXd>& phi, Eigen::VectorXd& x_new)
	{
		unsigned int n_ = phi.size();
		unsigned int np_ = x.size();

		std::vector<double> ratio_gradient_(np_ - 1);
		std::vector<double> ratio_curvature_(np_ - 1);
		std::vector<int>	variables_gradient_(np_ - 1);
		std::vector<int>	variables_curvature_(np_ - 1);
		std::vector<int>	variables_(np_ - 1);
		std::vector<bool>	mark_(np_ - 1);

		for (unsigned int k = 0; k < np_ - 1; k++)
		{
			ratio_gradient_[k] = 0.;
			ratio_curvature_[k] = 0.;
			variables_gradient_[k] = 0;
			variables_curvature_[k] = 0;
			variables_[k] = 0;
			mark_[k] = false;
		}

		unsigned int n_significant = 0;
		for (unsigned int j = 0; j < n_; j++)
		{
			const double lower = phi[j].minCoeff();
			const double upper = phi[j].maxCoeff();
			const double range = upper - lower;
			const double magnitude = std::max(std::fabs(lower), std::fabs(upper));

			const double range_level = threshold_* std::max(1., magnitude);

			const bool significant = (range>range_level) ? true : false;

			if (significant == true)
			{
				n_significant++;

				for (unsigned int k = 0; k < np_ - 1; k++)
				{
					const double dphi = std::fabs(phi[j](k + 1) - phi[j](k));

					ratio_gradient_[k] = std::max(ratio_gradient_[k], dphi / range);
					if (dphi > coeff_grad_*range)
						variables_gradient_[k]++;
				}

				double temp = (phi[j](1) - phi[j](0)) / (x(1) - x(0));
				double temp_lower = temp;
				double temp_upper = temp;

				for (unsigned int k = 1; k < np_ - 1; k++)
				{
					temp = (phi[j](k + 1) - phi[j](k)) / (x(k + 1) - x(k));
					temp_lower = std::min(temp_lower, temp);
					temp_upper = std::max(temp_upper, temp);
				}
				double temp_range = temp_upper - temp_lower;

				for (unsigned int k = 1; k < np_ - 1; k++)
				{
					const double left = (phi[j](k) - phi[j](k - 1)) / (x(k) - x(k - 1));
					const double right = (phi[j](k + 1) - phi[j](k)) / (x(k + 1) - x(k));
					const double differ = std::fabs(left - right);

					if (temp_range > 0.)
						ratio_curvature_[k] = std::max(ratio_curvature_[k], differ / temp_range);
					if (differ > coeff_curv_*temp_range)
						variables_curvature_[k]++;
				}

			}
		}

		variables_[0] = variables_gradient_[0];
		for (unsigned int k = 0; k < np_ - 2; k++)
			variables_[k] = variables_gradient_[k] + variables_curvature_[k] + variables_curvature_[k + 1];
		variables_[np_ - 2] = variables_gradient_[np_ - 2] + variables_curvature_[np_ - 2];

		unsigned int ideal = 0;
		for (unsigned int k = 0; k < np_ - 1; k++)
		if (variables_[k] > 0)	ideal++;

		if (ideal == 0)
		{
			std::cout << "No need to refine the grid: the refinement conditions are satisfied!" << std::endl;
			std::cout << " * max grad. ratio (cur./id.): " << std::fixed << std::setprecision(4) << *std::max_element(ratio_gradient_.begin(), ratio_gradient_.end()) << "/" << coeff_grad_ << std::endl;
			std::cout << " * max curv. ratio (cur./id.): " << std::fixed << std::setprecision(4) << *std::max_element(ratio_curvature_.begin(), ratio_curvature_.end()) << "/" << coeff_curv_ << std::endl;

			x_new = x;
			return NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED;
		}
		else
		{
			// SELECT THE INTERVALS TO HALVE.
			// HALVE INTERVALS ON WHICH COMPONENTS OR DERIVATIVES VARY TOO
			// GREATLY. IF IT IS NOT POSSIBLE TO HALVE ALL SUCH INTERVALS, THEN
			// FIND A LOWER BOUND FOR THE NUMBER OF VARYING COMPONENTS PER
			// INTERVAL SO THAT ALL INTERVALS WITH STRICTLY MORE CAN BE HALVED.

			// FIND THE BOUND

			bool found_bound = false;
			int bound = 0;
			unsigned int more = std::max( static_cast<unsigned int>(0), std::min(max_points_ - np_, max_points_to_be_added_));

			unsigned int equal;
			unsigned int strict;

			do
			{
				int higher = 0;
				equal = 0;
				strict = 0;

				for (unsigned int k = 0; k < np_ - 1; k++)
				{
					if (variables_[k] == bound)
						equal++;
					else if (variables_[k] > bound)
					{
						strict++;
						if (strict == 1)
							higher = variables_[k];
						else
							higher = std::min(higher, variables_[k]);
					}
				}

				if (strict > more)
					bound = higher;
				else
					found_bound = true;

			} while (found_bound == false);

			//  DETERMINE HOW MANY INTERVALS TO HALVE OF THOSE WHOSE NUMBER OF
			//  VARIATIONS EXACTLY EQUAL THE BOUND.

			if (bound > 0)
				equal = std::min(equal, more - strict);
			else
				equal = 0;

			//  MARK THE INTERVALS TO HALVE

			unsigned int count = 0;
			for (unsigned int k = 0; k < np_ - 1; k++)
			{
				if (bound < variables_[k])
					mark_[k] = true;
				else if (bound == variables_[k] && count < equal)
				{
					count++;
					mark_[k] = true;
				}
				else
					mark_[k] = false;
			}

			// HALVE THE INTERVALS, IF ANY

			if (equal + strict == 0)
			{
				x_new = x;
				//OpenSMOKE::FatalErrorMessage("Adapter_Grid1D Class: equal + strict = 0");
				return NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED;
			}
			else
			{
				std::vector<double> x_temp(np_);
				for (unsigned int k = 0; k < np_; k++)
					x_temp[k] = x(k);

				for (unsigned int k = 0; k < np_ - 1; k++)
				if (mark_[k] == true)
					x_temp.push_back(0.50*(x(k) + x(k + 1)));

				std::sort(x_temp.begin(), x_temp.end());

				x_new.resize(x_temp.size());
				for (unsigned int k = 0; k < x_temp.size(); k++)
					x_new(k) = x_temp[k];
			}

			// Refine grid
			std::cout << "Grid refinement" << std::endl;
			std::cout << " * old number of points:       " << np_ << std::endl;
			std::cout << " * new number of points:       " << x_new.size() << std::endl;
			std::cout << " * significant variables:      " << n_significant << "/" << n_ << std::endl;
			std::cout << " * added (curr./id.):          " << (x_new.size() - np_) << "/" << ideal << std::endl;
			std::cout << " * max grad. ratio (cur./id.): " << std::fixed << std::setprecision(4) << *std::max_element(ratio_gradient_.begin(), ratio_gradient_.end()) << "/" << coeff_grad_ << std::endl;
			std::cout << " * max curv. ratio (cur./id.): " << std::fixed << std::setprecision(4) << *std::max_element(ratio_curvature_.begin(), ratio_curvature_.end()) << "/" << coeff_curv_ << std::endl;

			if (x_new.size() < int(max_points_))
			{
				return NEW_POINTS_ARE_NEEDED;
			}
			else
			{
				return MAXIMUM_NUMBER_POINTS;
			}
		}
	}
}
