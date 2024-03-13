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
	class GaussianLegendre_BaseClass
	{
	public:

		unsigned int n() const { return n_; }

		double z(const unsigned int j) const  { return z_(j); }
		double w(const unsigned int j) const  { return w_(j); }

		const Eigen::VectorXd& z() const { return z_; }
		const Eigen::VectorXd& w() const { return w_; }

	protected:

		Eigen::VectorXd z_;
		Eigen::VectorXd w_;
		unsigned int n_;
	};

	class GaussianLegendre_Order10 : public GaussianLegendre_BaseClass
	{
	public:

		GaussianLegendre_Order10()
		{
			n_ = 10;
			z_.resize(n_);
			w_.resize(n_);

			w_(0) = 0.2955242247147529; z_(0) = -0.1488743389816312;
			w_(1) = 0.2955242247147529;	z_(1) = 0.1488743389816312;
			w_(2) = 0.2692667193099963; z_(2) = -0.4333953941292472;
			w_(3) = 0.2692667193099963;	z_(3) = 0.4333953941292472;
			w_(4) = 0.2190863625159820; z_(4) = -0.6794095682990244;
			w_(5) = 0.2190863625159820;	z_(5) = 0.6794095682990244;
			w_(6) = 0.1494513491505806; z_(6) = -0.8650633666889845;
			w_(7) = 0.1494513491505806;	z_(7) = 0.8650633666889845;
			w_(8) = 0.0666713443086881; z_(8) = -0.9739065285171717;
			w_(9) = 0.0666713443086881; z_(9) = 0.9739065285171717;
		}
	};

	void Analytical_1DSpherical::ErrorMessage(const std::string message)
	{
		std::cout << "Grid: " << message << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(0);
	}

	void Analytical_1DSpherical::Initialize(const unsigned int NPG_physical, const unsigned int NPG)
	{
		NPG_ = NPG;
		NPG_physical_ = NPG_physical;

		ChangeDimensions(NPG_, &temp_tilde_, true);
		ChangeDimensions(NPG_, &ro_tilde_, true);
		ChangeDimensions(NPG_, &akp_tilde_, true);

		ChangeDimensions(NPG_, &g2_, true);
		ChangeDimensions(NPG_, &h_, true);


		gaussian_ = new GaussianLegendre_Order10();
	}

	void Analytical_1DSpherical::Update(const OpenSMOKE::OpenSMOKEVectorDouble& akp, const OpenSMOKE::OpenSMOKEVectorDouble& TG, const OpenSMOKE::OpenSMOKEVectorDouble& rG)
	{
		// Normalized Planck mean absorption coefficient
		for (unsigned int j = 1; j <= NPG_physical_; j++)
			akp_tilde_[j] = akp[j] * rG[1];

		// Normalize temperature
		for (unsigned int j = 1; j <= NPG_physical_; j++)
			temp_tilde_[j] = TG[j] / TG[1];

		// Normalize radius
		for (unsigned int j = 1; j <= NPG_physical_; j++)
			ro_tilde_[j] = rG[j] / rG[1];

		// Fill 
		{
			for (unsigned int j = NPG_physical_ + 1; j <= NPG_; j++)
				akp_tilde_[j] = akp_tilde_[NPG_physical_];

			for (unsigned int j = NPG_physical_ + 1; j <= NPG_; j++)
				temp_tilde_[j] = temp_tilde_[NPG_physical_];

			const double last_delta_ro_tilde = ro_tilde_[NPG_physical_] - ro_tilde_[NPG_physical_ - 1];
			for (unsigned int j = NPG_physical_ + 1; j <= NPG_; j++)
				ro_tilde_[j] = ro_tilde_[j - 1] + last_delta_ro_tilde;
		}
	}

	void Analytical_1DSpherical::RadiationToInterface()
	{
		g3_ = G3();
	}


	double phi_function(const double r, const double rho, const double z)
	{
		return (r + rho - z)*(r + rho + z)*(r - rho + z)*(r - rho - z) / (4.*z*z);
	}

	double first_indefinite_integral(const double phi, const double u)
	{
		const double delta = u*u - phi;
		if (delta < 0.)
			OpenSMOKE::FatalErrorMessage("Analytical_1DSpherical: negative delta in first integral");

		const double sqrt_delta = std::sqrt(delta);
		return u*sqrt_delta / 2. - phi*std::log(u + sqrt_delta) / 2.;
	}

	double first_integral(const double phi, const double uA, const double uB)
	{
		return (first_indefinite_integral(phi, uB) - first_indefinite_integral(phi, uA));
	}

	void Analytical_1DSpherical::RadiationFromFlame()
	{
		for (unsigned int i = 1; i <= NPG_; i++)
			g2_[i] = G2(ro_tilde_[i]);
	}


	double second_integrand_function(const double ktilde, const double r, const double rho, const double z)
	{
		const double uA = (std::fabs(r*r - rho*rho) - z*z) / (2.*z);
		const double uB = (std::fabs(r*r - rho*rho) + z*z) / (2.*z);
		const double phi = phi_function(r, rho, z);

		const double f = first_integral(phi, uA, uB);
		return exp(-ktilde*f) / z;
	}

	double second_integrand_function_interface(const double ktilde, const double rho, const double z)
	{
		const double delta = z*z - rho*rho;
		const double uA = std::sqrt(delta + 1);
		const double uB = z;

		const double f = first_integral(delta, uA, uB);
		return exp(-ktilde*f);
	}

	double Analytical_1DSpherical::Kappa(const double r, const double rho, const double ktilde)
	{
		const double a = std::fabs(r - rho);
		const double b = std::sqrt(r*r - 1.) + std::sqrt(rho*rho - 1.);

		double sum = 0.;
		for (unsigned int i = 0; i < gaussian_->n(); i++)
		{
			const double x = 0.50*(a + b) + 0.50*(b - a)*gaussian_->z(i);
			sum += second_integrand_function(ktilde, r, rho, x)*gaussian_->w(i);
		}
		sum *= (b - a) / 2.;

		return sum;
	}

	double Analytical_1DSpherical::Psi(const double rho, const double ktilde)
	{
		const double a = std::sqrt(rho*rho - 1.);
		const double b = rho;

		double sum = 0.;
		for (unsigned int i = 0; i < gaussian_->n(); i++)
		{
			const double x = 0.50*(a + b) + 0.50*(b - a)*gaussian_->z(i);
			sum += second_integrand_function_interface(ktilde, rho, x)*gaussian_->w(i);
		}
		sum *= (b - a) / 2.;

		return sum;
	}

	double Analytical_1DSpherical::g2_integrand_function(const double rtilde, const unsigned int j)
	{
		return Kappa(rtilde, ro_tilde_[j], akp_tilde_[j])*ro_tilde_[j] * std::pow(temp_tilde_[j], 4.)*akp_tilde_[j];
	}

	double Analytical_1DSpherical::g2_integrand_function_interface(const unsigned int j)
	{
		return Psi(ro_tilde_[j], akp_tilde_[j])*ro_tilde_[j] * std::pow(temp_tilde_[j], 4.)*akp_tilde_[j];
	}

	double Analytical_1DSpherical::G2(const double rtilde)
	{
		for (unsigned int j = 1; j <= NPG_; j++)
			h_[j] = g2_integrand_function(rtilde, j);

		double sum = 0.;
		for (unsigned int j = 1; j < NPG_; j++)
			sum += 0.50*(ro_tilde_[j + 1] - ro_tilde_[j])*(h_[j] + h_[j + 1]);
		sum /= rtilde;

		return sum;
	}

	double Analytical_1DSpherical::G3()
	{
		for (unsigned int j = 1; j <= NPG_; j++)
			h_[j] = g2_integrand_function_interface(j);

		double sum = 0.;
		for (unsigned int j = 1; j < NPG_; j++)
			sum += 0.50*(ro_tilde_[j + 1] - ro_tilde_[j])*(h_[j] + h_[j + 1]);

		return sum;
	}
}