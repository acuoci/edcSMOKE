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

#ifndef OpenSMOKE_Analytical_1DSpherical_H
#define	OpenSMOKE_Analytical_1DSpherical_H

#include "math/OpenSMOKEVector.h"

namespace OpenSMOKE
{
	class GaussianLegendre_Order10;

	class Analytical_1DSpherical
	{
	public:

		void Initialize(const unsigned int NPG_physical_, const unsigned int NPG);
		void Update(const OpenSMOKE::OpenSMOKEVectorDouble& akp, const OpenSMOKE::OpenSMOKEVectorDouble& TG, const OpenSMOKE::OpenSMOKEVectorDouble& rG);

		void RadiationToInterface();
		void RadiationFromFlame();

		double g3() const { return g3_; }

		const OpenSMOKE::OpenSMOKEVectorDouble& temp_tilde() const { return temp_tilde_; }
		const OpenSMOKE::OpenSMOKEVectorDouble& g2() const { return g2_; }
		const OpenSMOKE::OpenSMOKEVectorDouble& akp() const { return akp_; }


	private:

		OpenSMOKE::OpenSMOKEVectorDouble ro_tilde_;
		OpenSMOKE::OpenSMOKEVectorDouble akp_tilde_;
		OpenSMOKE::OpenSMOKEVectorDouble temp_tilde_;

		OpenSMOKE::OpenSMOKEVectorDouble akp_;
		OpenSMOKE::OpenSMOKEVectorDouble g2_;
		OpenSMOKE::OpenSMOKEVectorDouble h_;

		unsigned int NPG_;
		unsigned int NPG_physical_;

		double g3_;

		double G2(const double rtilde);
		double G3();

		double g2_integrand_function(const double rtilde, const unsigned int j);
		double g2_integrand_function_interface(const unsigned int j);

		double Kappa(const double r, const double rho, const double ktilde);
		double Psi(const double rho, const double ktilde);

		GaussianLegendre_Order10* gaussian_;

		void ErrorMessage(const std::string message);
	};
}

#include "Analytical_1DSpherical.hpp"

#endif // OpenSMOKE_Analytical_1DSpherical_H