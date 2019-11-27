/*-----------------------------------------------------------------------*\
|                  _       _____ __  __  ____  _  ________                |
|                 | |     / ____|  \/  |/ __ \| |/ /  ____|               |
|          ___  __| | ___| (___ | \  / | |  | | ' /| |__                  |
|         / _ \/ _` |/ __|\___ \| |\/| | |  | |  < |  __|                 |
|        |  __/ (_| | (__ ____) | |  | | |__| | . \| |____                |
|         \___|\__,_|\___|_____/|_|  |_|\____/|_|\_\______|               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, M.R. Malik, Z. Li, A. Parente                      |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|   Contacts: Mohammad Rafi Malik, Zhiyi Li, Alessandro Parente           |
|   Aero-Thermo-Mechanical Department                                     |
|   Université Libre de Bruxelles                                         |
|   Avenue F. D. Roosevelt 50, 1050 Bruxelles (Belgium)                   |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of edcSMOKE solver.                                 |
|                                                                         |
|	License                                                           |
|                                                                         |
|   Copyright(C) 2017-2014 A. Cuoci, A. Parente                           |
|   edcSMOKE is free software: you can redistribute it and/or modify      |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   edcSMOKE is distributed in the hope that it will be useful,           |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with edcSMOKE. If not, see <http://www.gnu.org/licenses/>.      |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "math/OpenSMOKEVector.h"
#include "math/external-ode-solvers/OpenSMOKE_OdeSystemObject.h"
#include "ODESystemVirtualClassWithOpenSMOKEVectors.H"

#include "ODE_PFR_SPARC.H"

#ifndef ODE_PFR_SPARC_Interface_H
#define ODE_PFR_SPARC_Interface_H

namespace OpenSMOKE
{
	class ODE_PFR_SPARC_OpenSMOKE : public ODESystemVirtualClassWithOpenSMOKEVectors
	{
	public:

		ODE_PFR_SPARC_OpenSMOKE() {};

		void SetReactor(ODE_PFR_SPARC* reactor) { reactor_ = reactor; }

		virtual void GetEquations(const OpenSMOKE::OpenSMOKEVectorDouble& y, const double t, OpenSMOKE::OpenSMOKEVectorDouble& dy)
		{
			reactor_->Equations(t, y, dy);
		}
		virtual void PrintResults(const OpenSMOKE::OpenSMOKEVectorDouble& y, const double t)
		{
			reactor_->Print(t, y);
		}

	private:

		ODE_PFR_SPARC* reactor_;
	};
}

#endif // ODE_PFR_SPARC_Interface_H

