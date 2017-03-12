/*----------------------------------------------------------------------*\
|                  _       _____ __  __  ____  _  ________                |
|                 | |     / ____|  \/  |/ __ \| |/ /  ____|               |
|          ___  __| | ___| (___ | \  / | |  | | ' /| |__                  |
|         / _ \/ _` |/ __|\___ \| |\/| | |  | |  < |  __|                 |
|        |  __/ (_| | (__ ____) | |  | | |__| | . \| |____                |
|         \___|\__,_|\___|_____/|_|  |_|\____/|_|\_\______|               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, M.R. Malik, A. Parente                             |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|   Contacts: Mohammad Rafi Malik                                         |
|   Aero-Thermo-Mechanical Department                                     |
|   Université Libre de Bruxelles                                         |
|   Avenue F. D. Roosevelt 50, 1050 Bruxelles (Belgium)                   |
|                                                                         |
|   Contacts: Alessandro Parente                                          |
|   email: alessandro.parente@ulb.ac.be                                   |
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
|   Copyright(C) 2014 A. Cuoci, M.R. Malik, A. Parente                    |
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

#ifndef ODE_PSR_Interface_H
#define	ODE_PSR_Interface_H

#include "math/OpenSMOKEVector.h"
#include "math/external-ode-solvers/OpenSMOKE_OdeSystemObject.h"

#if OPENSMOKE_USE_DVODE == 1
#include "math/external-ode-solvers/OpenSMOKE_DVODE_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_DVODE.h"
#endif

#if OPENSMOKE_USE_ODEPACK == 1
#include "math/external-ode-solvers/OpenSMOKE_DLSODE_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_DLSODE.h"
#include "math/external-ode-solvers/OpenSMOKE_DLSODA_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_DLSODA.h"
#endif

#if OPENSMOKE_USE_DASPK == 1
#include "math/external-ode-solvers/OpenSMOKE_DASPK_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_DASPK.h"
#endif

#if OPENSMOKE_USE_RADAU == 1
#include "math/external-ode-solvers/OpenSMOKE_RADAU_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_RADAU.h"
#endif

#if OPENSMOKE_USE_SUNDIALS == 1
#include "math/external-ode-solvers/OpenSMOKE_CVODE_Sundials_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_CVODE_Sundials.h"
#endif

#if OPENSMOKE_USE_MEBDF == 1
#include "math/external-ode-solvers/OpenSMOKE_MEBDF_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_MEBDF.h"
#endif

#include "math/external-ode-solvers/OpenSMOKE_OdeSystemObject.h"
#include "ODESystemVirtualClassWithOpenSMOKEVectors.H"

#include "ODE_PSR.H"

namespace OpenSMOKE
{
	class ODE_PSR_OpenSMOKE : public ODESystemVirtualClassWithOpenSMOKEVectors
	{
	public:

		ODE_PSR_OpenSMOKE() {};

		void SetReactor(ODE_PSR* reactor) { reactor_ = reactor; }

		virtual void GetEquations(const OpenSMOKE::OpenSMOKEVectorDouble& y, const double t, OpenSMOKE::OpenSMOKEVectorDouble& dy)
		{
			reactor_->Equations(t, y, dy);
		}
		virtual void PrintResults(const OpenSMOKE::OpenSMOKEVectorDouble& y, const double t)
		{
			reactor_->Print(t, y);
		}

	private:

		ODE_PSR* reactor_;
	};

	#if OPENSMOKE_USE_DVODE == 1
	class ODE_PSR_DVODE : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_DVODE(ODE_PSR_DVODE)

		ODE_PSR* reactor_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetReactor(ODE_PSR* reactor)
		{
			reactor_ = reactor;
			ChangeDimensions(reactor_->NumberOfEquations(), &y_, true);
			ChangeDimensions(reactor_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DVODE(ODE_PSR_DVODE)
	#endif

	#if OPENSMOKE_USE_ODEPACK == 1
	class ODE_PSR_DLSODE : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_DLSODE(ODE_PSR_DLSODE)

		ODE_PSR* reactor_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetReactor(ODE_PSR* reactor)
		{
			reactor_ = reactor;
			ChangeDimensions(reactor_->NumberOfEquations(), &y_, true);
			ChangeDimensions(reactor_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DLSODE(ODE_PSR_DLSODE)

	class ODE_PSR_DLSODA : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_DLSODA(ODE_PSR_DLSODA)

		ODE_PSR* reactor_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetReactor(ODE_PSR* reactor)
		{
			reactor_ = reactor;
			ChangeDimensions(reactor_->NumberOfEquations(), &y_, true);
			ChangeDimensions(reactor_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DLSODA(ODE_PSR_DLSODA)
	#endif

	#if OPENSMOKE_USE_DASPK == 1
	class ODE_PSR_DASPK : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		ODE_PSR* reactor_;

		DEFINE_ODESOLVERINTERFACE_DASPK(ODE_PSR_DASPK)

		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetReactor(ODE_PSR* reactor)
		{
			reactor_ = reactor;
			ChangeDimensions(reactor_->NumberOfEquations(), &y_, true);
			ChangeDimensions(reactor_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DASPK(ODE_PSR_DASPK)
	#endif

	#if OPENSMOKE_USE_RADAU == 1
	class ODE_PSR_RADAU5 : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_RADAU(ODE_PSR_RADAU5)

		ODE_PSR* reactor_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetReactor(ODE_PSR* reactor)
		{
			reactor_ = reactor;
			ChangeDimensions(reactor_->NumberOfEquations(), &y_, true);
			ChangeDimensions(reactor_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_RADAU(ODE_PSR_RADAU5)
	#endif

	#if OPENSMOKE_USE_SUNDIALS == 1
	class ODE_PSR_CVODE : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_CVODE_Sundials(ODE_PSR_CVODE)

		ODE_PSR* reactor_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetReactor(ODE_PSR* reactor)
		{
			reactor_ = reactor;
			ChangeDimensions(reactor_->NumberOfEquations(), &y_, true);
			ChangeDimensions(reactor_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_CVODE_Sundials(ODE_PSR_CVODE)
	#endif

	#if OPENSMOKE_USE_MEBDF == 1
	class ODE_PSR_MEBDF : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_MEBDF(ODE_PSR_MEBDF)

		ODE_PSR* reactor_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetReactor(ODE_PSR* reactor)
		{
			reactor_ = reactor;
			ChangeDimensions(reactor_->NumberOfEquations(), &y_, true);
			ChangeDimensions(reactor_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = reactor_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_MEBDF(ODE_PSR_MEBDF)
	#endif
};

#endif // ODE_PSR_Interface_H

