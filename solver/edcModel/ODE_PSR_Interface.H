#include "math/OpenSMOKEVector.h"
#include "ode/OpenSMOKE_OdeSystemObject.h"

#if OPENSMOKE_USE_DVODE == 1
#include "ode/OpenSMOKE_DVODE_Interface.h"
#include "ode/OpenSMOKE_DVODE.h"
#endif

#if OPENSMOKE_USE_ODEPACK == 1
#include "ode/OpenSMOKE_DLSODE_Interface.h"
#include "ode/OpenSMOKE_DLSODE.h"
#include "ode/OpenSMOKE_DLSODA_Interface.h"
#include "ode/OpenSMOKE_DLSODA.h"
#endif

#if OPENSMOKE_USE_DASPK == 1
#include "ode/OpenSMOKE_DASPK_Interface.h"
#include "ode/OpenSMOKE_DASPK.h"
#endif

#if OPENSMOKE_USE_RADAU == 1
#include "ode/OpenSMOKE_RADAU_Interface.h"
#include "ode/OpenSMOKE_RADAU.h"
#endif

#if OPENSMOKE_USE_SUNDIALS == 1
#include "ode/OpenSMOKE_CVODE_Sundials_Interface.h"
#include "ode/OpenSMOKE_CVODE_Sundials.h"
#endif

#if OPENSMOKE_USE_MEBDF == 1
#include "ode/OpenSMOKE_MEBDF_Interface.h"
#include "ode/OpenSMOKE_MEBDF.h"
#endif

#include "ode/OpenSMOKE_OdeSystemObject.h"
#include "ODESystemVirtualClassWithOpenSMOKEVectors.H"

#include "ODE_PSR.H"

#ifndef ODE_PSR_Interface_H
#define ODE_PSR_Interface_H

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
}

#endif // ODE_PSR_Interface_H

