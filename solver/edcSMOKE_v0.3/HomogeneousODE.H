#include "math/OpenSMOKEVector.h"
#include "ode/OpenSMOKE_OdeSystemObject.h"
#include "math/stiff-ode-solvers/OdeSystemObject.h"
#include "math/stiff-ode-solvers/StiffOdeSolverObject_Dense.h"

#ifndef HomogeneousODE_H
#define HomogeneousODE_H

namespace OpenSMOKE
{
	template<typename type>
	class ThermodynamicsMap_CHEMKIN;

	template<typename type>
	class KineticsMap_CHEMKIN;
}


class HomogeneousODE
{
public:

	HomogeneousODE(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& 		thermodynamicsMapXML, 
				   OpenSMOKE::KineticsMap_CHEMKIN<double>& 				kineticsMapXML);

	void Set(const double T0, const double P0_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omegaMean, const double hMean, const double gammaStar, const double mDotStar, const double chi)
	{
		TStar_ 		= T0;
		P_Pa_ 		= P0_Pa;
		gammaStar_ 	= gammaStar;
		mDotStar_ 	= mDotStar;
		omegaMean_  = omegaMean;
		hMean_ 		= hMean;
		chi_ 		= chi;
	}

	unsigned int NumberOfEquations() const { return number_of_equations_; }

	int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);
	int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

private:

	unsigned int number_of_equations_;
	unsigned int number_of_gas_species_;

	double TStar_;
	double P_Pa_;
	double MWStar_;
	double cTotStar_;
	double rhoStar_;
	double mDotStar_;
	double gammaStar_;
	double cpStar_;
	double chi_;
	
	double hMean_;
	double hStar_;
	double hSurr_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& 	thermodynamicsMapXML_;
	OpenSMOKE::KineticsMap_CHEMKIN<double>& 		kineticsMapXML_;

	OpenSMOKE::OpenSMOKEVectorDouble omegaMean_;
	OpenSMOKE::OpenSMOKEVectorDouble omegaStar_;
	OpenSMOKE::OpenSMOKEVectorDouble omegaSurr_;
	OpenSMOKE::OpenSMOKEVectorDouble xStar_;
	OpenSMOKE::OpenSMOKEVectorDouble cStar_;
	OpenSMOKE::OpenSMOKEVectorDouble RStar_;
};


class ODESystem_OpenSMOKE_Homogeneous : public OpenSMOKE::UserOdeSystemObject 
{
public:

	ODESystem_OpenSMOKE_Homogeneous( HomogeneousODE& cell) :
	cell_(cell)
	{
		ChangeDimensions(cell_.NumberOfEquations(), &y_, true);
		ChangeDimensions(cell_.NumberOfEquations(), &dy_, false);
	}
	
	virtual void GetSystemFunctions(const Eigen::VectorXd &Y, const double t, Eigen::VectorXd &DY)
	{
		y_.CopyFrom(Y.data());
		cell_.Equations(t, y_, dy_);
		dy_.CopyTo(DY.data());
	}

	void MyPrint(const Eigen::VectorXd &Y, const double t)
	{
		y_.CopyFrom(Y.data());
		cell_.Print(t, y_);
	}

private:

	HomogeneousODE& cell_;
	OpenSMOKE::OpenSMOKEVectorDouble  y_;
	OpenSMOKE::OpenSMOKEVectorDouble dy_;
};

#include "HomogeneousODE.C"

#endif // HomogeneousODE_H

