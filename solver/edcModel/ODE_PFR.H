
#ifndef ODE_PFR_H
#define ODE_PFR_H

namespace OpenSMOKE
{
	template<typename type>
	class ThermodynamicsMap_CHEMKIN;

	template<typename type>
	class KineticsMap_CHEMKIN;
}

class ODE_PFR
{
public:

	ODE_PFR(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& 		thermodynamicsMapXML, 
			OpenSMOKE::KineticsMap_CHEMKIN<double>& 		kineticsMapXML);

	void Set(const double P0_Pa, const double gammaStar, const double mDotStar, const double chi)
	{
		P_Pa_ 		= P0_Pa;
		gammaStar_ 	= gammaStar;
		mDotStar_ 	= mDotStar;
		chi_ 		= chi;
	}

	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };
	unsigned int NumberOfEquations() const { return number_of_equations_; }

	int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);
	int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

private:

	unsigned int number_of_equations_;
	unsigned int number_of_gas_species_;

	double P_Pa_;
	double MWStar_;
	double cTotStar_;
	double rhoStar_;
	double mDotStar_;
	double gammaStar_;
	double cpStar_;
	double chi_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& 	thermodynamicsMapXML_;
	OpenSMOKE::KineticsMap_CHEMKIN<double>& 	kineticsMapXML_;

	OpenSMOKE::OpenSMOKEVectorDouble omegaStar_;
	OpenSMOKE::OpenSMOKEVectorDouble xStar_;
	OpenSMOKE::OpenSMOKEVectorDouble cStar_;
	OpenSMOKE::OpenSMOKEVectorDouble RStar_;

	bool checkMassFractions_;
};

#include "ODE_PFR.C"

#endif // ODE_PFR_H

