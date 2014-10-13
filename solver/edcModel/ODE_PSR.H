
#ifndef ODE_PSR_H
#define ODE_PSR_H

namespace OpenSMOKE
{
	template<typename type>
	class ThermodynamicsMap_CHEMKIN;

	template<typename type>
	class KineticsMap_CHEMKIN;
}

class ODE_PSR
{
public:

	ODE_PSR(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& 		thermodynamicsMapXML, 
			OpenSMOKE::KineticsMap_CHEMKIN<double>& 		kineticsMapXML);

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

	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };
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

	bool checkMassFractions_;
};

#include "ODE_PSR.C"

#endif // ODE_PSR_H

