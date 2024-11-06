// Default name of main dictionary
const Foam::word Foam::ChemistryLinearModel::dictName("ChemistryLinearModel");

// Constructor
Foam::ChemistryLinearModel::ChemistryLinearModel(const Foam::fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap)
:
	IOdictionary
	(
		IOobject
		(
			dictName,
			mesh.time().constant(),
			mesh,
			IOobject::MUST_READ_IF_MODIFIED,
			IOobject::NO_WRITE
		)
	),

	thermodynamicsMap_(thermodynamicsMap), 

	kineticsMap_(kineticsMap),

	sparseJacobian_(lookupOrDefault<Switch>("sparseJacobian", true)),

	jacobianUpdate_(lookupOrDefault<label>("jacobianUpdate", 1)),

	propertiesUpdate_(lookupOrDefault<label>("propertiesUpdate", 1)),

	implicitSourceTerm_(lookupOrDefault<Switch>("implicitSourceTerm", true))
{
	NC_ = thermodynamicsMap_.NumberOfSpecies();
	NE_ = NC_+1;
	
	omega_.resize(NC_);	omega_.setZero();
	x_.resize(NC_);		x_.setZero();
	R_.resize(NC_);		R_.setZero();
	c_.resize(NC_);		c_.setZero();

	y_plus_.resize(NE_);		y_plus_.setZero();
	dy_plus_.resize(NE_);		dy_plus_.setZero();
	dy_original_.resize(NE_);	dy_original_.setZero();
	Jdiagonal_.resize(NE_);		Jdiagonal_.setZero();

	word order_policy = lookupOrDefault<word>("orderSpecies", "constant");
	if (order_policy == "constant")		species_order_policy_ = SPECIES_ORDER_POLICY_CONSTANT;
	else if (order_policy == "sweep")	species_order_policy_ = SPECIES_ORDER_POLICY_SWEEP;
	else if (order_policy == "random")	species_order_policy_ = SPECIES_ORDER_POLICY_RANDOM_SHUFFLE;
	else if (order_policy == "rotate")	species_order_policy_ = SPECIES_ORDER_POLICY_ROTATE;
	else FatalError << "Wrong orderSpecies option. Allowed options: constant | sweep | random | rotate" << ::Foam::exit(FatalError);

	species_order_.resize(NC_);
	for(int i=0;i<NC_;i++)
	{
		species_order_[i] = i;
	}

	jacobianCounter_ = jacobianUpdate_;	
	propertiesCounter_ = propertiesUpdate_;
}

Foam::ChemistryLinearModel::~ChemistryLinearModel()
{}


void Foam::ChemistryLinearModel::reactionSourceTerms( const Eigen::VectorXd& y, const double P0, Eigen::VectorXd& S )
{
	for(unsigned int i=0;i<NC_;++i)
		omega_(i) = max(y(i), 0.);
	const double T = y(NC_);

	// Calculates the pressure and the concentrations of species
	double MW_ = 0.;
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_.data(), MW_, omega_.data());
	const double cTot_ = P0/PhysicalConstants::R_J_kmol/T;
    	const double rho_ = cTot_*MW_;
	c_ = cTot_*x_;

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T);
	thermodynamicsMap_.SetPressure(P0);

	// Calculates kinetics
	#if CHEMKINBASED == 1
	kineticsMap_.SetTemperature(T);
	kineticsMap_.SetPressure(P0);	
	ChemkinFormationRates(NC_, T, P0, omega_.data(), R_.data());
	#else
	// Calculates kinetics
	kineticsMap_.SetTemperature(T);
	kineticsMap_.SetPressure(P0);
	kineticsMap_.KineticConstants();
	kineticsMap_.ReactionRates(c_.data());
	kineticsMap_.FormationRates(R_.data());
	#endif

	// Species equations source term
	for (unsigned int i=0;i<NC_;++i)	
		S(i) = R_(i)*thermodynamicsMap_.MW(i);

	// Temperature equation source term
    	{
		const double CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.data());
		const double CpMixMass_ = CpMixMolar / MW_;
		const double QR_ = kineticsMap_.HeatRelease(R_.data());
		S(NC_)  = QR_;
	}
}

void Foam::ChemistryLinearModel::reactionSourceTerms( 	const Eigen::VectorXd& y, const double P0, Eigen::VectorXd& S,
							const std::vector<double>& cc, const std::vector<double>& cc_falloff, const std::vector<double>& cc_reverse)
{
	for(unsigned int i=0;i<NC_;++i)
		omega_(i) = max(y(i), 0.);
	const double T = y(NC_);

	// Calculates the pressure and the concentrations of species
	double MW_ = 0.;
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_.data(), MW_, omega_.data());
	const double cTot_ = P0/PhysicalConstants::R_J_kmol/T;
    	const double rho_ = cTot_*MW_;
	c_ = cTot_*x_;

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T);
	thermodynamicsMap_.SetPressure(P0);

	// Calculates kinetics
	kineticsMap_.SetTemperature(T);
	kineticsMap_.SetPressure(P0);
	kineticsMap_.KineticConstants(cc, cc_falloff, cc_reverse);
	kineticsMap_.ReactionRates(c_.data());
	kineticsMap_.FormationRates(R_.data());

	// Species equations source term
	for (unsigned int i=0;i<NC_;++i)	
		S(i) = R_(i)*thermodynamicsMap_.MW(i);

	// Temperature equation source term
    	{
		const double CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.data());
		const double CpMixMass_ = CpMixMolar / MW_;
		const double QR_ = kineticsMap_.HeatRelease(R_.data());
		S(NC_)  = QR_;
	}
}

void Foam::ChemistryLinearModel::reactionJacobian( const Eigen::VectorXd& y, const double P0, Eigen::VectorXd &J)
{
	// Calculated as suggested by Guido Buzzi-Ferraris (private communication)
	const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
	const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
	const double BETA = 1.e3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;
	const double TOLR = 1.e-7;
	const double TOLA = 1.e-12;

	// Start calculations
	y_plus_ = y;

	// Call equations
	reactionSourceTerms(y, P0, dy_original_);     

	// Derivatives with respect to y(kd)
	for(int kd=0;kd<NE_;kd++)
	{
		const double hf = 1.e0;
		const double error_weight = 1./(TOLA+TOLR*std::fabs(y(kd)));
		      double hJ = ETA2 * std::fabs(std::max(y(kd), 1./error_weight));
		const double hJf = hf/error_weight;
		hJ = std::max(hJ, hJf);
		hJ = std::max(hJ, ZERO_DER);

		// This is what is suggested by Guido Buzzi-Ferraris (private communication)
		const double dy = std::min(hJ, 1.e-3 + 1e-3*std::fabs(y(kd)));
		const double udy = 1. / dy;
		y_plus_(kd) += dy;

		reactionSourceTerms(y_plus_, P0, dy_plus_);

		J(kd) = (dy_plus_(kd)-dy_original_(kd)) * udy;

		y_plus_(kd) = y(kd);
	}
}

void Foam::ChemistryLinearModel::reactionJacobian( 	const Eigen::VectorXd& y, const double P0, Eigen::VectorXd &J,
							const std::vector<double>& cc, const std::vector<double>& cc_falloff, const std::vector<double>& cc_reverse)
{
	// Calculated as suggested by Guido Buzzi-Ferraris (private communication)
	const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
	const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
	const double BETA = 1.e3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;
	const double TOLR = 1.e-7;
	const double TOLA = 1.e-12;

	// Start calculations
	y_plus_ = y;

	// Call equations
	reactionSourceTerms(y, P0, dy_original_, cc, cc_falloff, cc_reverse);     

	// Derivatives with respect to y(kd)
	for(int kd=0;kd<NE_;kd++)
	{
		const double hf = 1.e0;
		const double error_weight = 1./(TOLA+TOLR*std::fabs(y(kd)));
		      double hJ = ETA2 * std::fabs(std::max(y(kd), 1./error_weight));
		const double hJf = hf/error_weight;
		hJ = std::max(hJ, hJf);
		hJ = std::max(hJ, ZERO_DER);

		// This is what is suggested by Guido Buzzi-Ferraris (private communication)
		const double dy = std::min(hJ, 1.e-3 + 1e-3*std::fabs(y(kd)));
		const double udy = 1. / dy;
		y_plus_(kd) += dy;

		reactionSourceTerms(y_plus_, P0, dy_plus_, cc, cc_falloff, cc_reverse);

		J(kd) = (dy_plus_(kd)-dy_original_(kd)) * udy;

		y_plus_(kd) = y(kd);
	}
}

void Foam::ChemistryLinearModel::reactionJacobianSparse( const Eigen::VectorXd& y, const double P0, Eigen::VectorXd &J, const bool energyEquation, const bool speciesEquations)
{
	J.setZero();

	// Sparse Jacobian
	if (speciesEquations == true)
	{
		// Recover mass fractions
		for(unsigned int i=0;i<NC_;++i)
			omega_(i) = max(y(i), 0.);
		const double T = y(NC_);

		kineticsMap_.jacobian_sparsity_pattern_map()->SetEpsilon(1.e-13); 
		kineticsMap_.jacobian_sparsity_pattern_map()->Jacobian(omega_.data(), T, P0, Jdiagonal_);
		for(unsigned int i=0;i<NC_;++i)
			J(i) = Jdiagonal_(i) * thermodynamicsMap_.MW(i);
	}

	// Only with respect to the temperature
	if (energyEquation == true)
	{
		// Calculated as suggested by Guido Buzzi-Ferraris (private communication)
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double BETA = 1.e3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;
		const double TOLR = 1.e-7;
		const double TOLA = 1.e-12;

		// Start calculations
		y_plus_ = y;

		// Call equations
		reactionSourceTerms(y, P0, dy_original_);     

		// Derivatives with respect to the temperature
		const int kd=NE_-1;
		{
			const double hf = 1.e0;
			const double error_weight = 1./(TOLA+TOLR*std::fabs(y(kd)));
			      double hJ = ETA2 * std::fabs(std::max(y(kd), 1./error_weight));
			const double hJf = hf/error_weight;
			hJ = std::max(hJ, hJf);
			hJ = std::max(hJ, ZERO_DER);

			// This is what is suggested by Guido Buzzi-Ferraris (private communication)
			const double dy = std::min(hJ, 1.e-3 + 1e-3*std::fabs(y(kd)));
			const double udy = 1. / dy;
			y_plus_(kd) += dy;

			reactionSourceTerms(y_plus_, P0, dy_plus_);

			J(kd) = (dy_plus_(kd)-dy_original_(kd)) * udy;

			y_plus_(kd) = y(kd);
		}
	}
}

void Foam::ChemistryLinearModel::reactionJacobianSparse( 	const Eigen::VectorXd& y, const double P0, Eigen::VectorXd &J, const bool energyEquation, const bool speciesEquations,
								const std::vector<double>& cc, const std::vector<double>& cc_falloff, const std::vector<double>& cc_reverse)
{
	J.setZero();

	// Sparse Jacobian
	if (speciesEquations == true)
	{
		// Recover mass fractions
		for(unsigned int i=0;i<NC_;++i)
			omega_(i) = max(y(i), 0.);
		const double T = y(NC_);

		kineticsMap_.jacobian_sparsity_pattern_map()->SetEpsilon(1.e-13); 
		kineticsMap_.jacobian_sparsity_pattern_map()->Jacobian(omega_.data(), T, P0, Jdiagonal_);
		for(unsigned int i=0;i<NC_;++i)
			J(i) = Jdiagonal_(i) * thermodynamicsMap_.MW(i);
	}

	// Only with respect to the temperature
	if (energyEquation == true)
	{
		// Calculated as suggested by Guido Buzzi-Ferraris (private communication)
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double BETA = 1.e3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;
		const double TOLR = 1.e-7;
		const double TOLA = 1.e-12;

		// Start calculations
		y_plus_ = y;

		// Call equations
		reactionSourceTerms(y, P0, dy_original_, cc, cc_falloff, cc_reverse);     

		// Derivatives with respect to the temperature
		const int kd=NE_-1;
		{
			const double hf = 1.e0;
			const double error_weight = 1./(TOLA+TOLR*std::fabs(y(kd)));
			      double hJ = ETA2 * std::fabs(std::max(y(kd), 1./error_weight));
			const double hJf = hf/error_weight;
			hJ = std::max(hJ, hJf);
			hJ = std::max(hJ, ZERO_DER);

			// This is what is suggested by Guido Buzzi-Ferraris (private communication)
			const double dy = std::min(hJ, 1.e-3 + 1e-3*std::fabs(y(kd)));
			const double udy = 1. / dy;
			y_plus_(kd) += dy;

			reactionSourceTerms(y_plus_, P0, dy_plus_, cc, cc_falloff, cc_reverse);

			J(kd) = (dy_plus_(kd)-dy_original_(kd)) * udy;

			y_plus_(kd) = y(kd);
		}
	}
}

void Foam::ChemistryLinearModel::updateSpeciesOrderPolicy()
{
	if (species_order_policy_ == SPECIES_ORDER_POLICY_SWEEP)
		std::reverse(species_order_.begin(),species_order_.end());
	else if (species_order_policy_ == SPECIES_ORDER_POLICY_RANDOM_SHUFFLE)
		std::random_shuffle ( species_order_.begin(), species_order_.end() );
	else if (species_order_policy_ == SPECIES_ORDER_POLICY_ROTATE)
		std::rotate(species_order_.begin(),species_order_.begin()+1,species_order_.end());
}

void Foam::ChemistryLinearModel::updateJacobianCounter()
{
	if (jacobianCounter_ == jacobianUpdate_)
		jacobianCounter_ = 0;
	jacobianCounter_++;
}

void Foam::ChemistryLinearModel::updatePropertiesCounter()
{
	if (propertiesCounter_ == propertiesUpdate_)
		propertiesCounter_ = 0;
	propertiesCounter_++;
}

