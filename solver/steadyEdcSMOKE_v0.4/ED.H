if (homogeneousReactions == true)
{
	const int unsigned ns = Y.size();

	//OpenSMOKE::OpenSMOKEVectorDouble yMean(ns);
	//OpenSMOKE::OpenSMOKEVectorDouble xMean(ns);
	//OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);
	
	//double hMean, MWMean;
	//const scalarField& pCells = thermo.p().internalField();
	const scalarField& rhoCells = rho.internalField();
	const scalarField& TCells = thermo.T().internalField();
	const scalarField& kCells = turbulence->k()->internalField();
	const scalarField& epsilonCells = turbulence->epsilon()->internalField();
	const scalarField& ShCells = Sh_.internalField();

	const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	forAll(TCells, celli)
	{	
		const double smallEpsilon = 1e-12;
		const double rho_epsilon_over_k = rhoCells[celli]*kCells[celli]/(epsilonCells[celli]+smallEpsilon);	
	
		for(unsigned int k=0;k<kineticsMapXML->NumberOfReactions();k++)
		{
			double sumOmegaProducts = 0.;
			for(unsigned int j=0;j<ed_product_index_.size();j++)
				sumOmegaProducts += YStar_[ed_product_index_[k][j]].internalField()[celli];

			double r = constA*constB*sumOmegaProducts/ed_denominator_products_[k];
			for(unsigned int j=0;j<ed_reactant_index_.size();j++)
				r = std::min(r, constA*YStar_[ed_reactant_index_[k][j]].internalField()[celli])/ed_denominator_reactants_[k][j];
		}
	
	//	for(unsigned int i=0;i<ns;i++)
	//		R_[i].internalField()[celli] = edc_coefficient[celli]*(YStar_[i].internalField()[celli]-Y[i].internalField()[celli]);

	//	double HeatRelease = 0.;
	//	for(unsigned int i=0;i<ns;i++)
	//		HeatRelease += -h_species[i+1]/thermodynamicsMapXML->MW()[i+1] * R_[i].internalField()[celli];
	//	ShCells[celli] = HeatRelease;
	}
	const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
	Info << tEnd - tStart << " s " << "(" << (tEnd - tStart)/double(TCells.size())*1e3 << " ms per cell)" << endl;
}
