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

{
	const int unsigned ns = Y.size();

	Eigen::VectorXd omega(ns);
	Eigen::VectorXd R(ns);
	Eigen::VectorXd r(kineticsMapXML->NumberOfReactions());
	OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);

	// These are used only for the EDFR model
	OpenSMOKE::OpenSMOKEVectorDouble omega_(ns);
	OpenSMOKE::OpenSMOKEVectorDouble c_(ns);
	
	// OpenFOAM scalar fields
	const scalarField& pCells = thermo.p().internalField();
	const scalarField& rhoCells = rho.internalField();
	const scalarField& TCells = thermo.T().internalField();
	      scalarField& ShCells = Sh_.internalField();

	const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	if (iLES == true)	// LES
	{
		// Nothing to do because uTau_ was already calculated
		//uTau_ = 1./U.db().objectRegistry::lookupObject<volScalarField>("tauMix").internalField();

		/*
		const double ksmall = 1e-5;
		const volScalarField& epsConst = turbulence->epsilon();
		const volScalarField& kConst   = turbulence->k();
		const scalarField& kCells      = kConst.internalField();
		const scalarField& epsCells    = epsConst.internalField();
		      scalarField& uTauCells   = uTau_.internalField();

		forAll(TCells, celli)
			uTauCells[celli] = epsCells[celli]/(kCells[celli]+ksmall);
		*/
	}
	else	// RANS
	{
		uTau_ = turbulence->epsilon()/(turbulence->k()+ksmall_);
	}

	forAll(TCells, celli)
	{	
		// Calculating constant coefficient
		const double rho_epsilon_over_k = rhoCells[celli]*uTau_.internalField()[celli];	

		// Vector containing mass fractions and mole formation rates
		for(unsigned int i=0;i<ns;i++)
			omega(i) = Y[i].internalField()[celli];
		R.setZero();

		// Set spark
		if( 	(spark == true) && 
			(runTime.value() >= ignitionTime)&&(runTime.value() <= (ignitionTime+ignitionDuration)) ) 
		{
			for(unsigned int i=0;i<ns;i++)
				omega(i) = (omega(i) + ignitionAdditionalMassFraction) / (1.+ns*ignitionAdditionalMassFraction);
		}

		// Calculating the reaction rates
		if (combustionModel == "ED")
		{
			for(unsigned int k=0;k<kineticsMapXML->NumberOfReactions();k++)
			{
				double sumOmegaProducts = 0.;
				for(unsigned int j=0;j<ed_product_index_[k].size();j++)
					sumOmegaProducts += omega(ed_product_index_[k][j]);

				r[k] = constA*constB*sumOmegaProducts/ed_denominator_products_[k];
				for(unsigned int j=0;j<ed_reactant_index_[k].size();j++)
					r[k] = std::min(r[k], constA*omega(ed_reactant_index_[k][j])/ed_denominator_reactants_[k][j]);
				r[k] *= rho_epsilon_over_k;
			}
		}
		else if (combustionModel == "EDFR")
		{
			double MW;
			for(unsigned int i=0;i<ns;i++)
				omega_[i+1] = omega(i);
			
			// Calculating the concentrations
			double cTot = pCells[celli]/(PhysicalConstants::R_J_kmol * TCells[celli]);
			thermodynamicsMapXML->MoleFractions_From_MassFractions(c_.GetHandle(), MW, omega_.GetHandle());
			Product(cTot, &c_);

			// Calculates kinetics
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			kineticsMapXML->SetTemperature(TCells[celli]);
			kineticsMapXML->SetPressure(pCells[celli]);
			kineticsMapXML->ReactionEnthalpiesAndEntropies();
			kineticsMapXML->KineticConstants();
			kineticsMapXML->ReactionRates(c_.GetHandle());

			for(unsigned int k=0;k<kineticsMapXML->NumberOfReactions();k++)
			{
				double sumOmegaProducts = 0.;
				for(unsigned int j=0;j<ed_product_index_[k].size();j++)
					sumOmegaProducts += omega(ed_product_index_[k][j]);

				r[k] = constA*constB*sumOmegaProducts/ed_denominator_products_[k];
				for(unsigned int j=0;j<ed_reactant_index_[k].size();j++)
					r[k] = std::min(r[k], constA*omega(ed_reactant_index_[k][j])/ed_denominator_reactants_[k][j]);
				
				// Constant coefficient
				r[k] *= rho_epsilon_over_k;

				// Comparison with chemistry
				r[k] = std::min(r[k], kineticsMapXML->GiveMeReactionRates()[k]);				
			}
		}		

		
		// Calculating the formation rates [kmol/m3/s]
		for(unsigned int k=0;k<kineticsMapXML->NumberOfReactions();k++)
		{	
			for(unsigned int j=0;j<ed_reactant_index_[k].size();j++)
				R[ed_reactant_index_[k][j]] -= ed_reactant_nu_[k][j]*r[k];
			for(unsigned int j=0;j<ed_product_index_[k].size();j++)
				R[ed_product_index_[k][j]]  += ed_product_nu_[k][j]*r[k];
		}

		// Formation rates [kg/m3/s]
		for(unsigned int j=0;j<ns;j++)
			R_[j].internalField()[celli] = R[j]*thermodynamicsMapXML->MW(j);

		// Calculating the heat release
		{
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->hMolar_Species(h_species.GetHandle());

			double HeatRelease = 0.;
			for(unsigned int i=0;i<ns;i++)
				HeatRelease += -h_species[i+1]*R[i];
			ShCells[celli] = HeatRelease;
		}
	}
	const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
	Info << "CPU time for chemistry (" << combustionModel << ") : " << tEnd - tStart << " s " << "(" << (tEnd - tStart)/double(TCells.size())*1e3 << " ms per cell)" << endl;
}
