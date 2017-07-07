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

{	 
	tmp<volScalarField> tepsilon(turbulence->epsilon());
        const volScalarField& epsilon = tepsilon();

        tmp<volScalarField> tk(turbulence->k());
        const volScalarField& k = tk();

        tmp<volScalarField> tmu(thermo.mu());
        const volScalarField& mu = tmu();

	tmp<volScalarField> trho(thermo.rho());
        const volScalarField& rho = trho();

	if (combustionModel == "ED" || combustionModel == "EDFR")
	{
		if (iLES == true)	// LES
		{
			tauStar = Ctau * sqrt( (mu+musmall_)/rho / (epsilon+epsilonsmall_) );
			one_over_tauStar = 1./tauStar;
		}
		else			// RANS
		{
			tauStar = k/(epsilon+epsilonsmall_);
			one_over_tauStar = 1./tauStar;
		}
		
		Info<< " * min/max(epsilon)   = " << min(epsilon).value() 	<< ", " << max(epsilon).value() << endl;
		Info<< " * min/max(k)         = " << min(k).value()		<< ", " << max(k).value() 	<< endl;
		Info<< " * min/max(tauStar)   = " << min(tauStar).value() 	<< ", " << max(tauStar).value() << endl;
	}

	if (combustionModel == "EDC")
	{
		if (edc_version == "v2016")
		{
			// Calculation of characteristic times
			{
				const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

				CharacteristicChemicalTimes* chemicalTimes = new CharacteristicChemicalTimes(*thermodynamicsMapXML, *kineticsMapXML);
				chemicalTimes->SetThresholdChemicalTime(threshold_chemical_time);
				chemicalTimes->SetType(characteristic_chemical_times_type);
				chemicalTimes->SetThresholdTemperatureChemicalTime(threshold_temperature_chemical_time);

				const unsigned int ns = thermodynamicsMapXML->NumberOfSpecies();
				OpenSMOKE::OpenSMOKEVectorDouble yStar(ns);

				const scalarField& TCells = T.internalField();
				const scalarField& pCells = p.internalField(); 
				#if OPENFOAM_VERSION == 4
				      scalarField& tauChemCells = tauChem.ref(); 
				      scalarField& DaCells = Da.ref(); 
				      scalarField& ReTCells = ReT.ref();
				#else
				      scalarField& tauChemCells = tauChem.internalField(); 
				      scalarField& DaCells = Da.internalField(); 
				      scalarField& ReTCells = ReT.internalField();
				#endif 

				forAll(tauStar, i)
				{
					const double Ti = TCells[i];
					const double pi = pCells[i];

					for(unsigned int j=0;j<ns;j++)
						yStar[j+1] = YStar_[j].internalField()[i];

					// Chemical time
					const double tauLocal = chemicalTimes->CalculateCharacteristicChemicalTime(Ti,pi,yStar);
					
					// Coefficients
					const scalar nu = (mu[i]+SMALL)/(rho[i] + SMALL);
					const scalar DaLocal = max(min(Foam::sqrt(nu/(epsilon[i] + SMALL))/tauLocal, maximum_Da), minimum_Da);
					const scalar ReTLocal = CRe_coefficient*sqr(k[i])/(nu*epsilon[i] + SMALL);

					const scalar C1 = Foam::sqrt(Ctau_coefficient/3.);
		        		const scalar CtauLocal = max( min(C1/(DaLocal*Foam::sqrt(ReTLocal + 1)), 5), 0.40825);
					const scalar CgammaLocal = max( min( Cgamma_coefficient *Foam::sqrt(DaLocal*(ReTLocal+1.)), 2.1377), 0.5);
					
					// Fill fields
					tauChemCells[i] = tauLocal;
					DaCells[i] = DaLocal;
					ReTCells[i] = ReTLocal;
				
					// Fine structure properties
					tauStar[i] = CtauLocal * Foam::sqrt( nu / (epsilon[i]+SMALL) );
					gammaL[i]  = CgammaLocal * pow025(nu*epsilon[i]/(sqr(k[i]) + SMALL));
				}

				const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

				// Correct boundary conditions
				tauChem.correctBoundaryConditions();
				Da.correctBoundaryConditions();
				ReT.correctBoundaryConditions();

				Info<< " * CPU time [s]   = " << tEnd - tStart		<< endl;
			}

			// Fine structures: characteristic time [s]
			mDotStar = 1./tauStar;

			// Fine structure properties: volumetric fraction
			gammaL = min(maxGammaL, gammaL);
			gammaStarNum = pow(gammaL,edc_exp_1);
			gammaStarDen = pow(gammaL,edc_exp_2);
		}
		else
		{

			// Fine structures: characteristic time [s]
			tauStar = Ctau * sqrt( (mu+musmall_)/rho / (epsilon+epsilonsmall_) );
			mDotStar = 1./tauStar;

			// Fine structure properties: volumetric fraction
			gammaL = Ccsi * pow025( (mu+musmall_)/rho * epsilon / sqr(k+ksmall_));	
			gammaL = min(maxGammaL, gammaL);
			gammaStarNum = pow(gammaL,edc_exp_1);
			gammaStarDen = pow(gammaL,edc_exp_2);
		}

		// EDC coefficient
		edc_coefficient = rho * mDotStar * gammaStarNum / (1.-gammaStarDen);
		tfStar = (1.-gammaStarDen)/(mDotStar);		

		// Info
		Info<< " * min/max(epsilon)   = " << min(epsilon).value() 	<< ", " << max(epsilon).value() 	<< endl;
		Info<< " * min/max(k)         = " << min(k).value()		<< ", " << max(k).value()		<< endl;
		Info<< " * min/max(gammaL)    = " << min(gammaL).value()	<< ", " << max(gammaL).value()		<< endl;
		Info<< " * min/max(gammaStar) = " << min(gammaStarDen)		<< ", " << max(gammaStarDen) 		<< endl;
		Info<< " * min/max(tauStar)   = " << min(tauStar).value()	<< ", " << max(tauStar).value() 	<< endl;
		Info<< " * min/max(tfStar)    = " << min(tfStar)	 	<< ", " << max(tfStar)	 		<< endl;

		if (edc_version == "v2016")
		{
			Info<< " * min/max(tc)   = " << min(tauChem).value() 	<< ", " << max(tauChem).value()	<< endl;
			Info<< " * min/max(Da)   = " << min(Da).value()		<< ", " << max(Da).value()	<< endl;
			Info<< " * min/max(ReT)  = " << min(ReT).value()	<< ", " << max(ReT).value()	<< endl;
		}
	}

	if (combustionModel == "PaSR")
	{
		#if STEADYSTATE == 0
		if (iLES == true)
		{
			#if OPENFOAM_VERSION == 4
			const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
			#else
			const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("LESProperties");
			#endif

		        epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
		        tauMix =  sqrt( (1.2247*lesModel.delta()/sqrt(k+ksmall_))*sqrt(mu/rho/epsilonSGS) ); //geometrical mean of Kolmogorov scale and filter size scale

			Info<< " * min/max(epsilon(ssg)) = " << min(epsilonSGS).value() 	<< ", " << max(epsilonSGS).value() 		<< endl;
			Info<< " * min/max(epsilon(les)) = " << min(lesModel.epsilon()).value() << ", " << max(lesModel.epsilon()).value() 	<< endl;
			Info<< " * min/max(mut)          = " << min(turbulence->mut()).value() 	<< ", " << max(turbulence->mut()).value() 	<< endl;
			Info<< " * min/max(taumix)       = " << min(tauMix).value() 		<< ", " << max(tauMix).value()	 		<< endl;

		}
		else
		#endif	// RANS
		{
			if ( tauMixType == SCALE_KOLMOGOROV)
			{
				Info << " * using Kolmogorov mixing time scale..." << endl;
				tauMix =  sqrt( (mu+musmall_)/rho/(epsilon+epsilonsmall_) );	// Kolmogorov mixing scale
			}
			else if ( tauMixType == SCALE_MEAN)
			{
				Info << " * using geometrical mean of integral scale and Kolmogorov scale as mixing time scale..." << endl;
				tauMix =  sqrt( ((k+ksmall_)/(epsilon+epsilonsmall_))*sqrt((mu+musmall_)/rho/(epsilon+epsilonsmall_)) ); //geometrical mean of integral scale and Kolmogorov scale
			}
			else if ( tauMixType == SCALE_GLOBAL)
			{
				Info << " * using global mixing time scale..." << endl;
				tauMix = Cmix*(k+ksmall_)/(epsilon+epsilonsmall_);  //global scale
			}

			Info<< " * min/max(epsilon)  = " << min(epsilon).value() 		<< ", " << max(epsilon).value() << endl;
			Info<< " * min/max(mut)      = " << min(turbulence->mut()).value()	<< ", " << max(turbulence->mut()).value() << endl;
			Info<< " * min/max(taumix)   = " << min(tauMix).value() 		<< ", " << max(tauMix).value() << endl;
		}

		// Calculation of characteristic times
		{
			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

			CharacteristicChemicalTimes* chemicalTimes = new CharacteristicChemicalTimes(*thermodynamicsMapXML, *kineticsMapXML);
			chemicalTimes->SetThresholdChemicalTime(threshold_chemical_time);
			chemicalTimes->SetType(characteristic_chemical_times_type);
			chemicalTimes->SetThresholdTemperatureChemicalTime(threshold_temperature_chemical_time);

			const unsigned int ns = thermodynamicsMapXML->NumberOfSpecies();
			OpenSMOKE::OpenSMOKEVectorDouble yStar(ns);

			const scalarField& TCells = T.internalField();
			const scalarField& pCells = p.internalField(); 
			const scalarField& tauMixCells = tauMix.internalField(); 
			#if OPENFOAM_VERSION == 4
			      scalarField& tauChemCells = tauChem.ref(); 
			      scalarField& DaCells = Da.ref(); 
			      scalarField& kappaCells = kappa.ref(); 
			#else
			      scalarField& tauChemCells = tauChem.internalField(); 
			      scalarField& DaCells = Da.internalField(); 
			      scalarField& kappaCells = kappa.internalField(); 
			#endif 
			      

			forAll(tauStar, i)
			{
				const double Ti = TCells[i];
				const double pi = pCells[i];

				for(unsigned int j=0;j<ns;j++)
					yStar[j+1] = YStar_[j].internalField()[i];

				// Chemical time
				const double tauLocal = chemicalTimes->CalculateCharacteristicChemicalTime(Ti,pi,yStar);

				// Mixing time				
				const scalar tauMixLocal = tauMixCells[i];

				// PASR coefficients
				const scalar kappaLocal = tauLocal/(tauLocal+tauMixLocal);
				const scalar DaLocal = max( min(tauMixLocal/tauLocal, maximum_Da), minimum_Da);

				// Fill fields
				tauChemCells[i] = tauLocal;
				DaCells[i] = DaLocal;
				kappaCells[i] = kappaLocal;
			}

			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			// Correct boundary conditions
			tauChem.correctBoundaryConditions();
			Da.correctBoundaryConditions();
			kappa.correctBoundaryConditions();

			Info<< " * min/max(tc)    = " << min(tauChem).value() 	<< ", " << max(tauChem).value() 	<< endl;
			Info<< " * min/max(Da)    = " << min(Da).value()	<< ", " << max(Da).value()		<< endl;
			Info<< " * min/max(kappa) = " << min(kappa).value()	<< ", " << max(kappa).value()		<< endl;
			Info<< " * CPU time [s]   = " << tEnd - tStart		<< endl;
		}
	}
}
    

