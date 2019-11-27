/*-----------------------------------------------------------------------*\
|                                                                         |
|                    ╔═══╦═╗╔═╦═══╦╗╔═╦═══╗                               |
|                    ║╔═╗║║╚╝║║╔═╗║║║╔╣╔══╝                               | 
|   ╔╗╔══╦╗╔╦╦═╗╔══╦═╣╚══╣╔╗╔╗║║ ║║╚╝╝║╚══╗                               |
|   ║║║╔╗║╚╝╠╣╔╗╣╔╗║╔╩══╗║║║║║║║ ║║╔╗║║╔══╝                               |
|   ║╚╣╔╗║║║║║║║║╔╗║║║╚═╝║║║║║║╚═╝║║║╚╣╚══╗                               |
|   ╚═╩╝╚╩╩╩╩╩╝╚╩╝╚╩╝╚═══╩╝╚╝╚╩═══╩╝╚═╩═══╝                               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, G. D'Alessio, A. Parente                           |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of laminarSMOKE solver.                             |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2018 A. Cuoci                                            |
|   laminarSMOKE is free software: you can redistribute it and/or modify  |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   laminarSMOKE is distributed in the hope that it will be useful,       |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with laminarSMOKE. If not, see <http://www.gnu.org/licenses/>.  |
|                                                                         |
\*-----------------------------------------------------------------------*/

#if SPARC==1

	typedef OdeSMOKE::KernelDense<OpenSMOKE::ODE_PFR_SPARC_OpenSMOKE> denseOdePFR_SPARC;
	typedef OdeSMOKE::MethodGear<denseOdePFR_SPARC> methodGearPFR_SPARC;

	PtrList< ODE_PFR_SPARC >    					vector_fine_structure_ode_pfr(nclusters);
	PtrList< OdeSMOKE::MultiValueSolver<methodGearPFR_SPARC> > 	vector_odeSolverPFR(nclusters);
	for(unsigned int i=0;i<nclusters;i++)
	{
		Info << "Creating PFR fine structure: " << i << " - "; 
		vector_fine_structure_ode_pfr.set(i, new ODE_PFR_SPARC(*vector_thermodynamicsMapXML[i], *vector_kineticsMapXML[i]));
		Info << "Number of species: " << vector_fine_structure_ode_pfr[i].NumberOfEquations() << endl; 

		Info << "Creating PFR fine structure ODE solver:   " << i << endl; 
		vector_odeSolverPFR.set(i, new OdeSMOKE::MultiValueSolver<methodGearPFR_SPARC>);
		vector_odeSolverPFR[i].SetReactor(&vector_fine_structure_ode_pfr[i]);
	}

#endif
