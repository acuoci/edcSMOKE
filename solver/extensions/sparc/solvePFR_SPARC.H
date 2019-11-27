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
	// Set initial conditions
	vector_odeSolverPFR[iCluster].SetInitialConditions(t0, yStar0);

	// Additional ODE solver options

	// Set linear algebra options
	vector_odeSolverPFR[iCluster].SetLinearAlgebraSolver(odeParametersFineStructures.linear_algebra());
	vector_odeSolverPFR[iCluster].SetFullPivoting(odeParametersFineStructures.full_pivoting());

	// Set relative and absolute tolerances
	vector_odeSolverPFR[iCluster].SetAbsoluteTolerances(odeParametersFineStructures.absolute_tolerance());
	vector_odeSolverPFR[iCluster].SetRelativeTolerances(odeParametersFineStructures.relative_tolerance());

	// Set minimum and maximum values
	vector_odeSolverPFR[iCluster].SetMinimumValues(yMin);
	vector_odeSolverPFR[iCluster].SetMaximumValues(yMax);
	
						
	// Solve
	OdeSMOKE::OdeStatus status = vector_odeSolverPFR[iCluster].Solve(tf);
	vector_odeSolverPFR[iCluster].Solution(yStarf);
}
