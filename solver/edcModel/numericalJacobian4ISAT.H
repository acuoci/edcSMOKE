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

template<typename ODESystem>
void numericalJacobian(const Eigen::VectorXd& y_eig, Eigen::MatrixXd &J, const Eigen::VectorXd& sF, ODESystem *os) 
{
     // Calculated as suggested by Buzzi (private communication)
     const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
     const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
     const double BETA = 1.e3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;
     const double TOLR = 1.e-7;
     const double TOLA = 1.e-12;

     const unsigned int NE = y_eig.size();

     OpenSMOKE::OpenSMOKEVectorDouble y(NE);
     OpenSMOKE::OpenSMOKEVectorDouble y_plus(NE);
     OpenSMOKE::OpenSMOKEVectorDouble dy_plus(NE);
     OpenSMOKE::OpenSMOKEVectorDouble dy_original(NE);

     for(unsigned int i=1;i<=NE;i++)
		y[i] = y_plus[i] = y_eig(i-1);

     // Call equations
     os->GetEquations(y, 0., dy_original);     

     // Derivatives with respect to y[kd]
     for(int kd=1;kd<=NE;kd++)
     {
         double hf = 1.e0;
         double error_weight = 1./(TOLA+TOLR*fabs(y[kd]));
         double hJ = ETA2 * fabs(std::max(y[kd], 1./error_weight));
         double hJf = hf/error_weight;
         hJ = std::max(hJ, hJf);
         hJ = std::max(hJ, ZERO_DER);

         // This is what is done by Buzzi
         double dy = std::min(hJ, 1.e-3 + 1e-3*fabs(y[kd]));
         double udy = 1. / dy;
         y_plus[kd] += dy;

         os->GetEquations(y_plus,0.,dy_plus);

         for(int j=1;j<=NE;j++)
             J(j-1,kd-1) = (dy_plus[j]-dy_original[j]) * udy;

         y_plus[kd] = y[kd];
     }

     for(unsigned int i=0;i<NE;i++)
         for(unsigned int j=0;j<NE;j++)
             J(i,j) *= sF(i);

//     for(unsigned int i=0;i<NE;i++)
//         for(unsigned int j=0;j<NE;j++)
//             J(i,j) *= (1./sF(j));

//     for(unsigned int i=0;i<NE;i++)
//         for(unsigned int j=0;j<NE;j++)
//             J(i,j) *= (sF(i)/sF(j));

 };	
