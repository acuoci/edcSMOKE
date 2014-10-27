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
