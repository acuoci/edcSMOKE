/******************************************************************************
LIBRARY MuTLib
Copyright (c) 2021 Universidad Nacional Española a Distancia UNED.
All rights reserved.
Authors:    Oscar Córdoba
            Manuel Arias-Zugasti
Transport.cpp implementation
******************************************************************************/

#include <math.h>
#include <iostream>
#include <algorithm>
#include <vector>

#if OPENSMOKE_USE_MKL == 1
#include "mkl.h"
#include "mkl_lapacke.h"
#elif OPENSMOKE_USE_OPENBLAS == 1
#include "cblas.h"
#include "lapacke.h"
#endif

namespace MuTLib
{
    /******************************************************************************
    Preliminary definitions
    *******************************************************************************/
    #define PI   3.14159265358979323846
    #define MAX(a, b) ((a) > (b) ? (a) : (b))

    #include "MuTLib_Omega.h"

    /******************************************************************************
    Class MuTLib_Transport constructor
    *******************************************************************************/
    MuTLib_Transport::MuTLib_Transport()
    {
        varValues.resize(TR_END);
    }

    /******************************************************************************
    Class Transport setup
    *******************************************************************************/
    void MuTLib_Transport::SetUp(MuTLib_Species* species_, MuTLib_KineticData* kindata_, MuTLib_ThermoData* thermdata_)
    {
        species=species_;
        kindata=kindata_;
        thermdata=thermdata_;
    }

    /******************************************************************************
    Kinetic parameters evaluation
        Molecular weight
        Lennard Jones Potential
        Collision diameter
        Rotational collision index
        Viscosity
        Integral collision A* B* C* Omega11 delta*
        deltap for resonant internal energy
    *******************************************************************************/
    void MuTLib_Transport::CalculateKinetic()
    {
        for(unsigned int i=0; i<TR_END; i++)
            varValues[i].resize(species->GetNum()*species->GetNum());

        mass.resize(species->GetNum());
        dens.resize(species->GetNum());
        cint.resize(species->GetNum());
        crot.resize(species->GetNum());
        deltap.resize(species->GetNum());
        chi.resize(species->GetNum());
        mu.resize(species->GetNum());

        double MASS=0;
        for (unsigned int i=0; i<species->GetNum();i++)
            MASS+=species->GetMweight(i)*species->GetMolFrac(i);
        for (unsigned int i=0; i<species->GetNum();i++)
            dens[i]=species->GetMweight(i)*species->GetMolFrac(i)/MASS;
        for (unsigned int i=0; i<species->GetNum();i++)
            mass[i]=species->GetMweight(i);

        for (unsigned int i=0; i<species->GetNum();i++)
        {
            switch(kindata->index[species->GetKPosition(i)])
            {
                case 0:
                    crot[i]=0.;
                    cint[i]=0;
                    break;
                case 1:
                    crot[i]=1.;
                    cint[i]=thermdata->GetCpR(species->GetTPosition(i),
                                              species->temp)-2.5;
                    break;
                case 2:
                    crot[i]=1.5;
                    cint[i]=thermdata->GetCpR(species->GetTPosition(i),
                                              species->temp)-2.5;
                    break;
            }
        }

        for (unsigned int i=0; i<species->GetNum();i++)
        {
            if(kindata->dipole[species->GetKPosition(i)])
                deltap[i]=2985./species->temp/sqrt(species->temp);
            else
                deltap[i]=.0;  
        }

        for (unsigned int i=0; i<species->GetNum();i++)
        {
            unsigned int posi=species->GetKPosition(i);
            chi[i]=MAX(kindata->collision[posi],1.)*
                    ZrotScale(i,298)/ZrotScale(i,species->temp);
        }

        for (unsigned int i=0; i<species->GetNum();i++)
            for (unsigned int j=0; j<species->GetNum();j++)
                CalculateMolecularWeight(i, j);

        for (unsigned int i=0; i<species->GetNum();i++)
            for (unsigned int j=0; j<species->GetNum();j++)
                CalculateZeta(i, j);

        for (unsigned int i=0; i<species->GetNum();i++)
            for (unsigned int j=0; j<species->GetNum();j++)
                CalculateLennardJones(i, j);

        for (unsigned int i=0; i<species->GetNum();i++)
            for (unsigned int j=0; j<species->GetNum();j++)
                CalculateCollisionDiameter(i, j);

        for (unsigned int i=0; i<species->GetNum();i++)
        {
            for (unsigned int j=0; j<species->GetNum();j++)
            {
                CalculateTstar(i, j);
                CalculateDipoleStar(i, j);
                CalculateAstar(i, j);
                CalculateBstar(i, j);
                CalculateCstar(i, j);
                CalculateEstar(i, j);
                CalculateFstar(i, j);
                CalculateOmega11(i, j);
                CalculateOmega22(i, j);
            }
        }

        for (unsigned int i=0; i<species->GetNum();i++)
        {
            double Sj=varValues[TR_COLLDIAM][i*species->GetNum()+i]*1e-10;
            mu[i]=5./16.*sqrt(PI*species->GetMweight(i)*GTOKG*BOLTZSI/AVNUMBER*
                              species->temp);
            mu[i]/=(PI*Sj*Sj*varValues[TR_OMEGA22][i*species->GetNum()+i]);
        }
    }

    /******************************************************************************
    Two species combined moleculat weight
    *******************************************************************************/
    void MuTLib_Transport::CalculateMolecularWeight(const unsigned int i, const unsigned int j)
    {
        varValues[TR_MOLWEIGHT][i*species->GetNum()+j]=
                mass[i]*mass[j]/(mass[i]+mass[j]);
    }

    /******************************************************************************
    Two species combined Lennard Jones potential
    *******************************************************************************/
    void MuTLib_Transport::CalculateLennardJones(const unsigned int i, const unsigned int j)
    {
        unsigned int posi=species->GetKPosition(i);
        unsigned int posj=species->GetKPosition(j);
        unsigned int N=species->GetNum();
        varValues[TR_LENNARD][i*N+j]=sqrt(kindata->potential[posi]*
            kindata->potential[posj])*varValues[TR_ZETA][i*N+j]*
                varValues[TR_ZETA][i*N+j];
    }

    /******************************************************************************
    Two species combined collision diameter
    *******************************************************************************/
    void MuTLib_Transport::CalculateCollisionDiameter(const unsigned int i, const unsigned int j)
    {
        unsigned int posi=species->GetKPosition(i);
        unsigned int posj=species->GetKPosition(j);
        unsigned int N=species->GetNum();
        varValues[TR_COLLDIAM][i*N+j]=(kindata->diameter[posi]+
                kindata->diameter[posj])/2./pow(varValues[TR_ZETA][i*N+j],1./6.);
    }

    /******************************************************************************
    Two species combined reduced temperature T*
    *******************************************************************************/
    void MuTLib_Transport::CalculateTstar(const unsigned int i, const unsigned int j)
    {
        varValues[TR_TSTAR][i*species->GetNum()+j]=species->temp/
            varValues[TR_LENNARD][i*species->GetNum()+j];
    }

    /******************************************************************************
    Two species combined delta*
    *******************************************************************************/
    void MuTLib_Transport::CalculateDipoleStar(const unsigned int i, const unsigned int j)
    {
        unsigned int N=species->GetNum();
        unsigned int posi=species->GetKPosition(i);
        unsigned int posj=species->GetKPosition(j);
        varValues[TR_DIPOLESTAR][i*N+j]=.5*
                kindata->dipole[posi]*kindata->dipole[posj]/(
                varValues[TR_COLLDIAM][i*N+j]*
                varValues[TR_COLLDIAM][i*N+j]*
                varValues[TR_COLLDIAM][i*N+j]*
                varValues[TR_LENNARD][i*N+j]);
    }

    /******************************************************************************
    Two species combined omega11
    *******************************************************************************/
    void MuTLib_Transport::CalculateOmega11(const unsigned int i, const unsigned int j)
    {
        //extern double omega11[37][8];
        unsigned int N=species->GetNum();

        double delta= varValues[TR_DIPOLESTAR][i*N+j];
        double tstar= species->temp/varValues[TR_LENNARD][i*N+j];
        varValues[TR_OMEGA11][i*N+j]=kindata->Interpol(delta, tstar, MuTLib_omega11);
    }

    /******************************************************************************
    Two species combined omega22
    *******************************************************************************/
    void MuTLib_Transport::CalculateOmega22(const unsigned int i, const unsigned int j)
    {
        //extern double omega22[37][8];
        unsigned int N=species->GetNum();

        double delta= varValues[TR_DIPOLESTAR][i*N+j];
        double tstar= species->temp/varValues[TR_LENNARD][i*N+j];
        varValues[TR_OMEGA22][i*N+j]=kindata->Interpol(delta, tstar, MuTLib_omega22);
    }

    /******************************************************************************
    Two species combined A*
    *******************************************************************************/
    void MuTLib_Transport::CalculateAstar(const unsigned int i, const unsigned int j)
    {
        //extern double astar[37][8];
        unsigned int N=species->GetNum();

        double delta= varValues[TR_DIPOLESTAR][i*N+j];
        double tstar= species->temp/varValues[TR_LENNARD][i*N+j];
        varValues[TR_ASTAR][i*N+j]=kindata->Interpol(delta, tstar, MuTLib_astar);
    }

    /******************************************************************************
    Two species combined B*
    *******************************************************************************/
    void MuTLib_Transport::CalculateBstar(const unsigned int i, const unsigned int j)
    {
        //extern double bstar[37][8];
        unsigned int N=species->GetNum();

        double delta= varValues[TR_DIPOLESTAR][i*N+j];
        double tstar= species->temp/varValues[TR_LENNARD][i*N+j];
        varValues[TR_BSTAR][i*N+j]=kindata->Interpol(delta, tstar, MuTLib_bstar);
    }

    /******************************************************************************
    Two species combined C*
    *******************************************************************************/
    void MuTLib_Transport::CalculateCstar(const unsigned int i, const unsigned int j)
    {
        //extern double cstar[37][8];
        unsigned int N=species->GetNum();

        double delta= varValues[TR_DIPOLESTAR][i*N+j];
        double tstar= species->temp/varValues[TR_LENNARD][i*N+j];
        varValues[TR_CSTAR][i*N+j]=kindata->Interpol(delta, tstar, MuTLib_cstar);
    }

    /******************************************************************************
    Two species combined E*
    *******************************************************************************/
    void MuTLib_Transport::CalculateEstar(const unsigned int i, const unsigned int j)
    {
        //extern double estar[37][8];
        unsigned int N=species->GetNum();

        double delta= varValues[TR_DIPOLESTAR][i*N+j];
        double tstar= species->temp/varValues[TR_LENNARD][i*N+j];
        varValues[TR_ESTAR][i*N+j]=kindata->Interpol(delta, tstar, MuTLib_estar);
    }

    /******************************************************************************
    Two species combined F*
    *******************************************************************************/
    void MuTLib_Transport::CalculateFstar(const unsigned int i, const unsigned int j)
    {
        //extern double fstar[37][8];
        unsigned int N=species->GetNum();

        double delta= varValues[TR_DIPOLESTAR][i*N+j];
        double tstar= species->temp/varValues[TR_LENNARD][i*N+j];
        varValues[TR_FSTAR][i*N+j]=kindata->Interpol(delta, tstar, MuTLib_fstar);
    }

    /******************************************************************************
    Two species combined polar-nonpolar zeta, necessary for collision dismeter and
    Lennard Jones potential
    *******************************************************************************/
    void MuTLib_Transport::CalculateZeta(const unsigned int i, const unsigned int j)
    {
        unsigned int posi=species->GetKPosition(i);
        unsigned int posj=species->GetKPosition(j);
        unsigned int N=species->GetNum();

        if ((kindata->dipole[posi] && kindata->dipole[posj]) ||
                (kindata->dipole[posi] && kindata->dipole[posj])){
            varValues[TR_ZETA][i*N+j]=1.;
            return;
        }
        unsigned int posp;
        unsigned int posn;
        if (kindata->dipole[posi]){
            posp=posi;
            posn=posj;
        }else{
            posp=posj;
            posn=posi;
        }

        double mupstar=kindata->dipole[posp]/sqrt(kindata->potential[posp]*
            kindata->diameter[posp]*kindata->diameter[posp]);
        double alphanstar=kindata->polarizability[posn]/kindata->diameter[posn]/
                kindata->diameter[posn]/kindata->diameter[posn];
        varValues[TR_ZETA][i*N+j]=1.+.25*alphanstar*mupstar*sqrt(
                    kindata->potential[posp]/kindata->potential[posn]);
    }

    /******************************************************************************
    Parker correction for rotational collision index
    *******************************************************************************/
    double MuTLib_Transport::ZrotScale(const unsigned int i, const double temp)
    {
        double scale;
        unsigned int posi=species->GetKPosition(i);
        double epst=kindata->potential[posi]/temp;
        scale=1.+.5*PI*sqrt(PI)*sqrt(epst);
        scale+=(PI*PI/4.+2.)*epst;
        scale+=PI*sqrt(PI)*epst*sqrt(epst);
        return scale;
    }

    /******************************************************************************
    Binary diffusion coefficient
    *******************************************************************************/
    double MuTLib_Transport::GetDij(const unsigned int i, const unsigned int j) const
    {
        double Dij;
        double Sij;
        unsigned int N=species->GetNum();
        Dij=sqrt(BOLTZSI*species->temp/2/PI/varValues[TR_MOLWEIGHT][i*N+j]/
                GTOKG*AVNUMBER);
        Sij=varValues[TR_COLLDIAM][i*N+j]*1e-10;
        Dij*=3./8./Sij/Sij/varValues[TR_OMEGA11][i*N+j];
        Dij/=(species->press/BOLTZSI/species->temp);

        return Dij;
    }

    /******************************************************************************
    Binary diffusion coefficients ratio. Involved the species N with higher molar
    fraction
    *******************************************************************************/
    double MuTLib_Transport::GetDjNDil(const unsigned int i, const unsigned int j, const unsigned int l, const unsigned int N) const
    {
        unsigned int indil;
        unsigned int indjN;
        indil=i*species->GetNum()+l;
        indjN=j*species->GetNum()+N;
        double DjNDil;
        DjNDil =varValues[TR_COLLDIAM][indil]*varValues[TR_COLLDIAM][indil];
        DjNDil/=(varValues[TR_COLLDIAM][indjN]*varValues[TR_COLLDIAM][indjN]);
        DjNDil*=sqrt(varValues[TR_MOLWEIGHT][indil]/varValues[TR_MOLWEIGHT][indjN]);
        DjNDil*=varValues[TR_OMEGA11][indil]/varValues[TR_OMEGA11][indjN];

        return DjNDil;
    }

    /******************************************************************************
    Binary diffusion coefficients ratio. Involved the diagonal species
    *******************************************************************************/
    double MuTLib_Transport::GetDjjDil(const unsigned int i, const unsigned int j, const unsigned int l) const
    {
        unsigned int indil;
        unsigned int indjj;
        indil=i*species->GetNum()+l;
        indjj=j*species->GetNum()+j;
        double DjjDil;
        DjjDil =varValues[TR_COLLDIAM][indil]*varValues[TR_COLLDIAM][indil];
        DjjDil/=(varValues[TR_COLLDIAM][indjj]*varValues[TR_COLLDIAM][indjj]);
        DjjDil*=sqrt(varValues[TR_MOLWEIGHT][indil]/
           varValues[TR_MOLWEIGHT][indjj]);
        DjjDil*=varValues[TR_OMEGA11][indil]/varValues[TR_OMEGA11][indjj];

        return DjjDil;
    }

    /******************************************************************************
    Matrix system elements calculation 00,00
    *******************************************************************************/
    double MuTLib_Transport::GetFijl_00_00(const unsigned int i, const unsigned int j, const unsigned int l) const
    {
        if (i==j)
        {
            if (l==j)
                return 0.;
            return 1.;
        }
        else
        {
            if (l!=j)
                return 0.;
            return -1.;
        }
    }

    /******************************************************************************
    Matrix system elements calculation 00,10
    *******************************************************************************/
    double MuTLib_Transport::GetFijl_00_10(const unsigned int i, const unsigned int j, const unsigned int l) const
    {
        double cstar_;

        double mi=species->GetMweight(i);
        if (i==j)
        {
            if (l==j)
                return 0.;

            double ml;
            ml=species->GetMweight(l);
            cstar_=varValues[TR_CSTAR][i*species->GetNum()+l];
            return -.5*ml/(mi+ml)*(6.*cstar_-5.);
        }
        else
        {
            if (l!=j)
                return 0.;

            double mj=species->GetMweight(j);
            cstar_=varValues[TR_CSTAR][i*species->GetNum()+l];
            return .5*mi/(mi+mj)*(6.*cstar_-5.);
        }
    }

    /******************************************************************************
    Matrix system elements calculation 10,00
    *******************************************************************************/
    double MuTLib_Transport::GetFijl_10_00(const unsigned int i, const unsigned int j, const unsigned int l) const
    {
        if (i==j)
            return GetFijl_00_10(j,i,l);

        if (l!=j)
            return 0.;

        double mi=species->GetMweight(i);
        double mj=species->GetMweight(j);
        double cstar_;

        cstar_=varValues[TR_CSTAR][i*species->GetNum()+j];

        return .5*mj/(mi+mj)*(6.*cstar_-5.);
    }

    /******************************************************************************
    Matrix system elements calculation 10,10
    *******************************************************************************/
    double MuTLib_Transport::GetFijl_10_10(const unsigned int i, const unsigned int j, const unsigned int l) const
    {
        unsigned int N=species->GetNum();

        double bstar_;
        double astar_;

        double polyatomic=0.;
        double rotation=0;
        double factor=0.;
        double mi;
        double ml;

        if (i==j)
            factor+=1.;
        if (j==l)
            factor+=1.;
        if (factor)
        {
            mi=species->GetMweight(i);
            ml=species->GetMweight(l);
            astar_=varValues[TR_ASTAR][i*N+l];
            polyatomic=20.*mi*ml/(mi+ml)/(mi+ml)/3./PI*astar_;
            rotation+=crot[i]/chi[i];
            rotation+=crot[l]/chi[l];
            polyatomic*=rotation*factor;
        }
        else
            polyatomic=.0;

        if (i==j)
        {
            astar_=varValues[TR_ASTAR][i*N+i];
            if (l==j)
                return factor*astar_+polyatomic;

            double mi;
            mi=species->GetMweight(i);
            double ml;
            ml=species->GetMweight(l);
            bstar_=varValues[TR_BSTAR][i*N+l];
            astar_=varValues[TR_ASTAR][i*N+l];
            return (4.*mi*ml*astar_+25./4.*ml*ml-3.*ml*ml*bstar_+15./2.*mi*mi)/
                    (mi+ml)/(mi+ml)+polyatomic;
        }
        else
        {
            if (l!=j)
                return 0.;
            double mi;
            double mj;
            mi=species->GetMweight(i);
            mj=species->GetMweight(j);
            bstar_=varValues[TR_BSTAR][i*N+j];
            astar_=varValues[TR_ASTAR][i*N+j];
            return mi*mj/(mi+mj)/(mi+mj)*(4.*astar_-55./4.+3.*bstar_)+polyatomic;
        }
    }

    /******************************************************************************
    Matrix system elements calculation 10,01
    *******************************************************************************/
    double MuTLib_Transport::GetFijl_10_01(const unsigned int i, const unsigned int j, const unsigned int l) const
    {

        double factor=0.;

        if (i==j)
            factor+=1.;
        if (j==l)
            factor+=1.;
        if (factor==0.)
            return 0.;

        unsigned int N=species->GetNum();
        double mi;
        double mj;
        double ml;
        double astar_;
        mi=species->GetMweight(i);
        mj=species->GetMweight(j);
        ml=species->GetMweight(l);
        astar_=varValues[TR_ASTAR][i*N+l];

        return -10.*mj*astar_/(mi+ml)/PI*crot[j]/
            chi[j]*factor;
    }

    /******************************************************************************
    Matrix system elements calculation 01,10
    *******************************************************************************/
    double MuTLib_Transport::GetFijl_01_10(const unsigned int i, const unsigned int j, const unsigned int l) const
    {
        double factor=0.;
        if (i==j)
            factor+=1.;
        if (j==l)
            factor+=1.;

        if (factor==0.)
            return 0.;

        unsigned int N=species->GetNum();
        double mi;
        double ml;
        double astar_;
        mi=species->GetMweight(i);
        ml=species->GetMweight(l);
        astar_=varValues[TR_ASTAR][i*N+l];

        return -10.*mi*astar_/(mi+ml)/PI*crot[i]/
            chi[i]*factor;
    }

    /******************************************************************************
    Matrix system elements calculation 01,01
    *******************************************************************************/
    double MuTLib_Transport::GetFijl_01_01(const unsigned int i, const unsigned int j, const unsigned int l) const
    {

        if (i!=j)
            return 0.;

        if (!kindata->index[species->GetKPosition(i)] && i==j)
                return 1.0;

        unsigned int N=species->GetNum();
        double mi;
        double ml;
        double astar_;
        mi=species->GetMweight(i);
        ml=species->GetMweight(l);
        astar_=varValues[TR_ASTAR][i*N+l];
        double polyatomic;
        polyatomic=12./5./PI/cint[i]*mi/ml*astar_*crot[i]/
                chi[i];

        polyatomic+=1.;
        if (i==l)
                polyatomic+=deltap[i];
        polyatomic*=25./4.*cint[j];

        return polyatomic;
    }

    /******************************************************************************
    Scaled Matrix system elements calculation 00,00
    *******************************************************************************/
    double MuTLib_Transport::GetFij_00_00(const unsigned int i, const unsigned int j, const unsigned int N) const
    {
        double sum=0;
        double Fij_00_00;
        for (unsigned int l=0; l<species->GetNum();l++)
        {
            Fij_00_00=GetFijl_00_00(i,j,l);
            if (Fij_00_00==0.)
                continue;
            if (i==j)
            {
                sum+=species->GetMolFrac(l)*GetDjNDil(i,j,l,N)*Fij_00_00;
            }
            else
            {
                if (j==l)
                    sum+=species->GetMolFrac(i)*GetDjNDil(i,j,l,N)*Fij_00_00;
                else
                    std::cout<<"Something wrong in GetFij_00_00"<<std::endl;
            }
        }
        return sum;
    }

    /******************************************************************************
    Scaled Matrix system elements calculation 00,10
    *******************************************************************************/
    double MuTLib_Transport::GetFij_00_10(const unsigned int i, const unsigned int j) const
    {
        double sum=0;
        double Fij_00_10;
        for (unsigned int l=0; l<species->GetNum();l++)
        {
            Fij_00_10=GetFijl_00_10(i,j,l);
            if (Fij_00_10==0.)
                continue;
            if (i==j)
                sum+=species->GetMolFrac(l)*GetDjjDil(i,j,l)*Fij_00_10;
            else
            {
                if (j==l)
                    sum+=species->GetMolFrac(i)*GetDjjDil(i,j,l)*Fij_00_10;
                else
                    std::cout<<"Something wrong in GetFij_00_10"<<std::endl;
            }
        }
        return sum;
    }

    /******************************************************************************
    Scaled Matrix system elements calculation 10,00
    *******************************************************************************/
    double MuTLib_Transport::GetFij_10_00(const unsigned int i, const unsigned int j, const unsigned int N) const
    {
        double num=0.;
        double Fij_10_00;

        for (unsigned int l=0; l<species->GetNum();l++){
            Fij_10_00=GetFijl_10_00(i,j,l);

            if (Fij_10_00==0.)
                continue;
            if (i==j)
                num+=species->GetMolFrac(l)*GetDjNDil(i,j,l,N)*Fij_10_00;
            else
            {
                if (j==l)
                    num+=species->GetMolFrac(i)*GetDjNDil(i,j,l,N)*Fij_10_00;
                else
                    std::cout << "Something wrong in GetFij_10_00 " << std::endl;
            }
        }

        double den=0.;
        double Fii_10_10;
        for (unsigned int l=0; l<species->GetNum();l++)
        {
            Fii_10_10=GetFijl_10_10(i,i,l);
            den+=species->GetMolFrac(l)*GetDjjDil(i,i,l)*Fii_10_10;
        }

        return num/den;
    }

    /******************************************************************************
    Scaled Matrix system elements calculation 10,10
    *******************************************************************************/
    double MuTLib_Transport::GetFij_10_10(const unsigned int i, const unsigned int j) const
    {
        double num=0;
        double Fij_10_10;
        if (i==j)
            return 1.;

        for (unsigned int l=0; l<species->GetNum();l++)
        {
            Fij_10_10=GetFijl_10_10(i,j,l);
            if (Fij_10_10==0.)
                continue;
            if (i==j)
                num+=species->GetMolFrac(l)*GetDjjDil(i,j,l)*Fij_10_10;
            else
            {
                if (j==l)
                    num+=species->GetMolFrac(i)*GetDjjDil(i,j,l)*Fij_10_10;
                else
                    std::cout<<"Something wrong in GetFij_10_10"<<std::endl;
            }
        }

        double den=0.;
        double Fii_10_10;
        for (unsigned int l=0; l<species->GetNum();l++)
        {
            Fii_10_10=GetFijl_10_10(i,i,l);
            den+=species->GetMolFrac(l)*GetDjjDil(i,i,l)*Fii_10_10;
        }

        return num/den;
    }

    /******************************************************************************
    Scaled Matrix system elements calculation 10,01
    *******************************************************************************/
    double MuTLib_Transport::GetFij_10_01(const unsigned int i, const unsigned int j) const
    {
        double num=0;
        double Fij_10_01;
        for (unsigned int l=0; l<species->GetNum();l++)
        {
            Fij_10_01=GetFijl_10_01(i,j,l);
            if (Fij_10_01==0.)
                continue;
            if (i==j)
                num+=species->GetMolFrac(l)*GetDjjDil(i,j,l)*Fij_10_01;
            else
            {
                if (j==l)
                    num+=species->GetMolFrac(i)*GetDjjDil(i,j,l)*Fij_10_01;
                else
                    std::cout<<"Something wrong in GetFij_10_01"<<std::endl;
            }
        }

        double den=0.;
        double Fii_10_10;
        for (unsigned int l=0; l<species->GetNum();l++)
        {
            Fii_10_10=GetFijl_10_10(i,i,l);
            den+=species->GetMolFrac(l)*GetDjjDil(i,i,l)*Fii_10_10;
        }

        return num/den;
    }

    /******************************************************************************
    Scaled Matrix system elements calculation 01,10
    *******************************************************************************/
    double MuTLib_Transport::GetFij_01_10(const unsigned int i, const unsigned int j) const
    {
        double num=0;
        double Fij_01_10;
        for (unsigned int l=0; l<species->GetNum();l++)
        {
            Fij_01_10=GetFijl_01_10(i,j,l);
            if (Fij_01_10==0.)
                continue;
            if (i==j)
                num+=species->GetMolFrac(l)*GetDjjDil(i,j,l)*Fij_01_10;
            else
            {
                if (j==l)
                    num+=species->GetMolFrac(i)*GetDjjDil(i,j,l)*Fij_01_10;
                else
                    std::cout<<"Something wrong in GetFij_01_10"<<std::endl;
            }
        }

        double den=0.;
        double Fii_01_01;
        for (unsigned int l=0; l<species->GetNum();l++){
            Fii_01_01=GetFijl_01_01(i,i,l);
            den+=species->GetMolFrac(l)*GetDjjDil(i,i,l)*Fii_01_01;
        }

        return num/den;

    }

    /******************************************************************************
    Scaled Matrix system elements calculation 01,01
    *******************************************************************************/
    double MuTLib_Transport::GetFij_01_01(const unsigned int i, const unsigned int j) const
    {
        if (i==j)
            return 1.;

        return 0.;
    }

    /******************************************************************************
    Scaled Matrix system elements calculation 00,00
    Elment N iwth higher molar mas substituted
    *******************************************************************************/
    double MuTLib_Transport::Get_Fij_00_00(const unsigned int i, const unsigned int j, const unsigned int N) const
    {
        return GetFij_00_00(i,j,N)-species->GetMweight(j)/species->GetMweight(N)*
                GetDjNDil(N,j,N,N)* GetFij_00_00(i,N,N);
    }

    /******************************************************************************
    Scaled Matrix system elements calculation 10,00
    Elment N iwth higher molar mas substituted
    *******************************************************************************/
    double MuTLib_Transport::Get_Fij_10_00(const unsigned int i, const unsigned int j, const unsigned int N) const
    {
        return GetFij_10_00(i,j,N)-species->GetMweight(j)/species->GetMweight(N)*
                GetDjNDil(N,j,N,N)* GetFij_10_00(i,N,N);
    }

    /******************************************************************************
    Neumann power series step
    *******************************************************************************/
    void MuTLib_Transport::NeumannStep(const unsigned int r, const unsigned int N)
    {
        std::vector<double> matrixB;
        matrixB.resize(N*N,0);

        for (unsigned int i=0; i<N && r; i++)
        {
            for (unsigned int j=0; j<N; j++)
            {
                for (unsigned int k=0; k<N; k++)
                {
                    matrixB[i*N+j]-=matrixAr[i*N+k]*matrixA[k*N+j];
                }
            }
        }

        for (unsigned int i=0; i<N; i++)
        {
            for (unsigned int j=0; j<N; j++)
            {
                if (r)
                    matrixAr[i*N+j]=matrixB[i*N+j];
                else
                {
                    if (i==j)
                        matrixAr[i*N+j]=1.;
                }
            }
        }

        for (unsigned int i=0; i<N; i++)
            for (unsigned int j=0; j<N; j++)
                matrixSAr[i*N+j]+=matrixAr[i*N+j];
    }

    /******************************************************************************
    PrintMatrix function, general method for debugging purposes
    *******************************************************************************/
    void MuTLib_Transport::PrintMatrix(const std::vector<double>& matrix, const unsigned int m, const unsigned int n)
    {
        for(unsigned int i=0; i<m; i++)
        {
            for (unsigned int j=0; j<n; j++)
            {
                std::cout << matrix[i*n+j] << " ";
            }

            std::cout << std::endl;
        }

        std::cout << std::endl;
    }

    /******************************************************************************
    Calculate, general driver to specific transport calculations
    *******************************************************************************/
    void MuTLib_Transport::Calculate(const bool exact, const unsigned int Fickr_, const int FickM_, const double FickThr_, const unsigned int Soretr_, const int SwapPolicy_)
    {
        CalculateKinetic();

        unsigned int N=species->GetNum();

        Fickr  = Fickr_;
        FickM  = FickM_;
        FickThr= FickThr_;
        Soretr = Soretr_;
        SwapPolicy = SwapPolicy_;

        // 1. Search for the swap species as the most abundant species
        if (SwapPolicy == -2)
        {
            double molarhigh = 0.;
            for (unsigned int i = 0; i < N; i++)
            {
                if (species->GetMolFrac(i) > molarhigh)
                {
                    molarhigh = species->GetMolFrac(i);
                    swap = i;
                }
            }
        }
        // 2. Set the swap species as the last species
        else if (SwapPolicy == -1)
        {
            swap = N - 1;
        }
        // 3. Set the swap species as the user-defined species (fixed)
        else
        {
            swap = SwapPolicy;
        }

        if (!exact)
        {
            matrixFick.clear();
            matrixFick.resize(N*N,0.);
            CalculateMatrixFick();
        }
        else
        {
            matrixFick.clear();
            matrixFick.resize(N*N,0.);
            CalculateMatrixFickExact();
        }

        CalculateSoretCoef();
        CalculateThermalConductivity();
    }

    /******************************************************************************
    Calculate, general driver to exact transport calculations
    *******************************************************************************/
    void MuTLib_Transport::CalculateExact()
    {
        CalculateKinetic();
        CalculateMatrixFickExact();
        CalculateSoretCoefExact();
        CalculateThermalConductivity();
    }


    /******************************************************************************
    Calculate, general driver to exact matrix Fick calculations
    *******************************************************************************/
    void MuTLib_Transport::CalculateMatrixFickExact()
    {
        unsigned int N=species->GetNum();
        matrixFick.clear();
        matrixFick.resize(N*N);
        
        // 1. Search for the swap species as the most abundant species
        if (SwapPolicy == -2)
        {
            double molarhigh = 0.;
            for (unsigned int i = 0; i < N; i++)
            {
                if (species->GetMolFrac(i) > molarhigh)
                {
                    molarhigh = species->GetMolFrac(i);
                    swap = i;
                }
            }
        }
        // 2. Set the swap species as the last species
        else if (SwapPolicy == -1)
        {
            swap = N - 1;
        }
        // 3. Set the swap species as the user-defined species (fixed)
        else
        {
            swap = SwapPolicy;
        }

        //MatrixFick should be called here
        unsigned int count=0;
        for(unsigned int i=0; i<N;i++)
        {
            if (i==swap)
                   continue;
            for(unsigned int j=0; j<N;j++)
            {
                if (j==swap)
                       continue;
                matrixFick[count++]= Get_Fij_00_00(i, j, swap);
            }
        }

        int NN=static_cast<int>(N-1);
        int LDA=NN;
        std::vector<int> IPIV(NN);
        
        int INFO1 = LAPACKE_dgetrf(CblasRowMajor, NN, NN, &matrixFick.front(), LDA, &IPIV.front());
        int INFO2 = LAPACKE_dgetri(CblasRowMajor, NN, &matrixFick.front(), LDA, &IPIV.front());
    }

    /******************************************************************************
    Calculates the scaled Fick diffusion coefficients matrix
    *******************************************************************************/
    void MuTLib_Transport::CalculateMatrixFick()
    {
        unsigned int N=species->GetNum()-1;
        unsigned int count;

        std::vector<double> molefractions(N+1);
        for (unsigned int i=0; i<N+1;i++)
            molefractions[i]=species->GetMolFrac(i);

        if(FickM>=0)
        {
            unsigned int uMuser=static_cast<unsigned int>(FickM);
            if (uMuser>N)
                uMuser=N;

            std::vector<double> sortedfractions;
            sortedfractions.reserve(N);
            for(unsigned int i=0;i<N+1; i++){
                if(i==swap)
                       continue;
                sortedfractions.push_back(molefractions[i]);
            }

            std::sort(sortedfractions.begin(),sortedfractions.end());

            if (uMuser==0)
                FickThr=1.;
            else if (uMuser==N)
                FickThr=.0;
            else
                FickThr=(sortedfractions[N-uMuser]+
                           sortedfractions[N-uMuser-1])/2;
            FickThr/=molefractions[swap];
        }

        matrixA.resize(N*N);
        matrixAr.resize(N*N);
        matrixSAr.resize(N*N);
        for (unsigned int i=0; i<N*N;i++)
        {
            matrixA[i]=0;
            matrixAr[i]=0;
            matrixSAr[i]=0;
        }

        //Sets M indexes and molar fractions to 0 for i=M+1,...,N
        std::vector<unsigned int> Mspecies;
        for (unsigned int i=0; i<N+1; i++)
        {
            if(i==swap)
                continue;
            if (molefractions[i]/molefractions[swap]>FickThr)
                Mspecies.push_back(i);
        }

        //Derives to Model 1 if M=0
        unsigned int M=Mspecies.size();
        if (M==0)
        {
            for (unsigned int i=0; i<N+1; i++)
                species->SetMolFrac(i,molefractions[i]);

            CalculateModel1Fick();

            return;
        }

        for (unsigned int i=0; i<N+1; i++)
        {
            if(i==swap)
            {
                species->SetMolFrac(i,molefractions[i]);
                continue;
            }
            if (molefractions[i]/molefractions[swap]>FickThr)
                species->SetMolFrac(i,molefractions[i]);
            else
                species->SetMolFrac(i,.0);
        }

        std::vector<double> matrixIA0;
        matrixIA0.resize(N*N);
        std::vector<double> matrixA11;
        matrixA11.resize(M*M,0.);
        std::vector<double> matrixA12;
        matrixA12.resize(M*(N-M),0.);
        std::vector<double> matrixAux;
        matrixAux.resize(M*(N-M),0.);
        std::vector<double> matrixA22;
        matrixA22.resize(N-M,0.);

        std::vector<bool> IsMspecies;
        IsMspecies.resize(N+1,false);

        //Calculates Matrix (I+A11)(0)
        for ( unsigned int i=0; i<M; i++)
        {
            IsMspecies[Mspecies[i]]=true;
            for ( unsigned int j=0; j<M; j++)
                matrixA11[i*M+j]=Get_Fij_00_00(Mspecies[i],Mspecies[j],swap);
        }

        //New indexes once the M species are known
        index.resize(N+1);
        unsigned int countM=0;
        unsigned int countN=0;
        for (unsigned int i=0; i<N+1; i++)
        {
            if (i==swap)
                continue;
            if (IsMspecies[i])
                index[i]=countM++;
            else
                index[i]=M+countN++;
        }

        //Calculates matrix (A12)(0)==(I+A12)(0)
        count=0;
        for(unsigned int i=0; i<M; i++)
        {
            for(unsigned int j=0; j<N+1; j++)
            {
                if (IsMspecies[j] || j==swap)
                    continue;
                matrixA12[count++]=Get_Fij_00_00(Mspecies[i],j,swap);
            }
        }

        //Calculates matrix (I+A22)(0)
        count=0;
        for ( unsigned int i=0; i<N+1; i++)
        {
            if (IsMspecies[i] || i==swap)
                continue;
            matrixA22[count++]=Get_Fij_00_00(i,i,swap);
        }

        for(unsigned int i=0; i<M;i++)
        {
            for(unsigned int j=0; j<M;j++)
            {
                matrixIA0[i*N+j]=matrixA11[i*M+j];
            }
        }
        for(unsigned int i=0; i<M;i++)
        {
            for(unsigned int j=0; j<(N-M);j++)
            {
                matrixIA0[i*N+(j+M)]=matrixA12[i*(N-M)+j];
            }
        }
        for(unsigned int i=0; i<(N-M);i++)
            matrixIA0[(i+M)*N+(i+M)]=matrixA22[i];

        //Inverting (I+A(0))==> (A12)(0)*inv((A22)(0))
        for ( unsigned int i=0; i<M; i++)
            for ( unsigned int j=0; j<N-M; j++)
                matrixA12[i*(N-M)+j]/=matrixA22[j];

        //Inverting (I+A(0))==> inv((A11)(0))
        int MM=static_cast<int>(M);
        
        //Calculates invA11
        int LDA=MM;
        std::vector<int> IPIV(M);

        int INFO1 = LAPACKE_dgetrf(CblasRowMajor, MM, MM, &matrixA11.front(), LDA, &IPIV.front());
        int INFO2 = LAPACKE_dgetri(CblasRowMajor, MM, &matrixA11.front(), LDA, &IPIV.front());

        //Inverting (I+A(0))==> inv((A11)(0))*(A12)(0)*inv((A22)(0))
        for(unsigned int i=0; i<M; i++)
        {
            for(unsigned int j=0; j<N-M; j++)
            {
                for (unsigned int k=0; k<M; k++)
                {
                    matrixAux[i*(N-M)+j]+=
                            matrixA11[i*M+k]*
                            matrixA12[k*(N-M)+j];
                }
            }
        }

        //Setting inv(I+A(0)) in a correct order==> inv(A11(0))
        std::vector<double> matrixInvIA0;
        std::vector<double> matrixA1;
        matrixInvIA0.resize(N*N);
        matrixA1.resize(N*N);
        for(unsigned int i=0; i<M;i++)
        {
            for(unsigned int j=0; j<M;j++)
            {
                matrixInvIA0[i*N+j]=matrixA11[i*M+j];
            }
        }

        //Setting inv(I+A(0)) in a correct order==>
        //inv((A11)(0))*(A12)(0)*inv((A22)(0))
        for(unsigned int i=0; i<M;i++)
        {
            for(unsigned int j=0; j<(N-M);j++){
                matrixInvIA0[i*N+(j+M)]=-matrixAux[i*(N-M)+j];
            }
        }

        //Setting inv(I+A(0)) in a correct order==> inv(A11(0))
        for(unsigned int i=0; i<(N-M);i++)
        {
            matrixInvIA0[(i+M)*N+(i+M)]=1/matrixA22[i];
        }

        for (unsigned int i=0; i<species->GetNum();i++)
               species->SetMolFrac(i,molefractions[i]);

        for(unsigned int i=0; i<N+1; i++)
        {
            if (i==swap)
                continue;
            for(unsigned int j=0; j<N+1; j++)
            {
                if (j==swap)
                    continue;
                matrixA1[index[i]*N+index[j]]=
                        Get_Fij_00_00(i,j,swap)-
                        matrixIA0[index[i]*N+index[j]];
            }
        }

        //MatrixA is B in the paper
        for(unsigned int i=0; i<N; i++)
        {
            for(unsigned int j=0; j<N; j++)
            {
                matrixA[i*N+j]=0;
                for (unsigned int k=0; k<N; k++)
                    matrixA[i*N+j]+=matrixInvIA0[i*N+k]*matrixA1[k*N+j];
            }
        }

        for (unsigned int i=0; i<=static_cast<unsigned int>(Fickr); i++)
            NeumannStep(i,N);

        //Final multiplication inv(I+B)*inv(I+A0)
        for(unsigned int i=0; i<N; i++)
        {
            for(unsigned int j=0; j<N; j++)
            {
                matrixA[i*N+j]=0;
                for (unsigned int k=0; k<N; k++)
                    matrixA[i*N+j]+=matrixSAr[i*N+k]*matrixInvIA0[k*N+j];
            }
        }

        unsigned int iminus=0;
        unsigned int jminus;
        for(unsigned int i=0; i<N+1; i++)
        {
            if (i==swap)
                continue;
            jminus=0;
            for(unsigned int j=0; j<N+1; j++)
            {
                if (j==swap)
                    continue;
                matrixFick[iminus*N+jminus]=matrixA[index[i]*N+index[j]];
                jminus++;
            }
            iminus++;
        }
    }

    /******************************************************************************
    Calculates the unscaled Fick diffusion coefficients matrix
    *******************************************************************************/
    void MuTLib_Transport::UnscaleMatrixFick()
    {
        unsigned int N=species->GetNum();

        FickCoef.clear();
        FickCoef.resize(N*N);

        for (unsigned int j=0; j<N; j++)
            FickCoef[swap*N+j]=0;

        for (unsigned int i=0; i<N; i++)
        {
            double sumj=0.;
            unsigned int jminus=0;
            for (unsigned int j=0; j<N; j++)
            {
                if (j==swap)
                    continue;
                double sumk=0;
                unsigned int kminus=0;
                for (unsigned int k=0; k<N; k++)
                {
                    if (k==swap)
                        continue;
                    double factor=dens[k];
                    if (i==k)
                        factor-=1.0;
                    sumk+=matrixFick[jminus*(N-1)+kminus]*factor;
                    kminus++;
                }
                sumj+=species->GetMweight(j)*GetDij(j,swap)*sumk;
                jminus++;
            }
            FickCoef[i*N+swap]=sumj/species->GetMweight(swap)/
                    species->GetMolFrac(swap);
        }

        unsigned int iminus=0;
        for (unsigned int i=0; i<N; i++)
        {
            if (i==swap)
                continue;
            unsigned int jminus=0;
            for (unsigned int j=0; j<N; j++)
            {
                if (j==swap)
                    continue;
                FickCoef[i*N+j]=species->GetMolFrac(i)*FickCoef[i*N+swap]+
                        matrixFick[iminus*(N-1)+jminus]*GetDij(i,swap);
                FickCoef[swap*N+j]-=FickCoef[i*N+j]*species->GetMweight(i);
                jminus++;
            }
            iminus++;
        }

        for (unsigned int i=0; i<N; i++)
            FickCoef[i*N+swap]*=species->GetMolFrac(i);

        for (unsigned int j=0; j<N; j++)
            if (j!=swap)
                FickCoef[swap*N+j]/=species->GetMweight(swap);
    }

    /******************************************************************************
    Model 1 for Fick diffusion coefficients
    *******************************************************************************/
    void MuTLib_Transport::CalculateModel1Fick()
    {

        unsigned int count=0;

        unsigned int N=species->GetNum()-1;
        matrixA.resize(N*N);
        matrixAr.resize(N*N);
        matrixSAr.resize(N*N);
        for (unsigned int i=0; i<N*N;i++)
        {
            matrixA[i]=0;
            matrixAr[i]=0;
            matrixSAr[i]=0;
        }

        for (unsigned int i=0; i<N+1;i++)
        {
            if (swap==i)
                continue;
            for (unsigned int j=0; j<N+1;j++)
            {
                if (swap==j)
                    continue;
                 matrixA[count++]=Get_Fij_00_00(i,j,swap);
                 if(i==j)
                     matrixA[count-1]-=1.;
            }
        }

        for (unsigned int i=0; i<=static_cast<unsigned int>(Fickr); i++)
            NeumannStep(i, N);

        //Restablishing (I+A)
        for (unsigned int i=0; i<N;i++)
             matrixA[i*N+i]+=1.;


        count=0;
        for (unsigned int i=0; i<N; i++)
        {
            for (unsigned int j=0; j<N; j++)
            {
                matrixFick[i*N+j]=matrixSAr[count++];
            }
        }
    }

    /******************************************************************************
    Returns the Soret coefficients under request
    Internallly the default storaged values are multiplied by the molar fraction to
    avoid singularities so this function does do the division
    *******************************************************************************/
    void MuTLib_Transport::GetSoretCoef(std::vector<double>& Coef) const
    {
        Coef = SoretCoef;
        for (unsigned int i=0; i<species->GetNum(); i++)
        {
            if (species->GetMolFrac(i))
                Coef[i]/=species->GetMolFrac(i);
            else
               Coef[i]=0;
        }
    }

    /******************************************************************************
    Fick fluxes calculation
    Returns the Fick Fluxes in [kg/m2/s] with drive in [1/m]
    *******************************************************************************/
    const Eigen::VectorXd& MuTLib_Transport::GetFickFluxes(const Eigen::VectorXd& drive)
    {
        unsigned int N=species->GetNum();
        FickFluxes.resize(N);
        unsigned int iminus=0;
        FickFluxes(swap)=0.;
        for (unsigned int i=0; i<N; i++)
        {
            if (i==swap)
                continue;
            FickFluxes(i)=0.;
            unsigned jminus=0;
            for(unsigned int j=0; j<N; j++)
            {
                if (j==swap)
                        continue;
                 FickFluxes(i) -= matrixFick[iminus*(N-1)+jminus]*(drive(j)*1e-3);
                 jminus++;
            }

            FickFluxes(i) *=    species->GetMweight(i)*GTOKG*species->press/
                                species->temp/BOLTZSI/AVNUMBER*GetDij(i,swap);
            FickFluxes(i) *= 1000.;         // (kg/m2/s)
            FickFluxes[swap]-=FickFluxes(i);
            iminus++;
        }

        return FickFluxes;
    }

    /******************************************************************************
    Thermal diffusion vector calculation
    *******************************************************************************/
    void MuTLib_Transport::CalculateSoretCoef()
    {
        unsigned int N=species->GetNum();
        SoretCoef.clear();
        SoretCoef.resize(N);

        std::vector<double> F_00_10;
        F_00_10.resize((N-1)*N,0);
        unsigned int count=0;
        for (unsigned int i=0; i<N;i++)
        {
            if (i==swap)
                continue;
            for (unsigned int j=0;j<N;j++)
                F_00_10[count++]=GetFij_00_10(i,j);
        }

        std::vector<double> F_10_00;
        F_10_00.resize(N*(N-1),0);
        count=0;
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int j=0;j<N;j++)
            {
                if (j==swap)
                    continue;
                F_10_00[count++]=Get_Fij_10_00(i,j, swap);
            }
        }

        std::vector<double> F_10_10;
        F_10_10.resize(N*N,0);
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int j=0;j<N;j++)
                F_10_10[i*N+j]=GetFij_10_10(i,j);
        }

        std::vector<double> F_10_01;
        F_10_01.resize(N*N);
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int j=0;j<N;j++)
                F_10_01[i*N+j]=GetFij_10_01(i,j);
        }

        std::vector<double> F_01_10;
        F_01_10.resize(N*N);
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int j=0;j<N;j++)
                F_01_10[i*N+j]=GetFij_01_10(i,j);
        }

        //Fii_10_10
        std::vector<double> Fii_10_10(N,0.);
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int l=0; l<N;l++)
            {
                Fii_10_10[i]+=species->GetMolFrac(l)*GetDjjDil(i,i,l)*
                        GetFijl_10_10(i,i,l);
            }
        }

        //Fii_01_01
        std::vector<double> Fii_01_01(N,0.);
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int l=0; l<N;l++)
            {
                Fii_01_01[i]+=species->GetMolFrac(l)*GetDjjDil(i,i,l)*
                        GetFijl_01_01(i,i,l);
            }
        }


        a_00.resize(N-1);
        a_10.resize(N);
        a_01.resize(N);

        std::vector<double> x_01;
        std::vector<double> w;
        std::vector<double> v1;
        std::vector<double> v2;
        std::vector<double> v3;
        std::vector<double> v4;
        std::vector<double> v5;

        x_01.resize(N);
        w.resize(N);
        v1.resize(N);
        v2.resize(N);
        v3.resize(N-1);
        v4.resize(N-1);
        v5.resize(N);

        for (unsigned int i=0; i<N; i++){
            w[i]=species->GetMolFrac(i)/Fii_10_10[i];
            x_01[i]=species->GetMolFrac(i)/Fii_01_01[i]*cint[i];
        }
        MatrixVectorMult(F_10_01, x_01, v1,N,N);
        for (unsigned int i=0; i<N; i++)
            w[i]-=v1[i];

        //inv(I+A_T) iterative loop
        a_10=w;
        for (unsigned int i=1; i<=static_cast<unsigned int>(Soretr);i++)
        {
            MatrixVectorMult(F_01_10, w, v1, N, N);
            MatrixVectorMult(F_10_01, v1, v2, N, N);
            MatrixVectorMult(F_00_10, w, v3, N-1, N);
            MatrixVectorMult(matrixFick, v3, v4, N-1, N-1);
            MatrixVectorMult(F_10_00, v4, v1, N, N-1);
            MatrixVectorMult(F_10_10, w, v5, N, N);
            for (unsigned j=0; j<N; j++){
                w[j]=v1[j]+v2[j]-v5[j]+w[j];
                a_10[j]+=w[j];
            }
        }

        MatrixVectorMult(F_01_10, a_10, a_01, N, N);
        for (unsigned int i=0; i<N; i++)
            a_01[i]=x_01[i]-a_01[i];

        MatrixVectorMult(F_00_10, a_10, v3, N-1, N);
        MatrixVectorMult(matrixFick, v3, a_00, N-1, N-1);
    
        for (unsigned int i=0; i<N; i++)
            a_00[i]=-a_00[i];
        
        unsigned int iminus=0;
        SoretCoef[swap]=0;   
        for (unsigned int i=0; i<N; i++)
        {
            if (i==swap)
                continue;
            SoretCoef[i]=-a_00[iminus]*5./2.*GetDij(i,swap);
            //Compatibility relation
            SoretCoef[swap]-=species->GetMweight(i)*SoretCoef[i];
            iminus++;
        }

        SoretCoef[swap]/=species->GetMweight(swap);
    }


    /******************************************************************************
    Exact thermal diffusion vector calculation
    *******************************************************************************/
    void MuTLib_Transport::CalculateSoretCoefExact()
    {
        unsigned int N=species->GetNum();
        SoretCoef.clear();
        SoretCoef.resize(N);

        std::vector<double> F_00_10;
        F_00_10.resize((N-1)*N,0);
        unsigned int count=0;
        for (unsigned int i=0; i<N;i++)
        {
            if (i==swap)
                continue;
            for (unsigned int j=0;j<N;j++)
                F_00_10[count++]=GetFij_00_10(i,j);
        }

        std::vector<double> F_10_00;
        F_10_00.resize(N*(N-1),0);
        count=0;
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int j=0;j<N;j++)
            {
                if (j==swap)
                    continue;
                F_10_00[count++]=Get_Fij_10_00(i,j, swap);
            }
        }

        std::vector<double> F_10_10;
        F_10_10.resize(N*N,0);
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int j=0;j<N;j++)
                F_10_10[i*N+j]=GetFij_10_10(i,j);
        }

        std::vector<double> F_10_01;
        F_10_01.resize(N*N);
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int j=0;j<N;j++)
                F_10_01[i*N+j]=GetFij_10_01(i,j);
        }

        std::vector<double> F_01_10;
        F_01_10.resize(N*N);
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int j=0;j<N;j++)
                F_01_10[i*N+j]=GetFij_01_10(i,j);
        }

        //Fii_10_10
        std::vector<double> Fii_10_10(N,0.);
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int l=0; l<N;l++)
            {
                Fii_10_10[i]+=species->GetMolFrac(l)*GetDjjDil(i,i,l)*
                        GetFijl_10_10(i,i,l);
            }
        }

        //Fii_01_01
        std::vector<double> Fii_01_01(N,0.);
        for (unsigned int i=0; i<N;i++)
        {
            for (unsigned int l=0; l<N;l++)
            {
                Fii_01_01[i]+=species->GetMolFrac(l)*GetDjjDil(i,i,l)*
                        GetFijl_01_01(i,i,l);
            }
        }

        a_00.resize(N-1);
        a_10.resize(N);
        a_01.resize(N);

        std::vector<double> x_01;
        std::vector<double> w;
        std::vector<double> v1;
        std::vector<double> v2;

        x_01.resize(N);
        w.resize(N);
        v1.resize(N);
        v2.resize(N-1);

        for (unsigned int i=0; i<N; i++){
            w[i]=species->GetMolFrac(i)/Fii_10_10[i];
            x_01[i]=species->GetMolFrac(i)/Fii_01_01[i]*cint[i];
        }

        MatrixVectorMult(F_10_01, x_01, v1, N, N);
        for (unsigned int i=0; i<N; i++)
            w[i]-=v1[i];

        std::vector<double> m1;
        m1.resize(N*N);
        MatrixMatrixMult(F_10_01, F_01_10, m1, N, N, N);

        std::vector<double> m2;
        m2.resize((N-1)*N);
        MatrixMatrixMult(matrixFick, F_00_10, m2, N-1, N-1,N);


        std::vector<double> m3;
        m3.resize(N*N);
        MatrixMatrixMult(F_10_00, m2, m3, N, N-1,N);
        for (unsigned int i=0; i<N*N;i++)
            m1[i]=F_10_10[i]-m3[i]-m1[i];

        int NN=static_cast<int>(N);
        int LDA=NN;
        std::vector<int> IPIV(N);

        int INFO1 = LAPACKE_dgetrf(CblasRowMajor, NN, NN, &m1.front(), LDA, &IPIV.front());
        int INFO2 = LAPACKE_dgetri(CblasRowMajor, NN, &m1.front(), LDA, &IPIV.front());

        std::vector<double> m4;
        m4.resize((N-1)*N);
        MatrixMatrixMult(F_00_10,m1,m4,N-1,N,N);

        std::vector<double> m5;
        m5.resize((N-1)*N);
        MatrixMatrixMult(matrixFick,m4,m5,N-1,N-1,N);

        MatrixVectorMult(m1, w, a_10, N, N);

        MatrixVectorMult(F_01_10, a_10, a_01, N, N);
        for (unsigned int i=0; i<N; i++)
            a_01[i]=x_01[i]-a_01[i];

        MatrixVectorMult(F_00_10, a_10, v2, N-1, N);
        MatrixVectorMult(matrixFick, v2, a_00, N-1, N-1);

        for (unsigned int i=0; i<N-1; i++)
            a_00[i]=-a_00[i];

        unsigned int iminus=0;
        SoretCoef[swap]=0;
        for (unsigned int i=0; i<N; i++)
        {
            if (i==swap)
                continue;
            SoretCoef[i]=-a_00[iminus]*5./2.*GetDij(i,swap);
            //Compatibility relation
            SoretCoef[swap]-=species->GetMweight(i)*SoretCoef[i];
            iminus++;
        }

        SoretCoef[swap]/=species->GetMweight(swap);
    }

    /******************************************************************************
    Soret fluxes calculation
    Returns the Fick Fluxes in [kg/m2/s] with gradlogT in [K/m]
    *******************************************************************************/
    const Eigen::VectorXd& MuTLib_Transport::GetSoretFluxes(const double gradLogT)
    {
        unsigned int N=species->GetNum();
        SoretFluxes.resize(N);
        for(unsigned int i=0; i<N; i++)
        {
            SoretFluxes(i)  =   -species->GetMweight(i)*GTOKG*species->press/BOLTZSI/
                                AVNUMBER/species->temp*SoretCoef[i]*(gradLogT*1.e-3);
            SoretFluxes(i) *= 1000.;    // (kg/m2/s)
        }
        return SoretFluxes;
    }

    /******************************************************************************
    Partial thermal conductivity calculation
    *******************************************************************************/
    void MuTLib_Transport::CalculateThermalConductivity()
    {
        const unsigned int N=species->GetNum();
        ThCond = 0.;
	
        //Unscaling coefficients and adding contributions
        for (unsigned int i=0; i<N; i++)
            ThCond += 25/4.*(a_10[i]+a_01[i]*cint[i])/species->temp*GetDij(i,i);
        ThCond *= species->press;
    }

    /******************************************************************************
    Generic matrix-vector multiplication
    *******************************************************************************/
    void MuTLib_Transport::MatrixVectorMult(const std::vector<double>& mat,
                                     const std::vector<double>& vec,
                                     std::vector<double>& out,
                                     const unsigned int M,
        const unsigned int N)
    {
        for (unsigned int i=0; i<M; i++)
        {
            out[i]=0;
            for (unsigned int j=0; j<N; j++)
                out[i]+=mat[i*N+j]*vec[j];
        }
    }

    /******************************************************************************
    Generic matrix-matrix multiplication
    *******************************************************************************/
    void MuTLib_Transport::MatrixMatrixMult(   const std::vector<double>& mat1,
                                        const std::vector<double>& mat2,
                                        std::vector<double>& out,
                                        const unsigned int M,
                                        const unsigned int N,
                                        const unsigned int L)
    {

        for (unsigned int i=0; i<M; i++)
        {
            for (unsigned int j=0; j<L; j++)
            {
                out[i*L+j]=0;
                for (unsigned int k=0; k<N; k++)
                    out[i*L+j]+=mat1[i*N+k]*mat2[k*L+j];
            }
        }
    }

}
