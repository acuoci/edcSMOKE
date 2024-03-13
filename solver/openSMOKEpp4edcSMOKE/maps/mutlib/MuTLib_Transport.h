/******************************************************************************
LIBRARY MuTLib
Copyright (c) 2021 Universidad Nacional Española a Distancia UNED.
All rights reserved.
Authors:    Oscar Córdoba
            Manuel Arias-Zugasti
transport.h header definitions
******************************************************************************/
#ifndef MUTLIB_TRANSPORT_H
#define MUTLIB_TRANSPORT_H


// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

#include <vector>
#include <string>

#include "MuTLib_Species.h"
#include "MuTLib_KineticData.h"
#include "MuTLib_ThermoData.h"

namespace MuTLib
{

    /******************************************************************************
    Enum to matricial parameters
    *******************************************************************************/
    enum TR_VAR 
    {
        TR_MOLWEIGHT,
        TR_LENNARD,
        TR_COLLDIAM,
        TR_TSTAR,
        TR_THCOEF,
        TR_DIPOLESTAR,
        TR_OMEGA11,
        TR_OMEGA22,
        TR_ASTAR,
        TR_BSTAR,
        TR_CSTAR,
        TR_ESTAR,
        TR_FSTAR,
        TR_ZETA,
        TR_END
    };


    /******************************************************************************
    Class transport. main container to coefficients, fluxes and
    conductivity members
    *******************************************************************************/
    class MuTLib_Transport 
    {

        //Pointers to main TransLib classes
        MuTLib_Species* species;
        MuTLib_KineticData* kindata;
        MuTLib_ThermoData* thermdata;

        //General matricial parameters
        std::vector<std::vector<double>>   varValues;

        //Vector parameters
        std::vector<double>                 mass;
        std::vector<double>                 dens;
        std::vector<double>                 cint;
        std::vector<double>                 crot;
        std::vector<double>                 deltap;
        std::vector<double>                 chi;
        std::vector<double>                 mu;

        //Transport output calculations
        std::vector<double>                 matrixFick;
        std::vector<double>                 FickCoef;
        std::vector<double>                 matrixSoret;
        std::vector<double>                 SoretCoef;
        std::vector<double>                 a_00;
        std::vector<double>                 a_10;
        std::vector<double>                 a_01;
        Eigen::VectorXd                     FickFluxes;
        Eigen::VectorXd                     SoretFluxes;

        double                              ThCond;

        //Algorithm flags
        unsigned int    swap;           // Index for the dominant species
        bool            exact;          // Flag if Fick diffusion coefficients are provided
        int             Fickr;          // Number of Neumann iterations in Fick  inversion
        int             FickM;          // M species for Fick diffusion 1+M model
        double          FickThr;        // Threshold to select M species in 1+M model
        int             Soretr;         // Number of Neumann iterations in Soret calculations
        int             SwapPolicy;     // Policy for selecting the swap species

        //Auxiliar members
        std::vector<double>                 matrixA;
        std::vector<double>                 matrixAr;
        std::vector<double>                 matrixSAr;
        std::vector<unsigned int>           index;

    public:

        MuTLib_Transport();

        void SetUp(MuTLib_Species* species_, MuTLib_KineticData* kindata_, MuTLib_ThermoData* thermdata);

        const std::vector<double>& GetVarValues(const int index) const
        {
            return (varValues[static_cast<unsigned int>(index)]);
        }

        //Methods for multi-species kinetic calculations
        void CalculateKinetic();

        void CalculateMolecularWeight(const unsigned int i, const unsigned  int j);

        void CalculateLennardJones(const unsigned int i, const unsigned  int j);

        void CalculateCollisionDiameter(const unsigned int i, const unsigned  int j);

        void CalculateTstar(const unsigned int i, const unsigned  int j);

        void CalculateDipoleStar(const unsigned int i, const unsigned  int j);

        void CalculateOmega11(const unsigned int i, const unsigned  int j);

        void CalculateOmega22(const unsigned int i, const unsigned  int j);

        void CalculateAstar(const unsigned int i, const unsigned  int j);

        void CalculateBstar(const unsigned int i, const unsigned  int j);

        void CalculateCstar(const unsigned int i, const unsigned  int j);

        void CalculateEstar(const unsigned int i, const unsigned  int j);

        void CalculateFstar(const unsigned int i, const unsigned  int j);

        void CalculateZeta(const unsigned int i, const unsigned  int j);

        double ZrotScale(const unsigned int i, const double temp);


        //Binary Diffusion coefficients
        double GetDij(const unsigned int i, const unsigned int j) const;

        //Binary Diffusion coefficients scaling ratios
        double GetDjNDil(const unsigned int i, const unsigned int j, const unsigned int l, const unsigned int N) const;

        double GetDjjDil(const unsigned int i, const unsigned int j, const unsigned int l) const;

        //System of equations assignation methods
        double GetFijl_00_00(const unsigned int i, const unsigned int j, const unsigned int l) const;

        double GetFijl_00_10(const unsigned int i, const unsigned int j, const unsigned int l) const;

        double GetFijl_10_00(const unsigned int i, const unsigned int j, const unsigned int l) const;

        double GetFijl_10_10(const unsigned int i, const unsigned int j, const unsigned int l) const;

        double GetFijl_10_01(const unsigned int i, const unsigned int j, const unsigned int l) const;

        double GetFijl_01_10(const unsigned int i, const unsigned int j, const unsigned int l) const;

        double GetFijl_01_01(const unsigned int i, const unsigned int j, const unsigned int l) const;

        //Scaled system of equations assignation methods
        double GetFij_00_00(const unsigned int i, const unsigned int j, const unsigned int N) const;

        double GetFij_00_10(const unsigned int i, const unsigned int j) const;

        double GetFij_10_00(const unsigned int i, const unsigned int j, const unsigned int N) const;

        double GetFij_10_10(const unsigned int i, const unsigned int j) const;
         
        double GetFij_10_01(const unsigned int i, const unsigned int j) const;

        double GetFij_01_10(const unsigned int i, const unsigned int j) const;

        double GetFij_01_01(const unsigned int i, const unsigned int j) const;

        //Compatibility equation substitution methods
        double Get_Fij_00_00(const unsigned int i, const unsigned int j, const unsigned int N) const;

        double Get_Fij_10_00(const unsigned int i, const unsigned int j, const unsigned int N) const;


        //Main driver
        void Calculate(
            const bool            exact = false,      // boolean to set exact Fick diffusion coefficients
            const unsigned int    Fickr = 1,          // Fick diffusion Neumann steps
            const int             FickM = 0,          // M number of species in model 1+M for Fick diffusion coefficients
            const double          FickThr = 0.1,      // Molar fraction threshold in model 1+M for Fick diffusion coefficients activated with FickM<0
            const unsigned int    Soretr = 1,         // Soret diffusion Neumann steps
            const int             SwapPolicy = -2);   // Policy for selecting the swap species

        void CalculateExact();

        //Transport calculation methods
        void                CalculateMatrixFick();

        void                CalculateMatrixFickExact();

        void                UnscaleMatrixFick();

        void                CalculateSoretCoef();

        void                CalculateSoretCoefExact();

        void                CalculateThermalConductivity();

        void                NeumannStep(unsigned int r, unsigned int N);

        const std::vector<double>& GetMatrixNeumann() const
        { 
            return matrixSAr;
        }

        const std::vector<double>& GetMatrixA() const
        { 
            return matrixA; 
        }

        void  CalculateModel1Fick();

        //Debugging method
        void PrintMatrix(const std::vector<double>& matrix, const unsigned int m, const unsigned int n);

        //Recovering methods
        const std::vector<double>& GetMatrixFick() const
        {
            return matrixFick;
        }
        const std::vector<double>& GetFickCoef()
        {
            UnscaleMatrixFick();

            return FickCoef;
        }

	void GetSoretCoef(std::vector<double>& Coef) const;

        const Eigen::VectorXd& GetFickFluxes(const Eigen::VectorXd& drive);

        const Eigen::VectorXd& GetSoretFluxes(const double gradLogT);

        double GetOmega11(const unsigned int i, const unsigned int j) const
        {
            return varValues[TR_OMEGA11][i * species->GetNum() + j];
        }

        double GetAstar(const unsigned int i, const unsigned int j) const
        {
            return varValues[TR_ASTAR][i * species->GetNum() + j];
        }

        double GetBstar(const unsigned int i, const unsigned int j) const
        {
            return varValues[TR_BSTAR][i * species->GetNum() + j];
        }

        double GetCstar(const unsigned int i, const unsigned int j) const
        {
            return varValues[TR_CSTAR][i * species->GetNum() + j];
        }

        double GetCollDiam(const unsigned int i, const unsigned int j) const
        {
            return varValues[TR_COLLDIAM][i * species->GetNum() + j];
        }

        double GetThCond() const
        {
            return ThCond;
        }

        void MatrixVectorMult(const std::vector<double>& mat, const std::vector<double>& vec,
            std::vector<double>& out, const unsigned int M, const unsigned int N);

        void MatrixMatrixMult(const std::vector<double>& mat1, const std::vector<double>& mat2,
            std::vector<double>& out, 
            const unsigned int M, const unsigned int N, const unsigned int L);
    };
}

#include "MuTLib_Transport.hpp"

#endif // MUTLIB_TRANSPORT_H
