/******************************************************************************
LIBRARY MuTLib
Copyright (c) 2021 Universidad Nacional Española a Distancia UNED.
All rights reserved.
Authors:    Oscar Córdoba
            Manuel Arias-Zugasti
kineticdata.h header definitions
******************************************************************************/

#ifndef MUTLIB_KINETICDATA_H
#define MUTLIB_KINETICDATA_H

#include <vector>
#include <string>

// Preliminary definitions
#define BOLTZ       1.38064852e-16    // CGS
#define BOLTZSI     1.38064852E-23    // SI
#define DEBYE       1e-18             // CGS
#define ARMSTRONG   1e-8              // cm
#define AVNUMBER    6.022140857e23    // Avogadro's number
#define GTOKG       0.001             // g to kg units conversion

namespace MuTLib
{
    class MuTLib_KineticData 
    {
        public:

            std::vector<std::string>    name;

            std::vector<int>            index;

            std::vector<double>         potential;

            std::vector<double>         diameter;

            std::vector<double>         dipole;

            std::vector<double>         polarizability;

            std::vector<double>         collision;

        public:

            void Clear();

            int CheckName(const std::string);

            double Interpol(const double, const double, const double[37][8]);

            bool ReadFile(const std::string name);

            bool ImportCoefficients(    const std::vector<std::string>& names, const std::vector<int>& shape_factor,
                                        const std::vector<double>& epsilon_over_kb, const std::vector<double>& sigma,
                                        const std::vector<double>& mu, const std::vector<double>& alfa,
                                        const std::vector<double>& zRot298); 
    };

}

#include "MuTLib_KineticData.hpp"

#endif // MUTLIB_KINETICDATA_H
