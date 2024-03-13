/******************************************************************************
LIBRARY MuTLib
Copyright (c) 2021 Universidad Nacional Española a Distancia UNED.
All rights reserved.
Authors:    Oscar Córdoba
            Manuel Arias-Zugasti
thermodata.h header definitions
******************************************************************************/
#ifndef MUTLIB_THERMODATA_H
#define MUTLIB_THERMODATA_H

#include <vector>
#include <string>

namespace MuTLib
{
    /******************************************************************************
    Class ThermoData
    Main container to JANAF thermodynamic properties interpolation
    ******************************************************************************/
    class MuTLib_ThermoData 
    {
        private:

            std::vector<std::string>           name;

            std::vector<std::vector<double>>   Trange;

            std::vector<std::vector<double>>   CoefL;

            std::vector<std::vector<double>>   CoefH;


        public:

            MuTLib_ThermoData();

            void    Clear();

            int     CheckName(const std::string);

            bool    ReadFile(const std::string name);

            bool    ImportCoefficients(const std::vector<std::string>& names, const std::vector<std::vector<double>>& CoeffLowT, const std::vector<std::vector<double>>& CoeffHighT, const std::vector<std::vector<double>>& T);

            double  GetCpR(const unsigned int index, const double temp);

            double  GetHRT(const unsigned int index, const double temp);

            double  GetSR(const unsigned int index, const double temp);
    };
}

#include "MuTLib_ThermoData.hpp"

#endif // MUTLIB_THERMODATA_H
