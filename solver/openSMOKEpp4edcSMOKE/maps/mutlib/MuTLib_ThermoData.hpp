/******************************************************************************
LIBRARY MuTLib
Copyright (c) 2021 Universidad Nacional Española a Distancia UNED.
All rights reserved.
Authors:    Oscar Córdoba
            Manuel Arias-Zugasti
thermodata.cpp implementation
JANAF thermodynamic data reading
******************************************************************************/

#include <vector>
#include <iostream>
#include <string>
#include <string.h>
#include <math.h>

namespace MuTLib
{
    /******************************************************************************
    Class constructor
    *******************************************************************************/
    MuTLib_ThermoData::MuTLib_ThermoData()
    {
    }

    /******************************************************************************
    Cleaning function
    *******************************************************************************/
    void MuTLib_ThermoData::Clear()
    {
        for (unsigned int i = 0; i < CoefL.size(); i++)
            CoefL[i].clear();
        CoefL.clear();
        for (unsigned int i = 0; i < CoefH.size(); i++)
            CoefH[i].clear();
        CoefH.clear();
        name.clear();
    }

    /******************************************************************************
    Returns the species index from the storage vector
    -1 is returned in case of name not found
    *******************************************************************************/
    int MuTLib_ThermoData::CheckName(const std::string specie)
    {
        for (unsigned int i = 0; i < name.size(); i++)
        {
            if (name[i] == specie)
                return static_cast<int>(i);
        }

        return -1;
    }

    /******************************************************************************
    Reads a file and stores thermodynamic data in  vectors
    *******************************************************************************/
    bool MuTLib_ThermoData::ReadFile(const std::string fileName)
    {
        const unsigned int MAXLINE = 256;
        char line[MAXLINE];
        char word[MAXLINE];

        std::vector<double> range(3);
        std::vector<double> lower(7);
        std::vector<double> higher(7);

        FILE* ptr;
        ptr = fopen(fileName.c_str(), "r");
        if (ptr == nullptr)
        {
            std::cout << "Error: File " << fileName << " not found..." << std::endl;
            return false;
        }

        Clear();


        if (!fgets(line, MAXLINE, ptr))
        {
            std::cout << "Error: File empty..." << std::endl;
            return false;
        }

        sscanf(line, "%s", word);

        while (strcmp(word, "THERMO")) 
        {
            if (fgets(line, MAXLINE, ptr) == NULL)
            {
                std::cout << "Error code 000..." << std::endl;
                return false;
            }
            sscanf(line, "%s", word);
        }

        if (!fgets(line, MAXLINE, ptr))
        {
            std::cout << "Error code 001..." << std::endl;
            return false;
        }

        while (fgets(line, MAXLINE, ptr) != NULL) 
        {
            sscanf(line, "%s", word);
            if (strcmp(word, "END") == 0)
                break;

            while (line[0] == '!') 
            {
                if (!fgets(line, MAXLINE, ptr))
                {
                    std::cout << "Error code 002..." << std::endl;
                    return false;
                }
            }
            sscanf(line, "%18s", word);
            sscanf(&(line[45]), "%lf%lf%lf", &range.front(), &range.front() + 1, &range.front() + 2);

            Trange.push_back(range);
            name.push_back(word);

            if (!fgets(line, MAXLINE, ptr))
            {
                std::cout << "Error code 003..." << std::endl;
                return false;
            }

            sscanf(line, "%lf%lf%lf%lf%lf", &higher.front() + 0, &higher.front() + 1, &higher.front() + 2, &higher.front() + 3, &higher.front() + 4);

            if (!fgets(line, MAXLINE, ptr))
            {
                std::cout << "Error code 004..." << std::endl;
                return false;
            }

            sscanf(line, "%lf%lf%lf%lf%lf", &higher.front() + 5, &higher.front() + 6, &lower.front() + 0, &lower.front() + 1, &lower.front() + 2);

            if (!fgets(line, MAXLINE, ptr))
            {
                std::cout << "Error code 001..." << std::endl;
                return false;
            }

            sscanf(line, "%lf%lf%lf%lf", &lower.front() + 3, &lower.front() + 4, &lower.front() + 5, &lower.front() + 6);

            CoefL.push_back(lower);
            CoefH.push_back(higher);
        }

        fclose(ptr);

        return true;
    }

    bool MuTLib_ThermoData::ImportCoefficients(const std::vector<std::string>& names, const std::vector<std::vector<double>>& CoeffLowT, const std::vector<std::vector<double>>& CoeffHighT, const std::vector<std::vector<double>>& T)
    {
        Clear();

        name = names;
        CoefL = CoeffLowT;
        CoefH = CoeffHighT;
        Trange = T;

        return true;
    }

    /******************************************************************************
    Returns specific heat Cp/R evaluating polynomial expressions
    *******************************************************************************/
    double MuTLib_ThermoData::GetCpR(const unsigned int index, const double temp)
    {
        bool low = true;
        if (temp > Trange[index][2])
            low = false;
        if (low) 
        {
            return CoefL[index][0] + CoefL[index][1] * temp + CoefL[index][2] * temp *
                temp + CoefL[index][3] * temp * temp * temp + CoefL[index][4] * temp *
                temp * temp * temp;
        }
        else 
        {
            return CoefH[index][0] + CoefH[index][1] * temp + CoefH[index][2] * temp *
                temp + CoefH[index][3] * temp * temp * temp + CoefH[index][4] * temp *
                temp * temp * temp;
        }
    }

    /******************************************************************************
    Returns enthalpy H/RT evaluating polynomial expressions
    *******************************************************************************/
    double MuTLib_ThermoData::GetHRT(const unsigned int index, const double temp)
    {
        bool low = true;
        if (temp > Trange[index][2])
            low = false;

        if (low) 
        {
            return CoefL[index][0] + CoefL[index][1] / 2. * temp + CoefL[index][2] / 3 *
                temp * temp + CoefL[index][3] / 4. * temp * temp * temp + CoefL[index][4] / 5. *
                temp * temp * temp * temp + CoefL[index][5] / temp;
        }
        else 
        {
            return CoefH[index][0] + CoefH[index][1] / 2. * temp + CoefH[index][2] / 3 *
                temp * temp + CoefH[index][3] / 4. * temp * temp * temp + CoefH[index][4] / 5. *
                temp * temp * temp * temp + CoefH[index][5] / temp;
        }
    }

    /******************************************************************************
    Returns entropy S/R evaluating polynomial expressions
    *******************************************************************************/
    double MuTLib_ThermoData::GetSR(const unsigned int index, const double temp)
    {
        bool low = true;
        if (temp > Trange[index][2])
            low = false;

        if (low) 
        {
            return CoefL[index][0] * log(temp) + CoefL[index][1] * temp + CoefL[index][2] /
                2. * temp * temp + CoefL[index][3] / 3. * temp * temp * temp + CoefL[index][4] /
                4. * temp * temp * temp * temp + CoefL[index][6];
        }
        else 
        {
            return CoefH[index][0] * log(temp) + CoefH[index][1] * temp + CoefH[index][2] /
                2. * temp * temp + CoefH[index][3] / 3. * temp * temp * temp + CoefH[index][4] /
                4. * temp * temp * temp * temp + CoefH[index][6];
        }
    }

}
