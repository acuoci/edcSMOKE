/******************************************************************************
LIBRARY MuTLib
Copyright (c) 2021 Universidad Nacional Española a Distancia UNED.
All rights reserved.
Authors:    Oscar Córdoba
            Manuel Arias-Zugasti
kineticdata.cpp implementation
Kinetic data file reading
******************************************************************************/

#include <iostream>

namespace MuTLib
{
    double delta[8] = { 0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5 };

    double Tstar[37] =
    {
        0.1,0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9,
        1.0,1.2,1.4,1.6,1.8,
        2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,9.0,
        10.0,12.0,14.0,16.0,18.0,20.0,25.0,30.0,35.0,40.0,50.0,75.0,100.0
    };

    bool MuTLib_KineticData::ReadFile(const std::string fileName) 
    {
        const unsigned int MAXLINE = 256;
        char line[MAXLINE];

        FILE* ptr;
        ptr = fopen(fileName.c_str(), "r");
        if (ptr == nullptr)
        {
            std::cout << "Error: File " << fileName << " not found..." << std::endl;
            return false;
        }

        Clear();

        int count = 0;
        while (fgets(line, MAXLINE, ptr)) 
        {
            if (line[0] == '!')
                continue;
            count++;
        }

        char    name_[MAXLINE];
        int     index_;
        double  potential_;
        double  diameter_;
        double  dipole_;
        double  polarizability_;
        double  collision_;

        rewind(ptr);
        while (fgets(line, MAXLINE, ptr)) 
        {
            if (line[0] == '!')
                continue;
            sscanf(line, "%s %i %lf %lf %lf %lf %lf", name_, &index_, &potential_, &diameter_, &dipole_, &polarizability_, &collision_);

            name.push_back(name_);
            index.push_back(index_);
            potential.push_back(potential_);
            diameter.push_back(diameter_);
            dipole.push_back(dipole_);
            polarizability.push_back(polarizability_);
            collision.push_back(collision_);
        }

        fclose(ptr);

        return true;
    }

    bool MuTLib_KineticData::ImportCoefficients(   const std::vector<std::string>& names, const std::vector<int>& shape_factor, 
                                            const std::vector<double>& epsilon_over_kb, const std::vector<double>& sigma, 
                                            const std::vector<double>& mu, const std::vector<double>& alfa, 
                                            const std::vector<double>& zRot298)
    {
        Clear();

        name = names;
        index = shape_factor;
        potential = epsilon_over_kb;
        diameter = sigma;
        dipole = mu;
        polarizability = alfa;
        collision = zRot298;

        return true;
    }

    /******************************************************************************
    Cleaning function
    *******************************************************************************/
    void MuTLib_KineticData::Clear() 
    {
        name.clear();
        index.clear();
        potential.clear();
        dipole.clear();
        polarizability.clear();
        collision.clear();
    }

    /******************************************************************************
    Returns the species index from the storage vector
    -1 is returned in case of name not found
    *******************************************************************************/
    int MuTLib_KineticData::CheckName(const std::string specie) 
    {
        for (unsigned int i = 0; i < name.size(); i++) 
        {
            if (name[i] == specie)
                return static_cast<int>(i);
        }

        return -1;
    }

    /******************************************************************************
    Interpolates in Tables
    Inputs:
        delta
        T*
        pointer to table
    Return the interpolated value
    *******************************************************************************/
    double MuTLib_KineticData::Interpol(const double delta_, const double Tstar_, const double value[37][8]) 
    {
        int     x;
        int     y;
        double   xx1;
        double   xx2;

        for (x = 0; x < 8; x++) 
        {
            if (delta_ < delta[x])
                break;
        }

        for (y = 0; y < 37; y++) 
        {
            if (Tstar_ < Tstar[y])
                break;
        }

        x--; 
        y--;

        double m;

        if (x == 7) 
        {
            xx1 = value[y][x];
            if (y == 36)
                return xx1;
            xx2 = value[y + 1][x];
        }
        else 
        {

            m = (value[y][x + 1] - value[y][x]) / (delta[x + 1] - delta[x]);
            xx1 = value[y][x] + m * (delta_ - delta[x]);

            if (y == 36)
                return xx1;
            m = (value[y + 1][x + 1] - value[y + 1][x]) / (delta[x + 1] - delta[x]);
            xx2 = value[y + 1][x] + m * (delta_ - delta[x]);
        }
        m = (xx2 - xx1) / (Tstar[y + 1] - Tstar[y]);

        return xx1 + m * (Tstar_ - Tstar[y]);
    }
}
