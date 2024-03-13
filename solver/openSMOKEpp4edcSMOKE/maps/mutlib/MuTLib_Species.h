/******************************************************************************
LIBRARY MuTLib
Copyright (c) 2021 Universidad Nacional Española a Distancia UNED.
All rights reserved.
Authors:    Oscar Córdoba
            Manuel Arias-Zugasti
species.h header definitions
******************************************************************************/
#ifndef MUTLIB_SPECIES_H
#define MUTLIB_SPECIES_H

#include <iostream>
#include <vector>
#include <string>
#include "MuTLib_KineticData.h"

namespace MuTLib
{

    /******************************************************************************
    Class Species
    Main container combustion state of the mixture
    ******************************************************************************/
    class MuTLib_Species 
    {

    private:

        std::vector<std::string>        name;

        std::vector<double>             mweight;    // Molar weight

        std::vector<double>             molfrac;    // Molar fraction

        std::vector<unsigned int>       Kposition;  // Kinetic data position

        std::vector<unsigned int>       Tposition;  // Thermo data position


    public:

        double                      temp;

        double                      press;

        void Clear();

        void Resize(const int size);

        void SetName(const unsigned int index, const std::string name_)
        {
            name[index] = name_;
        }

        void SetKIndex(const unsigned int index, const unsigned int Kposition_)
        {
            Kposition[index] = Kposition_;
        }

        void SetTIndex(const unsigned int index, const unsigned int Tposition_)
        {
            Tposition[index] = Tposition_;
            std::cout << index << " " << Tposition_ << std::endl;
        }

        void SetMweight(const unsigned int index, const double mweight_)
        {
            mweight[index] = mweight_;
        }

        void SetMolFrac(const unsigned int index, const double molfrac_)
        {
            molfrac[index] = molfrac_;
        }

        unsigned int GetNum() const
        {
            return name.size();
        }

        void AddSpecie(const std::string name, const double mweight, const double molfrac, const int indexK, const int indexT)
        {
            AddName(name);
            AddMweight(mweight);
            AddMolFrac(molfrac);
            AddKPosition(static_cast<unsigned int>(indexK));
            AddTPosition(static_cast<unsigned int>(indexT));
        }

        std::string GetName(const unsigned int index) const
        {
            return name[index];
        }
        
        void AddName(const std::string name_)
        {
            name.push_back(name_);
        }

        const std::vector<std::string>& GetNames() const
        {
            return name;
        }

        double GetMweight(const unsigned int index) const
        {
            return mweight[index];
        }

        void AddMweight(const double mweight_)
        {
            mweight.push_back(mweight_);
        }

        double GetMolFrac(const unsigned int index) const
        {
            return molfrac[index];
        }

        void AddMolFrac(const double molfrac_)
        {
            molfrac.push_back(molfrac_);
        }

        unsigned int GetKPosition(const unsigned int index) const
        {
            return Kposition[index];
        }

        void AddKPosition(const double position_)
        {
            Kposition.push_back(position_);
        }

        unsigned int GetTPosition(const unsigned int index) const
        {
            return Tposition[index];
        }

        void AddTPosition(const double position_)
        {
            Tposition.push_back(position_);
        }

    };
}

#include "MuTLib_Species.hpp"

#endif // MUTLIB_SPECIES_H
