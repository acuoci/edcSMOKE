/******************************************************************************
LIBRARY MuTLib
Copyright (c) 2021 Universidad Nacional Española a Distancia UNED.
All rights reserved.
Authors:    Oscar Córdoba
            Manuel Arias-Zugasti
species.cpp implementation
******************************************************************************/

#include <string.h>

namespace MuTLib
{

    /******************************************************************************
    Cleaning function
    *******************************************************************************/
    void MuTLib_Species::Clear() 
    {
        name.clear();
        mweight.clear();
        Kposition.clear();
        molfrac.clear();
    }

    /******************************************************************************
    Reallocating function
    *******************************************************************************/
    void MuTLib_Species::Resize(const int size) 
    {
        name.resize(static_cast<unsigned int>(size));
        mweight.resize(static_cast<unsigned int>(size));
        Kposition.resize(static_cast<unsigned int>(size));
        molfrac.resize(static_cast<unsigned int>(size));
    }
}

