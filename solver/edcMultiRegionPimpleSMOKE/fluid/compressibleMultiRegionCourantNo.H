    scalar CoNum = -great;

    // From reacting region
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                fluidRegions[0],
                runTime,
                rho,
                phi
            ),
            CoNum
        );
    }

    // From remaining (non-reacting) regions
    for(int regionI=1;regionI<fluidRegions.size();regionI++)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                fluidRegions[regionI],
                runTime,
                rhoFluid[regionI],
                phiFluid[regionI]
            ),
            CoNum
        );
    }
