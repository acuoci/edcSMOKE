    scalar CoNum = -great;

    forAll(fluidRegions, regionI)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                fluidRegions[regionI],
                runTime,
                rho,
                phi
            ),
            CoNum
        );
    }
