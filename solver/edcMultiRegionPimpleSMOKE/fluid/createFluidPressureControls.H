PtrList<pressureReference> pressureReferenceFluid(fluidRegions.size());

    pressureReferenceFluid.set
    (
        0,
        new pressureReference
        (
            thermo.p(),
            pimples.pimple(0).dict(),
            false
        )
    );

for(int i=1;i<fluidRegions.size();i++)
{
    pressureReferenceFluid.set
    (
        i,
        new pressureReference
        (
            thermoFluid[i].p(),
            pimples.pimple(i).dict(),
            false
        )
    );
}
