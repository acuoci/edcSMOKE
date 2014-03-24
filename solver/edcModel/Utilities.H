
void  WriteScalarField(fvMesh& mesh, scalarField& X, word fieldName)
{
    //output only if it is write time
    if ( mesh.time().outputTime() )
    {  
        // as paraView only undestands volScalarFields, make a volScalarField out of the scalarField and write it
        volScalarField volX
        (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName
        );

        volX.internalField() = X;
        volX.correctBoundaryConditions();
        volX.write();    
    }
}
