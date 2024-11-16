// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

UEqns[i] =
(
	fvm::ddt(rho, U) + fvm::div(phi, U)
	+ MRF.DDt(rho, U)
	+ turbulence.divDevTau(U)
	==
        rho*g + 
	fvModels.source(rho, U)
);
fvVectorMatrix& UEqnNonReacting = UEqns[i].ref();

UEqnNonReacting.relax();

fvConstraints.constrain(UEqnNonReacting);

if (pimple.momentumPredictor())
{
	solve(UEqnNonReacting == -fvc::grad(p));

    	fvConstraints.constrain(U);
    	K = 0.5*magSqr(U);
}

Info << " Umin: " << min(U).value() << " Umax: " << max(U).value() << endl;
