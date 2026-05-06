volScalarField rAU("rAU", 1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(HbyA)
     // + MRF.zeroFilter(rAUf*fvc::ddtCorr(rho, U, phi))
    )
  + phig
);

//MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
//constrainPressure(p_rgh, rho, U, phiHbyA, rAUf, MRF);
constrainPressure(p_rgh, rho, U, phiHbyA, rAUf);

{


    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        p_rghEqn.setReference
        (
            pRefCell,
            pRefValue
        );

        p_rghEqn.solve
        (
            p_rgh.select
            (
                (
                    oCorr == nOuterCorr-1
                 && corr == nCorr-1
                 && nonOrth == nNonOrthCorr
                )
            )
        );

        if (nonOrth == nNonOrthCorr)
        {
            phi = phiHbyA + p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            //fvOptions.correct(U);
            //K = 0.5*magSqr(U);
        }
    }

    p = p_rgh + rho*gh;

}

//pressureControl.limit(p);


