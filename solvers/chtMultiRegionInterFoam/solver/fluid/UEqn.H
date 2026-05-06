    // Solve the Momentum equation

    //MRF.correctBoundaryVelocity(U);

    UEqn =
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      //+ MRF.DDt(rho, U)
       + turbulence.divDevRhoReff(rho, U) 
      //==
      //  fvOptions(rho, U)
    );

    UEqn.relax();

    //fvOptions.constrain(UEqn);

    if (momentumPredictor)
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            ),
            U.select(finalIter)
        );

        //sfvOptions.correct(U);
        //K = 0.5*magSqr(U);
    }
