 volScalarField muEf //add
        (
               "muEf",
               rho*nu
             + rho*turbulence->nut()
        );


        fvVectorMatrix UEqn
        (
            fvm::ddt((1.0-alphap)*rho,U)
          + fvm::div(fvc::interpolate(rho)*phi, U)
          - fvm::laplacian(fvc::interpolate(muEf), U)
          + fvm::Sp(muEf*alphap/kpp,U) //add artificial viscosity like darcy term in porous media
            ==
            kinematicCloud.SU(U)
        );

        UEqn.relax();

        if (piso.momentumPredictor()) //change
        {
        	solve(UEqn == -fvc::grad(p));
        }
