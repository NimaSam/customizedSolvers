volScalarField rUA = 1.0/UEqn.A();
surfaceScalarField rUAf("rUAf", fvc::interpolate(rUA)); //add

            U = rUA*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf())
                  + fvc::interpolate(rho*rUA)*fvc::ddtCorr(U, phi); //TODO: check!
                //+ fvc::ddtPhiCorr(rUA,rho, U, phi); //replace

            adjustPhi(phi, U, p);

            //for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            while (piso.correctNonOrthogonal()) //replace
            {

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rUA, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                //if (nonOrth == nNonOrthCorr)
                if (piso.finalNonOrthogonalIter()) //replace
                {
                    phi -= pEqn.flux();
                }
            }

	    U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
