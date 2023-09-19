{
    dimensionedScalar vol = 4/3*Foam::constant::mathematical::pi*pow(rDET,3);
    dimensionedScalar ddtvol = 4*Foam::constant::mathematical::pi*pow(vcj,3)*pow(runTime,2);

    Sp=rhoExp*(QDET/vol)*ddtvol*pos(rDET-mag(mesh.C()-centre))*pos(-(runTime-t0)*(runTime-tDET));
	
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, e) + fvm::div(phi, e)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + fvc::div(fvc::absolute(phi/fvc::interpolate(rho), U), p, "div(phiv,p)")
      - fvm::laplacian(turbulence->alphaEff(), e)
     ==
        fvOptions(rho, e)
      + Sp
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(e);

    thermo.correct();
}
