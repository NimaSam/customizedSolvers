{

	fvScalarMatrix rhoEqn
	(
		fvm::ddt(rho)
		+ fvc::div(rhoPhi)
	);

	rhoEqn.relax();
	rhoEqn.solve();

	Info << "rhoAve = "
	<< rho.weightedAverage(mesh.V()).value()
	<< "  Min(rho) = " << min(rho).value()
	<< "  Max(rho) = " << max(rho).value()
	<< endl;

}

