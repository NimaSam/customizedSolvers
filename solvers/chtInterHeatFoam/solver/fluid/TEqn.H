{

	K = alpha1*K1 + (1.0 - alpha1)*K2;
	rhoC = alpha1*rho1*C1 + (1.0 - alpha1)*rho2*C2;


	fvScalarMatrix TEqn
	(
		  fvm::ddt(T)
		+ fvm::div(phi, T)
		- fvm::Sp(fvc::div(phi), T)//add
		- fvm::laplacian( K/rhoC , T,"laplacian(alphaEff,T)")
	);

	TEqn.relax();
    TEqn.solve();


	Info << "TAve = "
	<< T.weightedAverage(mesh.V()).value()
	<< "  Min(T) = " << min(T).value()
	<< "  Max(T) = " << max(T).value()
	<< endl;

}

