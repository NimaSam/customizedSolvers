if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        tmp<fvScalarMatrix> TEqn
        (
            fvm::ddt(T)
          - fvm::laplacian(K/(rho*Cp), T, "laplacian(K,T)")
        );

        TEqn().relax();

        TEqn().solve(mesh.solver(T.select(finalIter)));
    }
}


Info<< "Min/max T:" << min(T).value() << ' ' << max(T).value() << endl;

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
