tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

const volScalarField& m2= mesh.lookupObject<volScalarField>("moment.2");
Sg = growth->Kg(m2)*rho*ka*m2/2.0;
//Info << "Sg: " << Sg <<endl;

{
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(rho*(crystal.Dm() + turbulence->nut()/Sct), Yi)
              //+ fvOptions(rho, Yi)

            );

            if (Y[i].name()==growthSource){
            	Info <<"growth source for: " << Y[i].name() << endl;
            	YiEqn=YiEqn-Sg;
            }else if(Y[i].name()==growthSink) {
            	Info <<"growth sink for: " << Y[i].name() << endl;
            	YiEqn =YiEqn+Sg;
            }


            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));


            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
            Info<< Y[i].name() << " min/max:   " << min(Y[i]).value() << ", "
                   << max(Y[i]).value() << endl;
        }
        else
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
