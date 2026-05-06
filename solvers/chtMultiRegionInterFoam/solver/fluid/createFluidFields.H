// Initialise fluid field pointer lists
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> TFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size());

PtrList<thermalImmiscibleIncompressibleTwoPhaseMixture> mixtureFluid(fluidRegions.size());
PtrList<volScalarField> alpha1Fluid(fluidRegions.size());
PtrList<volScalarField> alpha2Fluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
//-thermal conductivity for fluid regions
PtrList<volScalarField> kappaFluid(fluidRegions.size());



PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<fvVectorMatrix> UEqFluid(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());


typedef incompressibleInterPhaseTransportModel
    <
        thermalImmiscibleIncompressibleTwoPhaseMixture
    > transportModelType;
PtrList<transportModelType> turbulenceFluid(fluidRegions.size());

//PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
//PtrList<fv::options> fluidFvOptions(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to p_rghFluid\n" << endl;
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to p_rghFluid\n" << endl;
    TFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(UFluid[i])
                & fluidRegions[i].Sf()
        )
    );


    Info<< "    Adding to mixtureFluid\n" << endl;
    mixtureFluid.set
    (
        i,
        new thermalImmiscibleIncompressibleTwoPhaseMixture
        (
            UFluid[i],
            phiFluid[i]
        )
    );

    Info<< "    Adding to alpha1Fluid\n" << endl;
    alpha1Fluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "alpha1",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mixtureFluid[i].alpha1()
        )
    );

    Info<< "    Adding to alpha2Fluid\n" << endl;
    alpha2Fluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "alpha2",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mixtureFluid[i].alpha2()
        )
    );

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mixtureFluid[i].alpha1()*mixtureFluid[i].rho1() + mixtureFluid[i].alpha2()*mixtureFluid[i].rho2()
        )
    );

    Info<< "    Adding to kappaFluid\n" << endl;
    kappaFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "kappa",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mixtureFluid[i].k() 
        )
    );

    Info<< "    Adding to rhoPhiFluid\n" << endl;
    rhoPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "rhoPhi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );
    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(word::null, dimLength, Zero)
        )
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        new transportModelType
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            rhoPhiFluid[i],
            mixtureFluid[i]
        )
    );

    pFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            //fluidRegions[i],
            p_rghFluid[i] + rhoFluid[i]*ghFluid[i]
        )
    );

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    

    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");
    // //pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);


    // /*Info<< "    Adding MRF\n" << endl;
    // MRFfluid.set
    // (
    //     i,
    //     new IOMRFZoneList(fluidRegions[i])
    // );

    // Info<< "    Adding fvOptions\n" << endl;
    // fluidFvOptions.set
    // (
    //     i,
    //     new fv::options(fluidRegions[i])
    // );*/

    UEqFluid.set
    (
        i,
        new fvVectorMatrix(UFluid[i], dimForce)
    );

    // //turbulenceFluid[i].validate();

    // all boudnay conditions be in  zero gradient form for p_rgh 
    pRefCellFluid[i] = -1;
    pRefValueFluid[i] = 0.0;

    if (p_rghFluid[i].needReference())
    {
        setRefCell
        (
            pFluid[i],
            p_rghFluid[i],
            pimpleDict,
            pRefCellFluid[i],
            pRefValueFluid[i]
        );
    }

}
