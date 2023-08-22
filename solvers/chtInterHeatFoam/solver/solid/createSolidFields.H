    // Initialise solid field pointer lists
    PtrList<volScalarField> rhoSolid(solidRegions.size());
    PtrList<volScalarField> CpSolid(solidRegions.size());
    PtrList<volScalarField> rhoCpSolid(solidRegions.size());
    PtrList<volScalarField> KSolid(solidRegions.size());
    PtrList<volScalarField> TSolid(solidRegions.size());

    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding to rhoSolid\n" << endl;
        rhoSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        Info<< "    Adding to CpSolid\n" << endl;
        CpSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Cp",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        rhoCpSolid.set
        (
            i,
            new volScalarField("rhoCp", rhoSolid[i]*CpSolid[i])
        );

        Info<< "    Adding to KSolid\n" << endl;
        KSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "K",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        Info<< "    Adding to TSolid\n" << endl;
        TSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );
    }
