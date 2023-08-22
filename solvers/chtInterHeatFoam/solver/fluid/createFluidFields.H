    // Initialise fluid field pointer lists
    PtrList<volScalarField> rhoFluid(fluidRegions.size());
    PtrList<volScalarField> p_rghFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
    PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size());
    PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
    PtrList<volScalarField> TFluid(fluidRegions.size());
    //PtrList<volScalarField> alpha1Fluid(fluidRegions.size());
    //PtrList<volScalarField> alpha2Fluid(fluidRegions.size());
    PtrList<volScalarField> KFluid(fluidRegions.size());
    PtrList<volScalarField> rhoCFluid(fluidRegions.size());
    PtrList<userIncompressibleTwoPhaseMixture> twoPhasePropertiesFluid(fluidRegions.size());
    PtrList<interfaceProperties> interfaceFluid(fluidRegions.size());
    PtrList<incompressible::turbulenceModel> turbulenceFluid(fluidRegions.size());
    PtrList<volScalarField> ghFluid(fluidRegions.size());
    PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
    PtrList<volScalarField> pFluid(fluidRegions.size());
    PtrList<surfaceScalarField> tphiAlphaCorr0Fluid(fluidRegions.size());


    List<label> pRefCellFluid(fluidRegions.size(),0);
    List<scalar> pRefValueFluid(fluidRegions.size(),0.0);

    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;
        Info<< "    Adding to pFluid\n" << endl;
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
        Info<< "    Adding to TFluid\n" << endl;
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

        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

       Info<< "    Adding twoPhasePropertiesFluid \n" << endl;
        twoPhasePropertiesFluid.set
        (
        	i,
        	new userIncompressibleTwoPhaseMixture(UFluid[i], phiFluid[i])
        );

       /* Info<< "    Adding to alpha1Fluid\n" << endl;
        alpha1Fluid.set
        (
        	i,
        	volScalarField& twoPhasePropertiesFluid[i].alpha1()
        );
        Info<< "    Adding to alpha1Fluid\n" << endl;
        alpha2Fluid.set
        (
        		i,
        		volScalarField& twoPhasePropertiesFluid[i].alpha2()
        )
        );
*/
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
        			IOobject::NO_READ,
        			IOobject::AUTO_WRITE
        		),
        		twoPhasePropertiesFluid[i].rho1()*twoPhasePropertiesFluid[i].alpha1()+
        		twoPhasePropertiesFluid[i].rho2()*(1.0-twoPhasePropertiesFluid[i].alpha1()),
        		twoPhasePropertiesFluid[i].alpha1().boundaryField().types()
        	)
        );
        Info<< " Adding to KFLuid \n" << endl;
        KFluid.set
        (
        	i,
        	new volScalarField
        	(
        		IOobject
        		(
        			"K",
        			runTime.timeName(),
        			fluidRegions[i],
        			IOobject::NO_READ,
        			IOobject::AUTO_WRITE
               	),
               	twoPhasePropertiesFluid[i].K1()*twoPhasePropertiesFluid[i].alpha1()+
               	twoPhasePropertiesFluid[i].K2()*(1.0-twoPhasePropertiesFluid[i].alpha1()),
               	zeroGradientFvPatchScalarField::typeName
               )
          );
        Info<< " Adding to CpFLuid \n" << endl;
        rhoCFluid.set
        (
        	i,
        	new volScalarField
        	(
        		IOobject
        		(
        			"rhoC",
        			runTime.timeName(),
        			fluidRegions[i],
        			IOobject::NO_READ,
        			IOobject::NO_WRITE
               	),
               	twoPhasePropertiesFluid[i].C1()*twoPhasePropertiesFluid[i].rho1()*twoPhasePropertiesFluid[i].alpha1()+
               	twoPhasePropertiesFluid[i].C2()*twoPhasePropertiesFluid[i].rho2()*(1.0 - twoPhasePropertiesFluid[i].alpha1()),
               	zeroGradientFvPatchScalarField::typeName
               )
          );

        Info<< "    Adding twoPhasePropertiesFluid \n" << endl;
        interfaceFluid.set
        (
        	i,
        	new interfaceProperties (twoPhasePropertiesFluid[i].alpha1(), UFluid[i], twoPhasePropertiesFluid[i])
        );

        Info<< "    Adding turbulence \n" << endl;
        turbulenceFluid.set
        (
            i,
           incompressible::turbulenceModel::New
            (
                UFluid[i],
                phiFluid[i],
                twoPhasePropertiesFluid[i]
            ).ptr()
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
        			IOobject::NO_READ,
        			IOobject::AUTO_WRITE
        		),
        		twoPhasePropertiesFluid[i].rho1()*phiFluid[i]
        )
        );

        Info<< "    Adding to ghFluid\n" << endl;
        ghFluid.set
        (
        	i,
        	new volScalarField("gh", gFluid[i] & fluidRegions[i].C())
        );

        Info<< "    Adding to ghfFluid\n" << endl;
        ghfFluid.set
        (
        	i,
        	new surfaceScalarField("ghf", gFluid[i] & fluidRegions[i].Cf())
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
        		 // Force p to be consistent with p_rgh
        		p_rghFluid[i]+rhoFluid[i]*ghFluid[i]
        	)
        );

        tphiAlphaCorr0Fluid.set
        (
               	i,
               	new surfaceScalarField("tphiAlphaCorr0", fluidRegions[i].magSf())//initial value is not important
        );

    }
