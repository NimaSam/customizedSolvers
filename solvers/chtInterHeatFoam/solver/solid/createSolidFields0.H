	    // Initialise fluid field pointer lists
    PtrList<IOdictionary> transportPropertiesDict(solidRegions.size());
	PtrList<volScalarField> TSolid(solidRegions.size());

	// Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

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

        Info<< "    Adding to solidTransportProperties\n" << endl;
        transportPropertiesDict.set
        (
        	i,
        	new Foam::IOdictionary
        	(
        		IOobject
        		(
        			"transportProperties",
        			runTime.constant(),
        			solidRegions[i],
        			IOobject::MUST_READ,
        			IOobject::NO_WRITE
        		)
        	)
        );

    }
