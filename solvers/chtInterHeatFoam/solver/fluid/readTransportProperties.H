 IOdictionary transportPropertiesDict
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word phaseChangeTwoPhaseMixtureTypeName
    (
        transportPropertiesDict.lookup("phaseChangeTwoPhaseMixture")
    );
