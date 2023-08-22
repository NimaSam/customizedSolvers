	fvMesh& mesh = solidRegions[i];
    IOdictionary& transportProperties=transportPropertiesDict[i];

    volScalarField& T=TSolid[i];
    const dimensionedScalar& rho = transportProperties.lookup("rho");
    const dimensionedScalar& Cp = transportProperties.lookup("C");
    const dimensionedScalar& K = transportProperties.lookup("K");
