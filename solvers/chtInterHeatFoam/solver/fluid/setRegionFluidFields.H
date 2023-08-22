    fvMesh& mesh = fluidRegions[i];


    volScalarField& alpha1 = twoPhasePropertiesFluid[i].alpha1();
    volScalarField& alpha2 = twoPhasePropertiesFluid[i].alpha2();
    volScalarField& p_rgh = p_rghFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];
    volScalarField& p = pFluid[i];


    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    volScalarField& T = TFluid[i];
    volScalarField& K = KFluid[i];
    volScalarField& rhoC = rhoCFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    surfaceScalarField& rhoPhi = rhoPhiFluid[i];
    uniformDimensionedVectorField& g=gFluid[i];
    tmp <surfaceScalarField> tphiAlphaCorr0 = tphiAlphaCorr0Fluid[i];



    userIncompressibleTwoPhaseMixture& twoPhaseProperties= twoPhasePropertiesFluid[i];
    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();
    const dimensionedScalar& C1 = twoPhaseProperties.C1();
    const dimensionedScalar& C2 = twoPhaseProperties.C2();
    const dimensionedScalar& K1 = twoPhaseProperties.K1();
    const dimensionedScalar& K2 = twoPhaseProperties.K2();


    interfaceProperties& interface =interfaceFluid[i];
    incompressible::turbulenceModel& turbulence=turbulenceFluid[i];

