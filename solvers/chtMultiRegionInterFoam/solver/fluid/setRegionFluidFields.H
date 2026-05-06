
    //rhoReactionThermo& thermo = reaction.thermo();
    //thermo.validate(args.executable(), "h", "e");

    
    volScalarField& rho = rhoFluid[i];
    volScalarField& kappa = kappaFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    surfaceScalarField& rhoPhi = rhoPhiFluid[i];

    transportModelType& turbulence = turbulenceFluid[i];
    
    
    
    
    
    
    //volScalarField& K = KFluid[i];
    //volScalarField& dpdt = dpdtFluid[i];

    volScalarField& p_rgh = p_rghFluid[i];
    volScalarField& p = pFluid[i];
    volScalarField& T = TFluid[i];

    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    //multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields =
    //    fieldsFluid[i];

    //volScalarField& Qdot = QdotFluid[i];

    
    //IOMRFZoneList& MRF = MRFfluid[i];
    //fv::options& fvOptions = fluidFvOptions[i];

    fvVectorMatrix& UEqn = UEqFluid[i];


    //bool frozenFlow = frozenFlowFluid[i];

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];

    
    //const pressureControl& pressureControl = pressureControls[i];

    thermalImmiscibleIncompressibleTwoPhaseMixture& mixture = mixtureFluid[i];
