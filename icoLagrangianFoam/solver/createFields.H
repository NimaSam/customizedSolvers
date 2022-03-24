    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


#   include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
    dimensionedScalar rhoVal
    (
        transportProperties.lookup("rho")
    );

    // dummy-fields to satisfy the requirements of the kinematicCloud
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        rhoVal,
        "zeroGradient"
    );

    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        rhoVal*nu,
        "zeroGradient"
    );

#   include "readGravitationalAcceleration.H"

    Info<< "Constructing kinematicCloud" << endl;
    basicNewKinematicCloud kinematicCloud
    //basicKinematicCloud kinematicCloud
    (
        "kinematicCloud",
        rho,
        U,
        mu,
        g
    );
    // Start with empty cloud
    /* passiveParticleCloud particles
     (
     		mesh,
     		"passiveParticleCloud",
     		IDLList<passiveParticle>()
     );*/
    basicNewKinematicCloud passiveCloud
     //basicKinematicCloud kinematicCloud
    (
    		mesh,
    		"passiveCloud",
    		kinematicCloud
    );

    singlePhaseTransportModel laminarTransport(U, phi);//add

    autoPtr<incompressible::turbulenceModel> turbulence //add
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    /*// define which cells passive particles are.
    volScalarField Marker
       (
           IOobject
           (
               "Marker",
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh,
           dimensionedScalar ("Marker", dimless,0)
       );
    */
    // define volume passive particle volume fraction
    volScalarField alphap
    (
    		IOobject
    		(
    			"alphap",
    			runTime.timeName(),
    			mesh,
    			IOobject::NO_READ,
    			IOobject::AUTO_WRITE
    		),
    		mesh,
    	    dimensionedScalar ("alphap", dimless,0)
    );
    volScalarField alphap0=alphap;

    #include "readTransportProperties.H"


