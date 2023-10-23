    Info<< "Reading field T\n" << endl;//add T
    volScalarField T
    (
    		IOobject
    		(
    			"T",
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
	#include "createPhi.H"
    
