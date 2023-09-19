  IOdictionary transportProperties
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

  // read diffusion coeffiecient from transportProperties
  dimensionedScalar t0
    (
        transportProperties.subDict("energyData").lookup("t0")
    );
  dimensionedScalar tDET
    (
          transportProperties.subDict("energyData").lookup("tDET")
    );

  dimensionedScalar QDET
    (
        transportProperties.subDict("energyData").lookup("QDET")
    );

  dimensionedScalar rDET
  (
	   transportProperties.subDict("energyData").lookup("rDET")
  );
  dimensionedScalar vcj
  (
	   transportProperties.subDict("energyData").lookup("vcj")
  );
  dimensionedScalar rhoExp
  (
  	   transportProperties.subDict("energyData").lookup("rhoExp")
  );
  dimensionedVector centre
  (
	   transportProperties.subDict("energyData").lookup("centre")
  );


