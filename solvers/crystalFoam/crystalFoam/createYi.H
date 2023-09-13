wordList speciesNames
    (
    		thermo.lookup("species")
    );

//word inertSpecie(thermo.lookup("inertSpecie"));
word inertSpecie(thermo.lookup("inertSpecie"));

//wordList growthSource(thermo.lookup("growthSource"));
word growthSource(thermo.lookup("growthSource"));
word growthSink(thermo.lookup("growthSink"));
scalar ka = readScalar(thermo.lookup("ka"));


PtrList<volScalarField> Y (speciesNames.size());
word phaseName ("mixture");

forAll(speciesNames, i)
    {
        IOobject header
        (
            //IOobject::groupName(speciesNames[i], phaseName),
        	speciesNames[i],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // Check if field exists and can be read
        if (header.headerOk())
        {
            Y.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        //IOobject::groupName(speciesNames[i], phaseName),
                    	speciesNames[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            volScalarField Ydefault
            (
                IOobject
                (
                    "Ydefault",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            Y.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                    	speciesNames[i],
                        //IOobject::groupName(speciesNames[i], phaseName),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Ydefault
                )
            );
        }
    }

    // Do not enforce constraint of sum of mass fractions to equal 1 here
    // - not applicable to all models

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}


crystallization crystal(thermo,mesh);
scalar Sct = readScalar(thermo.subDict("mixture").subDict("transport").lookup("Sct"));


