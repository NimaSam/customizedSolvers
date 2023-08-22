    const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

    const int nCorr =
        pimple.lookupOrDefault<int>("nCorrectors", 1);

    const int nNonOrthCorr =
        pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const bool momentumPredictor =
        pimple.lookupOrDefault("momentumPredictor", true);

    const bool evaporation =
            pimple.lookupOrDefault("evaporation", true);

    //check it ;) whether it makes solver to run when BCs are gradient?
       label pRefCell = 0;
       scalar pRefValue = 0.0;
       setRefCell(p_rgh, pimple, pRefCell, pRefValue);


