#   include "readTimeControls.H"
#   include "readPISOControls.H"

    bool correctPhi = true;
    if (pisoDict.found("correctPhi"))
    {
        correctPhi = Switch(pisoDict.lookup("correctPhi"));
    }

    bool checkMeshCourantNo = false;
    if (pisoDict.found("checkMeshCourantNo"))
    {
        checkMeshCourantNo = Switch(pisoDict.lookup("checkMeshCourantNo"));
    }
