/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    chtMultiRegionFoam

Description
    Combination of heatConductionFoam and buoyantFoam for conjugate heat
    transfer between a solid region and fluid region. It includes
    porous media in the primary fluid region treated explicitly.

    It handles secondary fluid or solid circuits which can be coupled
    thermally with the main fluid region. i.e radiators, etc.

    The secondary fluid region is

\*---------------------------------------------------------------------------*/
//combine headers :D
#include "fvCFD.H"
#include "CMULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "userIncompressibleTwoPhaseMixture.H"// add
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "regionProperties.H"
#include "solidRegionDiffNo.H"
#include "deltaTime/fluidRegionIncompressibleCoNo/fluidRegionIncompressibleCoNo.H"
#include "deltaTime/fluidRegionAlphaCoNo/fluidRegionAlphaCoNo.H"
#include "coordinateSystem.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    regionProperties rp(runTime);

    #include "createFluidMeshes.H"
    #include "createSolidMeshes.H"

    #include "createFluidFields.H"
    #include "createSolidFields.H"

    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "readSolidTimeControls.H"


    #include "fluidRegionIncompressibleCourantNo.H"
    #include "fluidRegionAlphaCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"
    Info <<"Enter solution loop:"<<endl;

    while (runTime.run())
    {

        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"


	#include "fluidRegionIncompressibleCourantNo.H"
	#include "fluidRegionAlphaCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "setRegionFluidFields.H"
                #include "storeOldFluidFields.H"
            }
        }

        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        {
            bool finalIter = oCorr == nOuterCorr-1;
            forAll(fluidRegions, i)
            {
            	Info<< "\nSolving for fluid region "
            			<< fluidRegions[i].name() << endl;
				#include "setRegionFluidFields.H"
				#include "readFluidMultiRegionPIMPLEControls.H"
				#include "solveFluid.H"
            }
            forAll(solidRegions, i)
            {
            	Info<< "\nSolving for solid region "
            	<< solidRegions[i].name() << endl;
				#include "setRegionSolidFields.H"
				#include "readSolidMultiRegionPIMPLEControls.H"
				#include "solveSolid.H"
            }

        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
