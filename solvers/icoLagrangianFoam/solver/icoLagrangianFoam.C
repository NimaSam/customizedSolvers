/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    icoFoam

Description
    Transient solver for incompressible, turbulent flow of fluids with particles

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicNewKinematicCloud.H"
#include "passiveParticleCloud.H"
//#include "basicKinematicCloud.H"
#include "turbulentTransportModel.H" //change
#include "singlePhaseTransportModel.H"//add
#include "pisoControl.H" //replace
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "createFields.H"
	#include "createControl.H" //add
	#include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

		//#include "readPISOControls.H"
		
		#include "CourantNo.H"

        Info<< "Evolving kinematicCloud" << endl;
        kinematicCloud.evolve();
        kinematicCloud.info();
	#include "filter.H"

	#include"UEqn.H"

        // --- PISO loop
       // for (int corr=0; corr<nCorr; corr++) replace
         while (piso.correct())
        {
			#include "pEqn.H"
        }

		#include "continuityErrs.H"

        turbulence->correct();



        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
