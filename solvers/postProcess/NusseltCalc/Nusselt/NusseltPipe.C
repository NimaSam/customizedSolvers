/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    wallHeatFluxIncompressible

Description
    Calculates and writes the heat flux for all patches as the boundary field
    of a volScalarField and also prints the integrated flux for all wall
    patches.
    Based on wallHeatFlux with changes to allow it on incompressible flows
    Also removed a bug at the typeid checkline
    Eelco van Vliet
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"
#include "scalarList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	argList::addOption
	(
		"direction",
		"vector",
		"define refrence direction <vector> - eg, '(1 0 0)'"
	);
	argList::addOption
	(
			"patchName",
			"word",
			"the name of patch"
	);
	argList::addOption
	(
			"Dh",
			"scalar",
			"hydraulic diameter of pipe"
	);

	argList::addOption
	(
			"pipeLength",
			"scalar",
			"pipe length"
	);
	argList::addOption
	(
			"cellNumber",
			"int",
			"number of cells in direction which nusselt is calculating"
	);

    timeSelector::addOptions();
    #include "setRootCase.H"
	#include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
	#include "createMesh.H"

    vector direction;
    word   patchName;
    scalar Dh;
    int nCellsx = 0;
    scalar pL;
    if (
    	args.optionReadIfPresent("direction", direction)
    	&& args.optionReadIfPresent("patchName", patchName)
    	&& args.optionReadIfPresent("Dh", Dh)
    	&& args.optionReadIfPresent("cellNumber", nCellsx)
    	&& args.optionReadIfPresent("pipeLength", pL)
    )
    {


	#include "bulkCellSelection.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate(); 
        #include "createFields.H"

        scalarField Tbulk (nCellsx, 0.0);
        scalarField vol(nCellsx, 0.0);
        Info <<"calculate local Tbulk:"<<endl;
        forAll(inBulkRegionList,listi)
        {
        	labelList cellGroup = inBulkRegionList[listi];
        	forAll(cellGroup,cellj)
        	{
        		label cellI = cellGroup[cellj];
        		//Tbulk[listi] += T[cellI]*rhoC[cellI]*mag(U[cellI])*mesh.V()[cellI];
        		//vol[listi] += rhoC[cellI]*mag(U[cellI])*mesh.V()[cellI];
        		vol[listi] += mag(U[cellI])*mesh.V()[cellI];
        		Tbulk[listi] += T[cellI]*mag(U[cellI])*mesh.V()[cellI];
        	}
        }

        forAll(inBulkRegionList,listi)
        {
        	Tbulk[listi] /= (vol[listi] + 1e-08);	
        }
        
        Info <<"************ Tbulk: "<< Tbulk <<"******************" <<endl;



        volScalarField Nusselt
        (
        		IOobject
        		(
        				"Nusselt",
        				runTime.timeName(),
        				mesh
        		),
        		mesh,
        		dimensionedScalar("Nusselt", dimless, 0.0)
        );
        Info <<"calculate grad T:"<<endl;
        surfaceScalarField gradT=fvc::snGrad(T);

        forAll(inPatchList,listi)
        {
        	labelList faceGroup = inPatchList[listi];
        	Info <<"Tbulk:"<<Tbulk[listi] <<endl;
        	forAll(faceGroup,facei)
        	{
        		label faceI = faceGroup[facei];
        		scalar Tface = T.boundaryField()[patchi][faceI];
        		Info <<"TFace:"<<Tface <<endl;
        		scalar gradTface=gradT.boundaryField()[patchi][faceI];
        		Info <<"gradT:"<<gradTface <<endl;
        		Nusselt.boundaryFieldRef()[patchi][faceI] =  gradTface*Dh / max(Tface - Tbulk[listi], 1e-3);
        		Info <<"Nusselt:"<< Nusselt.boundaryField()[patchi][faceI] <<endl;
        	}
		Info <<"******************************" <<endl;
        }
        Info <<"Calc Nusselt is complete!"<<endl;
        Nusselt.write();

    }

    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
