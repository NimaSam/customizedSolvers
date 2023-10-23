label patchi = mesh.boundaryMesh().findPatchID(patchName);
labelListList inBulkRegionList (nCellsx);
labelListList inPatchList(nCellsx);
volScalarField H = mesh.C() & direction;
scalar hRefMin = min(H.internalField()).value();
scalar hRefMax = max(H.internalField()).value();

Info <<"hRefMin: "<< hRefMin <<"	hRefMax: "<< hRefMax <<endl;
scalar li = mag(hRefMax-hRefMin)/(nCellsx-1.0);

for (int i = 0; i < nCellsx; ++i) {
  scalar h = hRefMin + i*li;
  Info <<"h: "<< h << endl;
  label nSelectedCell = 1;
  label nSelectedFace = 1;
  forAll(mesh.C(), celli)
  {
	scalar hCell = mesh.C()[celli] & direction;
	if ( hCell > 0.999*h && hCell < 1.001*h )
      {
		inBulkRegionList[i].setSize(nSelectedCell,celli);
		nSelectedCell ++;
      }
  }
  forAll(mesh.Cf().boundaryField()[patchi],facei){
	  scalar hFace = mesh.Cf().boundaryField()[patchi][facei] & direction;
	  if ( hFace > 0.999*h && hFace < 1.001*h )
	  {
		  Info <<"List: " << i <<" facei: "<< facei <<endl;
		  inPatchList[i].setSize(nSelectedFace,facei);
		  nSelectedFace ++;
	  }
  }
}


//Info <<"TBulkGroups: "<< inBulkRegionList <<endl;
//Info <<"TFaceGroup: "<< inPatchList <<endl;
