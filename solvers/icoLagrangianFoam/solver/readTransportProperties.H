const word filterPatchName (transportProperties.lookup("filterPatch"));

label filterIndex = mesh.boundaryMesh().findIndex(filterPatchName);
if ( filterIndex < 0 ){
FatalErrorInFunction
            << "Unable to find " << filterPatchName << " in patch list \n"
            << "Check the transportProperties filterPatch"
            << exit(FatalError);
}

dimensionedScalar kpp
   (
       transportProperties.lookup("kpp")
   );

dimensionedScalar alphaMax
   (
       transportProperties.lookup("alphaMax")
   );

scalar nPMax //maximum setteling particles number
   (
   		readScalar (transportProperties.lookup("nPMax"))
   );

