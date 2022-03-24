const word filterPatchName (transportProperties.lookup("filterPatch"));

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

