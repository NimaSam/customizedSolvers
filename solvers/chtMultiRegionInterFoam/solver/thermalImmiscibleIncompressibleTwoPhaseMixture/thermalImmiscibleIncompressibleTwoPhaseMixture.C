/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "thermalImmiscibleIncompressibleTwoPhaseMixture.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalImmiscibleIncompressibleTwoPhaseMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalImmiscibleIncompressibleTwoPhaseMixture::
thermalImmiscibleIncompressibleTwoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    incompressibleTwoPhaseMixture(U, phi),
    interfaceProperties(alpha1(), U, *this),
    cp1_("cp", dimensionSet(0, 2, -2, -1, 0), nuModel1().viscosityProperties()),
    cp2_("cp", dimensionSet(0, 2, -2, -1, 0), nuModel2().viscosityProperties()),
    k1_("k", dimensionSet(1, 1, -3, -1, 0), nuModel1().viscosityProperties()),
    k2_("k", dimensionSet(1, 1, -3, -1, 0), nuModel2().viscosityProperties())  
{
    Info << "Constructing thermalImmiscibleIncompressibleTwoPhaseMixture" << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::thermalImmiscibleIncompressibleTwoPhaseMixture::k() const
{
    const volScalarField limitedAlpha1
    (
        clamp(twoPhaseMixture::alpha1_, zero_one{})
    );

    return volScalarField::New
    (
        "k",
        IOobject::NO_REGISTER,
        limitedAlpha1*k1_
      + (scalar(1) - limitedAlpha1)*k2_
    );
}

Foam::tmp<Foam::volScalarField>
Foam::thermalImmiscibleIncompressibleTwoPhaseMixture::cp() const
{
    const volScalarField limitedAlpha1
    (
        clamp(twoPhaseMixture::alpha1_, zero_one{})
    );

    return volScalarField::New
    (
        "cp",
        IOobject::NO_REGISTER,
        limitedAlpha1*cp1_
      + (scalar(1) - limitedAlpha1)*cp2_
    );
}

Foam::tmp<Foam::volScalarField>
Foam::thermalImmiscibleIncompressibleTwoPhaseMixture::thermalDiffusivity() const
{
    const volScalarField limitedAlpha1
    (
        clamp(twoPhaseMixture::alpha1_, zero_one{})
    );

    return volScalarField::New
    (
        "thermalDiffusivity",
        IOobject::NO_REGISTER,
        limitedAlpha1*k1_/cp1_
      + (scalar(1) - limitedAlpha1)*k2_/cp2_
    );
}




bool Foam::thermalImmiscibleIncompressibleTwoPhaseMixture::read()
{
    return
        incompressibleTwoPhaseMixture::read()
     && interfaceProperties::read();
}


// ************************************************************************* //
