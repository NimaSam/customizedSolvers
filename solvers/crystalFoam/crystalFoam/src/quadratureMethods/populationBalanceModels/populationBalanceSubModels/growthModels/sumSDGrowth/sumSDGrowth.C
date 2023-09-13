/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Matteo Icardi
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "sumSDGrowth.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(sumSDGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        sumSDGrowth,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::sumSDGrowth
::sumSDGrowth
(
    const dictionary& dict
)
:
    growthModel(dict),
    BCF_(dict),
	PN_(dict),
	BS_(dict),
	diffusion_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::sumSDGrowth
::~sumSDGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::growthModels::sumSDGrowth::Kg
(
    const volScalarField& abscissa
) const
{   
	volScalarField kgSurf = BCF_.Kg(abscissa)+BS_.Kg(abscissa)+PN_.Kg(abscissa);
	dimensionedScalar kgMin("kgMin",kgSurf.dimensions(),1E-08);
    return kgSurf*diffusion_.Kg(abscissa)/max(kgSurf+diffusion_.Kg(abscissa),kgMin);
}

// ************************************************************************* //
