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

#include "BCFGrowth.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"
#include "../../../../crystallization/crystallization.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(BCFGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        BCFGrowth,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::BCFGrowth
::BCFGrowth
(
    const dictionary& dict
)
:
    growthModel(dict),
    minAbscissa_(dict.lookupOrDefault("minAbscissa", 1e-10)),
    maxAbscissa_(dict.lookupOrDefault("maxAbscissa", 1e10)),
    nus_(dict.lookupOrDefault("nus",1.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::BCFGrowth
::~BCFGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::growthModels::BCFGrowth::Kg
(
    const volScalarField& abscissa
) const
{   
    dimensionedScalar minAbs
    (
        "minAbs",
        abscissa.dimensions(),
        minAbscissa_.value()
    );

    dimensionedScalar maxAbs
    (
        "maxAbs",
        abscissa.dimensions(),
        maxAbscissa_.value()
    );

	IOdictionary thermophysicalProperties
	    (
	        IOobject
	        (
	            "thermophysicalProperties",
	            abscissa.mesh().time().constant(),
	            abscissa.mesh(),
	            IOobject::MUST_READ,
	            IOobject::NO_WRITE,
	            false
	        )
	    );

	crystallization crystal(thermophysicalProperties,abscissa.mesh());
	const volScalarField& Y =abscissa.mesh().lookupObject<volScalarField>(crystal.name());
	//molar concentration
	volScalarField C =crystal.Cc()*Y;
	//molar saturation
	volScalarField supSat = C/crystal.solubility();

    return 2.0*2.25E-03*sqr(nus_)*(crystal.Dm()/crystal.dm())*pow(crystal.solubility()/crystal.Cc(),4.0/3.0)
    		/log(crystal.solubility()/crystal.Cc())*sqr(1.0-supSat)
    		*pos(-abscissa+maxAbs)*pos(abscissa-minAbs);
}

// ************************************************************************* //
