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

#include "BSGrowth.H"
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
    defineTypeNameAndDebug(BSGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        BSGrowth,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::BSGrowth
::BSGrowth
(
    const dictionary& dict
)
:
    growthModel(dict),
    minAbscissa_(dict.lookupOrDefault("minAbscissa", 1e-10)),
    maxAbscissa_(dict.lookupOrDefault("maxAbscissa", 1e10)),
    nus_(dict.lookupOrDefault("nus",1.0)),
    K_(dict.lookupOrDefault("K",0.35)),
    kDs_(dict.lookupOrDefault("kDs_",0.1))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::BSGrowth
::~BSGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::growthModels::BSGrowth::Kg
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
	//the surface diffusion coefficient
	volScalarField Ds=kDs_*crystal.Dm();


    return 2.0*pow(2.0/Foam::constant::mathematical::pi,1.0/3.0)*(Ds/crystal.dm())*
    		pow(crystal.solubility()/crystal.Cc()*(supSat-1.0)*log(supSat),3.0/2.0)
    		*exp(-Foam::constant::mathematical::pi/3.0*sqr(K_*log(crystal.Cc()/crystal.solubility()))/(nus_*log(supSat)))
    		*pos(-abscissa+maxAbs)*pos(abscissa-minAbs);
}

// ************************************************************************* //
