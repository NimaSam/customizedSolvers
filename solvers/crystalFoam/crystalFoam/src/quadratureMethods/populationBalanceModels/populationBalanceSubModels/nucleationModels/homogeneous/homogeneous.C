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

#include "homogeneous.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"
#include "../../../../crystallization/crystallization.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(homogeneous, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        homogeneous,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::homogeneous
::homogeneous
(
    const dictionary& dict
)
:
    nucleationModel(dict),
    K_(dict.lookupOrDefault("K", 0.35)),
    nus_(dict.lookupOrDefault("nus", 1.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::homogeneous
::~homogeneous()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::homogeneous::Kn
(
		const volUnivariateMoment& moment
) const
{   
	IOdictionary thermophysicalProperties
		    (
		        IOobject
		        (
		            "thermophysicalProperties",
		            moment.mesh().time().constant(),
		            moment.mesh(),
		            IOobject::MUST_READ,
		            IOobject::NO_WRITE,
		            false
		        )
		    );

		crystallization crystal(thermophysicalProperties,moment.mesh());
		const volScalarField& Y =moment.mesh().lookupObject<volScalarField>(crystal.name());
		//molar concentration
		//dimensionedScalar CMin("CMin", dimensionSet(0, -3, 0, 0, 1), 1E-8);
		//volScalarField C =max(crystal.Cc()*Y,CMin);
		volScalarField C =crystal.Cc()*Y;
		//Info << "C" << C << endl;
		//molar saturation
		volScalarField supSat = C/crystal.solubility();
		//Info << "result" << 3.0/2.0*crystal.Dm()*sqr(crystal.dm())*pow(C*Foam::constant::physicoChemical::NA,7/3)
			//		*sqrt(K_*log(crystal.Cc()/crystal.solubility()))
				//	*exp(-16.0/3.0*Foam::constant::mathematical::pi*pow(K_*log(crystal.Cc()/crystal.solubility()),3)
				//	/sqr(nus_*log(supSat))) << endl;
    return 3.0/2.0*crystal.Dm()*sqr(crystal.dm())*pow(C*Foam::constant::physicoChemical::NA,7/3)
			*sqrt(K_*log(crystal.Cc()/crystal.solubility()))
			*exp(-16.0/3.0*Foam::constant::mathematical::pi*pow(K_*log(crystal.Cc()/crystal.solubility()),3)
			/sqr(nus_*log(supSat)));
}

// ************************************************************************* //
