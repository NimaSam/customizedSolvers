/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "particleDependentSurfaceTension.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(particleDependent, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionModel,
        particleDependent,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::particleDependent::particleDependent
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    surfaceTensionModel(mesh),
    cloudName_(dict.lookupOrDefault<word>("cloud", "kinematicCloud")),
    sigma_(Function1<scalar>::New("sigma", dict)),
    Gamma_(
    	    dict.lookupOrDefault<dimensionedScalar>
    	    (
    		"Gamma",
    		dimensionedScalar("Gamma", dimensionSet(0,2,0,0,0,0,0), 0)
    	    )
    	  ),
    interfaceArea_(
    		    IOobject
            	    (
                     "interfaceArea",
                      mesh.time().timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::NO_WRITE,
                      false
                    ),
                    mesh,
                    dimensionedScalar("interfaceArea", dimArea, 0.0),
                    zeroGradientFvPatchScalarField::typeName   
    		 )
    
{
interfaceArea_.ref() = pow(mesh.V(),2.0/3);
interfaceArea_.correctBoundaryConditions();

Info << interfaceArea_ << endl;

if ( Gamma_.value() == 0) {
  Warning<< "Gamma is zero!" << endl;
}

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::particleDependent::~particleDependent()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::particleDependent::sigma() const
{
    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimSigma
        )
    );
    volScalarField& sigma = tsigma.ref();
    const Cloud <basicKinematicMPPICParcel>  & cloud = mesh_.lookupObject< Cloud < basicKinematicMPPICParcel> >(cloudName_);
    
    tmp<volScalarField> tvar
    (
        new volScalarField
        (
            IOobject
            (
                "var",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimless // TODO: assign appropriate dimension
        )
    );
    volScalarField& var = tvar.ref();

    forAllConstIter(basicKinematicMPPICCloud, cloud, iter)
    {
        const basicKinematicMPPICParcel& p = iter();
        const label celli = p.cell();

        //var[celli] += p.nParticle()*p.volume(); // return the particle volume
        var[celli] += p.nParticle(); // return the number of particle in each cells
    };
    //Info << interfaceArea_ << endl;
    sigma.field() = sigma_->value(var.field())*(1.0 - 0.75*Gamma_*var/interfaceArea_);

    volScalarField::Boundary& sigmaBf = sigma.boundaryFieldRef();
    const volScalarField::Boundary& varBf = var.boundaryField();

    forAll(sigmaBf, patchi)
    {
        sigmaBf[patchi] = sigma_->value(varBf[patchi]);
    }

    return tsigma;
}


bool Foam::surfaceTensionModels::particleDependent::readDict
(
    const dictionary& dict
)
{
    const dictionary& sigmaDict = surfaceTensionModel::sigmaDict(dict);

    cloudName_ = sigmaDict.lookupOrDefault<word>("cloud", "basicKinematicCloud");
    sigma_ = Function1<scalar>::New("sigma", sigmaDict);

    return true;
}


bool Foam::surfaceTensionModels::particleDependent::writeData
(
    Ostream& os
) const
{
    if (surfaceTensionModel::writeData(os))
    {
        os  << sigma_() << token::END_STATEMENT << nl;
        return os.good();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
