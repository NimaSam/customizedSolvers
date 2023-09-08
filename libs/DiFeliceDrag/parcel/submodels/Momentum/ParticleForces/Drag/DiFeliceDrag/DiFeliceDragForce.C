/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "DiFeliceDragForce.H"

// * * * * * * * * * * * * *  Static Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::DiFeliceDragForce<CloudType>::CdRe(const scalar Re)
{

        return (23.04+6.048*pow(Re, 1/2.0)+0.3969*Re);

}
template<class CloudType>
Foam::scalar Foam::DiFeliceDragForce<CloudType>::X(const scalar Re)
{
        return 1.0-(3.7-0.65*exp(-0.5*pow(1.5-log10(Re+rootVSmall),2)));
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DiFeliceDragForce<CloudType>::DiFeliceDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    DenseDragForce<CloudType>(owner, mesh, dict, typeName)
{}


template<class CloudType>
Foam::DiFeliceDragForce<CloudType>::DiFeliceDragForce
(
    const DiFeliceDragForce<CloudType>& df
)
:
    DenseDragForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DiFeliceDragForce<CloudType>::~DiFeliceDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::DiFeliceDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{

    const scalar alphac =
        this->alphacInterp().interpolate
        (
            p.coordinates(),
            p.currentTetIndices()
        );
      
    return forceSuSp(Zero, (mass/p.rho())*0.75*muc*CdRe(alphac*Re)*pow(alphac,X(alphac*Re))/(alphac*sqr(p.d())));


}


// ************************************************************************* //
