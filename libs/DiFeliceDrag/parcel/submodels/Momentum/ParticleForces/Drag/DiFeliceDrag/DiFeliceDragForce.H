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

Class
    Foam::DiFeliceDragForce

Description
    Drag model for spheres.

    Reference:
    \verbatim
        Amsden, A. A., Butler, T. D., & O'Rourke, P. J. (1987).
        The KIVA-II computer program for transient multidimensional chemically
        reactive flows with sprays
        SAE Technical Paper.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef DiFeliceDragForce_H
#define DiFeliceDragForce_H

#include "DenseDragForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DiFeliceDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DiFeliceDragForce
:
    public DenseDragForce<CloudType>
{
public:

    // Static Member Functions

        //- Drag coefficient multiplied by Reynolds number
        static scalar CdRe(const scalar Re);
        
         //- limiter
        static scalar X(const scalar X);


    //- Runtime type information
    TypeName("DiFeliceDrag");


    // Constructors

        //- Construct from mesh
        DiFeliceDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        DiFeliceDragForce(const DiFeliceDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new DiFeliceDragForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DiFeliceDragForce();


    // Member Functions

        // Evaluation

            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DiFeliceDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
