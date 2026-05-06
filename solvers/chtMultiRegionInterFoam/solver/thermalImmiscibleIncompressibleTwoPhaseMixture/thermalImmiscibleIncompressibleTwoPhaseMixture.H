/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::thermalImmiscibleIncompressibleTwoPhaseMixture

Description
    An immiscible incompressible two-phase mixture transport model

SourceFiles
    thermalImmiscibleIncompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef thermalImmiscibleIncompressibleTwoPhaseMixture_H
#define thermalImmiscibleIncompressibleTwoPhaseMixture_H

#include "incompressibleTwoPhaseMixture.H"
#include "interfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class thermalImmiscibleIncompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class thermalImmiscibleIncompressibleTwoPhaseMixture
:
    public incompressibleTwoPhaseMixture,
    public interfaceProperties
{
    // Private Data

        //- Coefficients for the heat capacity of each phase
    dimensionedScalar cp1_;
    dimensionedScalar cp2_;

        //- Coefficients for the thermal conductivity of each phase
    dimensionedScalar k1_;
    dimensionedScalar k2_;

public:

    TypeName("thermalImmiscibleIncompressibleTwoPhaseMixture");
    // Constructors

        //- Construct from components
        thermalImmiscibleIncompressibleTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~thermalImmiscibleIncompressibleTwoPhaseMixture() = default;


    // Member Functions

        //- Return const-access to phase1 thermal conductivity
        const dimensionedScalar& cp1() const
        {
            return cp1_;
        }

        //- Return const-access to phase2 thermal conductivity
        const dimensionedScalar& cp2() const
        {
            return cp2_;
        }

        //- Return const-access to phase1 thermal conductivity
        const dimensionedScalar& k1() const
        {
            return k1_;
        }

        //- Return const-access to phase2 thermal conductivity
        const dimensionedScalar& k2() const
        {
            return k2_;
        }


        //- Return the mixture thermal conductivity
        tmp<volScalarField> k() const;
        

        //- Return the mixture heat capacity
        tmp<volScalarField> cp() const;


        //- Return the mixture thermal diffusivity
        tmp<volScalarField> thermalDiffusivity() const;


        //- Correct the transport and interface properties
        virtual void correct()
        {
            incompressibleTwoPhaseMixture::correct();
            interfaceProperties::correct();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
