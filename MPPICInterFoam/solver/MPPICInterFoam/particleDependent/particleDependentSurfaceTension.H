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

Class
    Foam::surfaceTensionModels::particleDependent

Description
    Temperature-dependent surface tension model.

    The surface tension is evaluated from the specified Foam::Function1 for the
    temperature field looked-up from the mesh database the name of which
    may optionally be provided.

Usage
    \table
        Property     | Description               | Required    | Default value
        T            | Temperature field name    | no          | T
        sigma        | Surface tension function  | yes         |
    \endtable

    Example of the surface tension specification:
    \verbatim
        sigma
        {
            type                particleDependent;
            sigma               constant 0.07;
        }
    \endverbatim

See also
    Foam::surfaceTensionModel
    Foam::Function1

SourceFiles
    particleDependentSurfaceTension.C

\*---------------------------------------------------------------------------*/

#ifndef particleDependentSurfaceTension_H
#define particleDependentSurfaceTension_H


#include "surfaceTensionModel.H"
#include "Function1.H"
#include "basicKinematicMPPICCloud.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    namespace surfaceTensionModels
    {

        /*---------------------------------------------------------------------------*\
                                   Class particleDependent Declaration
        \*---------------------------------------------------------------------------*/

        class particleDependent
            : public surfaceTensionModel
        {
            // Private data

            //- Name of cloud field, default = "kinematicCloud"
            word cloudName_;

            //- Surface-tension function
            autoPtr<Function1<scalar>> sigma_;
            
            dimensionedScalar Gamma_;
            
            volScalarField   interfaceArea_;

        public:
            
            //- Runtime type information
            TypeName("particleDependent");

            // Constructors

            //- Construct from dictionary and mesh
            particleDependent(
                const dictionary &dict,
                const fvMesh &mesh);

            //- Destructor
            virtual ~particleDependent();

            // Member Functions

            //- Surface tension coefficient
            virtual tmp<volScalarField> sigma() const;

            //- Update surface tension coefficient from given dictionary
            virtual bool readDict(const dictionary &dict);

            //- Write in dictionary format
            virtual bool writeData(Ostream &os) const;
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace surfaceTensionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
