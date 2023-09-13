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

Class
    Foam::populationBalanceSubModels::nucleationModels::inhomogeneous

Description
        based on
    "
    Investigation of Crystallization in a Jet Y-Mixer by a Hybrid Computational Fluid Dynamics and Process
    Simulation Approach

    CRYSTAL GROWTH & DESIGN, 2005 VOL. 5, NO. 3 959-968
    "
Autor: Nima Samkhaniani
Email:Nima.Samkhaniani@gmail.com

SourceFiles
    inhomogeneous.C

\*---------------------------------------------------------------------------*/

#ifndef inhomogeneous_H
#define inhomogeneous_H

#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                    Class inhomogeneous Declaration
\*---------------------------------------------------------------------------*/

class inhomogeneous
:
    public nucleationModel
{
private:
    
    // Private data
    	dimensionedScalar K_;
    	dimensionedScalar nus_;
    	dimensionedScalar f_;
    	dimensionedScalar het_;

public:

        //- Runtime type information
        TypeName("inhomogeneous");


    // Constructors

        //- Construct from components
        inhomogeneous(const dictionary& dict);


     //- Destructor
       virtual ~inhomogeneous();

    // Member Functions

        //- nucleation model
        virtual tmp<volScalarField> Kn
        (
        		const volUnivariateMoment& moment
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
