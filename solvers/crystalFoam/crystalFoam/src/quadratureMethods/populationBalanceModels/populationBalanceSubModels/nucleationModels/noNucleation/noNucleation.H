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
    Foam::populationBalanceSubModels::nucleationModels::noNucleation

Description
    Constant nucleation model. The value of the constant is read from dictionary
    or assumed equal to 1. The minimum and maximum value of the abscissa to which the
    model is applied can be defined by the user or it is assumed to be 1e-10 and 1e10.

SourceFiles
    noNucleation.C

\*---------------------------------------------------------------------------*/

#ifndef noNucleation_H
#define noNucleation_H

#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                    Class noNucleation Declaration
\*---------------------------------------------------------------------------*/

class noNucleation
:
    public nucleationModel
{
private:
    
    // Private data
    
public:

        //- Runtime type information
        TypeName("noNucleation");


    // Constructors

        //- Construct from components
        noNucleation(const dictionary& dict);


    //- Destructor
    virtual ~noNucleation();


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
