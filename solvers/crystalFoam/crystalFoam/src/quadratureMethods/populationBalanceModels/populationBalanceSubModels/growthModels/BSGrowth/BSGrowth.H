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
    Foam::populationBalanceSubModels::growthModels::BSGrowth

Description
    BS growth model.

        based on
    "
    Investigation of Crystallization in a Jet Y-Mixer by a Hybrid Computational Fluid Dynamics and Process
    Simulation Approach

    CRYSTAL GROWTH & DESIGN, 2005 VOL. 5, NO. 3 959-968
    "
Autor: Nima Samkhaniani
Email:Nima.Samkhaniani@gmail.com
SourceFiles
    BSGrowth.C

\*---------------------------------------------------------------------------*/

#ifndef BSGrowth_H
#define BSGrowth_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
                    Class BSGrowth Declaration
\*---------------------------------------------------------------------------*/

class BSGrowth
:
    public growthModel
{
private:
    
    // Private data

        //- Minimum value of the abscissa to apply the model
        dimensionedScalar minAbscissa_;

        //- Maximum value of the abscissa to apply the model
        dimensionedScalar maxAbscissa_;

        //- the stoichiometric coefficient
        dimensionedScalar nus_;

        //- the factor (0.310 < K < 0.414) related to the interfacial free energy of the crystal in the solution
        dimensionedScalar K_;

        // the surface diffusion coefficient, it is usually one order less than bulk diffusion
        dimensionedScalar kDs_;
    
public:

        //- Runtime type information
        TypeName("BS");


    // Constructors

        //- Construct from components
        BSGrowth(const dictionary& dict);


    //- Destructor
    virtual ~BSGrowth();


    // Member Functions

        //- Growth model
        virtual tmp<volScalarField> Kg
        (
            const volScalarField& abscissa
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
