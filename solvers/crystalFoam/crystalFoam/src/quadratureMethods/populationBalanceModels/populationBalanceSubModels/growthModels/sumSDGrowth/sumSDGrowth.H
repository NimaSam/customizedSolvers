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
    Foam::populationBalanceSubModels::growthModels::sumSDGrowth

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
    sumSDGrowth.C

\*---------------------------------------------------------------------------*/

#ifndef sumSDGrowth_H
#define sumSDGrowth_H

#include "growthModel.H"
#include "BSGrowth.H"
#include "PNGrowth.H"
#include "BCFGrowth.H"
#include "diffusionGrowth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
                    Class sumSDGrowth Declaration
\*---------------------------------------------------------------------------*/

class sumSDGrowth
:
    public growthModel
{
private:
    
    // Private data

        BCFGrowth BCF_;
        PNGrowth  PN_;
        BSGrowth  BS_;
        diffusionGrowth  diffusion_;
    
public:

        //- Runtime type information growth model due to surface integration and bulk diffusion
        TypeName("sumSD");


    // Constructors

        //- Construct from components
        sumSDGrowth(const dictionary& dict);


    //- Destructor
    virtual ~sumSDGrowth();


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
