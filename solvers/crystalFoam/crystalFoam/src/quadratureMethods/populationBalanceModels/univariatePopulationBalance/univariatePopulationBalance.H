/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
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
    Foam::univariatePopulationBalance

Description
    Solve a univariate population balance equation with the extended quadrature
    method of moments.

SourceFiles
    univariatePopulationBalance.C

\*---------------------------------------------------------------------------*/

#ifndef univariatePopulationBalance_H
#define univariatePopulationBalance_H

#include "populationBalanceModel.H"
#include "univariateQuadratureApproximation.H"
#include "aggregationKernel.H"
#include "breakupKernel.H"
#include "daughterDistribution.H"
#include "growthModel.H"
#include "diffusionModel.H"
#include "nucleationModel.H" //add

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceModels
{

/*---------------------------------------------------------------------------*\
                Class univariatePopulationBalance Declaration
\*---------------------------------------------------------------------------*/

class univariatePopulationBalance
:
    public populationBalanceModel
{
    // Private data

        //- Dictionary containing quadrature settings
        univariateQuadratureApproximation quadrature_;
        
        //- Switch indicating the status of the aggregation model
        Switch aggregation_;
        
        //- Switch indicating the status of the breakup model
        Switch breakup_;

        //- Switch indicating the status of the growth model
        Switch growth_;

        //- Switch indicating the status of the nucleation model
        Switch nucleation_; //add

        //- Aggregation kernel model
        autoPtr
        <
            Foam::populationBalanceSubModels::aggregationKernel
        > aggregationKernel_;
        
        //- Breakup kernel model
        autoPtr
        <
            Foam::populationBalanceSubModels::breakupKernel
        > breakupKernel_;
        
        //- Daughter distribution function
        autoPtr
        <
            Foam::populationBalanceSubModels::daughterDistribution
        > daughterDistribution_;
      
        //- Growth model function
        autoPtr
        <
            Foam::populationBalanceSubModels::growthModel
        > growthModel_;
      
        //- Diffusion model
        autoPtr
        <
            Foam::populationBalanceSubModels::diffusionModel
        > diffusionModel_;
        
        //- nucleation model function
        autoPtr
        <
        Foam::populationBalanceSubModels::nucleationModel
        > nucleationModel_; //add

    // Private member functions
        
        //- Update fields required to compute kinetic fluxes
        void updateAdvection
        (
            surfaceScalarField& phiOwn, 
            surfaceScalarField& phiNei
        );
        
        //- Calculate explicit advection term 
        tmp<volScalarField> advectMoment
        (
            const volUnivariateMoment& moment,
            const surfaceScalarField& phiOwn,
            const surfaceScalarField& phiNei
        );
        
        //- Calculate aggregation source terms
        tmp<volScalarField> aggregationSource
        (
            const volUnivariateMoment& moment
        );
        
        //- Calculate Breakup source term
        tmp<volScalarField> breakupSource(const volUnivariateMoment& moment);

        //- Calculate Growth source term
        tmp<volScalarField> growthSource(const volUnivariateMoment& moment);

        //- Calculate nucleation source term
        tmp<volScalarField> nucleationSource(const volUnivariateMoment& moment);//add


public:

    //- Runtime type information
    TypeName("univariatePopulationBalance");


    // Constructors

        //- Construct from components
        univariatePopulationBalance
        (
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~univariatePopulationBalance();


    // Member Functions

        //- Solve univariate population balance equation
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
