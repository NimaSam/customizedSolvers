/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::ParticleStatistics

Description
    Creates particle statistic to keep the total sum per cell and the average
    and rms dispersed phase velocity per cell. 

    Add to bubblePropertiesDict 

    cloudFunctions
    {
      particleStatistics
      {
         resetFields      false;          // on a restart, reset the last fields or continue
         prime2Mean       false;          // Also keep and write the ux'2, uy'2 and uz'2
         writeBubblePDF   false;          // next to the summation of the bubble, also write the bubble PDF
         timeStart        0.0;            // Start averaging time
         timeEnd          1.0;            // Stop averaging time
      }
   }

Eelco van Vliet
Tata Steel, Sep 2012
Nima Samkhaniani 2023


SourceFiles
    ParticleStatistics.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleStatistics_H
#define ParticleStatistics_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleStatistics Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleStatistics
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Void fraction field
        autoPtr<volScalarField> ParticlePDFPtr_;
        autoPtr<volScalarField> ParticleSumPtr_;
        autoPtr<volVectorField> ParticleUSumPtr_;
        autoPtr<volVectorField> ParticleUsqSumPtr_;


        // initialise fields to zero on start averaging
        Switch resetFields_; 

        // calculate uprime statistics
        Switch prime2Mean_; 

        // write the bubble PDF
        Switch writeBubblePDF_; 

        // start averaging at
        scalar timeStart_; 

        // end averaging at
        scalar timeEnd_; 

protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("particleStatistics");


    // Constructors

        //- Construct from dictionary
        ParticleStatistics(const dictionary& dict, CloudType& owner, const word& modelName);

        //- Construct copy
        ParticleStatistics(const ParticleStatistics<CloudType>& ps);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType> > clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType> >
            (
                new ParticleStatistics<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleStatistics();


    // Member Functions
    //
    // Access
      
        //- Return const access to the resetFields flag
        inline const Switch& resetFields() const;

        //- Return const access to the resetFields flag
        inline const Switch& prime2Mean() const;

        //- Return const access to the resetFields flag
        inline const Switch& writeBubblePDF() const;

        //- Return start avering time 
        inline label timeStart() const;

        //- Return end avering time 
        inline label timeEnd() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-evolve hook
            virtual void postEvolve();

            //- Post-move hook
            /*virtual void postMove
            (   
                parcelType& p,
                const label cellI,
                const scalar dt
            );*/
            virtual void postMove
            (
                parcelType& p,
                const scalar dt,
                const point& position0,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
#include "ParticleStatisticsI.H"

#ifdef NoRepository
#   include "ParticleStatistics.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
