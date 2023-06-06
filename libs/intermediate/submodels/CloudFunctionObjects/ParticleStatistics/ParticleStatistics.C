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

\*---------------------------------------------------------------------------*/

#include "ParticleStatistics.H"
#include "Pstream.H"
#include "ListListOps.H"
#include "surfaceWriter.H"
#include "Time.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleStatistics<CloudType>::write()
{
    // check time range
    const fvMesh& mesh = this->owner().mesh();
    if((mesh.time().value()<timeStart_) || mesh.time().value() > timeEnd_)
    {
         return;
    }

    volScalarField& ParticleSum        = ParticleSumPtr_();
    volVectorField& ParticleUSum       = ParticleUSumPtr_();
    volVectorField& ParticleUsqSum     = ParticleUsqSumPtr_();

    //const scalarField& nPart= ParticleSum.internalField();
     scalarField& nPart= ParticleSum.ref();	


    if (ParticleSumPtr_.valid())
    {
         
      // write the total number of particles counted so far
      ParticleSumPtr_->write();
      
      // create temporary field to store mean U and U2
      volVectorField tPartUMean
      (
         IOobject
         (
             this->owner().name() + "ParticleUMean",
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         ParticleUSum
      );
      vectorField& PUSum=tPartUMean.ref();

      volVectorField tPartUsqMean
        (
            IOobject
            (
                this->owner().name() + "ParticleUPrime2Mean",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ParticleUsqSum
      );
      vectorField& PUsqSum=tPartUsqMean.ref();


      volScalarField& ParticlePDF   = ParticlePDFPtr_();
      scalarField&    pdf           = ParticlePDF.ref();
      // loop over cells and calc mean and uprime if npart>0
      forAll (PUSum, cellI)
      {
         scalar np =nPart[cellI];

         if(np>0)
         {
           	
           vector UMean     = PUSum[cellI]/np;
           PUSum[cellI]     = UMean;
           if(prime2Mean_)
           {
             PUsqSum[cellI]   = PUsqSum[cellI]/np -
                                vector(UMean[0]*UMean[0],UMean[1]*UMean[1],UMean[2]*UMean[2]);
           }
         }
      }

      //calculate the total number of added particles. Use gSum to make sure to take all the processors 
      scalar nPartSum = gSum(ParticleSum);

      //scalar pdfSum=0;
      if(writeBubblePDF_)
      { 
         if(nPartSum>0)
         {
            //  calculate the pdf by deviding the number of particles by the total and the grid cell volumes
            pdf   = nPart/mesh.V()/nPartSum;
           Info << "pdf"<< pdf <<endl;  
         }
         ParticlePDF.write();
         // calculate the volume integral of the pdf to chekc if it is unity. For validation only
         //pdfSum=gSum(pdf*mesh.V());
      }

      tPartUMean.write();
      if(prime2Mean_)
      {
         tPartUsqMean.write();
      }

      Info << "Writing ParticleStatistics UMean" ;
      if(prime2Mean_)
      {
         Info << " and UPrime2Mean";
      }
      if(writeBubblePDF_)
      {
         Info << " and PDF";
      }
      Info << " with # Particles: " << nPartSum << endl ;
    }
    else
    {
        FatalErrorIn("void Foam::ParticleStatistics<CloudType>::write()")
            << "ParticleSumPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleStatistics<CloudType>::ParticleStatistics
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner,modelName, typeName),
    ParticlePDFPtr_(nullptr),
    ParticleSumPtr_(nullptr),
    ParticleUSumPtr_(nullptr),
    ParticleUsqSumPtr_(nullptr),
    resetFields_(this->coeffDict().lookup("resetFields")),
    prime2Mean_(this->coeffDict().lookup("prime2Mean")),
    writeBubblePDF_(this->coeffDict().lookup("writeBubblePDF")),
    timeStart_(readScalar(this->coeffDict().lookup("timeStart"))),
    timeEnd_(readScalar(this->coeffDict().lookup("timeEnd")))
{}


template<class CloudType>
Foam::ParticleStatistics<CloudType>::ParticleStatistics
(
    const ParticleStatistics<CloudType>& ps
)
:
    CloudFunctionObject<CloudType>(ps),
    ParticlePDFPtr_(nullptr),
    ParticleSumPtr_(nullptr),
    ParticleUSumPtr_(nullptr),
    ParticleUsqSumPtr_(nullptr),
    resetFields_(ps.resetFields_),
    prime2Mean_(ps.prime2Mean_),
    writeBubblePDF_(ps.writeBubblePDF_),
    timeStart_(ps.timeStart_),
    timeEnd_(ps.timeEnd_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleStatistics<CloudType>::~ParticleStatistics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleStatistics<CloudType>::preEvolve()
{
    const fvMesh& mesh = this->owner().mesh();

    if (!ParticleSumPtr_.valid())
    {
        ParticlePDFPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "ParticlePDF",
                    mesh.time().timeName(mesh.time().startTime().value()),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless/dimVolume, 0.0)
            )
        );

        ParticleSumPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "ParticleSum",
                    mesh.time().timeName(mesh.time().startTime().value()),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );

        Info << "Creating cloudFunction ParticleStatistics UMean summation array"  << nl;
        // since the cloudfunction is called in the first iteration, mesh.time() referces
        // not to the start time, but to the first next time step. You want to have the start time value in case
        // you want to continue with the previous run, therefore I enfore to take mesh.time.startTime()
        // This is not string (as timeName) but a dimensionedscalar. Therefor you take the value and pass
        // it to the Time.H function timeName to turn the scalar into a string with
        // mesh.time().timeName(mesh.time().startTime().value())
        //
        ParticleUSumPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    this->owner().name() + "ParticleUMean",
                    mesh.time().timeName(mesh.time().startTime().value()),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimVelocity, vector(0.0,0,0))
            )
        );

        if(prime2Mean_)
        {
            Info << "Creating cloudFunction ParticleStatistics UPrime2Mean summation array"  << nl;
        }
        const fvMesh& mesh = this->owner().mesh();

        ParticleUsqSumPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    this->owner().name() + "ParticleUMean2Prime",
                    mesh.time().timeName(mesh.time().startTime().value()),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimVelocity*dimVelocity, vector(0.0,0.0,0.0))
            )
        );
        volScalarField& ParticleSum     = ParticleSumPtr_();
        volVectorField& ParticleUSum    = ParticleUSumPtr_();
        volVectorField& ParticleUsqSum  = ParticleUsqSumPtr_();
        scalarField&    nPart=ParticleSum.ref();
        vectorField&    PUSum=ParticleUSum.ref();
        vectorField&    PUsqSum=ParticleUsqSum.ref();
        if(resetFields_)
        {
           Info << "Reseting ParticleStatistics averaging fields" << nl;
           nPart  = scalar(0);
           PUSum  = vector(0,0,0);
           PUsqSum= vector(0,0,0);
        }
        else
        {

          // the average fields are read, but the contain the averages <U>=Sum U_i/Sum n_i  and 
          // <u'2> = Sum_i U_i*U_i/Sum n_i - <U>^2 (where U=<U>+u')
          forAll (PUSum, cellI)
          {
            scalar np       = nPart[cellI];
            vector UMean    = PUSum[cellI];
            if(prime2Mean_)
            {
               PUsqSum[cellI]  = np*(PUsqSum[cellI]+vector(UMean[0]*UMean[0],UMean[1]*UMean[1],UMean[2]*UMean[2]));
            }
            PUSum[cellI]       = np*UMean;
          }
        }

   }
}


template<class CloudType>
void Foam::ParticleStatistics<CloudType>::postEvolve()
{
    CloudFunctionObject<CloudType>::postEvolve();
}


template<class CloudType>
void Foam::ParticleStatistics<CloudType>::postMove
(
    parcelType& p,
    const scalar dt,
    const point& position0,
    bool& keepParticle
)
{
    const fvMesh& mesh = this->owner().mesh();

    // question: is there at better way to filter on the time range: not in postMove but before
    // you call this function. Now you will always enter this function even if it is not necessary
    if((mesh.time().value()<timeStart_) || mesh.time().value() > timeEnd_)
    {
      // not in the time range: go back
      return;
    }

    volScalarField& ParticleSum     = ParticleSumPtr_();
    volVectorField& ParticleUSum    = ParticleUSumPtr_();
    volVectorField& ParticleUsqSum  = ParticleUsqSumPtr_();

    ParticleSum[p.cell()]  += p.nParticle();
    ParticleUSum[p.cell()] += p.nParticle()*p.U();
    
    if(prime2Mean_)
    {
      // keep the sum of the Ux^2 , Uy^2 and Uz^2 in a vector
      ParticleUsqSum[p.cell()] += p.nParticle()*
        vector(
               p.U()[0]*p.U()[0],
               p.U()[1]*p.U()[1],
               p.U()[2]*p.U()[2]
              );
    }
}


// ************************************************************************* //
