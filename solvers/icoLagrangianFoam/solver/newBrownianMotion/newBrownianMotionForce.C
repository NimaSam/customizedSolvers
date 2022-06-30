/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "newBrownianMotionForce.H"
#include "mathematicalConstants.H"
#include "demandDrivenData.H"
//#include "incompressible/turbulenceModel/turbulenceModel.H"
//#include "compressible/turbulenceModel/turbulenceModel.H"
 #include "turbulenceModel.H" //replaced
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::newBrownianMotionForce<CloudType>::erfInv(const scalar y) const
{
    const scalar a = 0.147;
    scalar k = 2.0/(mathematical::pi*a) +  0.5*log(1.0 - y*y);
    scalar h = log(1.0 - y*y)/a;
    scalar x = sqrt(-k + sqrt(k*k - h));

    if (y < 0.0)
    {
        return -x;
    }
    else
    {
        return x;
    }
}


template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::newBrownianMotionForce<CloudType>::kModel() const
{
    const objectRegistry& obr = this->owner().mesh();
    /*const word turbName = "turbulenceModel";

    if (obr.foundObject<compressible::turbulenceModel>(turbName))
    {
        const compressible::turbulenceModel& model =
            obr.lookupObject<compressible::turbulenceModel>(turbName);
        return model.k();
    }
    else if (obr.foundObject<incompressible::turbulenceModel>(turbName))
    {
        const incompressible::turbulenceModel& model =
            obr.lookupObject<incompressible::turbulenceModel>(turbName);
        return model.k();
    }*/
    //replaced 
    const word turbName =
         IOobject::groupName
         (
             turbulenceModel::propertiesName,
             this->owner().U().group()
         );
 
     if (obr.foundObject<turbulenceModel>(turbName))
     {
         const turbulenceModel& model =
             obr.lookupObject<turbulenceModel>(turbName);
         return model.k();
     }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::volScalarField>"
            "Foam::newBrownianMotionForce<CloudType>::kModel() const"
        )
            << "Turbulence model not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(NULL);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::newBrownianMotionForce<CloudType>::newBrownianMotionForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    rndGen_(owner.rndGen()),
    lambda_(readScalar(this->coeffs().lookup("lambda"))),
    Tc_(readScalar(this->coeffs().lookup("Tc"))),
    turbulence_(readBool(this->coeffs().lookup("turbulence"))),
    kPtr_(NULL),
    ownK_(false)
{}


template<class CloudType>
Foam::newBrownianMotionForce<CloudType>::newBrownianMotionForce
(
    const newBrownianMotionForce& bmf
)
:
    ParticleForce<CloudType>(bmf),
    rndGen_(bmf.rndGen_),
    lambda_(bmf.lambda_),
    Tc_(bmf.Tc_),
    turbulence_(bmf.turbulence_),
    kPtr_(NULL),
    ownK_(false)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::newBrownianMotionForce<CloudType>::~newBrownianMotionForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::newBrownianMotionForce<CloudType>::cacheFields(const bool store)
{
    if (turbulence_)
    {
        if (store)
        {
            tmp<volScalarField> tk = kModel();
            if (tk.isTmp())
            {
                kPtr_ = tk.ptr();
                ownK_ = true;
            }
            else
            {
                kPtr_ = tk.operator->();
                ownK_ = false;
            }
        }
        else
        {
            if (ownK_ && kPtr_)
            {
                deleteDemandDrivenData(kPtr_);
                ownK_ = false;
            }
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::newBrownianMotionForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(vector::zero, 0.0);

    const scalar dp = p.d();
    const scalar Tc = Tc_;//p.Tc();

//    const scalar eta = rndGen_.sample01<scalar>();
    const scalar alpha = 2.0*lambda_/dp;
    const scalar cc = 1.0 + alpha*(1.257 + 0.4*exp(-1.1/alpha));

    //const scalar sigma = physicoChemical::sigma.value();

    scalar f = 0.0;
  /*  if (turbulence_)
    {
        const label cellI = p.cell();
        const volScalarField& k = *kPtr_;
        const scalar kc = k[cellI];
        const scalar Dp = sigma*Tc*cc/(3*mathematical::pi*muc*dp);
        f = eta/mass*sqrt(2.0*sqr(kc)*sqr(Tc)/(Dp*dt));
    }
    else
    */
    {
        //const scalar rhoRatio = p.rho()/p.rhoc();
        const scalar rhoRatio = p.rho()/td.rhoc();
        //const scalar s0 =
           // 216*muc*sigma*Tc/(sqr(mathematical::pi)*pow5(dp)*(rhoRatio)*cc);
       	//	216*muc*1.3806488e-23*Tc/(sqr(mathematical::pi)*pow5(dp)*sqr(rhoRatio)*(p.rhoc())*cc); //add
       	const scalar s0 =
           // 216*muc*sigma*Tc/(sqr(mathematical::pi)*pow5(dp)*(rhoRatio)*cc);
       		216*muc*1.3806488e-23*Tc/(sqr(mathematical::pi)*pow5(dp)*sqr(rhoRatio)*(td.rhoc())*cc); //add

        f =sqrt(mathematical::pi*s0/dt);
/*
         Info<< "muc = " << muc << endl;
         Info<< "Tc = " << Tc << endl;
         Info<< "sigma = " << sigma << endl; 
         Info<< "dp = " << dp << endl;
         Info<< "f = " << f << endl;
        Info<< "dt = " << dt << endl;
         Info<< "mass = " << mass << endl;
*/
  
    }

    const scalar sqrt2 = sqrt(2.0);
    for (label i = 0; i < 2; i++)
    {
        const scalar x = rndGen_.sample01<scalar>();
        const scalar eta = sqrt2*erfInv(2*x - 1.0);
/*
      Info<< "eta = " << eta << endl;
      Info<< "eta = " << eta << endl;
      Info<< "f = " << f << endl;
      Info<< "mass = " << mass << endl;
      Info<< "mass*f*eta = " << mass*f*eta << endl;
*/
        value.Su()[i] = mass*f*eta; //add

    }

    return value;
}


// ************************************************************************* //
