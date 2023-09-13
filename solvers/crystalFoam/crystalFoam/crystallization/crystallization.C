/*
 * crystallization.C
 *
 *  Created on: Aug 14, 2017
 *      Author: nimasam
 */

#include "crystallization.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
Foam::crystallization::crystallization(Istream& is)
:
    name_(is),
    nMoles_(readScalar(is)),
    molWeight_(readScalar(is)),
    solubility_(readScalar(is))
{
    is.check("crystallization::crystallization(Istream& is)");
}
*/

Foam::crystallization::crystallization(const dictionary& dict, const fvMesh& mesh)
:
	mesh_(mesh),
    name_(dict.lookup("crystalSpecie")),
    nMoles_(readScalar(dict.subDict("mixture").subDict("specie").lookup("nMoles"))),
    molWeight_("molWeight", dimensionSet(1, 0, 0, 0, -1), dict.subDict("mixture").subDict("specie")),
    solubility_("solubility", dimensionSet(0, -3, 0, 0, 1), dict.subDict("mixture").subDict("specie")),
    Cc_
            (
                IOobject
                (
                    "Cc",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("Cc", dimensionSet(0, -3, 0, 0, 1), 0.0)
            ),
    Dm_
        (
            IOobject
            (
                "Dm",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("Dm", dimensionSet(0, 2, -1, 0, 0), 0.0)
        ),
   dm_
   (
		   IOobject
		   (
			  "dm",
			  mesh.time().timeName(),
			  mesh
		   ),
		   mesh,
		   dimensionedScalar("dm", dimLength, 0.0)
   )
{
correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::crystallization::correct()
{
	const fluidThermo& flThermo =mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
	Cc_=flThermo.rho()/molWeight_;
	dm_=pow(1.0/(Cc_*Foam::constant::physicoChemical::NA),1.0/3.0);
	Dm_=Foam::constant::physicoChemical::k*flThermo.T()
					/(2*Foam::constant::mathematical::pi*flThermo.mu()*dm_);
}
void Foam::crystallization::write(Ostream& os) const
{
    dictionary dict("crystallization");
    dict.add("nMoles", nMoles_);
    dict.add("molWeight", molWeight_);
    dict.add("solubility", solubility_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //
Foam::Ostream& Foam::operator<<(Ostream& os, const crystallization& st)
{
    os  << st.name_ << tab
        << st.nMoles_ << tab
        << st.molWeight_ << tab
        << st.solubility_ ;

    os.check("Ostream& operator<<(Ostream& os, const crystallization& st)");
    return os;
}

Foam::crystallization::~crystallization(){};
