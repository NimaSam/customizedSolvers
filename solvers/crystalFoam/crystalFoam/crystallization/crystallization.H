/*
 * crystallization.H
 *
 *  Created on: Aug 14, 2017
 *      Author: nimasam
 */

#ifndef crystallization_H_
#define crystallization_H_

#include "word.H"
#include "scalar.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam {

class Istream;
class Ostream;

class crystallization
{


	// Private data
	const fvMesh& mesh_;

	//- Name of specie
	word name_;

	//- Number of moles of this component in the mixture
	scalar nMoles_;

	//- Molecular weight of specie [kg/kmol]
	dimensionedScalar molWeight_;

	//- mass-based solubility
	dimensionedScalar solubility_;

	//- molar density
	volScalarField Cc_;

	//- volumetric diffusion coefficient
	volScalarField Dm_;

	//- molecular diameter
	volScalarField dm_;


public:
	 //- Runtime type information
	    ClassName("crystallization");


	    // Constructors

	    //- Construct from components without name
	    inline crystallization(const scalar nMoles, const scalar molWeight, const scalar solubility);

	    //- Construct from components with name
	    inline crystallization
	    (
	    	const word& name,
	    	const scalar nMoles,
	    	const scalar molWeight,
	    	const scalar solubility
	    );

	    //- Construct as copy
	    inline crystallization(const crystallization&);

	    //- Construct as named copy
	    inline crystallization(const word& name, const crystallization&);

	    //- Construct from Istream
	    //crystallization(Istream&);

	    //- Construct from dictionary
	    //crystallization(const dictionary& dict);

	    //- Construct from dictionary, volScalarField
	    crystallization(const dictionary& dict, const fvMesh& mesh);



	    // Member Functions
	    // Access molecular weight
	    inline const word& name() const {

	    	return name_;
	    };
	    // Access solubility
	    inline const dimensionedScalar& solubility() const {

	    return solubility_;

	    };
	    // Access molecular weight
	    inline const dimensionedScalar& W() const {

	    	return molWeight_;
	    };

	    // Access molecular diameter
	    inline const volScalarField& dm() const {

	    	return dm_;
	    };

	    // Access bulk diffusion
	    inline const volScalarField& Dm() const {

	    	return Dm_;
	    };
	    // Access molar density
	    inline const volScalarField& Cc() const {

	    	return Cc_;
	    };

	    virtual void correct();
	    // Ostream Operator
	    //- Write to Ostream
	    void write(Ostream& os) const;

	   friend Ostream& operator<<(Ostream&, const crystallization&);


	   virtual ~crystallization();
};

}

#endif /* crystallization_H_ */
