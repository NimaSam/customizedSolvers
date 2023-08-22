Info <<"read fluid Region Inerface courant Number"<<endl;
scalar alphaCoNum = 0;

    forAll(fluidRegions, i)
    {
    	//Info <<"calculate diffusion Number:" <<endl;
    	alphaCoNum = max
    		(
    			fluidRegionAlphaCoNo
    			(
    				fluidRegions[i],
    				runTime,
    				twoPhasePropertiesFluid[i].alpha1(),
    				phiFluid[i]
    	         ),
    	         alphaCoNum
    	    );
    }
