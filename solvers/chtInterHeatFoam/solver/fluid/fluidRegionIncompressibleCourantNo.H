//Info <<"read fluid Region Incompressible courant Number"<<endl;
scalar CoNum = 0;

    forAll(fluidRegions, i)
    {
    	CoNum = max
    		(
    			fluidRegionIncompressibleCourantNo
    			(
    				fluidRegions[i],
    				runTime,
    				phiFluid[i]
    	         ),
    	         CoNum
    	    );
    }
