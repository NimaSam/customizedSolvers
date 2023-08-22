Info <<"read solid Region diffusion Number"<<endl;
scalar DiNum = 0;

    forAll(solidRegions, i)
    {
    //	Info <<"calculate diffusion Number:" <<endl;
    	DiNum = max
    		(
    			solidRegionDiffNo
    			(
    				solidRegions[i],
    				runTime,
    				rhoSolid[i]*CpSolid[i],
    				KSolid[i]
    	         ),
    	         DiNum
    	    );
    }
