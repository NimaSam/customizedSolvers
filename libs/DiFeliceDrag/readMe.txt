readMe.txt
**********************************
1- The code is tested with OpenFOAM-9
2- Drag model is written for Di Felice Drag.
3- In order to use the code add the drag model in this address:

src/lagrangian/parcel/submodels/Momentum/ParticleForces/Drag/

4- Replace makeParcelForces.H with the original one in this address:

src/lagrangian/parcel/parcels/include/

5- then compile the whole lagrangian library again.

Or 

replace the zip file with the same directory in lagrangian directory in src/lagrangian
