label filterPatchLabel = mesh.boundaryMesh().findPatchID(filterPatchName);
vectorField filterFaces = mesh.Cf().boundaryField()[filterPatchLabel];

typedef typename Foam::newKinematicParcel<Foam::particle> parcelType ;

forAllIter( basicNewKinematicCloud , kinematicCloud, iter)
{
	parcelType& p = iter(); //paticle
	bool pdelete= false;
	if (p.active()){
	forAll(filterFaces, facei){
		vector posf = filterFaces[facei];
		scalar d = mag (posf - p.position());
		if (d < p.d()/2){
		//label celli = p.cell();
		//Marker[celli]=1.0;
		//alphap[celli]+=p.volume()/mesh.V()[celli];
		//particles.addParticle( new passiveParticle (mesh,p.position(),p.cell()));
		passiveCloud.addParticle( new parcelType (p));
		//Info << "Delete Particle on Patch Face" <<endl;
		//kinematicCloud.deleteParticle(p);
		pdelete=true;
		break;
		}
	}
	if (!pdelete)
	{
	/*forAll (Marker,celli){
	 if (Marker[celli] == 1.0){
		 vector posC = mesh.C()[celli];
		 scalar d = mag (posC - p.position());
		 if (d < p.d()){
			 label celli = p.cell();
			 Marker[celli]=1.0;
			 Info <<"Delete Particle on cells"<< endl;
			 kinematicCloud.deleteParticle(p);
			 break;
		 }
	 }
	}*/
	//forAllIter(passiveParticleCloud, particles, iter1)
	forAllIter(basicNewKinematicCloud, passiveCloud, iter1)
	{
		parcelType& passiveP = iter1();
		scalar d = mag (passiveP.position() - p.position());
		bool& active = p.active();
		vector& Up = p.U();
		if (d < ( p.d() + passiveP.d() )/2.0){
			//label celli = p.cell();
			//Marker[celli]=1.0;
			//alphap[celli]+=p.volume()/mesh.V()[celli];
			//particles.addParticle( new passiveParticle (mesh,p.position(),p.cell()));
			passiveCloud.addParticle( new parcelType (p));
			//Info <<"Delete Particle on cells"<< endl;
			//kinematicCloud.deleteParticle(p);
			active = false;
			Up = vector::zero;
			break;
		}
	  }
	}
  }
}

//reset when it is over loaded
if (nPMax < passiveCloud.size()){
	Info << "Warning: particles number:" << passiveCloud.size() << "is greater than maximum load"<< endl;
	forAllIter( basicNewKinematicCloud , kinematicCloud, iter)
	{
		parcelType& p = iter(); //paticle
		bool& active = p.active();
		active = true;
	}
	forAllIter( basicNewKinematicCloud , passiveCloud, iter)
	{
			parcelType& p = iter(); //paticle
			passiveCloud.deleteParticle(p);
	}
//	alphap=alphap0;

}

alphap = passiveCloud.theta();
alphap = min(alphap,alphaMax);




