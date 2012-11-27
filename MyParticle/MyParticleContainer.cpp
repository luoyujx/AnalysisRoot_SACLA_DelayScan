#include "MyParticleContainer.h"

//___________________________________________________________________________________________________________________________________________________________
void MyParticleContainer::Add(const char * name, double charge_au, double mass_amu,int kindParticle)
{
	//add particle info and set the properties//
	fPi.push_back(MyParticleInfo());
	fPi.back().SetProperties(name,charge_au,mass_amu,kindParticle);

	//add a particle//
	fP.push_back(MyParticle());

	//now add the pointer to both to the pointers vectors//
	fPip.push_back(&fPi.back());
	fPp.push_back(&fP.back());
}

