#include "AddParticles.h"

#include "./MyParticle/MyParticleContainer.h"
#include "./MyMomentaCalculator/MyMomentaCalculator.h"


//setup the particle Infos, this will also load the other settings from the ini files//
// kind of particle 0:atom 1:molecule 2:electron 3:Ion 4:Ion(simple) 

//_____SACLA 2012B iodomethane_____
void AddCH3I(MyParticleContainer &particles)
{
	particles.Add("H1p",1,MyMass::Hydrogen1(),0);//1

	particles.Add("C6p",6,MyMass::Carbon12(),1,0);//2
	particles.Add("C5p",5,MyMass::Carbon12(),1,0);//3
	particles.Add("C4p",4,MyMass::Carbon12(),1,0);
	particles.Add("C3p",3,MyMass::Carbon12(),1,0);
	particles.Add("C2p",2,MyMass::Carbon12(),1,0);
	particles.Add("C1p",1,MyMass::Carbon12(),1,0);
	
	particles.Add("I14p",14,MyMass::Iodine127(),1,1);
	particles.Add("I13p",13,MyMass::Iodine127(),1,1);
	particles.Add("I12p",12,MyMass::Iodine127(),1,1);
	particles.Add("I11p",11,MyMass::Iodine127(),1,1);
	particles.Add("I10p",10,MyMass::Iodine127(),1,1);
	particles.Add("I9p",9,MyMass::Iodine127(),1,1);
	particles.Add("I8p",8,MyMass::Iodine127(),1,1);
	particles.Add("I7p",7,MyMass::Iodine127(),1,1);
	particles.Add("I6p",6,MyMass::Iodine127(),1,1);
	particles.Add("I5p",5,MyMass::Iodine127(),1,1);
	particles.Add("I4p",4,MyMass::Iodine127(),1,1);
	particles.Add("I3p",3,MyMass::Iodine127(),1,1);
	particles.Add("I2p",2,MyMass::Iodine127(),1,1);
	particles.Add("I1p",1,MyMass::Iodine127(),1,1);
}
void AddNitrogen(MyParticleContainer &particles)
{
	particles.Add("N3P",3,MyMass::Nitrogen14(),1);
	particles.Add("N2P",2,MyMass::Nitrogen14(),1);
	particles.Add("N1P",1,MyMass::Nitrogen14(),1);
}
void AddArgon(MyParticleContainer &particles)
{
	//SACLA Argon TOF experiment
	particles.Add("Ar1P",1,MyMass::Argon40(),0);
	particles.Add("Ar2P",2,MyMass::Argon40(),0);
	particles.Add("Ar3P",3,MyMass::Argon40(),0);
	particles.Add("Ar4P",4,MyMass::Argon40(),0);
	particles.Add("Ar5P",5,MyMass::Argon40(),0);
	particles.Add("Ar6P",6,MyMass::Argon40(),0);
	particles.Add("Ar7P",7,MyMass::Argon40(),0);
	particles.Add("Ar8P",8,MyMass::Argon40(),0);
	particles.Add("Ar9P",9,MyMass::Argon40(),0);
	particles.Add("Ar10P",10,MyMass::Argon40(),0);
}
void AddXenon(MyParticleContainer &particles)
{
	//SACLA Xenon TOF experiment
	particles.Add("Xe1P",1,MyMass::Xenon(),0);
	particles.Add("Xe2P",2,MyMass::Xenon(),0);
	particles.Add("Xe3P",3,MyMass::Xenon(),0);
	particles.Add("Xe4P",4,MyMass::Xenon(),0);
	particles.Add("Xe5P",5,MyMass::Xenon(),0);
	particles.Add("Xe6P",6,MyMass::Xenon(),0);
	particles.Add("Xe7P",7,MyMass::Xenon(),0);
	particles.Add("Xe8P",8,MyMass::Xenon(),0);
	particles.Add("Xe9P",9,MyMass::Xenon(),0);
	particles.Add("Xe10P",10,MyMass::Xenon(),0);
	particles.Add("Xe11P",11,MyMass::Xenon(),0);
	particles.Add("Xe12P",12,MyMass::Xenon(),0);
	particles.Add("Xe13P",13,MyMass::Xenon(),0);

}
void AddXenon132(MyParticleContainer &particles)
{
	//SACLA Xenon TOF experiment
	particles.Add("132Xe1P",1,MyMass::Xenon132(),0);
	particles.Add("132Xe2P",2,MyMass::Xenon132(),0);
	particles.Add("132Xe3P",3,MyMass::Xenon132(),0);
	particles.Add("132Xe4P",4,MyMass::Xenon132(),0);
	particles.Add("132Xe5P",5,MyMass::Xenon132(),0);
	particles.Add("132Xe6P",6,MyMass::Xenon132(),0);
	particles.Add("132Xe7P",7,MyMass::Xenon132(),0);
	particles.Add("132Xe8P",8,MyMass::Xenon132(),0);
	particles.Add("132Xe9P",9,MyMass::Xenon132(),0);
	particles.Add("132Xe10P",10,MyMass::Xenon132(),0);
	particles.Add("132Xe11P",11,MyMass::Xenon132(),0);
	particles.Add("132Xe12P",12,MyMass::Xenon132(),0);
	particles.Add("132Xe13P",13,MyMass::Xenon132(),0);
	particles.Add("132Xe14P",14,MyMass::Xenon132(),0);
	particles.Add("132Xe15P",15,MyMass::Xenon132(),0);
	particles.Add("132Xe16P",16,MyMass::Xenon132(),0);
	particles.Add("132Xe17P",17,MyMass::Xenon132(),0);
	particles.Add("132Xe18P",18,MyMass::Xenon132(),0);
	particles.Add("132Xe19P",19,MyMass::Xenon132(),0);
	particles.Add("132Xe20P",20,MyMass::Xenon132(),0);
	particles.Add("132Xe21P",21,MyMass::Xenon132(),0);
	particles.Add("132Xe22P",22,MyMass::Xenon132(),0);
	particles.Add("132Xe23P",23,MyMass::Xenon132(),0);
	particles.Add("132Xe24P",24,MyMass::Xenon132(),0);
	particles.Add("132Xe25P",25,MyMass::Xenon132(),0);
	particles.Add("132Xe26P",26,MyMass::Xenon132(),0);

}