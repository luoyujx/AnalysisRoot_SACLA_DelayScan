#include "AddParticles.h"

#include "./MyParticle/MyParticleContainer.h"
#include "./MyMomentaCalculator/MyMomentaCalculator.h"


//setup the particle Infos, this will also load the other settings from the ini files//
// kind of particle 0:atom 1:molecule 2:electron 3:Ion 4:Ion(simple) 

//_____SACLA 2012B iodomethane_____
void AddCH3I(MyParticleContainer &particles)
{
	particles.Add("H1p",1,MyMass::Hydrogen1(),0,2);//1

	particles.Add("C4p",4,MyMass::Carbon12(),1,0);//2
	particles.Add("C3p",3,MyMass::Carbon12(),1,0);//3
	particles.Add("C2p",2,MyMass::Carbon12(),1,0);//4
	particles.Add("C1p",1,MyMass::Carbon12(),1,0);//5
	
	//particles.Add("I15p",15,MyMass::Iodine127(),1,1);
	//particles.Add("I14p",14,MyMass::Iodine127(),1,1);
	//particles.Add("I13p",13,MyMass::Iodine127(),1,1);
	//particles.Add("I12p",12,MyMass::Iodine127(),1,1);
	//particles.Add("I11p",11,MyMass::Iodine127(),1,1);
	//particles.Add("I10p",10,MyMass::Iodine127(),1,1);//7
	//particles.Add("I9p",9,MyMass::Iodine127(),1,1);//8
	//particles.Add("I8p",8,MyMass::Iodine127(),1,1);//9
	//particles.Add("I7p",7,MyMass::Iodine127(),1,1);//10
	particles.Add("I6p",6,MyMass::Iodine127(),1,1);//11
	particles.Add("I5p",5,MyMass::Iodine127(),1,1);//14
	particles.Add("I4p",4,MyMass::Iodine127(),1,1);//13
	particles.Add("I3p",3,MyMass::Iodine127(),1,1);//14
	particles.Add("I2p",2,MyMass::Iodine127(),1,1);//15
	particles.Add("I1p",1,MyMass::Iodine127(),1,1);//16
}

void AddIUracil(MyParticleContainer &particles)
{
	particles.Add("H1p",1,MyMass::Hydrogen1(),1,0);//1

	particles.Add("C3p",3,MyMass::Carbon12(),1,0);//3
	particles.Add("O3p",3,MyMass::Oxygen16(),1,0);//

	particles.Add("C2p",2,MyMass::Carbon12(),1,0);//4
	particles.Add("N2p",2,MyMass::Nitrogen14(),1,0);//
	particles.Add("O2p",2,MyMass::Oxygen16(),1,0);//
	//particles.Add("O2pL",2,MyMass::Oxygen16(),1,0);//
	//particles.Add("O2pH",2,MyMass::Oxygen16(),1,0);//

	particles.Add("C1p",1,MyMass::Carbon12(),1,0);//5
	particles.Add("N1p",1,MyMass::Nitrogen14(),1,0);//
	particles.Add("O1p",1,MyMass::Oxygen16(),1,0);//
	//particles.Add("O1pL",1,MyMass::Oxygen16(),1,2);//
	//particles.Add("O1pH",1,MyMass::Oxygen16(),1,2);//

	particles.Add("CC1p",1,MyMass::Carbon12()*2,1,0);//
	particles.Add("CN1p",1,MyMass::Carbon12() + MyMass::Nitrogen14(),1,0);//
	particles.Add("CO1p",1,MyMass::Carbon12() + MyMass::Oxygen16(),1,0);//
	
	particles.Add("I4p",4,MyMass::Iodine127(),1,1);//12
	particles.Add("I3p",3,MyMass::Iodine127(),1,1);//13
	particles.Add("I2p",2,MyMass::Iodine127(),1,1);//14
	particles.Add("I1p",1,MyMass::Iodine127(),1,1);//15
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
	//particles.Add("Ar9P",9,MyMass::Argon40(),0);
	//particles.Add("Ar10P",10,MyMass::Argon40(),0);
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

//SACLA 2014A Xe cluster pump-probe expriment
void AddXeCluster(MyParticleContainer &particles)
{
	particles.Add("Xe1P",1,MyMass::Xenon(),0);
	particles.Add("Xe2P",2,MyMass::Xenon(),0);
	particles.Add("Xe3P",3,MyMass::Xenon(),0);
	particles.Add("Xe4P",4,MyMass::Xenon(),0);
	particles.Add("XeXe1P",1,MyMass::Xenon()*2,0);
	particles.Add("XeXeXe1P",1,MyMass::Xenon()*3,0);
	//particles.Add("XeXeXeXe1P",1,MyMass::Xenon()*4,0);

}
//SACLA 2014A Ar cluster pump-probe expriment
void AddArCluster(MyParticleContainer &particles)
{
	particles.Add("Ar1P",1,MyMass::Argon40(),0);
	particles.Add("Ar2P",2,MyMass::Argon40(),0);
	particles.Add("Ar3P",3,MyMass::Argon40(),0);
	particles.Add("Ar4P",4,MyMass::Argon40(),0);
	particles.Add("Ar5P",5,MyMass::Argon40(),0);
	particles.Add("ArAr1P",1,MyMass::Argon40()*2,0);
	particles.Add("ArArAr1P",1,MyMass::Argon40()*3,0);
	//particles.Add("ArArArAr1P",1,MyMass::Argon40()*4,0);

}
//SACLA 2014A Xe cluster pump-probe expriment for isotope of Xe
void AddXeIsotope(MyParticleContainer &particles)
{	//particles.Add("128Xe1P",1,MyMass::Xenon128(),0);
	particles.Add("129Xe1P",1,MyMass::Xenon129(),0);
	//particles.Add("130Xe1P",1,MyMass::Xenon130(),0);
	particles.Add("131Xe1P",1,MyMass::Xenon131(),0);
	particles.Add("132Xe1P",1,MyMass::Xenon132(),0);
	particles.Add("134Xe1P",1,MyMass::Xenon134(),0);
	particles.Add("136Xe1P",1,MyMass::Xenon136(),0);
	//particles.Add("128Xe2P",2,MyMass::Xenon128(),0);
	particles.Add("129Xe2P",2,MyMass::Xenon129(),0);
	//particles.Add("130Xe2P",2,MyMass::Xenon130(),0);
	particles.Add("131Xe2P",2,MyMass::Xenon131(),0);
	particles.Add("132Xe2P",2,MyMass::Xenon132(),0);
	particles.Add("134Xe2P",2,MyMass::Xenon134(),0);
	particles.Add("136Xe2P",2,MyMass::Xenon136(),0);
	//particles.Add("128Xe3P",3,MyMass::Xenon128(),0);
	//particles.Add("129Xe3P",3,MyMass::Xenon129(),0);
	//particles.Add("130Xe3P",3,MyMass::Xenon130(),0);
	//particles.Add("131Xe3P",3,MyMass::Xenon131(),0);
	//particles.Add("132Xe3P",3,MyMass::Xenon132(),0);
	//particles.Add("134Xe3P",3,MyMass::Xenon134(),0);
	//particles.Add("136Xe3P",3,MyMass::Xenon136(),0);
}
//SACLA 2014B Kr cluster pump-probe expriment
void AddKrCluster(MyParticleContainer &particles)
{
	particles.Add("Kr1P",1,MyMass::Krypton(),0);
	particles.Add("Kr2P",2,MyMass::Krypton(),0);
	particles.Add("Kr3P",3,MyMass::Krypton(),0);
	particles.Add("Kr4P",4,MyMass::Krypton(),0);
	particles.Add("Kr5P",5,MyMass::Krypton(),0);
	particles.Add("KrKr1P",1,MyMass::Krypton()*2,0);
	particles.Add("KrKrKr1P",1,MyMass::Krypton()*3,0);

}
//SACLA 2014B KrAr mixture cluster pump-probe expriment
void AddKrArCluster(MyParticleContainer &particles)
{
	particles.Add("Kr1P",1,MyMass::Krypton(),0);
	particles.Add("Kr2P",2,MyMass::Krypton(),0);
	particles.Add("Kr3P",3,MyMass::Krypton(),0);
	particles.Add("Kr4P",4,MyMass::Krypton(),0);
	particles.Add("Kr5P",5,MyMass::Krypton(),0);
	particles.Add("KrKr1P",1,MyMass::Krypton()*2,0);
	particles.Add("KrKrKr1P",1,MyMass::Krypton()*3,0);
	particles.Add("Ar1P",1,MyMass::Argon40(),0);
	particles.Add("Ar2P",2,MyMass::Argon40(),0);
	particles.Add("Ar3P",3,MyMass::Argon40(),0);
	particles.Add("Ar4P",4,MyMass::Argon40(),0);
	particles.Add("Ar5P",5,MyMass::Argon40(),0);
	particles.Add("ArAr1P",1,MyMass::Argon40()*2,0);
	particles.Add("ArArAr1P",1,MyMass::Argon40()*3,0);
	particles.Add("KrAr1P",1,MyMass::Krypton()+MyMass::Argon40(),0);
	particles.Add("KrKrAr1P",1,MyMass::Krypton()*2+MyMass::Argon40(),0);
	particles.Add("KrArAr1P",1,MyMass::Krypton()+MyMass::Argon40()*2,0);
}
//SACLA 2014B Kr cluster pump-probe expriment for isotope of Kr
void AddKrIsotope(MyParticleContainer &particles)
{	
	particles.Add("78Kr1P",1,MyMass::Krypton78(),0);
	particles.Add("80Kr1P",1,MyMass::Krypton80(),0);
	particles.Add("82Kr1P",1,MyMass::Krypton82(),0);
	particles.Add("83Kr1P",1,MyMass::Krypton83(),0);
	particles.Add("84Kr1P",1,MyMass::Krypton84(),0);
	particles.Add("86Kr1P",1,MyMass::Krypton86(),0);
	particles.Add("78Kr2P",2,MyMass::Krypton78(),0);
	particles.Add("80Kr2P",2,MyMass::Krypton80(),0);
	particles.Add("82Kr2P",2,MyMass::Krypton82(),0);
	particles.Add("83Kr2P",2,MyMass::Krypton83(),0);
	particles.Add("84Kr2P",2,MyMass::Krypton84(),0);
	particles.Add("86Kr2P",2,MyMass::Krypton86(),0);
	particles.Add("78Kr3P",3,MyMass::Krypton78(),0);
	particles.Add("80Kr3P",3,MyMass::Krypton80(),0);
	particles.Add("82Kr3P",3,MyMass::Krypton82(),0);
	particles.Add("83Kr3P",3,MyMass::Krypton83(),0);
	particles.Add("84Kr3P",3,MyMass::Krypton84(),0);
	particles.Add("86Kr3P",3,MyMass::Krypton86(),0);
}

