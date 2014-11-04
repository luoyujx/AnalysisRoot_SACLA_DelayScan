#ifndef __AddParticles_h_
#define __AddParticles_h_

#include "./MyParticle/MyParticleContainer.h"

void AddArgon(MyParticleContainer &particles);
void AddXenon(MyParticleContainer &particles);
void AddXenon132(MyParticleContainer &particles);
void AddCH3I(MyParticleContainer &particles);
void AddNitrogen(MyParticleContainer &particles);
void AddIUracil(MyParticleContainer &particles);
void AddXeCluster(MyParticleContainer &particles);
void AddArCluster(MyParticleContainer &particles);
void AddXeIsotope(MyParticleContainer &particles);
void AddKrCluster(MyParticleContainer &particles);
void AddKrArCluster(MyParticleContainer &particles);
void AddKrIsotope(MyParticleContainer &particles);

namespace MyMass
{
	//source: http://physics.nist.gov/PhysRefData/Handbook/Tables/argontable1.htm (03.Feb.2009)//
	inline const double Argon36()		{return 35.967545;}		// 0.337%
	inline const double Argon38()		{return 37.962732;}		// 0.063%
	inline const double Argon40()		{return 39.962384;}		//99.600%

	inline const double Krypton78()		{return 77.920400;}		// 0.35%
	inline const double Krypton80()		{return 79.916380;}		// 2.25%
	inline const double Krypton82()		{return 81.913482;}		//11.6%
	inline const double Krypton83()		{return 82.914135;}		//11.5%
	inline const double Krypton84()		{return 83.911507;}		//57.0%
	inline const double Krypton86()		{return 85.910616;}		//17.3%
	inline const double Krypton()		{return 83.800025;}		// 100%

	inline const double Xenon128()		{return 127.903531;}	// 1.91%
	inline const double Xenon129()		{return 128.904780;}	//26.4%
	inline const double Xenon130()		{return 129.903509;}	// 4.1%
	inline const double Xenon131()		{return 130.905072;}	//21.2%
	inline const double Xenon132()		{return 131.904144;}	//26.9%
	inline const double Xenon134()		{return 133.905395;}	//10.4%
	inline const double Xenon136()		{return 135.907214;}	// 8.9%
	inline const double Xenon()			{return 131.293;}		// 100%

	inline const double Oxygen16()		{return 15.994915;}		//99.76%
	inline const double Oxygen17()		{return 16.999311;}		// 0.048%
	inline const double Oxygen18()		{return 17.999160;}		// 0.20%

	inline const double Iodine127()		{return 126.904473;}		// 100%
	inline const double Fluorine19()	{return 18.9984032;}
	inline const double Carbon12()		{return 12;}

	inline const double Nitrogen14()	{return 14.003074;}		//99.63%
	inline const double Nitrogen15()	{return 15.000108;}		// 0.37%

	inline const double Neon20()		{return 19.992435;}		//90.48%
	inline const double Neon21()		{return 20.993843;}		// 0.27%
	inline const double Neon22()		{return 21.991383;}		// 9.25%

	inline const double Helium4()		{return 4.00260;}		

	inline const double Hydrogen1()		{return 1.007825;}		// 99.985%

	//inline const double Electron()		{return 1.*MyUnitsConv::au2amu();}
}

#endif