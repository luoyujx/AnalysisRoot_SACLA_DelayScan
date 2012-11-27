#ifndef __MyParticleContainer_H__
#define __MyParticleContainer_H__

#include <vector>

#include "MyParticle.h"
#include "MyParticleInfo.h"


//----------------------------------the particle info class-------------------------------------------------------------------
typedef std::vector<MyParticleInfo> ParticleInfos;
typedef std::vector<MyParticleInfo*> ParticleInfoPointers;
typedef std::vector<MyParticle> Particles;
typedef std::vector<MyParticle*> ParticlePointers;
class MyParticleContainer
{
public:
	MyParticleContainer()						{}

public:
	void					 Add(const char * name, double charge_au, double mass_amu, int kindParticle = 0, int fCoinGroup = 100);
	void					 Init()						{for (size_t i=0; i<fPi.size();fP[i].ReadFromInfo(fPi[i++]));}
	void					 SaveParticleInfos()		{for (size_t i=0; i<fPi.size();fPi[i++].Save());}
	void					 ClearParticles()			{for (size_t i=0; i<fP.size() ;fP[i++].Clear());}

public:
	ParticleInfos			&GetParticleInfos()			{return fPi;}
	Particles				&GetParticles()				{return fP;}
	MyParticle				&GetParticle(int idx)		{return fP[idx];}
	const size_t			 GetNbrOfParticles()		{return fP.size();}

private:
	ParticleInfos			 fPi;
	ParticleInfoPointers	 fPip;
	Particles				 fP;
	ParticlePointers		 fPp;
};

#endif