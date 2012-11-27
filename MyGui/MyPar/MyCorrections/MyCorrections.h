#ifndef _MyCorrections_h_
#define _MyCorrections_h_

#include <TGFrame.h>

class MyParticleInfo;

class MyCorrections : public TGGroupFrame
{
public:
	MyCorrections(const TGWindow*, MyParticleInfo&);

public:
	void ValsChanged()		{ Emit("ValsChanged()"); }    // *SIGNAL*

	ClassDef(MyCorrections,0)
};
#endif
