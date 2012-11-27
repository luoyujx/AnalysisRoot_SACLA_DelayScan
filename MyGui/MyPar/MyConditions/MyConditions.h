#ifndef _MyConditions_h_
#define _MyConditions_h_

#include <TGFrame.h>

class MyParticleInfo;

class MyConditions : public TGGroupFrame
{
public:
	MyConditions(const TGWindow*, MyParticleInfo&);

public:
	void ValsChanged()		{ Emit("ValsChanged()"); }    // *SIGNAL*

	ClassDef(MyConditions,0)
};
#endif
