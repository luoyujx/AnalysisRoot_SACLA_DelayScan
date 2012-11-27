#ifndef _MyPar_h_
#define _MyPar_h_

#include <TGTab.h>

class MyParticleInfo;

class MyInput;
class MySpec;

class MyPar : public TGTab
{
public:
	MyPar(const TGWindow*, std::vector<MyParticleInfo>&);

public:
	void AddParticleTab(MyParticleInfo&);
	void ValsChanged()							{ Emit("ValsChanged()"); }    // *SIGNAL*

	ClassDef(MyPar,0)
};
#endif
