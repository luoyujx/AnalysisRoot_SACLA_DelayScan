#ifndef _MyGui_h_
#define _MyGui_h_

#include <TGFrame.h>

class MyPar;
class MyParticleInfo;

class MyGui : public TGMainFrame
{
public:
	MyGui(const TGWindow*, std::vector<MyParticleInfo>&);

	void	 ValsChanged()	{ Emit("ValsChanged()"); }    // *SIGNAL*

	ClassDef(MyGui,0);
};
#endif
