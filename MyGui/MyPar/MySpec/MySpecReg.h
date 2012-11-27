#ifndef _MySpecReg_h_
#define _MySpecReg_h_

#include <TGTab.h>

class MySpectrometer;
class MySpectrometerRegion;
class MyInput;

class MySpecReg : public TGTab
{
public:
	MySpecReg(const TGWindow*, std::vector<MySpectrometerRegion>&);

public:
	void AddSpecRegTab(const char * title, MySpectrometerRegion&);
	void ValsChanged()														{ Emit("ValsChanged()"); }    // *SIGNAL*

	ClassDef(MySpecReg,0)
};
#endif
