#ifndef _MySpec_h_
#define _MySpec_h_

#include <TGFrame.h>

class MySpectrometer;
class TGCheckButton;

class MySpec : public TGGroupFrame
{
public:
	MySpec(const TGWindow*, MySpectrometer&);

public:
	void ValsChanged()		{ Emit("ValsChanged()"); }    // *SIGNAL*
	void CheckButtonChanged();

private:
	TGCheckButton	*mfield;
	TGCheckButton	*rotcw;
	MySpectrometer	&fSpec;

	ClassDef(MySpec,0)
};
#endif
