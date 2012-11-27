#ifndef _MyInput_h_
#define _MyInput_h_

#include <TGFrame.h>

class TGNumberEntryField;
class TGHSlider;

class MyInput : public TGGroupFrame
{
public:
	MyInput(const TGWindow*, const char * title, double &val, double range);

public:
	void				 ValChanged()	{ Emit("ValChanged()"); }    // *SIGNAL*
	void				 ChangeVal();
	void				 SetNewCenterVal();
	void				 CalcValueFromPos(Int_t);

private:
	TGNumberEntryField	*fNumberEntry;
	TGHSlider			*fSlider;
	double				&fVal;
	double				 fRange;
	double				 fSliderCenterValue;

	ClassDef(MyInput,0);
};
#endif
