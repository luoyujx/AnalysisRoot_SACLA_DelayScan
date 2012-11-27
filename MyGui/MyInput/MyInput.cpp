#include <TGNumberEntry.h>
#include <TGSlider.h>

#include "MyInput.h"

//______________________________________________________________________________
MyInput::MyInput(const TGWindow *p, const char *title, double &val,double Range):
	TGGroupFrame(p,title),fVal(val),fSliderCenterValue(val),fRange(Range)
{
	//group frame
	SetTitlePos(TGGroupFrame::kLeft);

	//add number entry//
	fNumberEntry = new TGNumberEntryField(this);
	fNumberEntry->SetNumber(fVal);
	fNumberEntry->SetFormat((static_cast<int>(Range+0.1)==50) ? TGNumberFormat::kNESReal : TGNumberFormat::kNESRealThree);
	fNumberEntry->Resize(64,fNumberEntry->GetDefaultHeight());
	AddFrame(fNumberEntry, new TGLayoutHints(kLHintsLeft,18,0,0,0));

	//add slider
	fSlider = new TGHSlider(this,104);
	fSlider->SetRange(-50,50);
	fSlider->SetPosition(0);
	AddFrame(fSlider);

	//do something to this group frame
	SetLayoutManager(new TGVerticalLayout(this));
	Resize(136,78);

	//connections//
	//change the value in the number field while the slider is beeing moved//
	fSlider->Connect("PositionChanged(Int_t)","MyInput",this,"CalcValueFromPos(Int_t)");
	//Change the value, when either in the number field Return is pressed, or when the MouseButton over the Slider is released//
	fSlider->Connect("Released()","MyInput",this,"ChangeVal()");
	//whenn the number has changed int the number entry, this will also be the new center value for the slider
	fNumberEntry->Connect("ReturnPressed()","MyInput",this,"SetNewCenterVal()");

}
//______________________________________________________________________________
void MyInput::CalcValueFromPos(Int_t sliderPos)
{
	//calculate the acutal value from the position of the slider and the given range//
	double RealValue = sliderPos * fRange/50. + fSliderCenterValue;
	fNumberEntry->SetNumber(RealValue);
}
//______________________________________________________________________________
void MyInput::SetNewCenterVal()
{
	//change the center value and reset the slider to this postion//
	fSliderCenterValue = fNumberEntry->GetNumber();
	fSlider->SetPosition(0);
	//change the value and emit the changed signal//
	ChangeVal();
}
//______________________________________________________________________________
void MyInput::ChangeVal()
{ 
	//change the value and emit the changed signal//
	fVal = fNumberEntry->GetNumber();
	Emit("ValChanged()"); 
}
