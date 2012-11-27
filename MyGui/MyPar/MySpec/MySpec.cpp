#include <TGButton.h>

#include "MySpec.h"

#include "MySpecReg.h"
#include "../../../MyParticle/MySpectrometer/MySpectrometer.h"
#include "../../MyInput/MyInput.h"

//______________________________________________________________________________
MySpec::MySpec(const TGWindow *p, MySpectrometer &spec):
	TGGroupFrame(p,"Spectrometer"),fSpec(spec)
{
	//the check buttons//
	TGHorizontalFrame *fHorizontalFrame828 = new TGHorizontalFrame(this,236,21,kHorizontalFrame);

	//Magnetic field is on button//
	mfield = new TGCheckButton(fHorizontalFrame828,"M-Field is On");
	mfield->SetTextJustify(36);
	mfield->SetMargins(0,0,0,0);
	mfield->SetWrapLength(-1);
	mfield->SetState(spec.MagneticFieldIsOn() ? kButtonDown : kButtonUp);
	fHorizontalFrame828->AddFrame(mfield, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

	//rotation cw//
	rotcw = new TGCheckButton(fHorizontalFrame828,"Rotation Clockwise");
	rotcw->SetTextJustify(36);
	rotcw->SetMargins(0,0,0,0);
	rotcw->SetWrapLength(-1);
	rotcw->SetState(spec.RotationClockWise() ? kButtonDown : kButtonUp);
	fHorizontalFrame828->AddFrame(rotcw, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

	AddFrame(fHorizontalFrame828, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

	//add Cyclotron Frequency input//
	MyInput *cp = new MyInput(this,"Cyclotron Freq. [ns]",spec.CyclotronPeriod_ns(),50);
	AddFrame(cp,new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));


	//add Spectrometer Regions Tab//
	MySpecReg *sr = new MySpecReg(this,spec.GetSpectrometerRegions());
	AddFrame(sr,new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

	//set layout//
	SetLayoutManager(new TGVerticalLayout(this));
	Resize(310,258);

	//connections//
	mfield->Connect("Clicked()","MySpec",this,"CheckButtonChanged()");
	rotcw->Connect("Clicked()","MySpec",this,"CheckButtonChanged()");
	cp->Connect("ValChanged()","MySpec",this,"ValsChanged()");
	sr->Connect("ValsChanged()","MySpec",this,"ValsChanged()");
}

//______________________________________________________________________________
void MySpec::CheckButtonChanged()
{
	//get the state of the buttons
	if (mfield->GetState() == kButtonDown)
		fSpec.MagneticFieldIsOn() = true;
	else
		fSpec.MagneticFieldIsOn() = false;

	if (rotcw->GetState() == kButtonDown)
		fSpec.RotationClockWise() = true;
	else
		fSpec.RotationClockWise() = false;

	//emit the signal that something has changed//
	ValsChanged();
}
