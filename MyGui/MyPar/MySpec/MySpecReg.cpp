#include "MySpecReg.h"

#include "../../../MyParticle/MySpectrometer/MySpectrometer.h"
#include "../../../MyParticle/MySpectrometer/MySpectrometerRegion.h"
#include "../../MyInput/MyInput.h"

//______________________________________________________________________________
MySpecReg::MySpecReg(const TGWindow *p, SpecRegions &sr):
	TGTab(p,352,104)
{
	for (int i=0;i<sr.size();++i)
		AddSpecRegTab(Form("Region %d",i+1),sr[i]);

	SetTab(0);
	//Resize(this->GetDefaultSize());
}

//______________________________________________________________________________
void MySpecReg::AddSpecRegTab(const char *title, MySpectrometerRegion &sr)
{
	//Add a tab for a spectrometer region//
	TGCompositeFrame *fCompositeFrame936;
	fCompositeFrame936 = AddTab(title);
	fCompositeFrame936->SetLayoutManager(new TGHorizontalLayout(fCompositeFrame936));

	//add e-field input//
	MyInput *fField = new MyInput(fCompositeFrame936,"E-Field [V/cm]",sr.EField_Vpcm(),0.1);
	fCompositeFrame936->AddFrame(fField/*, new TGLayoutHints(kLHintsCenterX | kLHintsTop,0,0,0,0)*/);

	//add Length input//
	MyInput *fLength = new MyInput(fCompositeFrame936,"Length [mm]",sr.Length_mm(),1);
	fCompositeFrame936->AddFrame(fLength/*, new TGLayoutHints(kLHintsCenterX | kLHintsTop,0,0,0,0)*/);

	//connections//
	fField->Connect("ValChanged()","MySpecReg",this,"ValsChanged()");
	fLength->Connect("ValChanged()","MySpecReg",this,"ValsChanged()");
}