#include "MyPar.h"

#include "../../MyParticle/MySpectrometer/MySpectrometer.h"
#include "../../MyParticle/MySpectrometer/MySpectrometerRegion.h"
#include "../../MyParticle/MyParticleInfo.h"
#include "./MySpec/MySpec.h"
#include "./MyConditions/MyConditions.h"
#include "./MyCorrections/MyCorrections.h"

//______________________________________________________________________________
MyPar::MyPar(const TGWindow *p, std::vector<MyParticleInfo> &pi):
	TGTab(p,368,416)
{
	for (int i=0;i<pi.size();++i)
		AddParticleTab(pi[i]);

	SetTab(0);
	Resize(this->GetDefaultSize());
}

//______________________________________________________________________________
void MyPar::AddParticleTab(MyParticleInfo &pi)
{
	//Add a tab for a spectrometer region//
	TGCompositeFrame *fCompositeFrame1630;
	fCompositeFrame1630 = AddTab(pi.GetName());
	fCompositeFrame1630->SetLayoutManager(new TGHorizontalLayout(fCompositeFrame1630));

	//add conditions//
	MyConditions	*cond = new MyConditions(fCompositeFrame1630,pi);
	fCompositeFrame1630->AddFrame(cond, new TGLayoutHints(kLHintsCenterX));

	//add corrections//
	MyCorrections	*cor = new MyCorrections(fCompositeFrame1630,pi);
	fCompositeFrame1630->AddFrame(cor, new TGLayoutHints(kLHintsCenterX));

	//add spectrometer//
	MySpec *spec = new MySpec(fCompositeFrame1630,pi.GetSpectrometer());
	fCompositeFrame1630->AddFrame(spec, new TGLayoutHints(kLHintsCenterX));

	//connections//
	cond->Connect("ValsChanged()","MyPar",this,"ValsChanged()");
	cor->Connect("ValsChanged()","MyPar",this,"ValsChanged()");
	spec->Connect("ValsChanged()","MyPar",this,"ValsChanged()");
}