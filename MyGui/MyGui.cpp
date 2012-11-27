#include "MyGui.h"

#include "./MyPar/MyPar.h"

//______________________________________________________________________________
MyGui::MyGui(const TGWindow* p, std::vector<MyParticleInfo> &pi):
	TGMainFrame(p,10,10,kMainFrame | kVerticalFrame)
{
	//Particle tab//
	MyPar *fPar = new MyPar(this, pi);
	AddFrame(fPar, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
	
	//setup this window//
	SetMWMHints(kMWMDecorAll,kMWMFuncAll,kMWMInputModeless);
	SetWindowName("Calibrator");
	SetEditable(false);
	MapSubwindows();

	Resize(GetDefaultSize());
	MapWindow();

	//make this window unrezisable//
	SetWMSize(GetDefaultWidth(),GetDefaultHeight());
	SetWMSizeHints(GetDefaultWidth(),GetDefaultHeight(),
				   GetDefaultWidth(),GetDefaultHeight(), 0 ,0);

	//connections//
	fPar->Connect("ValsChanged()","MyGui",this,"ValsChanged()");
}
