#include "MyCorrections.h"

#include "../../../MyParticle/MyParticleInfo.h"
#include "../../MyInput/MyInput.h"

//______________________________________________________________________________
MyCorrections::MyCorrections(const TGWindow *p, MyParticleInfo &pi):
	TGGroupFrame(p,"Corrections")
{
	//add T0 input//
	MyInput *fT0 = new MyInput(this,"T0 [ns]",pi.GetT0(),0.1);
	AddFrame(fT0,new TGLayoutHints(kLHintsCenterX));

	//add Angle input//
	MyInput	*fAngle = new MyInput(this,"Angle [Deg]",pi.GetAngle(),0.1);
	AddFrame(fAngle,new TGLayoutHints(kLHintsCenterX));

	//add Posx input//
	MyInput *fPosx = new MyInput(this,"Pos X [mm]",pi.GetXcor(),0.1);
	AddFrame(fPosx,new TGLayoutHints(kLHintsCenterX));

	//add Posy input//
	MyInput	*fPosy = new MyInput(this,"Pos Y [mm]",pi.GetYcor(),0.1);
	AddFrame(fPosy, new TGLayoutHints(kLHintsCenterX));

	//add Sfx input//
	MyInput	*fSfx = new MyInput(this,"Scalefactor X",pi.GetSfx(),0.1);
	AddFrame(fSfx, new TGLayoutHints(kLHintsCenterX));

	//add Sfy input//
	MyInput	*fSfy = new MyInput(this,"Scalefactor Y",pi.GetSfy(),0.1);
	AddFrame(fSfy, new TGLayoutHints(kLHintsCenterX));

	//set layout//
	SetLayoutManager(new TGMatrixLayout(this,3,2,2,0));
	Resize(279,258);


	//connections//
	fT0->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	fAngle->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	fPosx->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	fPosy->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	fSfx->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	fSfy->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	
}
