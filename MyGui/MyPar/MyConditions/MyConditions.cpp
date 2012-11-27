#include "MyConditions.h"

#include "../../../MyParticle/MyParticleInfo.h"
#include "../../MyInput/MyInput.h"

//______________________________________________________________________________
MyConditions::MyConditions(const TGWindow *p, MyParticleInfo &pi):
	TGGroupFrame(p,"Conditions")
{
	//add Tof Range From input//
	MyInput	*toffr = new MyInput(this,"TRange From [ns]",pi.GetCondTofFr(),50);
	AddFrame(toffr, new TGLayoutHints(kLHintsCenterX));

	//add Tof Range To input//
	MyInput	*tofto = new MyInput(this,"TRange To [ns]",pi.GetCondTofTo(),50);
	AddFrame(tofto, new TGLayoutHints(kLHintsCenterX));


	MyInput *r;
	MyInput *xw;
	MyInput *yw;
	MyInput *rcx;
	MyInput *rcy;
	if (pi.GetPosFlag())
	{
		//add Center X input//
		rcx = new MyInput(this,"Center X [mm]",pi.GetCondRadX(),0.1);
		AddFrame(rcx,new TGLayoutHints(kLHintsCenterX));

		//add Center Y input//
		rcy = new MyInput(this,"Center Y [mm]",pi.GetCondRadY(),0.1);
		AddFrame(rcy,new TGLayoutHints(kLHintsCenterX));

		//add X Width input//
		xw = new MyInput(this,"Width X [mm]",pi.GetCondWidthX(),0.1);
		AddFrame(xw,new TGLayoutHints(kLHintsCenterX));
		
		//add Y Width input//
		yw = new MyInput(this,"Width Y [mm]",pi.GetCondWidthY(),0.1);
		AddFrame(yw,new TGLayoutHints(kLHintsCenterX));
		
	}
	else
	{
		//add Radius Center X input//
		rcx = new MyInput(this,"Rad Cent X [mm]",pi.GetCondRadX(),0.1);
		AddFrame(rcx,new TGLayoutHints(kLHintsCenterX));

		//add Radius Center Y input//
		rcy = new MyInput(this,"Rad Cent Y [mm]",pi.GetCondRadY(),0.1);
		AddFrame(rcy,new TGLayoutHints(kLHintsCenterX));

		//add Radius input//
		r = new MyInput(this,"Radius [mm]",pi.GetCondRad(),0.1);
		AddFrame(r,new TGLayoutHints(kLHintsCenterX));
	}

	//set layout//
	SetLayoutManager(new TGMatrixLayout(this,3,2,2,0));
	Resize(279,258);

	//connections//
	rcx->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	rcy->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	toffr->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	tofto->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	if (pi.GetPosFlag())
	{
		xw->Connect("ValChanged()","MyPar",this,"ValsChanged()");
		yw->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	}
	else
	{
		r->Connect("ValChanged()","MyPar",this,"ValsChanged()");
	}
	
}
