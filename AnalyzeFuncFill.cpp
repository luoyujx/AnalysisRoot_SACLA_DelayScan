#include "AnalyzeFuncFill.h"
#include "./MyAnalyzer/MyAnalyzer.h"
#include <TFile.h>

//-------------------------------------Fill Ion spectra---------------------------------------------------------------//
void fillSpectra(const MyParticle &p1, const MyParticle &p2, MyHistos &hi, int hiOff)
{
	TString Hname(p2.GetName());
	Hname += "GatedBy";
	Hname += p1.GetName();
	for (size_t i=0;i<p2.GetNbrOfParticleHits();++i)
	{
		hi.fill(hiOff+0,"Mass",p2[i].Mass(),"Mass/q",20000,0,200,Hname.Data());
		hi.fill(hiOff+1,"TofCor",p2[i].TofCor(),"tof [ns]",20000,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());
		hi.fill(hiOff+2,"Det",p2[i].X(),p2[i].Y(),"x [mm]","y [mm]",300,p2.GetCondRadX()-p2.GetCondRad()*1.3,p2.GetCondRadX()+p2.GetCondRad()*1.3,300,p2.GetCondRadY()-p2.GetCondRad()*1.3,p2.GetCondRadY()+p2.GetCondRad()*1.3,Hname.Data());
		hi.fill(hiOff+3,"Tof",p2[i].Tof(),"tof [ns]",10000,p2.GetCondTofFr()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()+p2.GetCondTofRange()*0.3,Hname.Data());
		hi.fill(hiOff+4,"XPosVsTof",p2[i].Tof(),p2[i].X(),"tof [ns]","x [mm]",5000,p2.GetCondTofFr()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()+p2.GetCondTofRange()*0.3,300,p2.GetCondRadX()-p2.GetCondRad()*1.3,p2.GetCondRadX()+p2.GetCondRad()*1.3,Hname.Data());
		hi.fill(hiOff+5,"YPosVsTof",p2[i].Tof(),p2[i].Y(),"tof [ns]","y [mm]",5000,p2.GetCondTofFr()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()+p2.GetCondTofRange()*0.3,300,p2.GetCondRadX()-p2.GetCondRad()*1.3,p2.GetCondRadX()+p2.GetCondRad()*1.3,Hname.Data());
	}
}
//--------------------------------------fill particle Histograms---------------------------------------------------//
void fillParticleHistograms(const MyParticle &p, const MyParticleHit &ph, std::vector<double>& intensity, MyHistos &hi, int hiOff, std::vector<double>& intPart )
{	
	double MomLim = 800;
	TString pname(p.GetName());
	if (pname == "H1p") MomLim = 100;
	const double SliceLim = 20;
	if (p.GetCharge_au()>8) MomLim = 1200;

	//Reconstruction Method for this particle Hit//
	hi.fill(hiOff++,"ReconstructionMethod",ph.RekMeth(),"Reconstruction Nbr",60,0,30,Form("%s/Raw",p.GetName()));
	if (intensity.size() && (intPart.size()>1)) hi.fill(hiOff++,Form("Intensity%s",p.GetName()),intensity[0],"Laser Power",intPart.size()-1,&intPart.front(),Form("%s",p.GetName()));

	//detektor pictures//
	hi.fill(hiOff++,"Det"			,ph.X()			,ph.Y()			,"x [mm]","y [mm]",300,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,300,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	hi.fill(hiOff++,"DetCor"		,ph.XCor()		,ph.YCor()		,"x [mm]","y [mm]",300,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,300,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	hi.fill(hiOff++,"DetCorRot"		,ph.XCorRot()	,ph.YCorRot()	,"x [mm]","y [mm]",300,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,300,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	hi.fill(hiOff++,"DetCorRotScale",ph.XCorRotScl(),ph.YCorRotScl(),"x [mm]","y [mm]",300,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,300,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));

	//tof//
	hi.fill(hiOff++,"Tof"   ,ph.Tof()   ,"tof [ns]",1000,p.GetCondTofFr()          -p.GetCondTofRange()*0.3,p.GetCondTofTo()          +p.GetCondTofRange()*0.3,Form("%s/Raw",p.GetName()));
	hi.fill(hiOff++,"TofCor",ph.TofCor(),"tof [ns]",1000,p.GetCondTofFr()-p.GetT0()-p.GetCondTofRange()*0.3,p.GetCondTofTo()-p.GetT0()+p.GetCondTofRange()*0.3,Form("%s/Raw",p.GetName()));

	//pos vs tof//
	hi.fill(hiOff++,"XPosVsTof",ph.TofCor(),ph.XCorRotScl(),"tof [ns]","x [mm]",300,p.GetCondTofFr()-p.GetT0()-p.GetCondTofRange()*0.3,p.GetCondTofTo()-p.GetT0()+p.GetCondTofRange()*0.3,300,p.GetCondRadX()-p.GetXcor()-p.GetCondRad()*1.3,p.GetCondRadX()-p.GetXcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	hi.fill(hiOff++,"YPosVsTof",ph.TofCor(),ph.YCorRotScl(),"tof [ns]","y [mm]",300,p.GetCondTofFr()-p.GetT0()-p.GetCondTofRange()*0.3,p.GetCondTofTo()-p.GetT0()+p.GetCondTofRange()*0.3,300,p.GetCondRadY()-p.GetYcor()-p.GetCondRad()*1.3,p.GetCondRadY()-p.GetYcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));

	//hi.fill(hiOff++,"XPosVsTofFine",ph.TofCor(),ph.XCorRotScl(),"tof [ns]","x [mm]",static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),300,p.GetCondRadX()-p.GetXcor()-p.GetCondRad()*1.3,p.GetCondRadX()-p.GetXcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	//hi.fill(hiOff++,"YPosVsTofFine",ph.TofCor(),ph.YCorRotScl(),"tof [ns]","y [mm]",static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),300,p.GetCondRadY()-p.GetYcor()-p.GetCondRad()*1.3,p.GetCondRadY()-p.GetYcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));

	//momenta//
	hi.fill(hiOff++,"PxPy",ph.Px(),ph.Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",p.GetName()));
	if (TMath::Abs(ph.Pz()) < SliceLim) 
		hi.fill(hiOff++,"PxPySlice",ph.Px(),ph.Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",p.GetName()));else hiOff++;

	hi.fill(hiOff++,"PxPz",ph.Pz(),ph.Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",p.GetName()));
	if (TMath::Abs(ph.Py()) < SliceLim) 
		hi.fill(hiOff++,"PxPzSlice",ph.Pz(),ph.Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",p.GetName()));else hiOff++;

	hi.fill(hiOff++,"PyPz",ph.Pz(),ph.Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",p.GetName()));
	if (TMath::Abs(ph.Px()) < SliceLim)
		hi.fill(hiOff++,"PyPzSlice",ph.Pz(),ph.Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",p.GetName()));else hiOff++;
	
	hi.fill(hiOff++,"TotalMomentum",ph.P(),"p [a.u.]",300,0,MomLim,Form("%s/Momenta",p.GetName()));

	//Energy//
	hi.fill(hiOff++,"Energy",ph.E(),"Energy [eV]",300,0,300,Form("%s/Energy",p.GetName()));
	//Angle//
	hi.fill(hiOff++,"ThetaX",ph.ThetaX(),"#theta [deg]",180,0,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"ThetaY",ph.ThetaY(),"#theta [deg]",180,0,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"ThetaZ",ph.ThetaZ(),"#theta [deg]",180,0,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiXY",ph.PhiXY(),"#phi [deg]",360,-180,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiYZ",ph.PhiYZ(),"#phi [deg]",360,-180,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiZX",ph.PhiZX(),"#phi [deg]",360,-180,180,Form("%s/Angle",p.GetName()));

	//hi.fill(hiOff++,"ThetaZvsEnergy",ph.ThetaZ(),ph.E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,100,Form("%s/Angle",p.GetName()));
	//hi.fill(hiOff++,"ThetaXvsEnergy",ph.ThetaX(),ph.E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,100,Form("%s/Angle",p.GetName()));
	//hi.fill(hiOff++,"ThetaYvsEnergy",ph.ThetaY(),ph.E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,100,Form("%s/Angle",p.GetName()));
	//hi.fill(hiOff++,"PhiXYvsEnergy",ph.PhiXY(),ph.E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,100,Form("%s/Angle",p.GetName()));
	//hi.fill(hiOff++,"PhiYZvsEnergy",ph.PhiYZ(),ph.E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,100,Form("%s/Angle",p.GetName()));
	//hi.fill(hiOff++,"PhiZXvsEnergy",ph.PhiZX(),ph.E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,100,Form("%s/Angle",p.GetName()));

	//hi.fill(hiOff,"ThetaY_double",ph.ThetaY(),"#theta [deg]",180,0,360,Form("%s/Angle",p.GetName()));
	//hi.fill(hiOff++,"ThetaY_double",360 - ph.ThetaY(),"#theta [deg]",180,0,360,Form("%s/Angle",p.GetName()));

	//double ThetaY = ph.ThetaY();
	//if (ph.Pz() < 0) ThetaY = 360 - ThetaY;
	//hi.fill(hiOff++,"ThetaY360",ThetaY,"#theta [deg]",360,0,360,Form("%s/Angle",p.GetName()));
	
	//Intensity//
	//for (size_t i=0; i< intensity.size();++i)
	//	hi.fill(hiOff++,Form("Int%d",i),intensity[i],Form("Int %d",i),1000,0,1000,Form("%s/Intensity",p.GetName()))
}
//-------------------Fill molecule histogram-----------------------------------------------------------------------------------------------------//
void fillMoleculeHistogram(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff, Molecule &mol, std::vector<double>& intPart)
{

	double MomLim = 800;
	double MomSumRotLim = 800;
	double MomScale = mol.momSumFactor;

	TString Hname(p1.GetName());
	if (Hname=="H1p") 
	{
		MomLim = 100;
		MomSumRotLim = 300;
		MomScale = mol.momSumFactor;
	}
	Hname += p2.GetName();
	
	const double pxSumWidth = mol.momSumWindowX;//10,9,7,5
	const double pySumWidth = mol.momSumWindowY;//10,8,5,4
	const double pzSumWidth = mol.momSumWindowZ;//5,6,4,2
	const double pxySumWidth = mol.momSumWindowX;//10,9,7,5
	const double pyzSumWidth = mol.momSumWindowY;//10,8,5,4
	const double pzxSumWidth = mol.momSumWindowZ;//5,6,4,2

	//std::cout<< Hname << ":" <<std::endl;
	//std::cout<<" pxSumWidth" << mol.momSumWindowX <<std::endl;
	//std::cout<<" pySumWidth" << mol.momSumWindowX <<std::endl;
	//std::cout<<" pzSumWidth" << mol.momSumWindowX <<std::endl;
	//std::cout<<" Factor: " <<MomScale <<std::endl;

	for (size_t i=0; i<p1.GetNbrOfParticleHits();++i)
	{
		for (size_t j=0;j<p2.GetNbrOfParticleHits();++j)//i=j
		{
			//skip if we are going to put in the same particlehit, for the case that we want coincidences for the same particle//
			if ((i>=j) && (p1==p2)) continue;//i==j

			//---p1="C" p2="I"---//
			//calcutrate momentum sum by Edwin's method
			const double ThetaXY = TMath::ATan2(p2[j].Py(),p2[j].Px());
			const double p1x_XY = p1[i].PxRotXY(ThetaXY);
			const double p1y_XY = p1[i].PyRotXY(ThetaXY);
			const double p2x_XY = p2[j].PxRotXY(ThetaXY);
			const double p2y_XY = p2[j].PyRotXY(ThetaXY);

			const double ThetaYZ = TMath::ATan2(p2[j].Pz(),p2[j].Py());
			const double p1y_YZ = p1[i].PyRotYZ(ThetaYZ);
			const double p1z_YZ = p1[i].PzRotYZ(ThetaYZ);
			const double p2y_YZ = p2[j].PyRotYZ(ThetaYZ);
			const double p2z_YZ = p2[j].PzRotYZ(ThetaYZ);

			const double ThetaZX = TMath::ATan2(p2[j].Px(),p2[j].Pz());
			const double p1z_ZX = p1[i].PzRotZX(ThetaZX);
			const double p1x_ZX = p1[i].PxRotZX(ThetaZX);
			const double p2z_ZX = p2[j].PzRotZX(ThetaZX);
			const double p2x_ZX = p2[j].PxRotZX(ThetaZX);

			const double p12xSumRot_XY = p1x_XY/MomScale+p2x_XY;
			const double p12ySumRot_XY = p1y_XY/MomScale+p2y_XY;
			const double p12ySumRot_YZ = p1y_YZ/MomScale+p2y_YZ;
			const double p12zSumRot_YZ = p1z_YZ/MomScale+p2z_YZ;
			const double p12zSumRot_ZX = p1z_ZX/MomScale+p2z_ZX;
			const double p12xSumRot_ZX = p1x_ZX/MomScale+p2x_ZX;

			hi.fill(hiOff+31,"PxPyRot",p1x_XY, p1y_XY,"px [a.u.]","py [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+31,"PxPyRot",p2x_XY, p2y_XY,"px [a.u.]","py [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+32,"PxPySumRot",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			
			hi.fill(hiOff+33,"PyPzRot",p1y_YZ, p1z_YZ,"py [a.u.]","pz [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+33,"PyPzRot",p2y_YZ, p2z_YZ,"py [a.u.]","pz [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+34,"PyPzSumRot",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			
			hi.fill(hiOff+35,"PzPxRot",p1z_ZX, p1x_ZX,"pz [a.u.]","px [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+35,"PzPxRot",p2z_ZX, p2x_ZX,"pz [a.u.]","px [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+36,"PzPxSumRot",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));

			//define the norm of momentum sum vector that rotated
			const double normOfSumRot_XY = TMath::Sqrt(p12xSumRot_XY*p12xSumRot_XY + p12ySumRot_XY*p12ySumRot_XY);
			const double normOfSumRot_YZ = TMath::Sqrt(p12ySumRot_YZ*p12ySumRot_YZ + p12zSumRot_YZ*p12zSumRot_YZ);
			const double normOfSumRot_ZX = TMath::Sqrt(p12zSumRot_ZX*p12zSumRot_ZX + p12xSumRot_ZX*p12xSumRot_ZX);

			//-----Gated stuff-----//
			//x momentum sum//
			hi.fill(hiOff+1,"PxSum",(p1[i].Px()/MomScale + p2[j].Px()),Form("Px_{%s} + Px_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-400,400,Form("%s/MomSums",Hname.Data()));
			if (normOfSumRot_YZ < pyzSumWidth)
				hi.fill(hiOff+2,"PxSumCondYZ",(p1[i].Px()/MomScale +p2[j].Px()),Form("Px_{%s} + Px_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-400,400, Form("%s/MomSums",Hname.Data()));
			//y momentum sum//
			hi.fill(hiOff+3,"PySum",(p1[i].Py()/MomScale + p2[j].Py()),Form("Py_{%s} + Py_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-400,400, Form("%s/MomSums",Hname.Data()));
			if (normOfSumRot_ZX < pzxSumWidth)
				hi.fill(hiOff+4,"PySumCondXZ",(p1[i].Py()/MomScale + p2[j].Py()),Form("Py_{%s} + Py_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-400,400, Form("%s/MomSums",Hname.Data()));
			//z momentum sum//
			hi.fill(hiOff+5,"PzSum",(p1[i].Pz()/MomScale + p2[j].Pz()),Form("Pz_{%s} + Pz_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-400,400, Form("%s/MomSums",Hname.Data()));
			if (normOfSumRot_XY < pxySumWidth)
				hi.fill(hiOff+6,"PzSumCondXY",(p1[i].Pz()/MomScale + p2[j].Pz()),Form("Pz_{%s} + Pz_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-400,400, Form("%s/MomSums",Hname.Data()));

			//PIPICO//
			hi.fill(hiOff+7,"PIPICO",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());
			hi.fill(hiOff+10,"PIPICOMom",p1[i].Pz(),p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			hi.fill(hiOff+11,"PIPICOMomFactor",p1[i].Pz()/MomScale,p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			if (normOfSumRot_XY < pxySumWidth)
			{
				hi.fill(hiOff+8,"PIPICOCondXY",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());
				hi.fill(hiOff+12,"PIPICOMomCondXY",p1[i].Pz(),p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
				hi.fill(hiOff+13,"PIPICOMomFactorCondXY",p1[i].Pz()/MomScale,p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			}

			//For confirmation. gated by only one plane
			if (normOfSumRot_XY < pxySumWidth)
			{
				hi.fill(hiOff+37,"PzPxSumRotCondXY",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
				hi.fill(hiOff+38,"PyPzSumRotCondXY",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			}
			if (normOfSumRot_YZ < pyzSumWidth)
			{
				hi.fill(hiOff+39,"PxPySumRotCondYZ",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
				hi.fill(hiOff+40,"PzPxSumRotCondYZ",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			}
			if (normOfSumRot_ZX < pzxSumWidth)
			{
				hi.fill(hiOff+41,"PyPzSumRotCondZX",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
				hi.fill(hiOff+42,"PxPySumRotCondZX",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			}
			///////////////////apply two sum condition//////////////////////
			if (normOfSumRot_XY < pxySumWidth)
			if (normOfSumRot_YZ < pyzSumWidth)
				hi.fill(hiOff+46,"PzPxSumRotCondXY-YZ",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));

			if (normOfSumRot_YZ < pyzSumWidth)
			if (normOfSumRot_ZX < pzxSumWidth)
				hi.fill(hiOff+47,"PxPySumRotCondYZ-ZX",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));

			if (normOfSumRot_ZX < pzxSumWidth)
			if (normOfSumRot_XY < pxySumWidth)
				hi.fill(hiOff+48,"PyPzSumRotCondZX-XY",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));

			///////////////////Momentum sum condition//////////////////////
			if (normOfSumRot_XY < pxySumWidth)
			if (normOfSumRot_YZ < pyzSumWidth)
			if (normOfSumRot_ZX < pzxSumWidth)
			{
				//Coincidence counter//
				mol.CoincidenceCount++;
				//PIPICO//
				hi.fill(hiOff+9,"PIPICOCondXYZ",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());

				//MomSumRot//
				hi.fill(hiOff+43,"PxPySumRotCondXYZ",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
				hi.fill(hiOff+44,"PyPzSumRotCondXYZ",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
				hi.fill(hiOff+45,"PzPxSumRotCondXYZ",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
				//Intensity//
				if (intensity.size() && (intPart.size()>1)) 
					hi.fill(hiOff+49,Form("intensity%s%s",p1.GetName(),p2.GetName()),intensity[0], "[arb. unit]",intPart.size()-1,&intPart.front(),"PowerDependence");

				//Momentum first Ion//
				int IDX = hiOff+50;
				hi.fill(IDX+0,Form("%sPxPy",p1.GetName()),p1[i].Px(),p1[i].Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p1[i].Pz()) < 30)
				{
					hi.fill(IDX+1,Form("%sPxPySlice",p1.GetName()),p1[i].Px(),p1[i].Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					//hi.fill(IDX+2,Form("%sPhiPxPySlice",p1.GetName()),TMath::ATan2(p1[i].Py(),p1[i].Px())*TMath::RadToDeg(),TMath::Sqrt(p1[i].Py()*p1[i].Py() + p1[i].Px()*p1[i].Px()),"#phi [deg]","#sqrt{px^{2} + py^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}

				hi.fill(IDX+3,Form("%sPxPz",p1.GetName()),p1[i].Pz(),p1[i].Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p1[i].Py()) < 30)
				{
					hi.fill(IDX+4,Form("%sPxPzSlice",p1.GetName()),p1[i].Pz(),p1[i].Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					//hi.fill(IDX+5,Form("%sPhiPxPzSlice",p1.GetName()),TMath::ATan2(p1[i].Px(),p1[i].Pz())*TMath::RadToDeg(),TMath::Sqrt(p1[i].Pz()*p1[i].Pz() + p1[i].Px()*p1[i].Px()),"#phi [deg]","#sqrt{px^{2} + pz^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}
		   
				hi.fill(IDX+6,Form("%sPyPz",p1.GetName()),p1[i].Pz(),p1[i].Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p1[i].Px()) < 30)
				{
					hi.fill(IDX+7,Form("%sPyPzSlice",p1.GetName()),p1[i].Pz(),p1[i].Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					//hi.fill(IDX+8,Form("%sPhiPyPzSlice",p1.GetName()),TMath::ATan2(p1[i].Py(),p1[i].Pz())*TMath::RadToDeg(),TMath::Sqrt(p1[i].Py()*p1[i].Py() + p1[i].Pz()*p1[i].Pz()),"#phi [deg]","#sqrt{pz^{2} + py^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}

				hi.fill(IDX+9,Form("%sTotalMomentum",p1.GetName()),p1[i].P(),"p [a.u.]",300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				//hi.fill(IDX+2,Form("%sPx",p1.GetName()),p1[i].Px(),"px [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				//hi.fill(IDX+5,Form("%sPy",p1.GetName()),p1[i].Py(),"py [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				//hi.fill(IDX+8,Form("%sPz",p1.GetName()),p1[i].Pz(),"pz [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));

				//Energy First Ion//
				hi.fill(IDX+10,Form("%sEnergy",p1.GetName()),p1[i].E(),"Energy [eV]",200,0,400,Form("%s/Energy",Hname.Data()));

				//Raw
				hi.fill(IDX+11,Form("%sTOF",p1.GetName()),p1[i].TofCor(),"Tof [ns]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,Form("%s/Raw",Hname.Data()));
				hi.fill(IDX+12,Form("%sDetCor",p1.GetName()),p1[i].XCor(),p1[i].YCor(),"x [mm]","y [mm]",300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				hi.fill(IDX+13,Form("%sXPosVsTof",p1.GetName()),p1[i].TofCor(),p1[i].XCorRotScl(),"tof [ns]","x [mm]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p1.GetXcor()-p1.GetCondRad()*1.3,p1.GetXcor()+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				hi.fill(IDX+14,Form("%sYPosVsTof",p1.GetName()),p1[i].TofCor(),p1[i].YCorRotScl(),"tof [ns]","y [mm]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p1.GetYcor()-p1.GetCondRad()*1.3,p1.GetYcor()+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));

				//Angular Distribution//added by motomura
				hi.fill(IDX+21,Form("%sThetaX",p1.GetName()),p1[i].ThetaX(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+22,Form("%sThetaY",p1.GetName()),p1[i].ThetaY(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+23,Form("%sThetaZ",p1.GetName()),p1[i].ThetaZ(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+24,Form("%sThetaYvsEnergy",p1.GetName()),p1[i].ThetaY(),p1[i].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,100,Form("%s/Angular",Hname.Data()));
			//	hi.fill(IDX+25,Form("%sThetaXvsEnergy",p1.GetName()),p1[i].ThetaX(),p1[i].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//hi.fill(IDX+26,Form("%sThetaZvsEnergy",p1.GetName()),p1[i].ThetaZ(),p1[i].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//hi.fill(IDX+27,Form("%sPhiXYvsEnergy",p1.GetName()),p1[i].PhiXY(),p1[i].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
			//	hi.fill(IDX+28,Form("%sPhiYZvsEnergy",p1.GetName()),p1[i].PhiYZ(),p1[i].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//hi.fill(IDX+29,Form("%sPhiZXvsEnergy",p1.GetName()),p1[i].PhiZX(),p1[i].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));

				//if (TMath::Abs(p1[i].Pz()) < 30) 
				//{
				//	hi.fill(IDX+30,Form("%sThetaYSlicePz",p1.GetName()),p1[i].ThetaY(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				//	hi.fill(IDX+31,Form("%sThetaYvsEnergySlicePz",p1.GetName()),p1[i].ThetaY(),p1[i].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//	hi.fill(IDX+32,Form("%sPhiXYvsEnergySlicePz",p1.GetName()),p1[i].PhiXY(),p1[i].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//	hi.fill(IDX+33,Form("%sPhiYZvsEnergySlicePz",p1.GetName()),p1[i].PhiYZ(),p1[i].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//	hi.fill(IDX+34,Form("%sPhiZXvsEnergySlicePz",p1.GetName()),p1[i].PhiZX(),p1[i].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//}

				//if((TMath::Abs(p1[i].PhiZX()) > 35)&&(TMath::Abs(p1[i].PhiZX()) < 150))
				//{
				//	hi.fill(IDX+61,Form("%sTOF",p1.GetName()),p1[i].TofCor(),"Tof [ns]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,Form("%s/KER_masked",Hname.Data()));
				//	hi.fill(IDX+62,Form("%sDetCor",p1.GetName()),p1[i].XCor(),p1[i].YCor(),"x [mm]","y [mm]",300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,Form("%s/KER_masked",Hname.Data()));
				//	hi.fill(IDX+63,Form("%sXPosVsTof",p1.GetName()),p1[i].TofCor(),p1[i].XCorRotScl(),"tof [ns]","x [mm]",500,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p1.GetXcor()-p1.GetCondRad()*1.3,p1.GetXcor()+p1.GetCondRad()*1.3,Form("%s/KER_masked",Hname.Data()));
				//	hi.fill(IDX+64,Form("%sYPosVsTof",p1.GetName()),p1[i].TofCor(),p1[i].YCorRotScl(),"tof [ns]","y [mm]",500,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p1.GetYcor()-p1.GetCondRad()*1.3,p1.GetYcor()+p1.GetCondRad()*1.3,Form("%s/KER_masked",Hname.Data()));
				//	hi.fill(IDX+65,Form("%sPhiZXvsEnergy",p1.GetName()),p1[i].PhiZX(),p1[i].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/KER_masked",Hname.Data()));
				//}

				//Momentum second Ion//
				IDX = (p1!=p2)?IDX+30:IDX;
				hi.fill(IDX+0,Form("%sPxPy",p2.GetName()),p2[j].Px(),p2[j].Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p2[j].Pz()) < 30)
				{
					hi.fill(IDX+1,Form("%sPxPySlice",p2.GetName()),p2[j].Px(),p2[j].Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					//hi.fill(IDX+2,Form("%sPhiPxPySlice",p2.GetName()),TMath::ATan2(p2[j].Py(),p2[j].Px())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Py()*p2[j].Py() + p2[j].Px()*p2[j].Px()),"#phi [deg]","#sqrt{px^{2} + py^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}
		   
				hi.fill(IDX+3,Form("%sPxPz",p2.GetName()),p2[j].Pz(),p2[j].Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p2[j].Py()) < 30)
				{
					hi.fill(IDX+4,Form("%sPxPzSlice",p2.GetName()),p2[j].Pz(),p2[j].Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					//hi.fill(IDX+5,Form("%sPhiPxPzSlice",p2.GetName()),TMath::ATan2(p2[j].Px(),p2[j].Pz())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Pz()*p2[j].Pz() + p2[j].Px()*p2[j].Px()),"#phi [deg]","#sqrt{px^{2} + pz^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}
		   
				hi.fill(IDX+6,Form("%sPyPz",p2.GetName()),p2[j].Pz(),p2[j].Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p2[j].Px()) < 30)
				{
					hi.fill(IDX+7,Form("%sPyPzSlice",p2.GetName()),p2[j].Pz(),p2[j].Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					//hi.fill(IDX+8,Form("%sPhiPyPzSlice",p2.GetName()),TMath::ATan2(p2[j].Py(),p2[j].Pz())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Py()*p2[j].Py() + p2[j].Pz()*p2[j].Pz()),"#phi [deg]","#sqrt{pz^{2} + py^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}

				hi.fill(IDX+9,Form("%sTotalMomentum",p2.GetName()),p2[j].P(),"p [a.u.]",300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				//hi.fill(IDX+2,Form("%sPx",p1.GetName()),p1[i].Px(),"px [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				//hi.fill(IDX+5,Form("%sPy",p1.GetName()),p1[i].Py(),"py [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				//hi.fill(IDX+8,Form("%sPz",p1.GetName()),p1[i].Pz(),"pz [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));

				//Energy second Ion//
				hi.fill(IDX+10,Form("%sEnergy",p2.GetName()),p2[j].E(),"Energy [eV]",200,0,200,Form("%s/Energy",Hname.Data()));

				//Raw
				hi.fill(IDX+11,Form("%sTOF",p2.GetName()),p2[j].TofCor(),"Tof [ns]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Form("%s/Raw",Hname.Data()));
				hi.fill(IDX+12,Form("%sDetCor",p2.GetName()),p2[j].XCor(),p2[j].YCor(),"x [mm]","y [mm]",300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				hi.fill(IDX+13,Form("%sXPosVsTof",p2.GetName()),p2[j].TofCor(),p2[j].XCorRotScl(),"tof [ns]","x [mm]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetXcor()-p2.GetCondRad()*1.3,p2.GetXcor()+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				hi.fill(IDX+14,Form("%sYPosVsTof",p2.GetName()),p2[j].TofCor(),p2[j].YCorRotScl(),"tof [ns]","y [mm]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetYcor()-p2.GetCondRad()*1.3,p2.GetYcor()+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));

				//Angular Distribution//added by motomura
				hi.fill(IDX+21,Form("%sThetaX",p2.GetName()),p2[j].ThetaX(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+22,Form("%sThetaY",p2.GetName()),p2[j].ThetaY(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+23,Form("%sThetaZ",p2.GetName()),p2[j].ThetaZ(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+24,Form("%sThetaYvsEnergy",p2.GetName()),p2[j].ThetaY(),p2[j].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,100,Form("%s/Angular",Hname.Data()));
			//	hi.fill(IDX+25,Form("%sThetaXvsEnergy",p2.GetName()),p2[j].ThetaX(),p2[j].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//hi.fill(IDX+26,Form("%sThetaZvsEnergy",p2.GetName()),p2[j].ThetaZ(),p2[j].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//hi.fill(IDX+27,Form("%sPhiXYvsEnergy",p2.GetName()),p2[j].PhiXY(),p2[j].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
			//	hi.fill(IDX+28,Form("%sPhiYZvsEnergy",p2.GetName()),p2[j].PhiYZ(),p2[j].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//hi.fill(IDX+29,Form("%sPhiZXvsEnergy",p2.GetName()),p2[j].PhiZX(),p2[j].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));

				//if (TMath::Abs(p2[j].Pz()) < 30) 
				//{
				//	hi.fill(IDX+30,Form("%sThetaYSlicePz",p2.GetName()),p2[j].ThetaY(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				//	hi.fill(IDX+31,Form("%sThetaYvsEnergySlicePz",p2.GetName()),p2[j].ThetaY(),p2[j].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//	hi.fill(IDX+32,Form("%sPhiXYvsEnergySlicePz",p2.GetName()),p2[j].PhiXY(),p2[j].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//	hi.fill(IDX+33,Form("%sPhiYZvsEnergySlicePz",p2.GetName()),p2[j].PhiYZ(),p2[j].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//	hi.fill(IDX+34,Form("%sPhiZXvsEnergySlicePz",p2.GetName()),p2[j].PhiZX(),p2[j].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
				//}

				//if((TMath::Abs(p2[j].PhiZX()) > 35)&&(TMath::Abs(p2[j].PhiZX()) < 150))
				//{
				//	hi.fill(IDX+61,Form("%sTOF",p2.GetName()),p2[j].TofCor(),"Tof [ns]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Form("%s/KER_masked",Hname.Data()));
				//	hi.fill(IDX+62,Form("%sDetCor",p2.GetName()),p2[j].XCor(),p2[j].YCor(),"x [mm]","y [mm]",300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,Form("%s/KER_masked",Hname.Data()));
				//	hi.fill(IDX+63,Form("%sXPosVsTof",p2.GetName()),p2[j].TofCor(),p2[j].XCorRotScl(),"tof [ns]","x [mm]",500,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetXcor()-p2.GetCondRad()*1.3,p2.GetXcor()+p2.GetCondRad()*1.3,Form("%s/KER_masked",Hname.Data()));
				//	hi.fill(IDX+64,Form("%sYPosVsTof",p2.GetName()),p2[j].TofCor(),p2[j].YCorRotScl(),"tof [ns]","y [mm]",500,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetYcor()-p2.GetCondRad()*1.3,p2.GetYcor()+p2.GetCondRad()*1.3,Form("%s/KER_masked",Hname.Data()));
				//	hi.fill(IDX+65,Form("%sPhiZXvsEnergy",p2.GetName()),p2[j].PhiZX(),p2[j].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/KER_masked",Hname.Data()));
				//}


				///////////////////KER//////////////////////////
				//hi.fill(IDX+50,Form("KER%s",Hname.Data()),p1[i].E()+p2[j].E(),"KER [eV]",200,0,100,Form("%s/KER",Hname.Data()));
				//hi.fill(IDX+51,Form("ThetaY1vsKER%s",Hname.Data()),p1[i].ThetaY(),p1[i].E()+p2[j].E(),"#theta [deg]","KER [eV]",360,0,180,200,0,100,Form("%s/KER",Hname.Data()));
				//hi.fill(IDX+52,Form("ThetaY2vsKER%s",Hname.Data()),p2[j].ThetaY(),p1[i].E()+p2[j].E(),"#theta [deg]","KER [eV]",360,0,180,200,0,100,Form("%s/KER",Hname.Data()));
				
				//if (((TMath::Abs(p1[i].PhiZX()) > 35)&&(TMath::Abs(p1[i].PhiZX()) < 150))
				//	&&((TMath::Abs(p2[j].PhiZX()) > 35)&&(TMath::Abs(p2[j].PhiZX()) < 150)))
				//{
				//	hi.fill(IDX+60,Form("KER%s_selectPhiXZ",Hname.Data()),p1[i].E()+p2[j].E(),"KER [eV]",200,0,100,Form("%s/KER_masked",Hname.Data()));
				//}

			}
		}
	}
	//std::cout << "\t" << hiOff;
}
void fillMoleculeHistogram2(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff)
{
	TString Hname(p2.GetName());
	Hname += "GatedBy";
	Hname += p1.GetName();

	const double MomLim = 800;

	for (size_t i=0;i<p2.GetNbrOfParticleHits();++i)//i=j
	{
		if (intensity.size())
			hi.fill(hiOff+0,Form("intensity%s%s",p1.GetName(),p2.GetName()),intensity[0], "[arb. unit]",300,0,1000,"Intensity");

		int IDX = hiOff+1;
		hi.fill(IDX+0,Form("%sPxPy",p2.GetName()),p2[i].Px(),p2[i].Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
		if (TMath::Abs(p2[i].Pz()) < 30)
		{
			hi.fill(IDX+1,Form("%sPxPySlice",p2.GetName()),p2[i].Px(),p2[i].Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
			//hi.fill(IDX+2,Form("%sPhiPxPySlice",p2.GetName()),TMath::ATan2(p2[i].Py(),p2[i].Px())*TMath::RadToDeg(),TMath::Sqrt(p2[i].Py()*p2[i].Py() + p2[i].Px()*p2[i].Px()),"#phi [deg]","#sqrt{px^{2} + py^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
		}

		hi.fill(IDX+3,Form("%sPxPz",p2.GetName()),p2[i].Pz(),p2[i].Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
		if (TMath::Abs(p2[i].Py()) < 30)
		{
			hi.fill(IDX+4,Form("%sPxPzSlice",p2.GetName()),p2[i].Pz(),p2[i].Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
			//hi.fill(IDX+5,Form("%sPhiPxPzSlice",p2.GetName()),TMath::ATan2(p2[i].Px(),p2[i].Pz())*TMath::RadToDeg(),TMath::Sqrt(p2[i].Pz()*p2[i].Pz() + p2[i].Px()*p2[i].Px()),"#phi [deg]","#sqrt{px^{2} + pz^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
		}

		hi.fill(IDX+6,Form("%sPyPz",p2.GetName()),p2[i].Pz(),p2[i].Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
		if (TMath::Abs(p2[i].Px()) < 30)
		{
			hi.fill(IDX+7,Form("%sPyPzSlice",p2.GetName()),p2[i].Pz(),p2[i].Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
			//hi.fill(IDX+8,Form("%sPhiPyPzSlice",p2.GetName()),TMath::ATan2(p2[i].Py(),p2[i].Pz())*TMath::RadToDeg(),TMath::Sqrt(p2[i].Py()*p2[i].Py() + p2[i].Pz()*p2[i].Pz()),"#phi [deg]","#sqrt{pz^{2} + py^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
		}

		hi.fill(IDX+9,Form("%sTotalMomentum",p2.GetName()),p2[i].P(),"p [a.u.]",300,0,MomLim,Form("%s/Momenta",Hname.Data()));
		//hi.fill(IDX+2,Form("%sPx",p2.GetName()),p2[i].Px(),"px [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
		//hi.fill(IDX+5,Form("%sPy",p2.GetName()),p2[i].Py(),"py [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
		//hi.fill(IDX+8,Form("%sPz",p2.GetName()),p2[i].Pz(),"pz [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));

		//Energy First Ion//
		hi.fill(IDX+10,Form("%sEnergy",p2.GetName()),p2[i].E(),"Energy [eV]",300,0,100,Form("%s/Energy",Hname.Data()));

		//Raw
		hi.fill(IDX+11,Form("%sTOF",p2.GetName()),p2[i].TofCor(),"Tof [ns]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Form("%s/Raw",Hname.Data()));
		hi.fill(IDX+12,Form("%sDetCor",p2.GetName()),p2[i].XCor(),p2[i].YCor(),"x [mm]","y [mm]",300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
		hi.fill(IDX+13,Form("%sXPosVsTof",p2.GetName()),p2[i].TofCor(),p2[i].XCorRotScl(),"tof [ns]","x [mm]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetXcor()-p2.GetCondRad()*1.3,p2.GetXcor()+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
		hi.fill(IDX+14,Form("%sYPosVsTof",p2.GetName()),p2[i].TofCor(),p2[i].YCorRotScl(),"tof [ns]","y [mm]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetYcor()-p2.GetCondRad()*1.3,p2.GetYcor()+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));

	}
}
//----------------------------------PIPICO ALL-----------------------------------------------------------------//
void fillPIPICO(const MyParticle &p,MyHistos &hi)
{	

	for (size_t i=0; i<p.GetNbrOfParticleHits();++i)
	{
		for (size_t j=i+1;j<p.GetNbrOfParticleHits();++j)
		{
			//PIPICO ALL//
			hi.fill(100,"PIPICO",p[i].TofCor(),p[j].TofCor(),"tof_{firstIon}","tof_{secondIon}",1000,p.GetCondTofFr(),p.GetCondTofTo(),1000,p.GetCondTofFr(),p.GetCondTofTo());
			//if (TMath::Abs(p[i].Px() + p[j].Px()) < pxSumWidth )
			//	if (TMath::Abs(p[i].Py() + p[j].Py()) < pySumWidth )
			//	{
			//		hi.fill(101,"PIPICOcondXY",p[i].TofCor(),p[j].TofCor(),"tof_{firstIon}","tof_{secondIon}",1000,1500,4000,1000,1500,4000);
			//		std::cout<<"sumZ"<<TMath::Abs(p[i].Pz() + p[j].Pz());

			//		if (TMath::Abs(p[i].Pz() + p[j].Pz()) < pzSumWidth )
			//		{
			//			hi.fill(102,"PIPICOcondXYZ",p[i].TofCor(),p[j].TofCor(),"tof_{firstIon}","tof_{secondIon}",1000,1500,4000,1000,1500,4000);
			//		}
			//	}
		}
	}
}

//----------------------------------fill pos condition histograms----------------------------------------------------------------//
void fillParticleConditionsPos(const MyOriginalEvent &oe, const MyDetektor &det, const MyParticle &p, const MyDetektorHit &dh, std::vector<double>& intensity, MyHistos &hi, int hiOff)
{
	//this detektor hit fits the position conditions, so we want to fill a tof histogram//
	const double maxTof	= oe.GetNbrSamples()*oe.GetSampleInterval()*1e9;
	hi.fill(hiOff++,"TofCond",dh.Tof(),"Tof [ns]",5000,0,maxTof,Form("%s/Raw",p.GetName()));
	//for (size_t i=0; i<intensity.size();++i)
	//	hi.fill(hiOff++,Form("Int%dCond",i),dh.Tof(),intensity[i],"Tof [ns]",Form("Int %d",i),300,0,maxTof,300,0,10000,Form("%s/Raw",p.GetName()));
}

//-----------------------------------fill tof condition histograms-----------------------------------------------------//
void fillParticleConditionsTof(const MyOriginalEvent &oe, const MyDetektor &det,const MyParticle &p, const MyDetektorHit &dh, std::vector<double>& intensity, MyHistos &hi, int hiOff)
{
	const double maxPos	= (det.GetRunTime()+30)*det.GetSfU();
	//this detektor hit fits the tof conditions, so we want to fill a detektor histogram//
	hi.fill(hiOff++,"DetCond",dh.X(),dh.Y(),"x [mm]","y [mm]",300,-maxPos,maxPos,300,-maxPos,maxPos,Form("%s/Raw",p.GetName()));
}

//----------------Calcutrate Factrial----------------------------------
double CalcFactorial(int n)
{
	if (n <= 0) return 1.;
	Double_t x = 1;
	Int_t b = 0;
	do {
		b++;
		x *= b;
	} while (b != n);
	return x;
}

void fillHistosAfterAnalyzis(const std::vector<MyParticle> &particles, MyHistos &hi,size_t nRegion)
{
	size_t idx = 9000;
	TH1D* nbrPartHitsNorm = dynamic_cast<TH1D*>(hi.create1d(idx,"NbrParticleHits_normalized","Number of Hits",particles.size()+1,0,particles.size()+1));
	TH1D* nbrPartHitsPois = dynamic_cast<TH1D*>(hi.create1d(idx+1,"NbrParticleHits_Poisson","Number of Hits",particles.size()+1,0,particles.size()+1));
	//TDirectory * rootDir = gDirectory;
	for(size_t i=1; i< particles.size();i++)
	{
		//create a histogram that will take the normalized Data//
		TH1D* nbrHitsNorm = dynamic_cast<TH1D*>(hi.create1d( (idx+20*i) + 1,"NbrOfHits_normalized","Number of Hits",100,0,100,Form("%s",particles[i].GetName())));
		TH1D* nbrHitsPoisson = dynamic_cast<TH1D*>(hi.create1d( (idx+20*i) + 2,"NbrOfHits_Poisson","Number of Hits",100,0,100,Form("%s",particles[i].GetName())));
		//get histogram of Hits
		TH1D* nbrHits = dynamic_cast<TH1D*>( gFile->GetDirectory(Form("%s",particles[i].GetName()))->FindObject("NumberOfHits") );
		TH1D* nbrParticleHits = dynamic_cast<TH1D*>(gFile->FindObject("NumberOfParticleHits"));

		if (nbrHits)
		{
			for (size_t x=0; x<nbrHits->GetXaxis()->GetNbins(); x++)
			{
				//nbrHitsNorm->SetBinContent(x,nbrHits->GetBinContent(x)/nbrHits->GetEntries());
				nbrHitsNorm->Fill(x,nbrHits->GetBinContent(x+1)/nbrHits->GetEntries());
			}
			const double n0 = nbrHitsNorm->GetBinContent(1);
			double mean = 0.;
			if (n0 > 0.) mean = TMath::Log(1/n0);
			else continue;

			for (Int_t j=0; j<100; j++)
			{
				double P = TMath::Exp(-mean)*TMath::Power(mean,j)/CalcFactorial(j);
				nbrHitsPoisson->Fill(j,P);
			}
		}
		nbrPartHitsNorm->Fill(i,nbrHitsNorm->GetMean());
		nbrPartHitsPois->Fill(i,nbrHitsPoisson->GetMean());
	}

}
