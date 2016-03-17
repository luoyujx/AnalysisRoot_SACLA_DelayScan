#include "AnalyzeFuncFill.h"
#include "AnalyzeFunktions.h"
#include "./MyAnalyzer/MyAnalyzer.h"
#include <TFile.h>
#include <TVector2.h>

//TRandom3 RandomGene;
//#define TEST_RANDOM2
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
//-------------------fill particle Histograms---------------------------------------------------//
void fillParticleHistograms(const MyParticle &p, const MyParticleHit &ph, std::vector<double>& intensity, std::vector<double>& delay, MyHistos &hi, int hiOff, std::vector<double>& intPart, int& delayBins, double& delayFrom, double& delayTo, bool& selectThetaZ, double& thetaZLowerLimit, double& thetaZUpperLimit)

{	
	//if (!(ph.E() > p.GetEnergyFrom() && ph.E() < p.GetEnergyTo())) return;
	double MomLim = 1000;
	const double SliceLim = 20;
	TString pname(p.GetName());

	if (pname == "H1p") MomLim = 200;
	if (p.GetCharge_au()>2) MomLim = 1500;
	if (p.GetCharge_au()>5) MomLim = 2000;
	if (p.GetCharge_au()>7) MomLim = 2500;

	double fdelay = 0.0;
	if (intensity.size()) fdelay = delay[2];
	
	double felIntensity = 0.0;
	if (intensity.size() > 1) felIntensity = intensity[0];

	//limit detection angle
	//if (!((ph.PhiZX() > 40 && ph.PhiZX() < 140) || (ph.PhiZX() > -140 && ph.PhiZX() < -40))) return;
	//if (ph.ThetaZ() > 90) return;
	if (selectThetaZ)
	{
		if (ph.ThetaZ()<thetaZLowerLimit) return;
		if (ph.ThetaZ()>thetaZUpperLimit) return;
	}

	//Reconstruction Method for this particle Hit//
	hi.fill(hiOff++,"ReconstructionMethod",ph.RekMeth(),"Reconstruction Nbr",60,0,30,Form("%s/Raw",p.GetName()));
	if (intensity.size() && (intPart.size()>1)) hi.fill(hiOff++,Form("Intensity%s",p.GetName()),intensity[0],"Laser Power",intPart.size()-1,&intPart.front(),Form("%s",p.GetName()));
	else hiOff++;
	if (intensity.size()) hi.fill(hiOff++,"Delay",fdelay,"Delay [fs]",delayBins,delayFrom,delayTo,Form("%s/Delay",p.GetName()));
	else hiOff++;

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

	//delay vs Tof
	hi.fill(hiOff++,"DelayVsTof",ph.TofCor(),fdelay,"tof [ns]","delay [fs]",300,p.GetCondTofFr()-p.GetT0()-p.GetCondTofRange()*0.3,p.GetCondTofTo()-p.GetT0()+p.GetCondTofRange()*0.3,delayBins,delayFrom,delayTo,Form("%s/Delay",p.GetName()));

	//intensity vs Tof
	//hi.fill(hiOff++,"FELIntensityVsTof",ph.TofCor(),felIntensity,"tof [ns]","FEL Intensity [ps]",300,p.GetCondTofFr()-p.GetT0()-p.GetCondTofRange()*0.3,p.GetCondTofTo()-p.GetT0()+p.GetCondTofRange()*0.3,1200,0,1200,Form("%s",p.GetName()));

	//hi.fill(hiOff++,"XPosVsTofFine",ph.TofCor(),ph.XCorRotScl(),"tof [ns]","x [mm]",static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),300,p.GetCondRadX()-p.GetXcor()-p.GetCondRad()*1.3,p.GetCondRadX()-p.GetXcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	//hi.fill(hiOff++,"YPosVsTofFine",ph.TofCor(),ph.YCorRotScl(),"tof [ns]","y [mm]",static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),300,p.GetCondRadY()-p.GetYcor()-p.GetCondRad()*1.3,p.GetCondRadY()-p.GetYcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));

	//momenta//
	hi.fill(hiOff++,"Px",ph.Px(),"px [a.u.]",300,-MomLim*3,MomLim*3,Form("%s/Momenta",p.GetName()));
	hi.fill(hiOff++,"Py",ph.Py(),"py [a.u.]",300,-MomLim*3,MomLim*3,Form("%s/Momenta",p.GetName()));
	hi.fill(hiOff++,"Pz",ph.Pz(),"pz [a.u.]",300,-MomLim*3,MomLim*3,Form("%s/Momenta",p.GetName()));
	//hi.fill(hiOff++,"PzVsBmPosX",ph.Pz(),BmPosX,"pz [a.u.]","Beam position [a.u.]",300,-MomLim,MomLim,300,-0.04,0.04,Form("%s/Momenta",p.GetName()));
	//hi.fill(hiOff++,"PzVsBmPosY",ph.Pz(),BmPosY,"pz [a.u.]","Beam position [a.u.]",300,-MomLim,MomLim,300,-0.04,0.04,Form("%s/Momenta",p.GetName()));
	//hi.fill(hiOff++,"PyVsBmPosX",ph.Py(),BmPosX,"pz [a.u.]","Beam position [a.u.]",300,-MomLim,MomLim,300,-0.04,0.04,Form("%s/Momenta",p.GetName()));
	//hi.fill(hiOff++,"PyVsBmPosY",ph.Py(),BmPosY,"pz [a.u.]","Beam position [a.u.]",300,-MomLim,MomLim,300,-0.04,0.04,Form("%s/Momenta",p.GetName()));

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
	hi.fill(hiOff++,"MomentumXY",TMath::Sqrt(ph.Px()*ph.Px() + ph.Py()*ph.Py()),"pXY [a.u.]",300,0,MomLim,Form("%s/Momenta",p.GetName()));
	
	//delay vs moumenta
	hi.fill(hiOff++, "DelayVsMomentum", ph.P(), fdelay, "p [a.u.]", "delay [fs]", 300, 0, MomLim, delayBins, delayFrom, delayTo, Form("%s/Delay", p.GetName()));
	
	//Energy//
	hi.fill(hiOff++,"Energy",ph.E(),"Energy [eV]",500,0,500,Form("%s/Energy",p.GetName()));
	//hi.fill(hiOff++,"EnergyVsPz",ph.Pz(),ph.E(),"pz [a.u.]","Energy [eV]",300,-MomLim,MomLim,500,0,500,Form("%s/Energy",p.GetName()));
	hi.fill(hiOff++,"DelayVsEnergy",ph.E(),fdelay,"Energy [eV]","delay [fs]",500,0,500,delayBins,delayFrom,delayTo,Form("%s/Delay",p.GetName()));
	//hi.fill(hiOff++,"XFELintensityVsEnergy",ph.E(),intensity[1],"Energy [eV]","XFEL intensity [arb. unit]",150,0,150,60,0,600,Form("%s/Energy",p.GetName()));
	//hi.fill(hiOff++,"DelayVsXFELintensityVsEnergy",ph.E(),delay,intensity[1],"Energy [eV]","delay [ps]","XFEL intensity [arb. unit]",150,0,150,delayBins,delayFrom,delayTo,60,0,600,Form("%s/Delay",p.GetName()));
	//hi.fill(hiOff++,"DelayVsXFELintensityVsEnergy",ph.E(),delay,intensity[1],"Energy [eV]","delay [ps]","XFEL intensity [arb. unit]",500,0,500,32,-5,11,30,0,600,Form("%s/Delay",p.GetName()));

	//Angle//
	hi.fill(hiOff++,"ThetaX",ph.ThetaX(),"#theta [deg]",180,0,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"ThetaY",ph.ThetaY(),"#theta [deg]",180,0,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"ThetaZ",ph.ThetaZ(),"#theta [deg]",180,0,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiXY",ph.PhiXY(),"#phi [deg]",360,-180,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiYZ",ph.PhiYZ(),"#phi [deg]",360,-180,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiZX",ph.PhiZX(),"#phi [deg]",360,-180,180,Form("%s/Angle",p.GetName()));

	hi.fill(hiOff++,"ThetaZvsEnergy",ph.ThetaZ(),ph.E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,200,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"ThetaXvsEnergy",ph.ThetaX(),ph.E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,200,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"ThetaYvsEnergy",ph.ThetaY(),ph.E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,200,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiXYvsEnergy",ph.PhiXY(),ph.E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,200,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiYZvsEnergy",ph.PhiYZ(),ph.E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,200,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiZXvsEnergy",ph.PhiZX(),ph.E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,200,Form("%s/Angle",p.GetName()));

	//hi.fill(hiOff++,"ThetaZvsEnergyNormSinThetaZ",ph.ThetaZ(),ph.E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,40,Form("%s/Angle",p.GetName()),ph.SinThetaZInv());
	//,ph.SinThetaZInv()
	//hi.fill(hiOff,"ThetaY_double",ph.ThetaY(),"#theta [deg]",180,0,360,Form("%s/Angle",p.GetName()));
	//hi.fill(hiOff++,"ThetaY_double",360 - ph.ThetaY(),"#theta [deg]",180,0,360,Form("%s/Angle",p.GetName()));

	//double ThetaY = ph.ThetaY();
	//if (ph.Pz() < 0) ThetaY = 360 - ThetaY;
	//hi.fill(hiOff++,"ThetaY360",ThetaY,"#theta [deg]",360,0,360,Form("%s/Angle",p.GetName()));

	//Intensity//
	//for (size_t i=0; i< intensity.size();++i)
	//	hi.fill(hiOff++,Form("Int%d",i),intensity[i],Form("Int %d",i),1000,0,1000,Form("%s/Intensity",p.GetName()))
}
//-------------------Fill molecule (2-body coincidence) histogram-----------------------------------------------------------------------------------------------------//
void fillMoleculeHistogramCH2I2(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff, Molecule &mol, std::vector<double>& delay, int& delayBins, double& delayFrom, double& delayTo)
{
	double MomLim = 1600;
	double MomSumRotLim = 1600;
	double MomScale = mol.momSumFactor;
	const int Mombins = 150;

	TString Hname(p1.GetName());
	Hname += p2.GetName();

	//if (mol.momSumWindowX * mol.momSumWindowY * mol.momSumWindowZ < 1e-20) return;
	const double pxSumWidth = mol.momSumWindowX;//10,9,7,5
	const double pySumWidth = mol.momSumWindowY;//10,8,5,4
	const double pzSumWidth = mol.momSumWindowZ;//5,6,4,2
	const double pxySumWidth = mol.momSumWindowX;//10,9,7,5
	const double pyzSumWidth = mol.momSumWindowY;//10,8,5,4
	const double pzxSumWidth = mol.momSumWindowZ;//5,6,4,2

	double fdelay = 0.0;
	if (delay.size()) fdelay = delay[2];

	for (size_t i = 0; i<p1.GetNbrOfParticleHits(); ++i)
	{
		for (size_t j = 0; j<p2.GetNbrOfParticleHits(); ++j)//i=j
		{
			//skip if we are going to put in the same particlehit, for the case that we want coincidences for the same particle//
			if ((i>=j) && (p1==p2)) continue;//i==j
			//---p1[i]="C" p2[j]="I"---//
			const double angleP1P2 = calcFormedAngle(p1[i], p2[j]);
			double weightPerSin = 1 / TMath::Sin(angleP1P2*TMath::DegToRad());
			if (weightPerSin > 100) weightPerSin = 100;
			const double ratioP1P2 = p1[i].P() / p2[j].P();
			const double angleP1P2XY = calcFormedAngleXY(p1[i], p2[j]);
			const double ratioP1P2XY = calcMagXY(p1[i]) / calcMagXY(p2[j]);

			hi.fill(hiOff + 17, "Delay", fdelay, "delay [fs]", delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
			hi.fill(hiOff + 18, "AngleVsRatio", angleP1P2, ratioP1P2, "angle", "Ratio", 180, 0, 180, 100, 0, 2, Form("%s/FormedAngle", Hname.Data()), weightPerSin);
			hi.fill(hiOff + 19, "AngleVsRatioXY", angleP1P2XY, ratioP1P2XY, "angle", "Ratio", 180, 0, 180, 100, 0, 2, Form("%s/FormedAngle", Hname.Data()));
			hi.fill(hiOff + 20, "Angle", angleP1P2, "angle", 180, 0, 180, Form("%s/FormedAngle", Hname.Data()), weightPerSin);
			hi.fill(hiOff + 21, "DelayVsAngle", angleP1P2, fdelay, "angle", "delay [fs]", 180, 0, 180, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()), weightPerSin);
			//1st condition
			if (angleP1P2 < mol.angleCondition) continue;
			if ((ratioP1P2 < mol.momSumFactorLow) || (ratioP1P2 > mol.momSumFactorUp)) continue;

			//Calculate momentum sum condition
			const double p12SumX = p1[i].Px() / MomScale + p2[j].Px();
			const double p12SumY = p1[i].Py() / MomScale + p2[j].Py();
			const double p12SumZ = p1[i].Pz() / MomScale + p2[j].Pz();


			//calcutrate momentum sum by Edwin's method
			const double ThetaXY = TMath::ATan2(p2[j].Py(), p2[j].Px());
			//const double ThetaXY = TMath::ATan2(p1[i].Py(),p1[i].Px());
			//Edwin's rotation
			//double ThetaXY = TMath::ASin(p2[j].Py()/calcMagXY(p2[j]));
			//if (p2[j].Px() < 0) ThetaXY = ThetaXY/TMath::Abs(ThetaXY)*(TMath::Pi()-TMath::Abs(ThetaXY));
			const double p1x_XY = p1[i].PxRotXY(ThetaXY);
			const double p1y_XY = p1[i].PyRotXY(ThetaXY);
			const double p2x_XY = p2[j].PxRotXY(ThetaXY);
			const double p2y_XY = p2[j].PyRotXY(ThetaXY);

			const double ThetaYZ = TMath::ATan2(p2[j].Pz(), p2[j].Py());
			//const double ThetaYZ = TMath::ATan2(p1[i].Pz(),p1[i].Py());
			const double p1y_YZ = p1[i].PyRotYZ(ThetaYZ);
			const double p1z_YZ = p1[i].PzRotYZ(ThetaYZ);
			const double p2y_YZ = p2[j].PyRotYZ(ThetaYZ);
			const double p2z_YZ = p2[j].PzRotYZ(ThetaYZ);

			const double ThetaZX = TMath::ATan2(p2[j].Px(), p2[j].Pz());
			//const double ThetaZX = TMath::ATan2(p1[i].Px(),p1[i].Pz());
			const double p1z_ZX = p1[i].PzRotZX(ThetaZX);
			const double p1x_ZX = p1[i].PxRotZX(ThetaZX);
			const double p2z_ZX = p2[j].PzRotZX(ThetaZX);
			const double p2x_ZX = p2[j].PxRotZX(ThetaZX);

			const double p12xSumRot_XY = p1x_XY / MomScale + p2x_XY;
			const double p12ySumRot_XY = p1y_XY / MomScale + p2y_XY;
			const double p12ySumRot_YZ = p1y_YZ / MomScale + p2y_YZ;
			const double p12zSumRot_YZ = p1z_YZ / MomScale + p2z_YZ;
			const double p12zSumRot_ZX = p1z_ZX / MomScale + p2z_ZX;
			const double p12xSumRot_ZX = p1x_ZX / MomScale + p2x_ZX;

			double normOfSumRot_XY = TMath::Sqrt(p12xSumRot_XY*p12xSumRot_XY + p12ySumRot_XY*p12ySumRot_XY);
			double normOfSumRot_YZ = TMath::Sqrt(p12ySumRot_YZ*p12ySumRot_YZ + p12zSumRot_YZ*p12zSumRot_YZ);
			double normOfSumRot_ZX = TMath::Sqrt(p12zSumRot_ZX*p12zSumRot_ZX + p12xSumRot_ZX*p12xSumRot_ZX);

			hi.fill(hiOff + 31, "PxPyRot", p1x_XY, p1y_XY, "px [a.u.]", "py [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));
			hi.fill(hiOff + 31, "PxPyRot", p2x_XY, p2y_XY, "px [a.u.]", "py [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));
			hi.fill(hiOff + 32, "PxPySumRot", p12xSumRot_XY, p12ySumRot_XY, "pxySum [a.u.]", "pxySum [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));

			hi.fill(hiOff+33,"PyPzRot",p1y_YZ, p1z_YZ,"py [a.u.]","pz [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+33,"PyPzRot",p2y_YZ, p2z_YZ,"py [a.u.]","pz [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff + 34, "PyPzSumRot", p12ySumRot_YZ, p12zSumRot_YZ, "pyzSum [a.u.]", "pyzSum [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));

			hi.fill(hiOff+35,"PzPxRot",p1z_ZX, p1x_ZX,"pz [a.u.]","px [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+35,"PzPxRot",p2z_ZX, p2x_ZX,"pz [a.u.]","px [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff + 36, "PzPxSumRot", p12zSumRot_ZX, p12xSumRot_ZX, "pzxSum [a.u.]", "pzxSum [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));

			//define the norm of momentum sum vector that rotated
			//const double normOfSumRot_XY = TMath::Sqrt(p12xSumRot_XY*p12xSumRot_XY + p12ySumRot_XY*p12ySumRot_XY);
			//const double normOfSumRot_YZ = TMath::Sqrt(p12ySumRot_YZ*p12ySumRot_YZ + p12zSumRot_YZ*p12zSumRot_YZ);
			//const double normOfSumRot_ZX = TMath::Sqrt(p12zSumRot_ZX*p12zSumRot_ZX + p12xSumRot_ZX*p12xSumRot_ZX);
			//divide by norm?

			//-----Gated stuff-----//
			// momentum sum//
			hi.fill(hiOff + 1, "PxSum", (p1[i].Px() / MomScale + p2[j].Px()), Form("Px_{%s} + Px_{%s} [a.u]", p1.GetName(), p2.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));
			if (normOfSumRot_YZ < pyzSumWidth)
				hi.fill(hiOff + 2, "PxSumCondYZ", (p1[i].Px() / MomScale + p2[j].Px()), Form("Px_{%s} + Px_{%s} [a.u]", p1.GetName(), p2.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));
			//y momentum sum//
			hi.fill(hiOff + 3, "PySum", (p1[i].Py() / MomScale + p2[j].Py()), Form("Py_{%s} + Py_{%s} [a.u]", p1.GetName(), p2.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));
			if (normOfSumRot_ZX < pzxSumWidth)
				hi.fill(hiOff + 4, "PySumCondXZ", (p1[i].Py() / MomScale + p2[j].Py()), Form("Py_{%s} + Py_{%s} [a.u]", p1.GetName(), p2.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));
			//z momentum sum//
			hi.fill(hiOff + 5, "PzSum", (p1[i].Pz() / MomScale + p2[j].Pz()), Form("Pz_{%s} + Pz_{%s} [a.u]", p1.GetName(), p2.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));
			if (normOfSumRot_XY < pxySumWidth)
				hi.fill(hiOff + 6, "PzSumCondXY", (p1[i].Pz() / MomScale + p2[j].Pz()), Form("Pz_{%s} + Pz_{%s} [a.u]", p1.GetName(), p2.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));

			const double pipicoFactor = 0.1;
			//PIPICO//
			hi.fill(hiOff + 7, "PIPICO", p1[i].TofCor(), p2[j].TofCor(), Form("tof_{%s} [ns]", p1.GetName()), Form("tof_{%s} [ns]", p2.GetName()), 300, p1.GetCondTofFr() - p1.GetT0() - p1.GetCondTofRange()*pipicoFactor, p1.GetCondTofTo() - p1.GetT0() + p1.GetCondTofRange()*pipicoFactor, 300, p2.GetCondTofFr() - p2.GetT0() - p2.GetCondTofRange()*pipicoFactor, p2.GetCondTofTo() - p2.GetT0() + p2.GetCondTofRange()*pipicoFactor, Hname.Data());
			//hi.fill(hiOff+10,"PIPICOMom",p1[i].Pz(),p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			//hi.fill(hiOff+11,"PIPICOMomFactor",p1[i].Pz()/MomScale,p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			if (normOfSumRot_XY < pxySumWidth)
			{
				hi.fill(hiOff + 8, "PIPICOCondXY", p1[i].TofCor(), p2[j].TofCor(), Form("tof_{%s} [ns]", p1.GetName()), Form("tof_{%s} [ns]", p2.GetName()), 300, p1.GetCondTofFr() - p1.GetT0() - p1.GetCondTofRange()*pipicoFactor, p1.GetCondTofTo() - p1.GetT0() + p1.GetCondTofRange()*pipicoFactor, 300, p2.GetCondTofFr() - p2.GetT0() - p2.GetCondTofRange()*pipicoFactor, p2.GetCondTofTo() - p2.GetT0() + p2.GetCondTofRange()*pipicoFactor, Hname.Data());
				hi.fill(hiOff + 12, "RatioPCPerPICondXY", ratioP1P2, "Ratio", 200, 0, 2, Form("%s/MomSums", Hname.Data()));
				hi.fill(hiOff + 13, "RatioXYPCPerPICondXY", ratioP1P2XY, "Ratio", 200, 0, 2, Form("%s/MomSums", Hname.Data()));
				hi.fill(5, "PIPICO_CondXY", p1[i].Tof(), p2[j].Tof(), "tof [ns]", "tof [ns]", 3000, 1000, 3500, 3000, 1000, 3500);
				//hi.fill(hiOff+12,"PIPICOMomCondXY",p1[i].Pz(),p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
				//hi.fill(hiOff+13,"PIPICOMomFactorCondXY",p1[i].Pz()/MomScale,p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			}

			//For confirmation. gated by only one plane
			if (normOfSumRot_XY < pxySumWidth)
			{
				hi.fill(hiOff + 37, "PzPxSumRotCondXY", p12zSumRot_ZX, p12xSumRot_ZX, "pzxSum [a.u.]", "pzxSum [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));
				hi.fill(hiOff + 38, "PyPzSumRotCondXY", p12ySumRot_YZ, p12zSumRot_YZ, "pyzSum [a.u.]", "pyzSum [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));

				//hi.fill(hiOff + 24, "AngleVsRatioXYCondXY", angleP1P2XY, ratioP1P2XY, "angle", "Ratio", 180, 0, 180, 100, 0, 2, Form("%s/MomSums", Hname.Data()));
				//hi.fill(hiOff+25,Form("%sXPosVsTofCondXY",p1.GetName()),p1[i].TofCor(),p1[i].XCorRotScl(),"tof [ns]","x [mm]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*pipicoFactor,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*pipicoFactor,300,p1.GetXcor()-p1.GetCondRad()*1.3,p1.GetXcor()+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				//hi.fill(hiOff+26,Form("%sYPosVsTofCondXY",p1.GetName()),p1[i].TofCor(),p1[i].YCorRotScl(),"tof [ns]","y [mm]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*pipicoFactor,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*pipicoFactor,300,p1.GetYcor()-p1.GetCondRad()*1.3,p1.GetYcor()+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				//hi.fill(hiOff+27,Form("%sXPosVsTofCondXY",p2.GetName()),p2[j].TofCor(),p2[j].XCorRotScl(),"tof [ns]","x [mm]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*pipicoFactor,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*pipicoFactor,300,p2.GetXcor()-p2.GetCondRad()*1.3,p2.GetXcor()+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				//hi.fill(hiOff+28,Form("%sYPosVsTofCondXY",p2.GetName()),p2[j].TofCor(),p2[j].YCorRotScl(),"tof [ns]","y [mm]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*pipicoFactor,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*pipicoFactor,300,p2.GetYcor()-p2.GetCondRad()*1.3,p2.GetYcor()+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				//hi.fill(hiOff+29,Form("%sDetCorCondXY",p1.GetName()),p1[i].XCor(),p1[i].YCor(),"x [mm]","y [mm]",300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				//hi.fill(hiOff+30,Form("%sDetCorCondXY",p2.GetName()),p2[j].XCor(),p2[j].YCor(),"x [mm]","y [mm]",300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
			}
			//if (normOfSumRot_YZ < pyzSumWidth)
			//{
			//	hi.fill(hiOff+39,"PxPySumRotCondYZ",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			//	hi.fill(hiOff+40,"PzPxSumRotCondYZ",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			//}
			//if (normOfSumRot_ZX < pzxSumWidth)
			//{
			//	hi.fill(hiOff+41,"PyPzSumRotCondZX",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			//	hi.fill(hiOff+42,"PxPySumRotCondZX",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			//}
			///////////////////apply two sum condition//////////////////////
			if (normOfSumRot_XY < pxySumWidth)
				if (normOfSumRot_YZ < pyzSumWidth)
				{
					hi.fill(hiOff + 39, "PzPxSumRotCondXY-YZ", p12zSumRot_ZX, p12xSumRot_ZX, "pzxSum [a.u.]", "pzxSum [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));
					//hi.fill(hiOff + 40, "NormOfSumRot_ZX_Gated", normOfSumRot_ZX, "norm [a.u.]", 1000, 0, 0.001, Form("%s/MomSums", Hname.Data()), 1 / (2 * TMath::Pi()*normOfSumRot_ZX));
				}
			if (normOfSumRot_YZ < pyzSumWidth)
				if (normOfSumRot_ZX < pzxSumWidth)
				{
					hi.fill(hiOff + 41, "PxPySumRotCondYZ-ZX", p12xSumRot_XY, p12ySumRot_XY, "pxySum [a.u.]", "pxySum [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));
					//hi.fill(hiOff + 42, "NormOfSumRot_XY_Gated", normOfSumRot_XY, "norm [a.u.]", 1000, 0, 0.001, Form("%s/MomSums", Hname.Data()), 1 / (2 * TMath::Pi()*normOfSumRot_XY));
				}
			if (normOfSumRot_ZX < pzxSumWidth)
				if (normOfSumRot_XY < pxySumWidth)
				{
					hi.fill(hiOff + 43, "PyPzSumRotCondZX-XY", p12ySumRot_YZ, p12zSumRot_YZ, "pyzSum [a.u.]", "pyzSum [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));
					//hi.fill(hiOff + 45, "NormOfSumRot_YZ_Gated", normOfSumRot_YZ, "norm [a.u.]", 1000, 0, 0.001, Form("%s/MomSums", Hname.Data()), 1 / (2 * TMath::Pi()*normOfSumRot_YZ));
				}

			///////////////////3-plane Momentum sum condition//////////////////////
			if (normOfSumRot_XY < pxySumWidth)
				if (normOfSumRot_YZ < pyzSumWidth)
					if (normOfSumRot_ZX < pzxSumWidth)
					{
						//Coincidence counter//
						mol.CoincidenceCount++;
						mol.CoinHitNbrC.push_back(i);
						mol.CoinHitNbrI.push_back(j);

						//Coincidence Tof//
						hi.fill(1, "TOF_Coincidence", p1[i].Tof(), "Tof [ns]", 4000, 0, 4000);
						hi.fill(1, "TOF_Coincidence", p2[j].Tof(), "Tof [ns]", 4000, 0, 4000);
						hi.fill(2, "XTOF_Coincidence", p1[i].Tof(), p1[i].X(), "tof [ns]", "x [mm]", 4000, 0, 4000, 300, -50, 50);
						hi.fill(2, "XTOF_Coincidence", p2[j].Tof(), p2[j].X(), "tof [ns]", "x [mm]", 4000, 0, 4000, 300, -50, 50);
						hi.fill(3, "YTOF_Coincidence", p1[i].Tof(), p1[i].Y(), "tof [ns]", "y [mm]", 4000, 0, 4000, 300, -50, 50);
						hi.fill(3, "YTOF_Coincidence", p2[j].Tof(), p2[j].Y(), "tof [ns]", "y [mm]", 4000, 0, 4000, 300, -50, 50);
						hi.fill(4, "PIPICO_CondXYZ", p1[i].Tof(), p2[j].Tof(), "tof [ns]", "tof [ns]", 3000, 1000, 3500, 3000, 1000, 3500);

						//PIPICO//
						hi.fill(hiOff + 9, "PIPICOCondXYZ", p1[i].TofCor(), p2[j].TofCor(), Form("tof_{%s} [ns]", p1.GetName()), Form("tof_{%s} [ns]", p2.GetName()), 300, p1.GetCondTofFr() - p1.GetT0() - p1.GetCondTofRange()*pipicoFactor, p1.GetCondTofTo() - p1.GetT0() + p1.GetCondTofRange()*pipicoFactor, 300, p2.GetCondTofFr() - p2.GetT0() - p2.GetCondTofRange()*pipicoFactor, p2.GetCondTofTo() - p2.GetT0() + p2.GetCondTofRange()*pipicoFactor, Hname.Data());
						hi.fill(hiOff + 10, "PzSumCondXYZ", (p1[i].Pz() / MomScale + p2[j].Pz()), Form("Pz_{%s} + Pz_{%s} [a.u]", p1.GetName(), p2.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));

						//MomSumRot//
						hi.fill(hiOff + 46, "PxPySumRotCondXYZ", p12xSumRot_XY, p12ySumRot_XY, "pxySum [a.u.]", "pxySum [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));
						hi.fill(hiOff + 47, "PyPzSumRotCondXYZ", p12ySumRot_YZ, p12zSumRot_YZ, "pyzSum [a.u.]", "pyzSum [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));
						hi.fill(hiOff + 48, "PzPxSumRotCondXYZ", p12zSumRot_ZX, p12xSumRot_ZX, "pzxSum [a.u.]", "pzxSum [a.u.]", 200, -MomSumRotLim, MomSumRotLim, 200, -MomSumRotLim, MomSumRotLim, Form("%s/MomSums", Hname.Data()));
						//Intensity//
						//hi.fill(hiOff + 49, Form("intensity%s%s", p1.GetName(), p2.GetName()), intensity[0], "[arb. unit]", intPart.size() - 1, &intPart.front(), "PowerDependence");
						//Ratio PI/PC
						hi.fill(hiOff + 15, "RatioPCPerPICondXYZ", ratioP1P2, "Ratio", 200, 0, 2, Form("%s/MomSums", Hname.Data()));
						hi.fill(hiOff + 16, "RatioXYPCPerPICondXYZ", ratioP1P2XY, "Ratio", 200, 0, 2, Form("%s/MomSums", Hname.Data()));
						//Formed angle VS Ratio
						hi.fill(hiOff + 24, "AngleVsRatioCondXYZ", angleP1P2, ratioP1P2, "angle", "Ratio", 180, 0, 180, 100, 0, 2, Form("%s/FormedAngle", Hname.Data()), weightPerSin);
						//hi.fill(hiOff + 23, "AngleVsRatioXYCondXYZ", angleP1P2XY, ratioP1P2XY, "angle", "Ratio", 180, 0, 180, 100, 0, 2, Form("%s/FormedAngle", Hname.Data()));
						hi.fill(hiOff + 22, "AngleCondXYZ", angleP1P2, "angle", 180, 0, 180, Form("%s/FormedAngle", Hname.Data()), weightPerSin);
						hi.fill(hiOff + 23, "DelayVsAngleCondXYZ", angleP1P2, fdelay, "angle", "delay [fs]", 180, 0, 180, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()), weightPerSin);

						//Momentum first Ion//
					
						int IDX = hiOff + 50;
						hi.fill(IDX + 0, Form("%sPxPy", p1.GetName()), p1[i].Px(), p1[i].Py(), "px [a.u.]", "py [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
						if (TMath::Abs(p1[i].Pz()) < 30)
						{
							hi.fill(IDX + 1, Form("%sPxPySlice", p1.GetName()), p1[i].Px(), p1[i].Py(), "px [a.u.]", "py [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
						}

						hi.fill(IDX + 3, Form("%sPzPx", p1.GetName()), p1[i].Pz(), p1[i].Px(), "pz [a.u.]", "px [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
						if (TMath::Abs(p1[i].Py()) < 30)
						{
							hi.fill(IDX + 4, Form("%sPxPzSlice", p1.GetName()), p1[i].Pz(), p1[i].Px(), "pz [a.u.]", "px [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
						}

						hi.fill(IDX + 6, Form("%sPyPz", p1.GetName()), p1[i].Pz(), p1[i].Py(), "pz [a.u.]", "py [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
						if (TMath::Abs(p1[i].Px()) < 30)
						{
							hi.fill(IDX + 7, Form("%sPyPzSlice", p1.GetName()), p1[i].Pz(), p1[i].Py(), "pz [a.u.]", "py [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
						}
						//Momentum2D
						hi.fill(IDX + 2, Form("%sMomentumXY", p1.GetName()), TMath::Sqrt(p1[i].Px()*p1[i].Px() + p1[i].Py()*p1[i].Py()), "pXY [a.u.]", Mombins, 0, MomLim, Form("%s/Momenta", Hname.Data()));
						hi.fill(IDX + 5, Form("%sMomentumXZ", p1.GetName()), TMath::Sqrt(p1[i].Px()*p1[i].Px() + p1[i].Pz()*p1[i].Pz()), "pXY [a.u.]", Mombins, 0, MomLim, Form("%s/Momenta", Hname.Data()));
						hi.fill(IDX + 8, Form("%sMomentumYZ", p1.GetName()), TMath::Sqrt(p1[i].Py()*p1[i].Py() + p1[i].Pz()*p1[i].Pz()), "pXY [a.u.]", Mombins, 0, MomLim, Form("%s/Momenta", Hname.Data()));
						//hi.fill(IDX + 21, Form("%sMomentumXY_k", p1.GetName()), TMath::Sqrt(p1[i].Px()*p1[i].Px() + p1[i].Py()*p1[i].Py()) / MomScale, "pXY [a.u.]", Mombins, 0, MomLim, Form("%s/Momenta", Hname.Data()));

						hi.fill(IDX + 9, Form("%sTotalMomentum", p1.GetName()), p1[i].P(), "p [a.u.]", Mombins, 0, MomLim, Form("%s/Momenta", Hname.Data()));
						//hi.fill(IDX + 22, Form("%sTotalMomentum_k", p1.GetName()), p1[i].P() / MomScale, "p [a.u.]", Mombins, 0, MomLim, Form("%s/Momenta", Hname.Data()));
						//hi.fill(IDX+2,Form("%sPx",p1.GetName()),p1[i].Px(),"px [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+5,Form("%sPy",p1.GetName()),p1[i].Py(),"py [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+8,Form("%sPz",p1.GetName()),p1[i].Pz(),"pz [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));

						//Energy First Ion//
						hi.fill(IDX + 10, Form("%sEnergy", p1.GetName()), p1[i].E(), "Energy [eV]", 200, 0, 400, Form("%s/Energy", Hname.Data()));

						//Raw
						hi.fill(IDX + 11, Form("%sTOF", p1.GetName()), p1[i].TofCor(), "Tof [ns]", 300, p1.GetCondTofFr() - p1.GetT0() - p1.GetCondTofRange()*0.3, p1.GetCondTofTo() - p1.GetT0() + p1.GetCondTofRange()*0.3, Form("%s/Raw", Hname.Data()));
						hi.fill(IDX + 12, Form("%sDetCor", p1.GetName()), p1[i].XCor(), p1[i].YCor(), "x [mm]", "y [mm]", 300, 0 - p1.GetCondRad()*1.3, 0 + p1.GetCondRad()*1.3, 300, 0 - p1.GetCondRad()*1.3, 0 + p1.GetCondRad()*1.3, Form("%s/Raw", Hname.Data()));
						hi.fill(IDX + 13, Form("%sXPosVsTof", p1.GetName()), p1[i].TofCor(), p1[i].XCorRotScl(), "tof [ns]", "x [mm]", 300, p1.GetCondTofFr() - p1.GetT0() - p1.GetCondTofRange()*0.3, p1.GetCondTofTo() - p1.GetT0() + p1.GetCondTofRange()*0.3, 300, p1.GetXcor() - p1.GetCondRad()*1.3, p1.GetXcor() + p1.GetCondRad()*1.3, Form("%s/Raw", Hname.Data()));
						hi.fill(IDX + 14, Form("%sYPosVsTof", p1.GetName()), p1[i].TofCor(), p1[i].YCorRotScl(), "tof [ns]", "y [mm]", 300, p1.GetCondTofFr() - p1.GetT0() - p1.GetCondTofRange()*0.3, p1.GetCondTofTo() - p1.GetT0() + p1.GetCondTofRange()*0.3, 300, p1.GetYcor() - p1.GetCondRad()*1.3, p1.GetYcor() + p1.GetCondRad()*1.3, Form("%s/Raw", Hname.Data()));

						//Angular Distribution
						//hi.fill(IDX+21,Form("%sThetaX",p1.GetName()),p1[i].ThetaX(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						//hi.fill(IDX+22,Form("%sThetaY",p1.GetName()),p1[i].ThetaY(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						//hi.fill(IDX+23,Form("%sThetaZ",p1.GetName()),p1[i].ThetaZ(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX + 24, Form("%sThetaYvsEnergy", p1.GetName()), p1[i].ThetaY(), p1[i].E(), "#theta [deg]", "Energy [eV]", 180, 0, 180, 200, 0, 200, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 25, Form("%sThetaXvsEnergy", p1.GetName()), p1[i].ThetaX(), p1[i].E(), "#theta [deg]", "Energy [eV]", 180, 0, 180, 200, 0, 200, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 26, Form("%sThetaZvsEnergy", p1.GetName()), p1[i].ThetaZ(), p1[i].E(), "#theta [deg]", "Energy [eV]", 180, 0, 180, 200, 0, 200, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 27, Form("%sPhiXYvsEnergy", p1.GetName()), p1[i].PhiXY(), p1[i].E(), "#phi [deg]", "Energy [eV]", 180, -180, 180, 200, 0, 200, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 28, Form("%sPhiYZvsEnergy", p1.GetName()), p1[i].PhiYZ(), p1[i].E(), "#phi [deg]", "Energy [eV]", 180, -180, 180, 200, 0, 200, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 29, Form("%sPhiZXvsEnergy", p1.GetName()), p1[i].PhiZX(), p1[i].E(), "#phi [deg]", "Energy [eV]", 180, -180, 180, 200, 0, 200, Form("%s/Angular", Hname.Data()));

						//Momentum second Ion//
						IDX = (p1 != p2) ? IDX + 30 : IDX;
						hi.fill(IDX + 0, Form("%sPxPy", p2.GetName()), p2[j].Px(), p2[j].Py(), "px [a.u.]", "py [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
						if (TMath::Abs(p2[j].Pz()) < 30)
						{
							hi.fill(IDX + 1, Form("%sPxPySlice", p2.GetName()), p2[j].Px(), p2[j].Py(), "px [a.u.]", "py [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
							//hi.fill(IDX+2,Form("%sPhiPxPySlice",p2.GetName()),TMath::ATan2(p2[j].Py(),p2[j].Px())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Py()*p2[j].Py() + p2[j].Px()*p2[j].Px()),"#phi [deg]","#sqrt{px^{2} + py^{2}} [a.u.]",360,-180,180,Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						}

						hi.fill(IDX + 3, Form("%sPzPx", p2.GetName()), p2[j].Pz(), p2[j].Px(), "pz [a.u.]", "px [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
						if (TMath::Abs(p2[j].Py()) < 30)
						{
							hi.fill(IDX + 4, Form("%sPxPzSlice", p2.GetName()), p2[j].Pz(), p2[j].Px(), "pz [a.u.]", "px [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
							//hi.fill(IDX+5,Form("%sPhiPxPzSlice",p2.GetName()),TMath::ATan2(p2[j].Px(),p2[j].Pz())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Pz()*p2[j].Pz() + p2[j].Px()*p2[j].Px()),"#phi [deg]","#sqrt{px^{2} + pz^{2}} [a.u.]",360,-180,180,Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						}

						hi.fill(IDX + 6, Form("%sPyPz", p2.GetName()), p2[j].Pz(), p2[j].Py(), "pz [a.u.]", "py [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
						if (TMath::Abs(p2[j].Px()) < 30)
						{
							hi.fill(IDX + 7, Form("%sPyPzSlice", p2.GetName()), p2[j].Pz(), p2[j].Py(), "pz [a.u.]", "py [a.u.]", Mombins, -MomLim, MomLim, Mombins, -MomLim, MomLim, Form("%s/Momenta", Hname.Data()));
							//hi.fill(IDX+8,Form("%sPhiPyPzSlice",p2.GetName()),TMath::ATan2(p2[j].Py(),p2[j].Pz())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Py()*p2[j].Py() + p2[j].Pz()*p2[j].Pz()),"#phi [deg]","#sqrt{pz^{2} + py^{2}} [a.u.]",360,-180,180,Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						}
						//Momentum2D
						hi.fill(IDX + 2, Form("%sMomentumXY", p2.GetName()), TMath::Sqrt(p2[j].Px()*p2[j].Px() + p2[j].Py()*p2[j].Py()), "pXY [a.u.]", Mombins, 0, MomLim, Form("%s/Momenta", Hname.Data()));
						hi.fill(IDX + 5, Form("%sMomentumXZ", p2.GetName()), TMath::Sqrt(p2[j].Px()*p2[j].Px() + p2[j].Pz()*p2[j].Pz()), "pXY [a.u.]", Mombins, 0, MomLim, Form("%s/Momenta", Hname.Data()));
						hi.fill(IDX + 8, Form("%sMomentumYZ", p2.GetName()), TMath::Sqrt(p2[j].Py()*p2[j].Py() + p2[j].Pz()*p2[j].Pz()), "pXY [a.u.]", Mombins, 0, MomLim, Form("%s/Momenta", Hname.Data()));

						hi.fill(IDX + 9, Form("%sTotalMomentum", p2.GetName()), p2[j].P(), "p [a.u.]", 300, 0, MomLim, Form("%s/Momenta", Hname.Data()));
						//hi.fill(IDX+2,Form("%sPx",p1.GetName()),p1[i].Px(),"px [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+5,Form("%sPy",p1.GetName()),p1[i].Py(),"py [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+8,Form("%sPz",p1.GetName()),p1[i].Pz(),"pz [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));

						//Energy second Ion//
						hi.fill(IDX + 10, Form("%sEnergy", p2.GetName()), p2[j].E(), "Energy [eV]", 200, 0, 200, Form("%s/Energy", Hname.Data()));

						//Raw
						hi.fill(IDX + 11, Form("%sTOF", p2.GetName()), p2[j].TofCor(), "Tof [ns]", 300, p2.GetCondTofFr() - p2.GetT0() - p2.GetCondTofRange()*0.3, p2.GetCondTofTo() - p2.GetT0() + p2.GetCondTofRange()*0.3, Form("%s/Raw", Hname.Data()));
						hi.fill(IDX + 12, Form("%sDetCor", p2.GetName()), p2[j].XCor(), p2[j].YCor(), "x [mm]", "y [mm]", 300, 0 - p2.GetCondRad()*1.3, 0 + p2.GetCondRad()*1.3, 300, 0 - p2.GetCondRad()*1.3, 0 + p2.GetCondRad()*1.3, Form("%s/Raw", Hname.Data()));
						hi.fill(IDX + 13, Form("%sXPosVsTof", p2.GetName()), p2[j].TofCor(), p2[j].XCorRotScl(), "tof [ns]", "x [mm]", 300, p2.GetCondTofFr() - p2.GetT0() - p2.GetCondTofRange()*0.3, p2.GetCondTofTo() - p2.GetT0() + p2.GetCondTofRange()*0.3, 300, p2.GetXcor() - p2.GetCondRad()*1.3, p2.GetXcor() + p2.GetCondRad()*1.3, Form("%s/Raw", Hname.Data()));
						hi.fill(IDX + 14, Form("%sYPosVsTof", p2.GetName()), p2[j].TofCor(), p2[j].YCorRotScl(), "tof [ns]", "y [mm]", 300, p2.GetCondTofFr() - p2.GetT0() - p2.GetCondTofRange()*0.3, p2.GetCondTofTo() - p2.GetT0() + p2.GetCondTofRange()*0.3, 300, p2.GetYcor() - p2.GetCondRad()*1.3, p2.GetYcor() + p2.GetCondRad()*1.3, Form("%s/Raw", Hname.Data()));

						//Angular Distribution
						hi.fill(IDX + 21, Form("%sThetaX", p2.GetName()), p2[j].ThetaX(), "#theta [deg]", 36, 0, 180, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 22, Form("%sThetaY", p2.GetName()), p2[j].ThetaY(), "#theta [deg]", 36, 0, 180, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 23, Form("%sThetaZ", p2.GetName()), p2[j].ThetaZ(), "#theta [deg]", 36, 0, 180, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 24, Form("%sThetaYvsEnergy", p2.GetName()), p2[j].ThetaY(), p2[j].E(), "#theta [deg]", "Energy [eV]", 180, 0, 180, 200, 0, 100, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 25, Form("%sThetaXvsEnergy", p2.GetName()), p2[j].ThetaX(), p2[j].E(), "#theta [deg]", "Energy [eV]", 180, 0, 180, 200, 0, 100, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 26, Form("%sThetaZvsEnergy", p2.GetName()), p2[j].ThetaZ(), p2[j].E(), "#theta [deg]", "Energy [eV]", 180, 0, 180, 200, 0, 100, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 27, Form("%sPhiXYvsEnergy", p2.GetName()), p2[j].PhiXY(), p2[j].E(), "#phi [deg]", "Energy [eV]", 180, -180, 180, 200, 0, 100, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 28, Form("%sPhiYZvsEnergy", p2.GetName()), p2[j].PhiYZ(), p2[j].E(), "#phi [deg]", "Energy [eV]", 180, -180, 180, 200, 0, 100, Form("%s/Angular", Hname.Data()));
						hi.fill(IDX + 29, Form("%sPhiZXvsEnergy", p2.GetName()), p2[j].PhiZX(), p2[j].E(), "#phi [deg]", "Energy [eV]", 180, -180, 180, 200, 0, 100, Form("%s/Angular", Hname.Data()));

						///////////////////KER//////////////////////////
						hi.fill(IDX + 30, Form("KESum_%s", Hname.Data()), p1[i].E() + p2[j].E(), "KE [eV]", 200, 0, 250, Form("%s/KESum", Hname.Data()));
						hi.fill(IDX + 31, "DelayVsKESum", p1[i].E() + p2[j].E(), fdelay, "KE [eV]", "delay [fs]", 200, 0, 250, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
						hi.fill(IDX + 32, "DelayCondXYZ", fdelay, "delay [fs]", delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
						//hi.fill(IDX+31,Form("KER_IperQ_%s",Hname.Data()),
						//	p2[j].E()*9.447078522/(p1.GetCharge_au()*(p2.GetCharge_au()+3)),
						//	"KER [eV]",200,0,20,Form("%s/KER",Hname.Data()));
						//hi.fill(IDX+32,Form("KEIplusKEC_%s",Hname.Data()),p1[i].E()+p2[j].E(),"KE [eV]",200,0,500,Form("%s/KER",Hname.Data()));
						//hi.fill(IDX+33,Form("KEIplusKECperQ_%s",Hname.Data()),
						//	(p1[i].E()+p2[j].E())/(p1.GetCharge_au()*p2.GetCharge_au()),
						//	"KE [eV]",200,0,20,Form("%s/KER",Hname.Data()));

						//momentum of H
						//TVector3 h3Pvec = p1[i].Pvec() + p2[j].Pvec();
						//hi.fill(IDX+30,"MomentumPC_PI",h3Pvec.Mag(),"p [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+31,"MomentumPC_PI-H",h3Pvec.Mag()/3/TMath::Cos(70*TMath::DegToRad()),"p [a.u.]",Mombins,0,MomLim/2,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+32,"PxPyPC_PI",h3Pvec.X(),h3Pvec.Y(),"px [a.u.]","py [a.u.]",Mombins,-MomLim/2,MomLim/2,Mombins,-MomLim/2,MomLim/2,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+33,"PyPzPC_PI",h3Pvec.Y(),h3Pvec.Z(),"py [a.u.]","pz [a.u.]",Mombins,-MomLim/2,MomLim/2,Mombins,-MomLim/2,MomLim/2,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+34,"PzPxPC_PI",h3Pvec.Z(),h3Pvec.X(),"pz [a.u.]","px [a.u.]",Mombins,-MomLim/2,MomLim/2,Mombins,-MomLim/2,MomLim/2,Form("%s/Momenta",Hname.Data()));
					}
		}
	}
	//std::cout << "\t" << hiOff;

}


//-------------------Fill molecule CH2I2(3-body coincidence) histogram-----------------------------------------------------------------------------------------------------//
void fillMoleculeHistogramCH2I2_3body(const MyParticle &p1, const MyParticle &p2, const MyParticle &p3, std::vector<double>& intensity, MyHistos &hi, int hiOff, std::map< std::string, Molecule > molecule3, std::vector<double>& delay, int& delayBins, double& delayFrom, double& delayTo)
{
	double MomLim = 1600;
	double MomSumRotLim = 1600;
//	double MomScale = mol.momSumFactor;
	const int Mombins = 150;

	TString Hname(p1.GetName());
	Hname += p2.GetName();
	Hname += p3.GetName();

	Molecule mol = molecule3[Hname.Data()];
	if (mol.momSumWindowX * mol.momSumWindowY * mol.momSumWindowZ < 1e-20) return;
	const double pxSumWidth = mol.momSumWindowX;//10,9,7,5
	const double pySumWidth = mol.momSumWindowY;//10,8,5,4
	const double pzSumWidth = mol.momSumWindowZ;//5,6,4,2
	//const double pxySumWidth = mol.momSumWindowX;//10,9,7,5
	//const double pyzSumWidth = mol.momSumWindowY;//10,8,5,4
	//const double pzxSumWidth = mol.momSumWindowZ;//5,6,4,2
	//std::cout << Hname << std::endl;
	const double angleCondition = mol.angleCondition;
	const double momSumFactor = 1.0 / mol.momSumFactor;
	const double momSumFactorLow = mol.momSumFactorLow;
	const double momSumFactorUp = mol.momSumFactorUp;
	//std::cout << Hname << momSumFactor << " " << pxSumWidth << ":" << pySumWidth << ":" << pzSumWidth << std::endl;

	double fdelay = 0.0;
	if (delay.size()) fdelay = delay[2];


	for (size_t i = 0; i < p1.GetNbrOfParticleHits(); ++i)
	{
		for (size_t j = 0; j < p2.GetNbrOfParticleHits(); ++j)//i=j
		{
			for (size_t k = 0; k < p3.GetNbrOfParticleHits(); ++k)//i=j
			{
				//skip if we are going to put in the same particlehit, for the case that we want coincidences for the same particle//
				if ((j >= k) && (p2 == p3)) continue;//i==j
				//---p1[i]="C" p2[j]="I" p3[k]="I"---//
				const TVector3 &pvecC = p1[i].Pvec();
				const TVector3 &pvecI1 = p2[j].Pvec();
				const TVector3 &pvecI2 = p3[k].Pvec();
				const TVector3 pvecSumII = pvecI1 + pvecI2;//molecule axis
				
				const double angleC_II = pvecC.Angle(pvecSumII) * TMath::RadToDeg();
				double weightPerSin = 1 / TMath::Sin(angleC_II * TMath::DegToRad());
				if (weightPerSin > 100) weightPerSin = 100;
				const double ratioC_II = pvecC.Mag() / pvecSumII.Mag();

				hi.fill(hiOff + 0, "Delay", fdelay, "delay [fs]", delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
				hi.fill(hiOff + 1, "AngleVsRatio", angleC_II, ratioC_II, "angle", "Ratio", 180, 0, 180, 100, 0, 2, Form("%s/FormedAngle", Hname.Data()), weightPerSin);
				hi.fill(hiOff + 2, "Angle", angleC_II, "angle", 180, 0, 180, Form("%s/FormedAngle", Hname.Data()), weightPerSin);
				hi.fill(hiOff + 3, "DelayVsAngle", angleC_II, fdelay, "angle", "delay [fs]", 180, 0, 180, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()), weightPerSin);
				hi.fill(hiOff + 4, "MomRatioCII", ratioC_II, "Ratio", 200, 0, 2, Form("%s/FormedAngle", Hname.Data()));

				//1st condition
				if (angleC_II < angleCondition) continue;
				if ((ratioC_II < momSumFactorLow) || (ratioC_II > momSumFactorUp)) continue;

				int hiOffMom = hiOff + 5;
				const TVector3 pvecSumCII = momSumFactor * pvecC + pvecSumII;
				hi.fill(hiOffMom + 0, "PSum", pvecSumCII.Mag(), Form("P_{%s} + P_{%s} + P_{%s} [a.u]", p1.GetName(), p2.GetName(), p3.GetName()), 200, 0, 800, Form("%s/MomSums", Hname.Data()));
				hi.fill(hiOffMom + 1, "PxSum", pvecSumCII.X(), Form("Px_{%s} + Px_{%s} + Px_{%s} [a.u]", p1.GetName(), p2.GetName(), p3.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));
				hi.fill(hiOffMom + 2, "PySum", pvecSumCII.Y(), Form("Py_{%s} + Py_{%s} + Py_{%s}[a.u]", p1.GetName(), p2.GetName(), p3.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));
				hi.fill(hiOffMom + 3, "PzSum", pvecSumCII.Z(), Form("Pz_{%s} + Pz_{%s} + Pz_{%s}[a.u]", p1.GetName(), p2.GetName(), p3.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));
				if ((pvecSumCII.Y() < pySumWidth) && (pvecSumCII.Z() < pzSumWidth))
					hi.fill(hiOffMom + 4, "PxSumCondYZ", pvecSumCII.X(), Form("Px_{%s} + Px_{%s} + Px_{%s} [a.u]", p1.GetName(), p2.GetName(), p3.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));
				if ((pvecSumCII.X() < pxSumWidth) && (pvecSumCII.Z() < pzSumWidth))
					hi.fill(hiOffMom + 5, "PySumCondXZ", pvecSumCII.Y(), Form("Py_{%s} + Py_{%s} + Py_{%s}[a.u]", p1.GetName(), p2.GetName(), p3.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));
				if ((pvecSumCII.X() < pxSumWidth) && (pvecSumCII.Y() < pySumWidth))
					hi.fill(hiOffMom + 6, "PzSumCondXY", pvecSumCII.Z(), Form("Pz_{%s} + Pz_{%s} + Pz_{%s}[a.u]", p1.GetName(), p2.GetName(), p3.GetName()), 200, -400, 400, Form("%s/MomSums", Hname.Data()));

				int hiOffResult = hiOffMom + 7;
				if ((pvecSumCII.X() < pxSumWidth) && (pvecSumCII.Y() < pySumWidth) && (pvecSumCII.Z() < pzSumWidth))
				{
					const double angleII = pvecI1.Angle(pvecI2)*TMath::RadToDeg();
					const double angleCI1 = pvecC.Angle(pvecI1)*TMath::RadToDeg();
					const double angleCI2 = pvecC.Angle(pvecI2)*TMath::RadToDeg();

					hi.fill(hiOffResult + 0, "DelayCondXYZ", fdelay, "delay [fs]", delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					//Ratio Pc/(Pi1+Pi2)
					hi.fill(hiOffResult + 1, "MomRatioCIICondXYZ", ratioC_II, "Ratio", 200, 0, 2, Form("%s/FormedAngle", Hname.Data()));
					hi.fill(hiOffResult + 2, "MomRatioIICondXYZ", pvecI1.Mag() / pvecI2.Mag(), "Ratio", 200, 0, 2, Form("%s/FormedAngle", Hname.Data()));
					hi.fill(hiOffResult + 3, "MomRatioCI1CondXYZ", pvecC.Mag() / pvecI1.Mag(), "Ratio", 200, 0, 2, Form("%s/FormedAngle", Hname.Data()));
					hi.fill(hiOffResult + 4, "MomRatioCI2CondXYZ", pvecC.Mag() / pvecI2.Mag(), "Ratio", 200, 0, 2, Form("%s/FormedAngle", Hname.Data()));

					// Formed angle
					hi.fill(hiOffResult + 5, "AngleCondXYZ", angleC_II, "angle", 180, 0, 180, Form("%s/FormedAngle", Hname.Data()));
					hi.fill(hiOffResult + 6, "AngleIICondXYZ", angleII, "angle", 180, 0, 180, Form("%s/FormedAngle", Hname.Data()));
					hi.fill(hiOffResult + 7, "AngleCI1CondXYZ", angleCI1, "angle", 180, 0, 180, Form("%s/FormedAngle", Hname.Data()));
					hi.fill(hiOffResult + 8, "AngleCI2CondXYZ", angleCI2, "angle", 180, 0, 180, Form("%s/FormedAngle", Hname.Data()));

					//Formed angle VS Ratio
					hi.fill(hiOffResult + 9, "AngleVsRatioCondXYZ", angleC_II, ratioC_II, "angle", "Ratio", 180, 0, 180, 100, 0, 2, Form("%s/FormedAngle", Hname.Data()), weightPerSin);
					hi.fill(hiOffResult + 10, "AngleVsIIvsRatioCIICondXYZ", angleII, ratioC_II, "angle", "Ratio", 180, 0, 180, 100, 0, 2, Form("%s/FormedAngle", Hname.Data()));
					hi.fill(hiOffResult + 11, "AngleCI1vsRatioCIICondXYZ", angleCI1, ratioC_II, "angle", "Ratio", 180, 0, 180, 100, 0, 2, Form("%s/FormedAngle", Hname.Data()));
					hi.fill(hiOffResult + 12, "AngleCI2vsRatioCIICondXYZ", angleCI2, ratioC_II, "angle", "Ratio", 180, 0, 180, 100, 0, 2, Form("%s/FormedAngle", Hname.Data()));
					
					// delay vs angle
					hi.fill(hiOffResult + 13, "DelayVsAngleCondXYZ", angleC_II, fdelay, "angle", "delay [fs]", 180, 0, 180, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()), weightPerSin);
					hi.fill(hiOffResult + 14, "DelayVsAngleIICondXYZ", angleII, fdelay, "angle", "delay [fs]", 180, 0, 180, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					hi.fill(hiOffResult + 15, "DelayVsAngleCI1CondXYZ", angleCI1, fdelay, "angle", "delay [fs]", 180, 0, 180, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					hi.fill(hiOffResult + 16, "DelayVsAngleCI2CondXYZ", angleCI2, fdelay, "angle", "delay [fs]", 180, 0, 180, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					hi.fill(hiOffResult + 17, "DelayVsAngleCICondXYZ", angleCI1, fdelay, "angle", "delay [fs]", 180, 0, 180, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					hi.fill(hiOffResult + 17, "DelayVsAngleCICondXYZ", angleCI2, fdelay, "angle", "delay [fs]", 180, 0, 180, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));

					//Sum of kinetic energy
					hi.fill(hiOffResult + 20, "DelayVsKESumCII", p1[i].E() + p2[j].E() + p3[k].E(), fdelay, "KE [eV]", "delay [fs]", 200, 0, 250, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					hi.fill(hiOffResult + 21, "DelayVsKESumII", p2[j].E() + p3[k].E(), fdelay, "KE [eV]", "delay [fs]", 200, 0, 250, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					hi.fill(hiOffResult + 22, "DelayVsKESumCI1", p1[i].E() + p2[j].E(), fdelay, "KE [eV]", "delay [fs]", 200, 0, 250, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					hi.fill(hiOffResult + 23, "DelayVsKESumCI2", p1[i].E() + p3[k].E(), fdelay, "KE [eV]", "delay [fs]", 200, 0, 250, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					hi.fill(hiOffResult + 24, "DelayVsKESumCI", p1[i].E() + p2[j].E(), fdelay, "KE [eV]", "delay [fs]", 200, 0, 250, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					hi.fill(hiOffResult + 24, "DelayVsKESumCI", p1[i].E() + p3[k].E(), fdelay, "KE [eV]", "delay [fs]", 200, 0, 250, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					
					// KE
					hi.fill(hiOffResult + 25, "DelayVsKEC", p1[i].E(), fdelay, "KE [eV]", "delay [fs]", 200, 0, 250, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					hi.fill(hiOffResult + 26, "DelayVsKEI1",p2[j].E(), fdelay, "KE [eV]", "delay [fs]", 200, 0, 250, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));
					hi.fill(hiOffResult + 27, "DelayVsKEI2",p3[k].E(), fdelay, "KE [eV]", "delay [fs]", 200, 0, 250, delayBins, delayFrom, delayTo, Form("%s/DelayDep", Hname.Data()));

					double angleTheta = pvecI1.Angle(pvecI2);
					double angle3Body = pvecI1.Cross(pvecI2).Dot(pvecC) / (pvecI1.Mag()*pvecI2.Mag()*pvecC.Mag()*TMath::Sin(angleTheta));
					hi.fill(hiOffResult + 28, Form("3BodyAngleCos-%s-%s-%s", p1.GetName(), p2.GetName(), p3.GetName()), angle3Body, "cos(Phi)", 100, -1, 1, "3BodyAngularCorr");
					// Newton diagram
					double I1 = pvecI1.Mag();
					TVector3 vI1(1, 0, 0);
					TVector3 vI2(pvecI2.Mag() / I1, 0, 0);
					TVector3 vC(pvecC.Mag() / I1, 0, 0);
					vI2.RotateZ(pvecI1.Angle(pvecI2));
					vC.RotateZ(-pvecI1.Angle(pvecC));
					//hi.fill(hiOffResult + 29, Form("NewtonDia-%s-%s-%s", p1.GetName(), p2.GetName(), p3.GetName()),
					//	vI1.X(), vI1.Y(), "Normalized Momentum", "Normalized Momentum", 200, -1.5, 1.5, 200, -1.5, 1.5, "NewtonDaiagram");
					hi.fill(hiOffResult + 29, Form("NewtonDia-%s-%s-%s", p1.GetName(), p2.GetName(), p3.GetName()), 
						vI2.X(), vI2.Y(), "Normalized Momentum", "Normalized Momentum", 200, -1.5, 1.5, 200, -1.5, 1.5, "NewtonDaiagram");
					hi.fill(hiOffResult + 29, Form("NewtonDia-%s-%s-%s", p1.GetName(), p2.GetName(), p3.GetName()),
						vC.X(), vC.Y(), "Normalized Momentum", "Normalized Momentum", 200, -1.5, 1.5, 200, -1.5, 1.5, "NewtonDaiagram");

					hi.fill(hiOffResult + 30, "KE_C", p1[i].E(), "KE [eV]", 200, 0, 200, Form("%s/KE", Hname.Data()));
					hi.fill(hiOffResult + 31, "KE_I1", p2[j].E(), "KE [eV]", 200, 0, 200, Form("%s/KE", Hname.Data()));
					hi.fill(hiOffResult + 32, "KE_I2", p3[k].E(), "KE [eV]", 200, 0, 200, Form("%s/KE", Hname.Data()));

					//double I1 = pvecI1.Mag();
					//TVector3 pvecI1N = (1.0 / I1) * pvecI1;
					//TVector3 pvecI2N = (1.0 / I1) * pvecI2;
					//TVector3 pvecCN = (1.0 / I1) * pvecC;
					//const TVector3 crossI1I2 = pvecI1.Cross(pvecI2);
					//const TVector3 zAxis(0, 0, 1);
					//const TVector3 xAxis(1, 0, 0);
					//double angleZaxis = crossI1I2.Angle(zAxis);
					//pvecI1N.Rotate(-angleZaxis, crossI1I2.Cross(zAxis));
					//pvecI2N.Rotate(-angleZaxis, crossI1I2.Cross(zAxis));
					//pvecCN.Rotate(-angleZaxis, crossI1I2.Cross(zAxis));
					//double angleXaxis = pvecI1N.Angle(xAxis);
					//pvecI1N.RotateZ(-angleXaxis);
					//pvecI2N.RotateZ(-angleXaxis);
					//pvecCN.RotateZ(-angleXaxis);
					//hi.fill(hiOffResult + 29, Form("NewtonDia-%s-%s-%s", p1.GetName(), p2.GetName(), p3.GetName()),
					//	pvecI1N.X(), pvecI1N.Y(), "Momenta", "Momenta", 200, -1.5, 1.5, 200, -1.5, 1.5, "NewtonDaiagram");
					//hi.fill(hiOffResult + 29, Form("NewtonDia-%s-%s-%s", p1.GetName(), p2.GetName(), p3.GetName()), 
					//	pvecI2N.X(), pvecI2N.Y(), "Momenta", "Momenta", 200, -1.5, 1.5, 200, -1.5, 1.5, "NewtonDaiagram");
					//hi.fill(hiOffResult + 29, Form("NewtonDia-%s-%s-%s", p1.GetName(), p2.GetName(), p3.GetName()),
					//	pvecCN.X(), pvecCN.Y(), "Momenta", "Momenta", 200, -1.5, 1.5, 200, -1.5, 1.5, "NewtonDaiagram");

				}
			}

		}
	}
}

//-------------------Fill Analog histogram-----------------------------------------------------------------------------------------------------//
void fillAnalogHistogram(const MyOriginalEvent &oe, MyHistos &hi, const MyParticle &p, double delay, int delayBins, double delayFrom, double delayTo)
{

}
//-------------------Fill molecule (2-body coincidence) histogram-----------------------------------------------------------------------------------------------------//
void fillMoleculeHistogram(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff, Molecule &mol, std::vector<double>& intPart)
{

	double MomLim = 800;
	double MomSumRotLim = 800;
	double MomScale = mol.momSumFactor;
	const int Mombins = 150;

	TString Hname(p1.GetName());
	if (Hname=="H1p") 
	{
		MomLim = 200;
		MomSumRotLim = 300;
		MomScale = mol.momSumFactor;
	}
	Hname += p2.GetName();
	
	if (mol.momSumWindowX * mol.momSumWindowY * mol.momSumWindowZ < 1e-20) return;
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
			//if ((i>=j) && (p1==p2)) continue;//i==j

			//---p1[i]="C" p2[j]="I"---//
			const double angleP1P2 = calcFormedAngle(p1[i],p2[j]);
			double weightPerSin = 1/TMath::Sin(angleP1P2*TMath::DegToRad());
			if (weightPerSin > 100) weightPerSin = 100;
			const double ratioP1P2 = p1[i].P()/p2[j].P();
			const double angleP1P2XY = calcFormedAngleXY(p1[i],p2[j]);
			const double ratioP1P2XY = calcMagXY(p1[i])/calcMagXY(p2[j]);


			hi.fill(hiOff+20,"AngleVsRatio",angleP1P2, ratioP1P2,"angle","Ratio",180,0,180,100,0,2, Form("%s/MomSums",Hname.Data()),weightPerSin);
			hi.fill(hiOff+21,"AngleVsRatioXY",angleP1P2XY, ratioP1P2XY,"angle","Ratio",180,0,180,100,0,2, Form("%s/MomSums",Hname.Data()));

			//1st condition
			if (angleP1P2 < mol.angleCondition) continue;
			if ((ratioP1P2 < mol.momSumFactorLow) || (ratioP1P2 > mol.momSumFactorUp)) continue;

			//Calculate momentum sum condirion (divide by norm)
			const double p12SumX = p1[i].Px()/MomScale + p2[j].Px();
			const double p12SumY = p1[i].Py()/MomScale + p2[j].Py();
			const double p12SumZ = p1[i].Pz()/MomScale + p2[j].Pz();

			const double pminLength = 5e-22 / 1.992851E-24;
			const double normXY = (p1[i].Px()/MomScale)*(p1[i].Px()/MomScale) + p2[j].Px()*p2[j].Px() + (p1[i].Py()/MomScale)*(p1[i].Py()/MomScale) + p2[j].Py()*p2[j].Py() + 2*pminLength*pminLength;
			const double normYZ = (p1[i].Py()/MomScale)*(p1[i].Py()/MomScale) + p2[j].Py()*p2[j].Py() + (p1[i].Pz()/MomScale)*(p1[i].Pz()/MomScale) + p2[j].Pz()*p2[j].Pz() + 2*pminLength*pminLength;
			const double normZX = (p1[i].Pz()/MomScale)*(p1[i].Pz()/MomScale) + p2[j].Pz()*p2[j].Pz() + (p1[i].Px()/MomScale)*(p1[i].Px()/MomScale) + p2[j].Px()*p2[j].Px() + 2*pminLength*pminLength;

			const double normOfSumRot_XY = (p12SumX*p12SumX + p12SumY*p12SumY)/normXY;
			const double normOfSumRot_YZ = (p12SumY*p12SumY + p12SumZ*p12SumZ)/normYZ;
			const double normOfSumRot_ZX = (p12SumZ*p12SumZ + p12SumX*p12SumX)/normZX;

			hi.fill(hiOff+25,"NormOfSumRot_ZX",normOfSumRot_ZX,"norm [a.u.]",100,0,0.001, Form("%s/MomSums",Hname.Data()),1/(2*TMath::Pi()*normOfSumRot_ZX));
			hi.fill(hiOff+26,"NormOfSumRot_XY",normOfSumRot_XY,"norm [a.u.]",100,0,0.001, Form("%s/MomSums",Hname.Data()),1/(2*TMath::Pi()*normOfSumRot_XY));
			hi.fill(hiOff+27,"NormOfSumRot_YZ",normOfSumRot_YZ,"norm [a.u.]",100,0,0.001, Form("%s/MomSums",Hname.Data()),1/(2*TMath::Pi()*normOfSumRot_YZ));

			//calcutrate momentum sum by Edwin's method
			//MyParticleHit pC(p1[i]);
			//pC.MultiplyP(1.0/MomScale);
			const double ThetaXY = TMath::ATan2(p2[j].Py(),p2[j].Px());
			//const double ThetaXY = TMath::ATan2(p1[i].Py(),p1[i].Px());
			//Edwin's rotation
			//double ThetaXY = TMath::ASin(p2[j].Py()/calcMagXY(p2[j]));
			//if (p2[j].Px() < 0) ThetaXY = ThetaXY/TMath::Abs(ThetaXY)*(TMath::Pi()-TMath::Abs(ThetaXY));
			const double p1x_XY = p1[i].PxRotXY(ThetaXY);
			const double p1y_XY = p1[i].PyRotXY(ThetaXY);
			//const double p1x_XY = pC.PxRotXY(ThetaXY);
			//const double p1y_XY = pC.PyRotXY(ThetaXY);
			const double p2x_XY = p2[j].PxRotXY(ThetaXY);
			const double p2y_XY = p2[j].PyRotXY(ThetaXY);

			const double ThetaYZ = TMath::ATan2(p2[j].Pz(),p2[j].Py());
			//const double ThetaYZ = TMath::ATan2(p1[i].Pz(),p1[i].Py());
			const double p1y_YZ = p1[i].PyRotYZ(ThetaYZ);
			const double p1z_YZ = p1[i].PzRotYZ(ThetaYZ);
			const double p2y_YZ = p2[j].PyRotYZ(ThetaYZ);
			const double p2z_YZ = p2[j].PzRotYZ(ThetaYZ);

			const double ThetaZX = TMath::ATan2(p2[j].Px(),p2[j].Pz());
			//const double ThetaZX = TMath::ATan2(p1[i].Px(),p1[i].Pz());
			const double p1z_ZX = p1[i].PzRotZX(ThetaZX);
			const double p1x_ZX = p1[i].PxRotZX(ThetaZX);
			const double p2z_ZX = p2[j].PzRotZX(ThetaZX);
			const double p2x_ZX = p2[j].PxRotZX(ThetaZX);

			//const double p12xSumRot_XY = p1x_XY+p2x_XY;
			//const double p12ySumRot_XY = p1y_XY+p2y_XY;

			const double p12xSumRot_XY = p1x_XY/MomScale+p2x_XY;
			const double p12ySumRot_XY = p1y_XY/MomScale+p2y_XY;
			const double p12ySumRot_YZ = p1y_YZ/MomScale+p2y_YZ;
			const double p12zSumRot_YZ = p1z_YZ/MomScale+p2z_YZ;
			const double p12zSumRot_ZX = p1z_ZX/MomScale+p2z_ZX;
			const double p12xSumRot_ZX = p1x_ZX/MomScale+p2x_ZX;

			//double normOfSumRot_XY = TMath::Sqrt(p12xSumRot_XY*p12xSumRot_XY + p12ySumRot_XY*p12ySumRot_XY);
			//double normOfSumRot_YZ = TMath::Sqrt(p12ySumRot_YZ*p12ySumRot_YZ + p12zSumRot_YZ*p12zSumRot_YZ);
			//double normOfSumRot_ZX = TMath::Sqrt(p12zSumRot_ZX*p12zSumRot_ZX + p12xSumRot_ZX*p12xSumRot_ZX);

			//const double normXY = TMath::Sqrt((p1x_XY/MomScale)*(p1x_XY/MomScale) + p2x_XY*p2x_XY + (p1y_XY/MomScale)*(p1y_XY/MomScale) + p2y_XY*p2y_XY);
			//const double normYZ = TMath::Sqrt((p1y_YZ/MomScale)*(p1y_YZ/MomScale) + p2y_YZ*p2y_YZ + (p1z_YZ/MomScale)*(p1z_YZ/MomScale) + p2z_YZ*p2z_YZ);
			//const double normZX = TMath::Sqrt((p1z_ZX/MomScale)*(p1z_ZX/MomScale) + p2z_ZX*p2z_ZX + (p1x_ZX/MomScale)*(p1x_ZX/MomScale) + p2x_ZX*p2x_ZX);

			//double normOfSumRot_XY = p12xSumRot_XY*p12xSumRot_XY + p12ySumRot_XY*p12ySumRot_XY;
			//double normOfSumRot_YZ = p12ySumRot_YZ*p12ySumRot_YZ + p12zSumRot_YZ*p12zSumRot_YZ;
			//double normOfSumRot_ZX = p12zSumRot_ZX*p12zSumRot_ZX + p12xSumRot_ZX*p12xSumRot_ZX;

			//const double normXY = (p1x_XY/MomScale)*(p1x_XY/MomScale) + p2x_XY*p2x_XY + (p1y_XY/MomScale)*(p1y_XY/MomScale) + p2y_XY*p2y_XY;
			//const double normYZ = (p1y_YZ/MomScale)*(p1y_YZ/MomScale) + p2y_YZ*p2y_YZ + (p1z_YZ/MomScale)*(p1z_YZ/MomScale) + p2z_YZ*p2z_YZ;
			//const double normZX = (p1z_ZX/MomScale)*(p1z_ZX/MomScale) + p2z_ZX*p2z_ZX + (p1x_ZX/MomScale)*(p1x_ZX/MomScale) + p2x_ZX*p2x_ZX;

			//const double normXY = (p1[i].Px()/MomScale)*(p1[i].Px()/MomScale) + (p1[i].Py()/MomScale)*(p1[i].Py()/MomScale)

			hi.fill(hiOff+31,"PxPyRot",p1x_XY, p1y_XY,"px [a.u.]","py [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+31,"PxPyRot",p2x_XY, p2y_XY,"px [a.u.]","py [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+32,"PxPySumRot",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",400,-MomSumRotLim,MomSumRotLim,400,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));

			//hi.fill(hiOff+33,"PyPzRot",p1y_YZ, p1z_YZ,"py [a.u.]","pz [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			//hi.fill(hiOff+33,"PyPzRot",p2y_YZ, p2z_YZ,"py [a.u.]","pz [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+34,"PyPzSumRot",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",400,-MomSumRotLim,MomSumRotLim,400,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));

			//hi.fill(hiOff+35,"PzPxRot",p1z_ZX, p1x_ZX,"pz [a.u.]","px [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			//hi.fill(hiOff+35,"PzPxRot",p2z_ZX, p2x_ZX,"pz [a.u.]","px [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			hi.fill(hiOff+36,"PzPxSumRot",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",400,-MomSumRotLim,MomSumRotLim,400,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			
			//define the norm of momentum sum vector that rotated
			//const double normOfSumRot_XY = TMath::Sqrt(p12xSumRot_XY*p12xSumRot_XY + p12ySumRot_XY*p12ySumRot_XY);
			//const double normOfSumRot_YZ = TMath::Sqrt(p12ySumRot_YZ*p12ySumRot_YZ + p12zSumRot_YZ*p12zSumRot_YZ);
			//const double normOfSumRot_ZX = TMath::Sqrt(p12zSumRot_ZX*p12zSumRot_ZX + p12xSumRot_ZX*p12xSumRot_ZX);
			//divide by norm?

			//-----Gated stuff-----//
			// momentum sum//
			hi.fill(hiOff+1,"PxSum",(p1[i].Px()/MomScale + p2[j].Px()),Form("Px_{%s} + Px_{%s} [a.u]",p1.GetName(),p2.GetName()),200,-400,400,Form("%s/MomSums",Hname.Data()));
			if (normOfSumRot_YZ < pyzSumWidth)
				hi.fill(hiOff+2,"PxSumCondYZ",(p1[i].Px()/MomScale +p2[j].Px()),Form("Px_{%s} + Px_{%s} [a.u]",p1.GetName(),p2.GetName()),200,-400,400, Form("%s/MomSums",Hname.Data()));
			//y momentum sum//
			hi.fill(hiOff+3,"PySum",(p1[i].Py()/MomScale + p2[j].Py()),Form("Py_{%s} + Py_{%s} [a.u]",p1.GetName(),p2.GetName()),200,-400,400, Form("%s/MomSums",Hname.Data()));
			if (normOfSumRot_ZX < pzxSumWidth)
				hi.fill(hiOff+4,"PySumCondXZ",(p1[i].Py()/MomScale + p2[j].Py()),Form("Py_{%s} + Py_{%s} [a.u]",p1.GetName(),p2.GetName()),200,-400,400, Form("%s/MomSums",Hname.Data()));
			//z momentum sum//
			hi.fill(hiOff+5,"PzSum",(p1[i].Pz()/MomScale + p2[j].Pz()),Form("Pz_{%s} + Pz_{%s} [a.u]",p1.GetName(),p2.GetName()),200,-400,400, Form("%s/MomSums",Hname.Data()));
			if (normOfSumRot_XY < pxySumWidth)
				hi.fill(hiOff+6,"PzSumCondXY",(p1[i].Pz()/MomScale + p2[j].Pz()),Form("Pz_{%s} + Pz_{%s} [a.u]",p1.GetName(),p2.GetName()),200,-400,400, Form("%s/MomSums",Hname.Data()));

			//PIPICO//
			hi.fill(hiOff+7,"PIPICO",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());
			//hi.fill(hiOff+10,"PIPICOMom",p1[i].Pz(),p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			//hi.fill(hiOff+11,"PIPICOMomFactor",p1[i].Pz()/MomScale,p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			if (normOfSumRot_XY < pxySumWidth)
			{
				hi.fill(hiOff+8,"PIPICOCondXY",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());
				hi.fill(hiOff+12,"RatioPCPerPICondXY",ratioP1P2,"Ratio",200,0,2,Form("%s/MomSums",Hname.Data()));
				hi.fill(hiOff+13,"RatioXYPCPerPICondXY",ratioP1P2XY,"Ratio",200,0,2,Form("%s/MomSums",Hname.Data()));
				//hi.fill(hiOff+12,"PIPICOMomCondXY",p1[i].Pz(),p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
				//hi.fill(hiOff+13,"PIPICOMomFactorCondXY",p1[i].Pz()/MomScale,p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			}

			//For confirmation. gated by only one plane
			if (normOfSumRot_XY < pxySumWidth)
			{
				hi.fill(hiOff+37,"PzPxSumRotCondXY",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
				hi.fill(hiOff+38,"PyPzSumRotCondXY",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
				
				hi.fill(hiOff+24,"AngleVsRatioXYCondXY",angleP1P2XY, ratioP1P2XY,"angle","Ratio",180,0,180,100,0,2, Form("%s/MomSums",Hname.Data()));				
				//hi.fill(hiOff+25,Form("%sXPosVsTofCondXY",p1.GetName()),p1[i].TofCor(),p1[i].XCorRotScl(),"tof [ns]","x [mm]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p1.GetXcor()-p1.GetCondRad()*1.3,p1.GetXcor()+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				//hi.fill(hiOff+26,Form("%sYPosVsTofCondXY",p1.GetName()),p1[i].TofCor(),p1[i].YCorRotScl(),"tof [ns]","y [mm]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p1.GetYcor()-p1.GetCondRad()*1.3,p1.GetYcor()+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				//hi.fill(hiOff+27,Form("%sXPosVsTofCondXY",p2.GetName()),p2[j].TofCor(),p2[j].XCorRotScl(),"tof [ns]","x [mm]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetXcor()-p2.GetCondRad()*1.3,p2.GetXcor()+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				//hi.fill(hiOff+28,Form("%sYPosVsTofCondXY",p2.GetName()),p2[j].TofCor(),p2[j].YCorRotScl(),"tof [ns]","y [mm]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetYcor()-p2.GetCondRad()*1.3,p2.GetYcor()+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				//hi.fill(hiOff+29,Form("%sDetCorCondXY",p1.GetName()),p1[i].XCor(),p1[i].YCor(),"x [mm]","y [mm]",300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
				//hi.fill(hiOff+30,Form("%sDetCorCondXY",p2.GetName()),p2[j].XCor(),p2[j].YCor(),"x [mm]","y [mm]",300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
			}
			//if (normOfSumRot_YZ < pyzSumWidth)
			//{
			//	hi.fill(hiOff+39,"PxPySumRotCondYZ",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			//	hi.fill(hiOff+40,"PzPxSumRotCondYZ",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			//}
			//if (normOfSumRot_ZX < pzxSumWidth)
			//{
			//	hi.fill(hiOff+41,"PyPzSumRotCondZX",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			//	hi.fill(hiOff+42,"PxPySumRotCondZX",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
			//}
			///////////////////apply two sum condition//////////////////////
			if (normOfSumRot_XY < pxySumWidth)
				if (normOfSumRot_YZ < pyzSumWidth)
				{
					hi.fill(hiOff+39,"PzPxSumRotCondXY-YZ",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
					hi.fill(hiOff+40,"NormOfSumRot_ZX_Gated",normOfSumRot_ZX,"norm [a.u.]",1000,0,0.001, Form("%s/MomSums",Hname.Data()),1/(2*TMath::Pi()*normOfSumRot_ZX));
				}
			if (normOfSumRot_YZ < pyzSumWidth)
				if (normOfSumRot_ZX < pzxSumWidth)
				{
					hi.fill(hiOff+41,"PxPySumRotCondYZ-ZX",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
					hi.fill(hiOff+42,"NormOfSumRot_XY_Gated",normOfSumRot_XY,"norm [a.u.]",1000,0,0.001, Form("%s/MomSums",Hname.Data()),1/(2*TMath::Pi()*normOfSumRot_XY));
				}
			if (normOfSumRot_ZX < pzxSumWidth)
				if (normOfSumRot_XY < pxySumWidth)
				{
					hi.fill(hiOff+43,"PyPzSumRotCondZX-XY",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
					hi.fill(hiOff+45,"NormOfSumRot_YZ_Gated",normOfSumRot_YZ,"norm [a.u.]",1000,0,0.001, Form("%s/MomSums",Hname.Data()),1/(2*TMath::Pi()*normOfSumRot_YZ));
				}

			///////////////////3-plane Momentum sum condition//////////////////////
			if (normOfSumRot_XY < pxySumWidth)
				if (normOfSumRot_YZ < pyzSumWidth)
					if (normOfSumRot_ZX < pzxSumWidth)
					{
						//Coincidence counter//
						mol.CoincidenceCount++;
						mol.CoinHitNbrC.push_back(i);
						mol.CoinHitNbrI.push_back(j);

						//Coincidence Tof//
						hi.fill(1,"TOF_Coincidence",p1[i].Tof(),"Tof [ns]",4000,0,4000);
						hi.fill(1,"TOF_Coincidence",p2[j].Tof(),"Tof [ns]",4000,0,4000);
						hi.fill(2,"XTOF_Coincidence",p1[i].Tof(),p1[i].X(),"tof [ns]","x [mm]",4000,0,4000,300,-50,50);
						hi.fill(2,"XTOF_Coincidence",p2[j].Tof(),p2[j].X(),"tof [ns]","x [mm]",4000,0,4000,300,-50,50);
						hi.fill(3,"YTOF_Coincidence",p1[i].Tof(),p1[i].Y(),"tof [ns]","y [mm]",4000,0,4000,300,-50,50);
						hi.fill(3,"YTOF_Coincidence",p2[j].Tof(),p2[j].Y(),"tof [ns]","y [mm]",4000,0,4000,300,-50,50);
						hi.fill(4,"PIPICO_Coincidence",p1[i].Tof(),p2[j].Tof(),"tof [ns]","tof [ns]",3000,1000,3000,3000,1000,3000);

						//PIPICO//
						hi.fill(hiOff+9,"PIPICOCondXYZ",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());
						hi.fill(hiOff+10,"PzSumCondXYZ",(p1[i].Pz()/MomScale + p2[j].Pz()),Form("Pz_{%s} + Pz_{%s} [a.u]",p1.GetName(),p2.GetName()),200,-400,400, Form("%s/MomSums",Hname.Data()));

						//MomSumRot//
						hi.fill(hiOff+46,"PxPySumRotCondXYZ",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
						hi.fill(hiOff+47,"PyPzSumRotCondXYZ",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
						hi.fill(hiOff+48,"PzPxSumRotCondXYZ",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
						//Intensity//
						if (intensity.size() && (intPart.size()>1)) 
							hi.fill(hiOff+49,Form("intensity%s%s",p1.GetName(),p2.GetName()),intensity[0], "[arb. unit]",intPart.size()-1,&intPart.front(),"PowerDependence");
						//Ratio PI/PC
						hi.fill(hiOff+15,"RatioPCPerPICondXYZ",ratioP1P2,"Ratio",200,0,2,Form("%s/MomSums",Hname.Data()));
						hi.fill(hiOff+16,"RatioXYPCPerPICondXYZ",ratioP1P2XY,"Ratio",200,0,2,Form("%s/MomSums",Hname.Data()));
						//Formed angle VS Ratio
						hi.fill(hiOff+22,"AngleVsRatioCondXYZ",angleP1P2, ratioP1P2,"angle","Ratio",180,0,180,100,0,2, Form("%s/MomSums",Hname.Data()),weightPerSin);
						hi.fill(hiOff+23,"AngleVsRatioXYCondXYZ",angleP1P2XY, ratioP1P2XY,"angle","Ratio",180,0,180,100,0,2, Form("%s/MomSums",Hname.Data()));

						//Momentum first Ion//
						int IDX = hiOff+50;
						hi.fill(IDX+0,Form("%sPxPy",p1.GetName()),p1[i].Px(),p1[i].Py(),"px [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						if (TMath::Abs(p1[i].Pz()) < 30)
						{
							hi.fill(IDX+1,Form("%sPxPySlice",p1.GetName()),p1[i].Px(),p1[i].Py(),"px [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						}

						hi.fill(IDX+3,Form("%sPxPz",p1.GetName()),p1[i].Pz(),p1[i].Px(),"pz [a.u.]","px [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						if (TMath::Abs(p1[i].Py()) < 30)
						{
							hi.fill(IDX+4,Form("%sPxPzSlice",p1.GetName()),p1[i].Pz(),p1[i].Px(),"pz [a.u.]","px [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						}

						hi.fill(IDX+6,Form("%sPyPz",p1.GetName()),p1[i].Pz(),p1[i].Py(),"pz [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						if (TMath::Abs(p1[i].Px()) < 30)
						{
							hi.fill(IDX+7,Form("%sPyPzSlice",p1.GetName()),p1[i].Pz(),p1[i].Py(),"pz [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						}
						//Momentum2D
						hi.fill(IDX+2,Form("%sMomentumXY",p1.GetName()),TMath::Sqrt(p1[i].Px()*p1[i].Px() + p1[i].Py()*p1[i].Py()),"pXY [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						hi.fill(IDX+5,Form("%sMomentumXZ",p1.GetName()),TMath::Sqrt(p1[i].Px()*p1[i].Px() + p1[i].Pz()*p1[i].Pz()),"pXY [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						hi.fill(IDX+8,Form("%sMomentumYZ",p1.GetName()),TMath::Sqrt(p1[i].Py()*p1[i].Py() + p1[i].Pz()*p1[i].Pz()),"pXY [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						hi.fill(IDX+21,Form("%sMomentumXY_k",p1.GetName()),TMath::Sqrt(p1[i].Px()*p1[i].Px() + p1[i].Py()*p1[i].Py())/MomScale,"pXY [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));

						hi.fill(IDX+9,Form("%sTotalMomentum",p1.GetName()),p1[i].P(),"p [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						hi.fill(IDX+22,Form("%sTotalMomentum_k",p1.GetName()),p1[i].P()/MomScale,"p [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
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

						//Angular Distribution
						//hi.fill(IDX+21,Form("%sThetaX",p1.GetName()),p1[i].ThetaX(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						//hi.fill(IDX+22,Form("%sThetaY",p1.GetName()),p1[i].ThetaY(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						//hi.fill(IDX+23,Form("%sThetaZ",p1.GetName()),p1[i].ThetaZ(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+24,Form("%sThetaYvsEnergy",p1.GetName()),p1[i].ThetaY(),p1[i].E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,200,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+25,Form("%sThetaXvsEnergy",p1.GetName()),p1[i].ThetaX(),p1[i].E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,200,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+26,Form("%sThetaZvsEnergy",p1.GetName()),p1[i].ThetaZ(),p1[i].E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,200,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+27,Form("%sPhiXYvsEnergy",p1.GetName()),p1[i].PhiXY(),p1[i].E(),"#phi [deg]","Energy [eV]",180,-180,180,200,0,200,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+28,Form("%sPhiYZvsEnergy",p1.GetName()),p1[i].PhiYZ(),p1[i].E(),"#phi [deg]","Energy [eV]",180,-180,180,200,0,200,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+29,Form("%sPhiZXvsEnergy",p1.GetName()),p1[i].PhiZX(),p1[i].E(),"#phi [deg]","Energy [eV]",180,-180,180,200,0,200,Form("%s/Angular",Hname.Data()));

						//Momentum second Ion//
						IDX = (p1!=p2)?IDX+30:IDX;
						hi.fill(IDX+0,Form("%sPxPy",p2.GetName()),p2[j].Px(),p2[j].Py(),"px [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						if (TMath::Abs(p2[j].Pz()) < 30)
						{
							hi.fill(IDX+1,Form("%sPxPySlice",p2.GetName()),p2[j].Px(),p2[j].Py(),"px [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
							//hi.fill(IDX+2,Form("%sPhiPxPySlice",p2.GetName()),TMath::ATan2(p2[j].Py(),p2[j].Px())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Py()*p2[j].Py() + p2[j].Px()*p2[j].Px()),"#phi [deg]","#sqrt{px^{2} + py^{2}} [a.u.]",360,-180,180,Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						}

						hi.fill(IDX+3,Form("%sPzPx",p2.GetName()),p2[j].Pz(),p2[j].Px(),"pz [a.u.]","px [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						if (TMath::Abs(p2[j].Py()) < 30)
						{
							hi.fill(IDX+4,Form("%sPxPzSlice",p2.GetName()),p2[j].Pz(),p2[j].Px(),"pz [a.u.]","px [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
							//hi.fill(IDX+5,Form("%sPhiPxPzSlice",p2.GetName()),TMath::ATan2(p2[j].Px(),p2[j].Pz())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Pz()*p2[j].Pz() + p2[j].Px()*p2[j].Px()),"#phi [deg]","#sqrt{px^{2} + pz^{2}} [a.u.]",360,-180,180,Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						}

						hi.fill(IDX+6,Form("%sPyPz",p2.GetName()),p2[j].Pz(),p2[j].Py(),"pz [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						if (TMath::Abs(p2[j].Px()) < 30)
						{
							hi.fill(IDX+7,Form("%sPyPzSlice",p2.GetName()),p2[j].Pz(),p2[j].Py(),"pz [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
							//hi.fill(IDX+8,Form("%sPhiPyPzSlice",p2.GetName()),TMath::ATan2(p2[j].Py(),p2[j].Pz())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Py()*p2[j].Py() + p2[j].Pz()*p2[j].Pz()),"#phi [deg]","#sqrt{pz^{2} + py^{2}} [a.u.]",360,-180,180,Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						}
						//Momentum2D
						hi.fill(IDX+2,Form("%sMomentumXY",p2.GetName()),TMath::Sqrt(p2[j].Px()*p2[j].Px() + p2[j].Py()*p2[j].Py()),"pXY [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						hi.fill(IDX+5,Form("%sMomentumXZ",p2.GetName()),TMath::Sqrt(p2[j].Px()*p2[j].Px() + p2[j].Pz()*p2[j].Pz()),"pXY [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						hi.fill(IDX+8,Form("%sMomentumYZ",p2.GetName()),TMath::Sqrt(p2[j].Py()*p2[j].Py() + p2[j].Pz()*p2[j].Pz()),"pXY [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));

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

						//Angular Distribution
						hi.fill(IDX+21,Form("%sThetaX",p2.GetName()),p2[j].ThetaX(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+22,Form("%sThetaY",p2.GetName()),p2[j].ThetaY(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+23,Form("%sThetaZ",p2.GetName()),p2[j].ThetaZ(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+24,Form("%sThetaYvsEnergy",p2.GetName()),p2[j].ThetaY(),p2[j].E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,100,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+25,Form("%sThetaXvsEnergy",p2.GetName()),p2[j].ThetaX(),p2[j].E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,100,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+26,Form("%sThetaZvsEnergy",p2.GetName()),p2[j].ThetaZ(),p2[j].E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,100,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+27,Form("%sPhiXYvsEnergy",p2.GetName()),p2[j].PhiXY(),p2[j].E(),"#phi [deg]","Energy [eV]",180,-180,180,200,0,100,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+28,Form("%sPhiYZvsEnergy",p2.GetName()),p2[j].PhiYZ(),p2[j].E(),"#phi [deg]","Energy [eV]",180,-180,180,200,0,100,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+29,Form("%sPhiZXvsEnergy",p2.GetName()),p2[j].PhiZX(),p2[j].E(),"#phi [deg]","Energy [eV]",180,-180,180,200,0,100,Form("%s/Angular",Hname.Data()));

						///////////////////KER//////////////////////////
						//hi.fill(IDX+30,Form("KER_I_%s",Hname.Data()),p2[j].E()*9.447078522,"KER [eV]",200,0,500,Form("%s/KER",Hname.Data()));
						//hi.fill(IDX+31,Form("KER_IperQ_%s",Hname.Data()),
						//	p2[j].E()*9.447078522/(p1.GetCharge_au()*(p2.GetCharge_au()+3)),
						//	"KER [eV]",200,0,20,Form("%s/KER",Hname.Data()));
						//hi.fill(IDX+32,Form("KEIplusKEC_%s",Hname.Data()),p1[i].E()+p2[j].E(),"KE [eV]",200,0,500,Form("%s/KER",Hname.Data()));
						//hi.fill(IDX+33,Form("KEIplusKECperQ_%s",Hname.Data()),
						//	(p1[i].E()+p2[j].E())/(p1.GetCharge_au()*p2.GetCharge_au()),
						//	"KE [eV]",200,0,20,Form("%s/KER",Hname.Data()));

						//momentum of H
						//TVector3 h3Pvec = p1[i].Pvec() + p2[j].Pvec();
						//hi.fill(IDX+30,"MomentumPC_PI",h3Pvec.Mag(),"p [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+31,"MomentumPC_PI-H",h3Pvec.Mag()/3/TMath::Cos(70*TMath::DegToRad()),"p [a.u.]",Mombins,0,MomLim/2,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+32,"PxPyPC_PI",h3Pvec.X(),h3Pvec.Y(),"px [a.u.]","py [a.u.]",Mombins,-MomLim/2,MomLim/2,Mombins,-MomLim/2,MomLim/2,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+33,"PyPzPC_PI",h3Pvec.Y(),h3Pvec.Z(),"py [a.u.]","pz [a.u.]",Mombins,-MomLim/2,MomLim/2,Mombins,-MomLim/2,MomLim/2,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+34,"PzPxPC_PI",h3Pvec.Z(),h3Pvec.X(),"pz [a.u.]","px [a.u.]",Mombins,-MomLim/2,MomLim/2,Mombins,-MomLim/2,MomLim/2,Form("%s/Momenta",Hname.Data()));
					}
		}
	}
	//std::cout << "\t" << hiOff;
}
//-----
void fillHydrogenHistogram(const MyParticle &p1, const MyParticle &p2, const MyParticle &p3, MyHistos &hi, int hiOff, Molecule &mol)
{
	double MomLim = 300;
	double MomSumRotLim = 800;
	double MomScale = 1;
	const int Mombins = 150;
	//const double rotAngle = 70 * TMath::DegToRad();

	TString Hname(p1.GetName());
	Hname += p2.GetName();
	Hname += p3.GetName();

	const double pxSumWidth = 40;
	const double pySumWidth = 40;
	const double pzSumWidth = 40;
	const double pxyzSumWidth = 40;

	hi.fill(hiOff+0,"NbrOfHitsH1p",p1.GetNbrOfParticleHits(),"Number of hits",20,0,20, Form("%s",Hname.Data()));
	for (size_t i=0; i<p1.GetNbrOfParticleHits();++i)
	{
		//if (p1[i].E() > 80) continue;
		//for (size_t j=0;j<p2.GetNbrOfParticleHits();++j)
		//for (size_t k=0;k<p3.GetNbrOfParticleHits();++k)
		for (size_t coin = 0; coin < mol.CoincidenceCount; ++coin) 
		{
			//if (j != mol.CoinHitNbrC[0]) continue;
			//if (k != mol.CoinHitNbrI[0]) continue;
			size_t j = mol.CoinHitNbrC[coin];
			size_t k = mol.CoinHitNbrI[coin];
			//---p1[i]="H" p2[j]="C" p3[k]="I"---//
			const TVector3 &pvecH = p1[i].Pvec();
			const TVector3 &pvecC = p2[j].Pvec();
			const TVector3 &pvecI = p3[k].Pvec();
			const TVector3 pvecIC = pvecI - pvecC;//molecule axis
			const TVector3 pvecSumIC = pvecI + pvecC;//3*H+

			
#ifdef TEST_RANDOM1
			//random test Sphere
			double rndPx, rndPy, rndPz;
			RandomGene.Sphere(rndPx, rndPy, rndPz, p1[i].P());
			const TVector3 pvecH(rndPx, rndPy, rndPz);
#endif // TEST_RANDOM

			//if (pvecIC.Z() > 0) continue;

			//hi.fill(hiOff+29,"Momentum_H",pvecH.Mag(),"p [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
			//hi.fill(hiOff+30,"MomentumPC_PI",pvecSumIC.Mag(),"p [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
			//hi.fill(hiOff+31,"MomentumPC_PI-H",pvecSumIC.Mag()/3/TMath::Cos(70*TMath::DegToRad()),"p [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
			//hi.fill(hiOff+32,"PxPyPC_PI",pvecSumIC.X(),pvecSumIC.Y(),"px [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
			//hi.fill(hiOff+33,"PyPzPC_PI",pvecSumIC.Y(),pvecSumIC.Z(),"py [a.u.]","pz [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
			//hi.fill(hiOff+34,"PzPxPC_PI",pvecSumIC.Z(),pvecSumIC.X(),"pz [a.u.]","px [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));

			const double angleHC = pvecH.Angle(-pvecC)*TMath::RadToDeg();
			const double angleHI = pvecH.Angle(pvecI)*TMath::RadToDeg();
			const double angleHCI = pvecH.Angle(pvecIC)*TMath::RadToDeg();
			double weightPerSin = 1/TMath::Sin(angleHCI*TMath::DegToRad());
			if (weightPerSin > 100) weightPerSin = 100;
			//if ((angleHCI < 90) || (angleHCI > 130)) continue;

			//const double ratioHC = p1[i].P()/p2[j].P();
			const double ratioHCI = pvecH.Mag()/(pvecSumIC.Mag());
			const double ratioHCIcos = pvecH.Mag()*TMath::Cos((180-angleHCI)*TMath::DegToRad())*3/(pvecSumIC.Mag());
			//const TVector3 crossHC = pvecC.Cross(pvecH);
			const TVector3 crossHI = pvecI.Cross(pvecH);
			//const TVector3 crossHCI = pvecIC.Cross(pvecH);
			const double planeHCI = pvecC * crossHI/(pvecH.Mag()*pvecC.Mag()*pvecI.Mag());
			//test condition
			//if ((planeHCI < -0.05) || (planeHCI > 0.05)) continue;

			TVector3 pvecHRot1(pvecH);
			TVector3 pvecHRot2(pvecH);
			TVector3 pvecHRot3(pvecH);
			pvecHRot2.Rotate(120*TMath::DegToRad(), pvecIC);//pvecIC
			pvecHRot3.Rotate(-120*TMath::DegToRad(), pvecIC);//pvecIC
			TVector3 pvec3Hp(pvecHRot1+pvecHRot2+pvecHRot3);
			const double ratioH3CI = pvec3Hp.Mag()/pvecSumIC.Mag();
			//Total momentum sum Vector
			TVector3 pvecMomSum(pvec3Hp+pvecI+pvecC);

			//TVector3 pvecHRot(pvecH);
			//pvecHRot.Rotate((180-angleHCI)*TMath::DegToRad() , crossHCI);

			//TVector3 pvecHRot(-pvecIC);
			//pvecHRot = pvecHRot.Unit() * 3 * TMath::Cos((180-angleHCI)*TMath::DegToRad()) * pvecH.Mag();

			hi.fill(hiOff+17,"AngleVsMomSumH3CI",angleHCI, pvecMomSum.Mag(),"angle [deg]","MomSum",90,0,180,100,0,200, Form("%s/Angular",Hname.Data()),weightPerSin);
			hi.fill(hiOff+19,"AngleVsRatioH3_sumCI",angleHCI, ratioH3CI,"angle [deg]","Ratio",90,0,180,50,0,4, Form("%s/Angular",Hname.Data()),weightPerSin);
			hi.fill(hiOff+20,"AngleVsRatioHCI",angleHCI, ratioHCI,"angle [deg]","Ratio",90,0,180,50,0,4, Form("%s/Angular",Hname.Data()),weightPerSin);
			hi.fill(hiOff+21,"PlaneHCI",planeHCI,"cos",100,-1,1, Form("%s/Angular",Hname.Data()));
			//hi.fill(hiOff+22,"AngleVsRatioRotHCI",pvecHRot.Angle(pvecIC)*TMath::RadToDeg(), ratioHCI,"angle","Ratio",90,0,180,50,0,4, Form("%s/Angular",Hname.Data()));
			hi.fill(hiOff+23,"AngleVsRatioRotHI",pvec3Hp.Angle(pvecI)*TMath::RadToDeg(), ratioHCI,"angle [deg]","Ratio",90,0,180,50,0,4, Form("%s/Angular",Hname.Data()));
			hi.fill(hiOff+24,"AngleHCI",angleHCI,"angle [deg]",60,0,180, Form("%s/Angular",Hname.Data()),weightPerSin);
			hi.fill(hiOff+33,"AngleHI",angleHI,"angle [deg]",60,0,180, Form("%s/Angular",Hname.Data()),weightPerSin);
			hi.fill(hiOff+34,"AngleHC",angleHC,"angle [deg]",60,0,180, Form("%s/Angular",Hname.Data()),weightPerSin);
			//hi.fill(hiOff+18,"AngleH-CIvsMomH",angleHCI, p1[i].P(),"angle [deg]","H Momentum [a.u.]",90,0,180,50,0,MomLim, Form("%s/Angular",Hname.Data()),weightPerSin);
			//hi.fill(hiOff+35,"AngleH-CvsMomH",angleHC, p1[i].P(),"angle [deg]","H Momentum [a.u.]",90,0,180,50,0,MomLim, Form("%s/Angular",Hname.Data()),weightPerSin);
			//hi.fill(hiOff+36,"AngleH-IvsMomH",angleHI, p1[i].P(),"angle [deg]","H Momentum [a.u.]",90,0,180,50,0,MomLim, Form("%s/Angular",Hname.Data()),weightPerSin);
			hi.fill(hiOff+18,"AngleH-CIvsMomH",angleHCI, pvecH.Mag(),"angle [deg]","H Momentum [a.u.]",90,0,180,50,0,MomLim, Form("%s/Angular",Hname.Data()),weightPerSin);
			hi.fill(hiOff+35,"AngleH-CvsMomH",angleHC, pvecH.Mag(),"angle [deg]","H Momentum [a.u.]",90,0,180,50,0,MomLim, Form("%s/Angular",Hname.Data()),weightPerSin);
			hi.fill(hiOff+36,"AngleH-IvsMomH",angleHI, pvecH.Mag(),"angle [deg]","H Momentum [a.u.]",90,0,180,50,0,MomLim, Form("%s/Angular",Hname.Data()),weightPerSin);


			//1st condition
			//if ((angleHC < 20) || (angleHC > 120)) continue;
			//if ((ratioHC < 0.1) || (ratioHC > 0.4)) continue;
			//hi.fill(hiOff+21,"AngleVsRatioCond",angleP1P2, ratioP1P2,"angle","Ratio",180,0,180,100,0,2, Form("%s/MomSums",Hname.Data()));
			
			double momSumX = pvecMomSum.X();
			double momSumY = pvecMomSum.Y();
			double momSumZ = pvecMomSum.Z();

			//-----Gated stuff-----//
			// momentum sum//
			hi.fill(hiOff+1,"PxSum",(momSumX),"PxSum [a.u]",1000,-400,400,Form("%s/MomSums",Hname.Data()));
			if ((momSumY < pySumWidth)&&(momSumZ < pzSumWidth))
				hi.fill(hiOff+2,"PxSumCondYZ", momSumX,"PxSum [a.u]",1000,-400,400, Form("%s/MomSums",Hname.Data()));
			//y momentum sum//
			hi.fill(hiOff+3,"PySum",(momSumY),"PxSum [a.u]",1000,-400,400, Form("%s/MomSums",Hname.Data()));
			if ((momSumZ < pzSumWidth)&&(momSumX < pxSumWidth))
				hi.fill(hiOff+4,"PySumCondZX", momSumY,"PxSum [a.u]",1000,-400,400, Form("%s/MomSums",Hname.Data()));
			//z momentum sum//
			hi.fill(hiOff+5,"PzSum",(momSumZ),"PxSum [a.u]",1000,-400,400, Form("%s/MomSums",Hname.Data()));
			if ((momSumX < pxSumWidth)&&(momSumY < pySumWidth))
				hi.fill(hiOff+6,"PzSumCondXY", momSumZ,"PxSum [a.u]",1000,-400,400, Form("%s/MomSums",Hname.Data()));

			//PIPICO//
			//hi.fill(hiOff+7,"PIPICO",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());
			//hi.fill(hiOff+10,"PIPICOMom",p1[i].Pz(),p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			//hi.fill(hiOff+11,"PIPICOMomFactor",p1[i].Pz()/MomScale,p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			//if ((momSumX < pxSumWidth)&&(momSumY < pySumWidth))
			//{
			//	hi.fill(hiOff+8,"PIPICOCondXY",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());
			//	hi.fill(hiOff+12,"PIPICOMomCondXY",p1[i].Pz(),p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			//	hi.fill(hiOff+13,"PIPICOMomFactorCondXY",p1[i].Pz()/MomScale,p2[j].Pz(),Form("pz_{%s} [a.u]",p1.GetName()),Form("pz_{%s} [a.u]",p2.GetName()),200,-1000,1000,200,-1000,1000,Hname.Data());
			//}

			///////////////////Momentum sum condition//////////////////////
			//if ((planeHCI > -0.1)&&(planeHCI < 0.1))
			if (pvecMomSum.Mag() < (pxyzSumWidth))
			//if (momSumX < pxSumWidth)
			//	if (momSumY < pySumWidth)
			//		if (momSumZ < pzSumWidth)
					{
						//PIPICO//
						//hi.fill(hiOff+9,"PIPICOCondXYZ",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());

						////MomSumRot//
						//hi.fill(hiOff+43,"PxPySumRotCondXYZ",p12xSumRot_XY, p12ySumRot_XY,"pxySum [a.u.]","pxySum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
						//hi.fill(hiOff+44,"PyPzSumRotCondXYZ",p12ySumRot_YZ, p12zSumRot_YZ,"pyzSum [a.u.]","pyzSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
						//hi.fill(hiOff+45,"PzPxSumRotCondXYZ",p12zSumRot_ZX, p12xSumRot_ZX,"pzxSum [a.u.]","pzxSum [a.u.]",200,-MomSumRotLim,MomSumRotLim,200,-MomSumRotLim,MomSumRotLim, Form("%s/MomSums",Hname.Data()));
						//Intensity//
						//Ratio PI/PC
						//hi.fill(hiOff+15,"RatioPHPerPCPICondXYZ",p1[i].P()/(p3[k].P()-p2[j].P()),"Ratio",200,0,2,Form("%s/MomSums",Hname.Data()));
						//Formed angle VS Ratio
						hi.fill(hiOff+25,"AngleVsRatioH3-CICondXYZ",pvec3Hp.Angle(pvecIC)*TMath::RadToDeg(), ratioHCI,"angle [deg]","Ratio",90,0,180,50,0,4, Form("%s/Angular",Hname.Data()));
						hi.fill(hiOff+26,"AngleVsRatioHCICondXYZ",angleHCI, ratioHCI,"angle [deg]","Ratio",90,0,180,50,0,4, Form("%s/Angular",Hname.Data()),weightPerSin);
						//hi.fill(hiOff+27,"AngleVsRatioHCIcosCondXYZ",angleHCI, ratioHCIcos,"angle [deg]","Ratio",90,0,180,50,-4,4, Form("%s/Angular",Hname.Data()),weightPerSin);
						hi.fill(hiOff+14,"AngleHCIVsMomHCondXYZ",angleHCI, pvecH.Mag(),"angle [deg]","H Momentum [a.u.]",90,0,180,50,0,MomLim, Form("%s/Angular",Hname.Data()),weightPerSin);
						hi.fill(hiOff+37,"AngleHCVsMomHCondXYZ",angleHC, pvecH.Mag(),"angle [deg]","H Momentum [a.u.]",90,0,180,50,0,MomLim, Form("%s/Angular",Hname.Data()),weightPerSin);
						hi.fill(hiOff+38,"AngleHIVsMomHCondXYZ",angleHI, pvecH.Mag(),"angle [deg]","H Momentum [a.u.]",90,0,180,50,0,MomLim, Form("%s/Angular",Hname.Data()),weightPerSin);
						hi.fill(hiOff+27,"AngleVsMomSumH3CICondXYZ",angleHCI, pvecMomSum.Mag(),"angle [deg]","MomSum",90,0,180,100,0,200, Form("%s/Angular",Hname.Data()),weightPerSin);
						
						hi.fill(hiOff+10,"AngleHCICondXYZ",angleHCI,"angle [deg]",60,0,180, Form("%s/Angular",Hname.Data()),weightPerSin);
						hi.fill(hiOff+11,"AngleHCCondXYZ",angleHC,"angle [deg]",60,0,180, Form("%s/Angular",Hname.Data()),weightPerSin);
						hi.fill(hiOff+12,"AngleHICondXYZ",angleHI,"angle [deg]",60,0,180, Form("%s/Angular",Hname.Data()),weightPerSin);

						hi.fill(hiOff+30,"PxSumCondXYZ",momSumX,"PxSum [a.u]",1000,-400,400, Form("%s/MomSums",Hname.Data()));
						hi.fill(hiOff+31,"PySumCondXYZ",momSumY,"PxSum [a.u]",1000,-400,400, Form("%s/MomSums",Hname.Data()));
						hi.fill(hiOff+32,"PzSumCondXYZ",momSumZ,"PxSum [a.u]",1000,-400,400, Form("%s/MomSums",Hname.Data()));

						
						//Momentum Proton//
						int IDX = hiOff+40;
						hi.fill(IDX+0,Form("%sPxPy",p1.GetName()),p1[i].Px(),p1[i].Py(),"px [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						hi.fill(IDX+2,Form("%s-%s_PxPy",p2.GetName(),p3.GetName()),p2[j].Px()+p3[k].Px(),p2[j].Py()+p3[k].Py(),"px [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						if (TMath::Abs(p1[i].Pz()) < 30)
						{
							hi.fill(IDX+1,Form("%sPxPySlice",p1.GetName()),p1[i].Px(),p1[i].Py(),"px [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						}

						hi.fill(IDX+3,Form("%sPzPx",p1.GetName()),p1[i].Pz(),p1[i].Px(),"pz [a.u.]","px [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						hi.fill(IDX+5,Form("%s-%s_PzPx",p2.GetName(),p3.GetName()),p2[j].Pz()+p3[k].Pz(),p2[j].Px()+p3[k].Px(),"pz [a.u.]","px [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						if (TMath::Abs(p1[i].Py()) < 30)
						{
							hi.fill(IDX+4,Form("%sPzPxSlice",p1.GetName()),p1[i].Pz(),p1[i].Px(),"pz [a.u.]","px [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						}

						hi.fill(IDX+6,Form("%sPyPz",p1.GetName()),p1[i].Pz(),p1[i].Py(),"py [a.u.]","pz [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						hi.fill(IDX+8,Form("%s-%s_PyPz",p2.GetName(),p3.GetName()),p2[j].Py()+p3[k].Py(),p2[j].Pz()+p3[k].Pz(),"py [a.u.]","pz [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						if (TMath::Abs(p1[i].Px()) < 30)
						{
							hi.fill(IDX+7,Form("%sPyPzSlice",p1.GetName()),p1[i].Pz(),p1[i].Py(),"pz [a.u.]","py [a.u.]",Mombins,-MomLim,MomLim,Mombins,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						}

						hi.fill(IDX+9,Form("%sTotalMomentum",p1.GetName()),p1[i].P(),"Momentum [a.u.]",Mombins,0,MomLim,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+2,Form("%sPx",p1.GetName()),p1[i].Px(),"px [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+5,Form("%sPy",p1.GetName()),p1[i].Py(),"py [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						//hi.fill(IDX+8,Form("%sPz",p1.GetName()),p1[i].Pz(),"pz [a.u.]",300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
						
						//Energy Ion//
						hi.fill(IDX+10,Form("%sEnergy",p1.GetName()),p1[i].E(),"Energy [eV]",200,0,400,Form("%s/Energy",Hname.Data()));
						hi.fill(IDX+11,Form("%sEnergy",p2.GetName()),p2[j].E(),"Energy [eV]",200,0,400,Form("%s/Energy",Hname.Data()));
						hi.fill(IDX+12,Form("%sEnergy",p3.GetName()),p3[k].E(),"Energy [eV]",200,0,400,Form("%s/Energy",Hname.Data()));
						//Raw
						//hi.fill(IDX+11,Form("%sTOF",p1.GetName()),p1[i].TofCor(),"Tof [ns]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,Form("%s/Raw",Hname.Data()));
						//hi.fill(IDX+12,Form("%sDetCor",p1.GetName()),p1[i].XCor(),p1[i].YCor(),"x [mm]","y [mm]",300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
						hi.fill(IDX+13,Form("%sXPosVsTof",p1.GetName()),p1[i].TofCor(),p1[i].XCorRotScl(),"tof [ns]","x [mm]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p1.GetXcor()-p1.GetCondRad()*1.3,p1.GetXcor()+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
						hi.fill(IDX+14,Form("%sYPosVsTof",p1.GetName()),p1[i].TofCor(),p1[i].YCorRotScl(),"tof [ns]","y [mm]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p1.GetYcor()-p1.GetCondRad()*1.3,p1.GetYcor()+p1.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
						
						//Angular Distribution
						hi.fill(IDX+15,Form("%sThetaX",p1.GetName()),p1[i].ThetaX(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+16,Form("%sThetaY",p1.GetName()),p1[i].ThetaY(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+17,Form("%sThetaZ",p1.GetName()),p1[i].ThetaZ(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+18,Form("%sThetaYvsEnergy",p1.GetName()),p1[i].ThetaY(),p1[i].E(),"#theta [deg]","Energy [eV]",90,0,180,100,0,400,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+19,Form("%sThetaXvsEnergy",p1.GetName()),p1[i].ThetaX(),p1[i].E(),"#theta [deg]","Energy [eV]",90,0,180,100,0,400,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+20,Form("%sThetaZvsEnergy",p1.GetName()),p1[i].ThetaZ(),p1[i].E(),"#theta [deg]","Energy [eV]",90,0,180,100,0,400,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+21,Form("%sPhiXYvsEnergy",p1.GetName()),p1[i].PhiXY(),p1[i].E(),"#phi [deg]","Energy [eV]",90,-180,180,100,0,400,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+22,Form("%sPhiYZvsEnergy",p1.GetName()),p1[i].PhiYZ(),p1[i].E(),"#phi [deg]","Energy [eV]",90,-180,180,100,0,400,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+23,Form("%sPhiZXvsEnergy",p1.GetName()),p1[i].PhiZX(),p1[i].E(),"#phi [deg]","Energy [eV]",90,-180,180,100,0,400,Form("%s/Angular",Hname.Data()));
						hi.fill(IDX+24,Form("%sThetaZvs%sEnergy",p3.GetName(), p1.GetName()),p3[k].ThetaZ(),p1[i].E(),"#theta [deg]","Energy [eV]",90,0,180,100,0,400,Form("%s/Angular",Hname.Data()));

						///////////////////KER//////////////////////////
						//hi.fill(IDX+29,Form("KER",Hname.Data()),p1[i].E()*3 + p2[j].E() + p3[k].E(),"KER [eV]",200,0,2000,Form("%s/KER",Hname.Data()));
						//hi.fill(IDX+31,Form("KER_IperQ_%s",Hname.Data()),
						//	p2[j].E()*9.447078522/(p1.GetCharge_au()*(p2.GetCharge_au()+3)),
						//	"KER [eV]",200,0,20,Form("%s/KER",Hname.Data()));
						//hi.fill(IDX+32,Form("KEIplusKEC_%s",Hname.Data()),p1[i].E()+p2[j].E(),"KE [eV]",200,0,500,Form("%s/KER",Hname.Data()));
						//hi.fill(IDX+33,Form("KEIplusKECperQ_%s",Hname.Data()),
						//	(p1[i].E()+p2[j].E())/(p1.GetCharge_au()*p2.GetCharge_au()),
						//	"KE [eV]",200,0,20,Form("%s/KER",Hname.Data()));
					
					}
						
		}
	}

}

//----For gated particles
void fillMoleculeHistogram2(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff)
{
	TString Hname(p2.GetName());
	Hname += "GatedBy";
	Hname += p1.GetName();

	const double MomLim = 800;

	for (size_t i=0;i<p2.GetNbrOfParticleHits();++i)
	{
		if (!(p1[0].E() > p1.GetEnergyFrom() && p1[0].E() < p1.GetEnergyTo())) return;
		//limit detection angle
		//if (!( (p1[0].PhiZX() > 40 && p1[0].PhiZX() < 140) || (p1[0].PhiZX() > -140 && p1[0].PhiZX() < -40) )) return;
		//if (!( (p2[i].PhiZX() > 40 && p2[i].PhiZX() < 140) || (p2[i].PhiZX() > -140 && p2[i].PhiZX() < -40) )) continue;


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
		hi.fill(IDX+10,Form("%sEnergy",p2.GetName()),p2[i].E(),"Energy [eV]",300,0,200,Form("%s/Energy",Hname.Data()));

		//Raw
		hi.fill(IDX+11,Form("%sTOF",p2.GetName()),p2[i].TofCor(),"Tof [ns]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Form("%s/Raw",Hname.Data()));
		hi.fill(IDX+12,Form("%sDetCor",p2.GetName()),p2[i].XCor(),p2[i].YCor(),"x [mm]","y [mm]",300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
		hi.fill(IDX+13,Form("%sXPosVsTof",p2.GetName()),p2[i].TofCor(),p2[i].XCorRotScl(),"tof [ns]","x [mm]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetXcor()-p2.GetCondRad()*1.3,p2.GetXcor()+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
		hi.fill(IDX+14,Form("%sYPosVsTof",p2.GetName()),p2[i].TofCor(),p2[i].YCorRotScl(),"tof [ns]","y [mm]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetYcor()-p2.GetCondRad()*1.3,p2.GetYcor()+p2.GetCondRad()*1.3,Form("%s/Raw",Hname.Data()));
		
		//Angular Distribution
		hi.fill(IDX+15,Form("%sThetaX",p2.GetName()),p2[i].ThetaX(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
		hi.fill(IDX+16,Form("%sThetaY",p2.GetName()),p2[i].ThetaY(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
		hi.fill(IDX+17,Form("%sThetaZ",p2.GetName()),p2[i].ThetaZ(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
		hi.fill(IDX+18,Form("%sThetaYvsEnergy",p2.GetName()),p2[i].ThetaY(),p2[i].E(),"#theta [deg]","Energy [eV]",90,0,180,100,0,400,Form("%s/Angular",Hname.Data()));
		hi.fill(IDX+19,Form("%sThetaXvsEnergy",p2.GetName()),p2[i].ThetaX(),p2[i].E(),"#theta [deg]","Energy [eV]",90,0,180,100,0,400,Form("%s/Angular",Hname.Data()));
		hi.fill(IDX+20,Form("%sThetaZvsEnergy",p2.GetName()),p2[i].ThetaZ(),p2[i].E(),"#theta [deg]","Energy [eV]",90,0,180,100,0,400,Form("%s/Angular",Hname.Data()));
		hi.fill(IDX+21,Form("%sPhiXYvsEnergy",p2.GetName()),p2[i].PhiXY(),p2[i].E(),"#phi [deg]","Energy [eV]",90,-180,180,100,0,400,Form("%s/Angular",Hname.Data()));
		hi.fill(IDX+22,Form("%sPhiYZvsEnergy",p2.GetName()),p2[i].PhiYZ(),p2[i].E(),"#phi [deg]","Energy [eV]",90,-180,180,100,0,400,Form("%s/Angular",Hname.Data()));
		hi.fill(IDX+23,Form("%sPhiZXvsEnergy",p2.GetName()),p2[i].PhiZX(),p2[i].E(),"#phi [deg]","Energy [eV]",90,-180,180,100,0,400,Form("%s/Angular",Hname.Data()));

		//anglar
		//for (int j = 0; j < p1.GetNbrOfParticleHits(); j++)
		//{
		//	if ((i>=j) && (p1==p2)) continue;
		//	const TVector3 &PIodine = p1[j].Pvec();
		//	const TVector3 &PTarget = p2[i].Pvec();
		//	double formedAngle =  PIodine.Angle(PTarget) * TMath::RadToDeg();			
		//	double weightPerSin = 1/TMath::Sin(formedAngle*TMath::DegToRad());
		//	if (weightPerSin > 100) weightPerSin = 100;

		//	hi.fill(IDX+15,Form("Angle%s-%s",p1.GetName(),p2.GetName()),formedAngle,"Formed angle [deg]",180,0,180,Form("%s/Angle",Hname.Data()),weightPerSin);
		//}
	}
}
void fillAngleHistogram(const MyParticle &p1, const MyParticle &p2, MyHistos &hi, int hiOff)
{
	//TString Hname(p2.GetName());
	//Hname += "GatedBy";
	//Hname += p1.GetName();

	for (size_t i = 0; i < p1.GetNbrOfParticleHits(); ++i)
	{
		for (int j = 0; j < p2.GetNbrOfParticleHits(); ++j)
		{
			if ((i>=j) && (p1==p2)) continue;
			if (!(p1[i].E() > p1.GetEnergyFrom() && p1[i].E() < p1.GetEnergyTo())) return;//energy in range
			if (!(p2[j].E() > p2.GetEnergyFrom() && p2[j].E() < p2.GetEnergyTo())) return;
			//limit detection angle
			//if (!((p1[i].PhiZX() > 40 && p1[i].PhiZX() < 140) || (p1[i].PhiZX() > -140 && p1[i].PhiZX() < -40))) continue;
			//if (!((p2[j].PhiZX() > 40 && p2[j].PhiZX() < 140) || (p2[j].PhiZX() > -140 && p2[j].PhiZX() < -40))) continue;

			const TVector3 &PIodine = p1[i].Pvec();
			const TVector3 &PTarget = p2[j].Pvec();
			double formedAngle =  PIodine.Angle(PTarget) * TMath::RadToDeg();			
			double weightPerSin = 1/TMath::Sin(formedAngle*TMath::DegToRad());
			if (weightPerSin > 100) weightPerSin = 100;

			hi.fill(hiOff+0,Form("Angle%s-%s",p1.GetName(),p2.GetName()),formedAngle,"Formed angle [deg]",180,0,180,"AngularDist",weightPerSin);
			hi.fill(hiOff+1,Form("Angle%s-%sVsMomRatio",p1.GetName(),p2.GetName()),formedAngle, p2[j].P()/p1[i].P(),"Formed angle [deg]","Energy [eV]",180,0,180,100,0,3,"AngleVsMomRatio",weightPerSin);
			hi.fill(hiOff+2,Form("Energy%sVsEnergy%s",p2.GetName(),p1.GetName()),p1[i].E(),p2[j].E(),Form("Energy%s [eV]",p1.GetName()),Form("Energy%s [eV]",p2.GetName()),100,0,300,100,0,300,"KEcorrelation");
		}
	}
}
//________________3-body angular distribution___________________________
void fill3BodyHistogram(const MyParticle &p1, const MyParticle &p2,const MyParticle &p3, MyHistos &hi, int hiOff)
{
	for (size_t i = 0; i < p1.GetNbrOfParticleHits(); ++i)
	{
		for (int j = 0; j < p2.GetNbrOfParticleHits(); ++j)
		{
			for (int k = 0; k < p3.GetNbrOfParticleHits(); ++k)
			{
			if ((p1==p2) && (i==j)) continue;
			if ((p2==p3) && (j==k)) continue;
			if ((p3==p1) && (k==i)) continue;
			//if (!((p1[i].PhiZX() > 40 && p1[i].PhiZX() < 140) || (p1[i].PhiZX() > -140 && p1[i].PhiZX() < -40))) continue;
			//if (!((p2[j].PhiZX() > 40 && p2[j].PhiZX() < 140) || (p2[j].PhiZX() > -140 && p2[j].PhiZX() < -40))) continue;
			//if (!((p3[k].PhiZX() > 40 && p3[k].PhiZX() < 140) || (p3[k].PhiZX() > -140 && p3[k].PhiZX() < -40))) continue;


			const TVector3 &pvec1 = p1[i].Pvec();
			const TVector3 &pvec2 = p2[j].Pvec();
#ifndef TEST_RANDOM2
			const TVector3 &pvec3 = p3[k].Pvec();
#endif
#ifdef TEST_RANDOM2
			//random test Sphere
			double rndPx, rndPy, rndPz;
			//RandomGene.Sphere(rndPx, rndPy, rndPz, p1[i].P());
			//const TVector3 pvec1(rndPx, rndPy, rndPz);
			//RandomGene.Sphere(rndPx, rndPy, rndPz, p2[j].P());
			//const TVector3 pvec2(rndPx, rndPy, rndPz);
			RandomGene.Sphere(rndPx, rndPy, rndPz, p3[k].P());
			const TVector3 pvec3(rndPx, rndPy, rndPz);
#endif // TEST_RANDOM
			//cos(phi) = (P1xP2).P3/|P1||P2||P3|
			//calc A
			//double angle3Body =  pvec1.Cross(pvec2).Dot(pvec3)/(pvec1.Mag()*pvec2.Mag()*pvec3.Mag());
			//calc B
			//double angle3Body =  pvec1.Cross(pvec2).Dot(pvec3)/(pvec1.Cross(pvec2).Mag() *pvec3.Mag());
			//calc C
			double angleTheta = pvec1.Angle(pvec2);
			double angle3Body =  pvec1.Cross(pvec2).Dot(pvec3)/(pvec1.Mag()*pvec2.Mag()*pvec3.Mag()*TMath::Sin(angleTheta));

			double formedAngle = acos(angle3Body) * TMath::RadToDeg();
			//double weightPerSin = 1/TMath::Sin(formedAngle*TMath::DegToRad());
			//if (weightPerSin > 100) weightPerSin = 100;

			hi.fill(hiOff+0,Form("3BodyAngleCos-%s-%s-%s",p1.GetName(),p2.GetName(),p3.GetName()),angle3Body,"cos(Phi)",100,-1,1,"3BodyAngularCorr");
			//hi.fill(hiOff+1,Form("3BodyAngle-%s-%s-%s",p1.GetName(),p2.GetName(),p3.GetName()),formedAngle,"angle [deg]",90,0,180,"3BodyAngularCorr",weightPerSin);
			//hi.fill(hiOff+1,Form("3BodyAngle%s-%sVsEnegy%s",p1.GetName(),p2.GetName(),p1.GetName()),formedAngle,p1[i].E(),"Formed angle [deg]","Energy [eV]",180,0,180,100,0,100,"AngularEnergyDist",weightPerSin);
			//hi.fill(hiOff+2,Form("Angle%s-%sVsEnegy%s",p1.GetName(),p2.GetName(),p2.GetName()),formedAngle,p2[j].E(),"Formed angle [deg]","Energy [eV]",180,0,180,100,0,300,"AngularEnergyDist",weightPerSin);
		
			}
		}
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
			hi.fill(100,"PIPICO",p[i].TofCor(),p[j].TofCor(),"tof","tof",2000,p.GetCondTofFr(),p.GetCondTofTo(),2000,p.GetCondTofFr(),p.GetCondTofTo());
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

void fillMCPToFHistograms(const MyOriginalEvent &oe, MyHistos &hi, std::vector<MCPToFRegion> &mcpTofRegion)
{

}
void fillHistosAfterAnalyzis(const std::vector<MyParticle> &particles, MyHistos &hi,size_t nRegion,int delayBins,double delayFrom,double delayTo )
{
	std::cout << std::endl << "Running post analysis." << std::endl;
	size_t idx = 9000;

	TH1D* delayVsShots = dynamic_cast<TH1D*>(gFile->FindObject("DelayVsShots"));

	TH2D* delayVsTOF = dynamic_cast<TH2D*>( gFile->GetDirectory("Ion")->FindObject("DelayVsTOF") );
	DivideHisto2Dby1D(delayVsTOF,delayVsShots);
	//TH2D* delayVsTOFCor = dynamic_cast<TH2D*>( gFile->GetDirectory("Ion")->FindObject("DelayVsTOFCor") );
	//DivideHisto2Dby1D(delayVsTOFCor,delayVsShots);
	//TH2D* FELintensityVsTOFCor = dynamic_cast<TH2D*>(gFile->GetDirectory("Ion")->FindObject("FELIntensityVsTOFCor"));
	//DivideHisto2Dby1D(FELintensityVsTOFCor, FELIntensityVsShots);
	//TH2D* delayVsMCPSignal = dynamic_cast<TH2D*>( gFile->FindObject("DelayVsMCPSignal") );
	//DivideHisto2Dby1D(delayVsMCPSignal, delayVsShots);

	//Normalize histogram in particle directory
	TDirectory * rootDir = gDirectory;
	for(size_t i=1; i< particles.size(); i++)
	{
		//create a histogram that will take the normalized Data//
		TH1D* delayVsCountNorm = dynamic_cast<TH1D*>(hi.create1d( (idx+10*i) + 1,"DelayVsCount_Norm","delay [fs]", delayBins,delayFrom,delayTo,Form("%s/Delay",particles[i].GetName())));
		//TH2D* DealyVsEnergyNorm = dynamic_cast<TH2D*>(hi.create2d( (idx+10*i) + 2,"DealyVsEnergy_Norm","Number of Hits", delayBins,delayFrom,delayTo,Form("%s",particles[i].GetName())));
		
		//get histogram
		TH1D* delayVsCount = dynamic_cast<TH1D*>( gFile->GetDirectory(Form("%s/Delay",particles[i].GetName()))->FindObject("Delay") );

		//TH1D* nbrParticleHits = dynamic_cast<TH1D*>(gFile->FindObject("NumberOfParticleHits"));
		if (delayVsCount) 
		{
			delayVsCountNorm->Divide(delayVsCount,delayVsShots);
		}
		
		TH2D* delayVsMomentum = dynamic_cast<TH2D*>(gFile->GetDirectory(Form("%s/Delay", particles[i].GetName()))->FindObject("DelayVsMomentum"));
		if (delayVsMomentum)
		{
			DivideHisto2Dby1D(delayVsMomentum, delayVsShots);
		}
		TH2D* delayVsEnergy = dynamic_cast<TH2D*>(gFile->GetDirectory(Form("%s/Delay", particles[i].GetName()))->FindObject("DelayVsEnergy"));
		if (delayVsEnergy)
		{
			DivideHisto2Dby1D(delayVsEnergy, delayVsShots);
		}

		//TH2D* delayVsEnergy = dynamic_cast<TH2D*>( gFile->GetDirectory(Form("%s/Delay",particles[i].GetName()))->FindObject("DelayVsEnergy") );
		//TH2D* DealyVsEnergyNorm = dynamic_cast<TH2D*>(delayVsEnergy ->Clone("DealyVsEnergy"));
		//DealyVsEnergyNorm -> SetUniqueID((idx+20*i));
		//DealyVsEnergyNorm -> SetTitle("DealyVsEnergy_Norm");
		//DivideHisto2Dby1D(DealyVsEnergyNorm,delayVsShots);
	}


	////-------------------------------------------------------------------------------------------------------------------------------
	//TH1D* nbrPartHitsNorm = dynamic_cast<TH1D*>(hi.create1d(idx,"NbrParticleHits_normalized","Number of Hits",particles.size()+1,0,particles.size()+1));
	//TH1D* nbrPartHitsPois = dynamic_cast<TH1D*>(hi.create1d(idx+1,"NbrParticleHits_Poisson","Number of Hits",particles.size()+1,0,particles.size()+1));
	////TDirectory * rootDir = gDirectory;
	//for(size_t i=1; i< particles.size();i++)
	//{
	//	//create a histogram that will take the normalized Data//
	//	TH1D* nbrHitsNorm = dynamic_cast<TH1D*>(hi.create1d( (idx+20*i) + 1,"NbrOfHits_normalized","Number of Hits",100,0,100,Form("%s",particles[i].GetName())));
	//	TH1D* nbrHitsPoisson = dynamic_cast<TH1D*>(hi.create1d( (idx+20*i) + 2,"NbrOfHits_Poisson","Number of Hits",100,0,100,Form("%s",particles[i].GetName())));
	//	//get histogram of Hits
	//	TH1D* nbrHits = dynamic_cast<TH1D*>( gFile->GetDirectory(Form("%s",particles[i].GetName()))->FindObject("NumberOfHits") );
	//	TH1D* nbrParticleHits = dynamic_cast<TH1D*>(gFile->FindObject("NumberOfParticleHits"));

	//	if (nbrHits)
	//	{
	//		for (size_t x=0; x<nbrHits->GetXaxis()->GetNbins(); x++)
	//		{
	//			//nbrHitsNorm->SetBinContent(x,nbrHits->GetBinContent(x)/nbrHits->GetEntries());
	//			nbrHitsNorm->Fill(x,nbrHits->GetBinContent(x+1)/nbrHits->GetEntries());
	//		}
	//		const double n0 = nbrHitsNorm->GetBinContent(1);
	//		double mean = 0.;
	//		if (n0 > 0.) mean = TMath::Log(1/n0);
	//		else continue;

	//		for (Int_t j=0; j<100; j++)
	//		{
	//			double P = TMath::Exp(-mean)*TMath::Power(mean,j)/CalcFactorial(j);
	//			nbrHitsPoisson->Fill(j,P);
	//		}
	//	}
	//	nbrPartHitsNorm->Fill(i,nbrHitsNorm->GetMean());
	//	nbrPartHitsPois->Fill(i,nbrHitsPoisson->GetMean());
	//}

}

