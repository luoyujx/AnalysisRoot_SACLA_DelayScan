#include "AnalyzeFunktions.h"
#include "AddParticles.h"

#include "./MyMomentaCalculator/MyMomentaCalculator.h"
#include "./MyParticle/MyParticleContainer.h"
#include "./MyAnalyzer/MyAnalyzer.h"
#include "FilesFromLma2Root/MyRootManager/MyHistos.h"
#include "FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalChannel/MyOriginalChannel.h"
#include "FilesFromLma2Root/MyEvent/MySortedEvent/MyDetektor/MyDetektor.h"
#include "FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalEventInfo.h"
#include "FilesFromLma2Root/MyEvent/MySortedEvent/MySortedEventInfo.h"
#include "FilesFromLma2Root/MyEvent/MySignalAnalyzedEvent/MySignalAnalyzedEventInfo.h"

#include <iostream>
#include <math.h>
#include <TTree.h>
#include <TFile.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TGClient.h>

//#define XVelocity 0.0003704 //(mm/ns)
//#define YVelocity -0.00037775
#define XVelocity 0. //(mm/ns)
#define YVelocity 0.

//_____________________________functions______________________________________________________________________________________________________________________________
//momentum calculation
double calcPx(const MyParticle &p, const MyParticleHit &ph)
{
	return MyMomentaCalculator::px(ph.XCorRotScl(),ph.YCorRotScl(),ph.TofCor(),p.GetMass_au(),p.GetCharge_au(),p.GetSpectrometer());
}
double calcPy(const MyParticle &p, const MyParticleHit &ph)
{
	return MyMomentaCalculator::py(ph.XCorRotScl(),ph.YCorRotScl(),ph.TofCor(),p.GetMass_au(),p.GetCharge_au(),p.GetSpectrometer());
}
double calcPz(const MyParticle &p, const MyParticleHit &ph)
{
	return MyMomentaCalculator::pz(ph.TofCor(),p.GetMass_au(),p.GetCharge_au(),p.GetSpectrometer());
}
double calcMass(const MyParticle &p, const MyParticleHit &ph)
{
	return MyMomentaCalculator::mass(ph.TofCor(),p.GetSpectrometer());
}



//___________________________________________________________________________________________________________________________________________________________
void DefineParticlesAndRootFile(MyParticleContainer &particles, MyHistos &hi)
{
	//setup the particle Infos, this will also load the other settings from the ini files//
	// kind of particle 0:atom 1:molecule 2:electron 3:Ion 4:Ion(simple) 

	particles.Add("Ion",1,1,0);//------------particle 0 --- Do not comment out!

	//---SACLA Ar atom
	//AddArgon(particles);

	//---SACLA Xe atom
	//AddXenon(particles);
	//AddXenon132(particles);

	//---SACLA CH3I molecule
	//AddCH3I(particles);

	//---Test N2 molecule
	AddNitrogen(particles);
}
//_______SACLA 2012A______________________________________________________________________________
void TofCorrection(MyDetektorHit &dh, const double alpha, const double k2, const double k4, const double xc, const double t0)
{
	const double t = dh.Time()-t0;
	const double x = dh.X()-xc;
	dh.SetTof( t - alpha*t*(k2*x*x+k4*x*x*x*x) );
}

//________Analysis main loop___________________________________________________________________________________________________________________________________________________
void MyAnalyzer::Analyze()
{
	MyDetektor &rd = fSE.GetDetektor(0);
	
	int startIdx=0;
	
	//clear the containers//
	fIntensities.clear();
	fParticles.ClearParticles();

	//Get Tag number on this event
	const unsigned int TagNumber = fOE.GetEventID();

	//Get intensity from loaded data map
	if (intFileName != "")
	{
		std::map<unsigned int, double>::iterator itTagInt;//BM1
		std::map<unsigned int, double>::iterator itTagInt2;//PD
		itTagInt = tagIntensity.find(TagNumber);
		itTagInt2 = tagIntensity2.find(TagNumber);

		if ((itTagInt == tagIntensity.end())||(itTagInt2 == tagIntensity2.end()))
		{
			std::cout <<TagNumber<< " is not found!!" << std::endl;
			return;
		}

		fIntensities.push_back(itTagInt2->second * 1000);//[0]
		fIntensities.push_back(itTagInt->second);//[1]

		////------------------------------------------SACLA 2012A
		////PhotoDiode intensity
		//if ( itTagInt->first < 156088080 )
		//	//before changing Al plate
		//	//fIntensities.push_back((tagIntensity2.find(TagNumber)->second)*1000*36.15);
		//	fIntensities.push_back((tagIntensity2.find(TagNumber)->second)/1.926E-5);//5.5keV PD->microJ(weak setting)
		//else
		//	fIntensities.push_back((tagIntensity2.find(TagNumber)->second)/6.963E-4);//5.5keV PD->microJ
		//	//fIntensities.push_back((tagIntensity2.find(TagNumber)->second)/1.17E-3);//5.0keV PD->microJ

		////FEL is attenuated, set tag region
		//const unsigned int tagFrom = 167982762;	//att0p4(5.5keV):159612312		Al25micro(5.0keV):167982762
		//const unsigned int tagTo = 170567544;	//att0p4(5.5keV):161807544		Al25micro(5.0keV):170567544
		//const double attFactor = 0.381;//5.5keV:0.381  5.5keV:0.52(test)  5.0keV:0.271
		//if (( tagFrom <= itTagInt->first ) && ( itTagInt->first <= tagTo ))
		//	fIntensities.push_back(itTagInt->second * attFactor);
		//else 
		//	fIntensities.push_back(itTagInt->second);
		////----------------------------------------------------

		for (size_t i=0;i<fIntensities.size();++i)
		{
			fHi.fill(startIdx++,Form("Int%d",i),fIntensities[i],Form("Int%d",i),1000,0,1000,"Intensity");
			for (size_t j=i+1; j<fIntensities.size();++j)
			{
				fHi.fill(startIdx++,Form("IntDep(%d)(%d)",i,j),fIntensities[i], fIntensities[j],Form("Int %d",i),Form("Int %d",j),1000,0,1000,1000,0,1000,"Intensity");
			}
		}
		if (fIntensities.size() > 1)
			if ( fIntensities[0]*fIntensities[1]>0)
				fHi.fill(startIdx+1,"IntensityDivide0by1",fIntensities[0]/fIntensities[1],"IntensityDivide0by1",1000,0,5,"Intensity");
		startIdx+=2;
	}



	//-------------------------get the raw mcp events called times//
	const MySignalAnalyzedChannel &sac = fSAE.GetChannel(7-1);
	for (size_t i=0; i<sac.GetNbrPeaks();++i)
	{
		const MyPeak &p = sac.GetPeak(i);
		fHi.fill(startIdx,"Times_MCP",p.GetTime(),"tof [ns]",5000,0,fOE.GetNbrSamples()*fOE.GetSampleInterval()*1e9);
	}
	startIdx++;

	//if(fSAE.GetChannel(7-1).GetNbrPeaks())
	//{
	//	const MyPeak &p_t0 = fSAE.GetChannel(7-1).GetPeak(0);
	//	fHi.fill(startIdx,"Times_MCP_first",p_t0.GetTime(),"tof [ns]",10000,0,100000);
	//}
	//startIdx++;

	//-------------------------------------------------------------------------------------------------------------------------------
	if (fIntensities.size()) fHi.fill(startIdx,"IntensityBM1",fIntensities[1],"microJ",1000,0,1000);//added by moto 2009/05/21(50,0,2500)
	startIdx++;
	if (fIntensities.size()) fHi.fill(startIdx,"IntensityPD",fIntensities[0],"arb. unit",1000,0,1000);//added by moto 2009/05/21(50,0,2500)
	startIdx++;

	//---skip when FEL is stopped
	if (fIntensities.size())
		if (fIntensities[0]<1) 
		{
			//std::cout << "skip this event" << std::endl;
			return;
		}

	if (fIntensities.size()) fHi.fill(startIdx,"IntensitySelectPD",fIntensities[0],"arb unit",1000,0,1000);//added by moto 2009/05/21
	startIdx++;

	fHi.fill(startIdx,"NumberOfHits",rd.GetNbrOfHits(),"Number of Hits",100,0,100);
	startIdx++;
	int secondStartIdx=startIdx+50;

	//-----------------------Set jet speed
	//for (size_t i=0;i<fParticles.GetNbrOfParticles();++i)
	//{
	//	MyParticle &p = fParticles.GetParticle(i);
	//	p.SetXVelocity(XVelocity);
	//	p.SetYVelocity(YVelocity);
	//}

	//go through all resorted detektorhits//
	for (size_t i=0; i<rd.GetNbrOfHits();++i)
	{
		MyDetektorHit &dh = rd.GetHit(i);

		//the tof is just the timing of the mcp signal//
		dh.SetTof(dh.Time());

		////rotate det image
		//-----------Position correction(SACLA 2012A Spectromertor D" 520V)----------//
		//const double angle = -8*TMath::DegToRad();
		//dh.SetXmm(TMath::Cos(angle)*dh.X() + TMath::Sin(angle)*dh.Y());
		//dh.SetYmm(-TMath::Sin(angle)*dh.X() + TMath::Cos(angle)*dh.Y());

		//detektor & tof for all Hits//
		const double maxPos	= (rd.GetRunTime()+30)*rd.GetSfU();
		MyParticle &LastParticle = fParticles.GetParticle(fParticles.GetNbrOfParticles()-1);
		//const double maxTof = LastParticle.GetCondTofTo() + LastParticle.GetCondTofRange()*0.3;
		const double maxTof	= fOE.GetNbrSamples()*fOE.GetSampleInterval()*1e9;

		fHi.fill(startIdx+1,"DetAll",dh.X(),dh.Y(),"x [mm]","y [mm]",300,-maxPos,maxPos,300,-maxPos,maxPos);
		fHi.fill(startIdx+2,"TofAll",dh.Tof(),"tof [ns]",10000,0,maxTof);
		fHi.fill(startIdx+3,"XPosVsTofAll",dh.Tof(),dh.X(),"tof [ns]","x [mm]",5000,0,maxTof,300,-maxPos,maxPos);
		fHi.fill(startIdx+4,"YPosVsTofAll",dh.Tof(),dh.Y(),"tof [ns]","y [mm]",5000,0,maxTof,300,-maxPos,maxPos);
		fHi.fill(startIdx+5,"TofFine",dh.Tof(),"tof [ns]",static_cast<int>(maxTof),0,maxTof);
		fHi.fill(startIdx+6,"XPosVsTofFine",dh.Tof(),dh.X(),"tof [ns]","x [mm]",static_cast<int>(maxTof/5),0,maxTof,300,-maxPos,maxPos);
		fHi.fill(startIdx+7,"YPosVsTofFine",dh.Tof(),dh.Y(),"tof [ns]","y [mm]",static_cast<int>(maxTof/5),0,maxTof,300,-maxPos,maxPos);

		
		//-----Particle(0)---Ion---
		//get the particle from the vector//
		MyParticle &p = fParticles.GetParticle(0);
		//if this hit fits both conditions then add the Hit to this Particle and fill the histo for this hit//
		if (p.CheckTofAndPos(dh))
		//select hit by reconstruction method//
		if (dh.RekMeth() < rekmeth)//added by motomura
		{
			const MyParticleHit &ph = p.AddHit(dh);

			fHi.fill(startIdx+8,"Mass",ph.Mass(),"Mass/q",20000,0,200,"Ion");
			fHi.fill(startIdx+9,"TofCor",ph.TofCor(),"tof [ns]",20000,p.GetCondTofFr()-p.GetT0()-p.GetCondTofRange()*0.3,p.GetCondTofTo()-p.GetT0()+p.GetCondTofRange()*0.3,"Ion");
			fHi.fill(startIdx+10,"Det",ph.X(),ph.Y(),"x [mm]","y [mm]",300,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,300,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,"Ion");
			fHi.fill(startIdx+11,"Tof",ph.Tof(),"tof [ns]",10000,p.GetCondTofFr()-p.GetCondTofRange()*0.3,p.GetCondTofTo()+p.GetCondTofRange()*0.3,"Ion");
			fHi.fill(startIdx+12,"XPosVsTof",ph.Tof(),ph.X(),"tof [ns]","x [mm]",5000,p.GetCondTofFr()-p.GetCondTofRange()*0.3,p.GetCondTofTo()+p.GetCondTofRange()*0.3,300,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,"Ion");
			fHi.fill(startIdx+13,"YPosVsTof",ph.Tof(),ph.Y(),"tof [ns]","y [mm]",5000,p.GetCondTofFr()-p.GetCondTofRange()*0.3,p.GetCondTofTo()+p.GetCondTofRange()*0.3,300,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,"Ion");
		}

		//-----------Tof correction by position (SACLA 2012A Spectromertor D" 520V)----------//
		//if (extraCondition) TofCorrection(dh,  1.25, 2.6852e-5, -1.1172e-8, 0., 770);

		secondStartIdx=startIdx+20;
		//check wether this hit fits the conditions of the particles we created//
		//if so, the add this hit to the particle and fill the particle histograms//
		//go through all particles//
		for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)
		{
			//get the particle from the vector//
			MyParticle &p = fParticles.GetParticle(j);
			//if this hit fits the tof conditions of this particle fill histos for the tof condition//
			//if (p.CheckTof(dh)) fillParticleConditionsTof(fOE,rd,p,dh,fIntensities,fHi,secondStartIdx + 0);
			//if this hit fits the pos conditions of this particle fill histos for the pos condition//
			//if (p.CheckPos(dh)) fillParticleConditionsPos(fOE,rd,p,dh,fIntensities,fHi,secondStartIdx + 5);
			//if this hit fits both conditions then add the Hit to this Particle and fill the histo for this hit//
			if (p.CheckTofAndPos(dh))
			//select hit by reconstruction method//
			if (dh.RekMeth() < rekmeth)
			{
				const MyParticleHit &ph = p.AddHit(dh);
				fillParticleHistograms(p,ph,fIntensities,fHi,secondStartIdx);
			}
			//we reserve 100 histograms for one particle//
			secondStartIdx +=100;
			//std::cout <<j<<" "<< secondStartIdx<<std::endl;
		}
	}//------------------------------------------------------------------------------------------------------//

	startIdx += (fParticles.GetNbrOfParticles()*100);

	for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)//particle index j=0 is "ion" 
	{
		fHi.fill(startIdx+j,"NumberOfHits",fParticles.GetParticle(j).GetNbrOfParticleHits(),"Number of Hits",100,0,100,Form("%s",fParticles.GetParticle(j).GetName()));

		for (size_t k=0;k<fParticles.GetParticle(j).GetNbrOfParticleHits();++k)
			fHi.fill(startIdx+fParticles.GetNbrOfParticles(),"NumberOfParticleHits",j,"Particle Number",fParticles.GetNbrOfParticles()+1,0,fParticles.GetNbrOfParticles()+1);
	}

	startIdx += (fParticles.GetNbrOfParticles())+1;

	//now you have found the particles//
	//we can look for coincidences//
	if (MoleculeAnalysis)
		for (size_t i=1;i<fParticles.GetNbrOfParticles();++i)//from particle No.1
		{
			const MyParticle &ip = fParticles.GetParticle(i);
			for (size_t j=i;j<fParticles.GetNbrOfParticles();++j)
			{
				const MyParticle &jp = fParticles.GetParticle(j);
				if ((ip.GetKindParticle() == 1)&&(jp.GetKindParticle() == 1))
				if (
					((ip.GetCoinGroup() != jp.GetCoinGroup()))
					||((ip.GetCoinGroup()==100)&&(jp.GetCoinGroup()==100))
					)
				{
					fillMoleculeHistogram(ip,jp,fIntensities,fHi,startIdx);
					startIdx += 300;
				}
			}
		}
	if (MoleculeAnalysis) fillPIPICO(fParticles.GetParticle(0),fHi);

}

//--------------------------------------fill particle Histograms---------------------------------------------------//
void fillParticleHistograms(const MyParticle &p, const MyParticleHit &ph, std::vector<double>& intensity, MyHistos &hi, int hiOff )
{	
	const double MomLim = 400;
	const double SliceLim = 20;

	//Reconstruction Method for this particle Hit//
	hi.fill(hiOff++,"ReconstructionMethod",ph.RekMeth(),"Reconstruction Nbr",60,0,30,Form("%s/Raw",p.GetName()));
	if (intensity.size()) hi.fill(hiOff++,Form("Intensity%s",p.GetName()),intensity[0],"Laser Power",1000,0,1000,Form("%s",p.GetName()));

	//detektor pictures//
	hi.fill(hiOff++,"Det"			,ph.X()			,ph.Y()			,"x [mm]","y [mm]",1000,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,1000,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	hi.fill(hiOff++,"DetCor"		,ph.XCor()		,ph.YCor()		,"x [mm]","y [mm]",1000,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,1000,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	hi.fill(hiOff++,"DetCorRot"		,ph.XCorRot()	,ph.YCorRot()	,"x [mm]","y [mm]",1000,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,1000,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	hi.fill(hiOff++,"DetCorRotScale",ph.XCorRotScl(),ph.YCorRotScl(),"x [mm]","y [mm]",1000,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,1000,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));

	//tof//
	hi.fill(hiOff++,"Tof"   ,ph.Tof()   ,"tof [ns]",1000,p.GetCondTofFr()          -p.GetCondTofRange()*0.3,p.GetCondTofTo()          +p.GetCondTofRange()*0.3,Form("%s/Raw",p.GetName()));
	hi.fill(hiOff++,"TofCor",ph.TofCor(),"tof [ns]",1000,p.GetCondTofFr()-p.GetT0()-p.GetCondTofRange()*0.3,p.GetCondTofTo()-p.GetT0()+p.GetCondTofRange()*0.3,Form("%s/Raw",p.GetName()));

	//pos vs tof//
	hi.fill(hiOff++,"XPosVsTof",ph.TofCor(),ph.XCorRotScl(),"tof [ns]","x [mm]",500,p.GetCondTofFr()-p.GetT0()-p.GetCondTofRange()*0.3,p.GetCondTofTo()-p.GetT0()+p.GetCondTofRange()*0.3,1000,p.GetCondRadX()-p.GetXcor()-p.GetCondRad()*1.3,p.GetCondRadX()-p.GetXcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	hi.fill(hiOff++,"YPosVsTof",ph.TofCor(),ph.YCorRotScl(),"tof [ns]","y [mm]",500,p.GetCondTofFr()-p.GetT0()-p.GetCondTofRange()*0.3,p.GetCondTofTo()-p.GetT0()+p.GetCondTofRange()*0.3,1000,p.GetCondRadY()-p.GetYcor()-p.GetCondRad()*1.3,p.GetCondRadY()-p.GetYcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));

	hi.fill(hiOff++,"XPosVsTofFine",ph.TofCor(),ph.XCorRotScl(),"tof [ns]","x [mm]",static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),300,p.GetCondRadX()-p.GetXcor()-p.GetCondRad()*1.3,p.GetCondRadX()-p.GetXcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	hi.fill(hiOff++,"YPosVsTofFine",ph.TofCor(),ph.YCorRotScl(),"tof [ns]","y [mm]",static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),300,p.GetCondRadY()-p.GetYcor()-p.GetCondRad()*1.3,p.GetCondRadY()-p.GetYcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));

	//const double tofPosCorr= ph.TofCor();
	//hi.fill(hiOff++,"TofPosCorr"   ,tofPosCorr   ,"tof [ns]",1000,p.GetCondTofFr()          -p.GetCondTofRange()*0.3,p.GetCondTofTo()          +p.GetCondTofRange()*0.3,Form("%s/Raw",p.GetName()));
	//hi.fill(hiOff++,"XPosVsTofFinePosCorr",tofPosCorr,ph.XCorRotScl(),"tof [ns]","x [mm]",static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),100,p.GetCondRadX()-p.GetXcor()-p.GetCondRad()*1.3,p.GetCondRadX()-p.GetXcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	//hi.fill(hiOff++,"YPosVsTofFinePosCorr",tofPosCorr,ph.YCorRotScl(),"tof [ns]","y [mm]",static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),100,p.GetCondRadY()-p.GetYcor()-p.GetCondRad()*1.3,p.GetCondRadY()-p.GetYcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	//hi.fill(hiOff++,"YPosVsTofFinePosCorr0.5ns",tofPosCorr,ph.YCorRotScl(),"tof [ns]","y [mm]",2*static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),100,p.GetCondRadY()-p.GetYcor()-p.GetCondRad()*1.3,p.GetCondRadY()-p.GetYcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	////hi.fill(hiOff++,"YPosVsTofFinePosCorr0.25ns",tofPosCorr,ph.YCorRotScl(),"tof [ns]","y [mm]",4*static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),100,p.GetCondRadY()-p.GetYcor()-p.GetCondRad()*1.3,p.GetCondRadY()-p.GetYcor()+p.GetCondRad()*1.3,Form("%s/Raw",p.GetName()));
	//hi.fill(hiOff++,"TofFinePosCorr",tofPosCorr,"tof [ns]",static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),Form("%s/Raw",p.GetName()));
	//hi.fill(hiOff++,"TofFinePosCorr0.5ns",tofPosCorr,"tof [ns]",2*static_cast<int>(p.GetCondTofRange()),p.GetCondTofFr()-p.GetT0(),p.GetCondTofTo()-p.GetT0(),Form("%s/Raw",p.GetName()));

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
	hi.fill(hiOff++,"Energy",ph.E(),"Energy [eV]",300,0,40,Form("%s/Energy",p.GetName()));
	//Angle//
	hi.fill(hiOff++,"ThetaX",ph.ThetaX(),"#theta [deg]",180,0,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"ThetaY",ph.ThetaY(),"#theta [deg]",180,0,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"ThetaZ",ph.ThetaZ(),"#theta [deg]",180,0,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiXY",ph.PhiXY(),"#phi [deg]",360,-180,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiYZ",ph.PhiYZ(),"#phi [deg]",360,-180,180,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiZX",ph.PhiZX(),"#phi [deg]",360,-180,180,Form("%s/Angle",p.GetName()));

	hi.fill(hiOff++,"ThetaZvsEnergy",ph.ThetaZ(),ph.E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,100,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"ThetaXvsEnergy",ph.ThetaX(),ph.E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,100,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"ThetaYvsEnergy",ph.ThetaY(),ph.E(),"#theta [deg]","Energy [eV]",180,0,180,200,0,100,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiXYvsEnergy",ph.PhiXY(),ph.E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,100,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiYZvsEnergy",ph.PhiYZ(),ph.E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,100,Form("%s/Angle",p.GetName()));
	hi.fill(hiOff++,"PhiZXvsEnergy",ph.PhiZX(),ph.E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,100,Form("%s/Angle",p.GetName()));

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
void fillMoleculeHistogram(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff)
{
	TString Hname(p1.GetName());
	Hname += p2.GetName();
	
	const double pxSumWidth = 20;//10,9,7,5
	const double pySumWidth = 10;//10,8,5,4
	const double pzSumWidth = 3;//5,6,4,2

	const double MomLim = 400;

	for (size_t i=0; i<p1.GetNbrOfParticleHits();++i)
	{
		for (size_t j=0;j<p2.GetNbrOfParticleHits();++j)//i=j
		{
			//skip if we are going to put in the same particlehit, for the case that we want coincidences for the same particle//
			if ((i>=j) && (p1==p2)) continue;//i==j

			//x momentum sum//
			hi.fill(hiOff+1,"PxSum",(p1[i].Px() + p2[j].Px()),Form("Px_{%s} + Px_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-100,100,Form("%s/MomSums",Hname.Data()));
			if (TMath::Abs(p1[i].Py() + p2[j].Py()) < pySumWidth )
			if (TMath::Abs(p1[i].Pz() + p2[j].Pz()) < pzSumWidth )
				hi.fill(hiOff+2,"PxSumCond",(p1[i].Px() + p2[j].Px()),Form("Px_{%s} + Px_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-100,100, Form("%s/MomSums",Hname.Data()));
			
			//y momentum sum//
			hi.fill(hiOff+3,"PySum",(p1[i].Py() + p2[j].Py()),Form("Py_{%s} + Py_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-100,100, Form("%s/MomSums",Hname.Data()));
			if (TMath::Abs(p1[i].Px() + p2[j].Px()) < pxSumWidth )
			if (TMath::Abs(p1[i].Pz() + p2[j].Pz()) < pzSumWidth )
				hi.fill(hiOff+4,"PySumCond",(p1[i].Py() + p2[j].Py()),Form("Py_{%s} + Py_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-100,100, Form("%s/MomSums",Hname.Data()));

			//z momentum sum//
			hi.fill(hiOff+5,"PzSum",(p1[i].Pz() + p2[j].Pz()),Form("Pz_{%s} + Pz_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-100,100, Form("%s/MomSums",Hname.Data()));
			if (TMath::Abs(p1[i].Py() + p2[j].Py()) < pySumWidth )
			if (TMath::Abs(p1[i].Px() + p2[j].Px()) < pxSumWidth )
				hi.fill(hiOff+6,"PzSumCond",(p1[i].Pz() + p2[j].Pz()),Form("Pz_{%s} + Pz_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-100,100, Form("%s/MomSums",Hname.Data()));

			//PIPICO//
			hi.fill(hiOff+7,"PIPICO",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());
			if (TMath::Abs(p1[i].Px() + p2[j].Px()) < pxSumWidth )
			if (TMath::Abs(p1[i].Py() + p2[j].Py()) < pySumWidth )
				hi.fill(hiOff+8,"PIPICOCondXY",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());

			{
				//x momentum sum//
				if (TMath::Abs(p1[i].Py() + p2[j].Py()) < pySumWidth*3 )
				if (TMath::Abs(p1[i].Pz() + p2[j].Pz()) < pzSumWidth*3 )
					hi.fill(hiOff+12,"PxSumCond_Width3",(p1[i].Px() + p2[j].Px()),Form("Px_{%s} + Px_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-100,100, Form("%s/MomSums",Hname.Data()));
				
				//y momentum sum//
				if (TMath::Abs(p1[i].Px() + p2[j].Px()) < pxSumWidth*3 )
				if (TMath::Abs(p1[i].Pz() + p2[j].Pz()) < pzSumWidth*3 )
					hi.fill(hiOff+14,"PySumCond_Width3",(p1[i].Py() + p2[j].Py()),Form("Py_{%s} + Py_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-100,100, Form("%s/MomSums",Hname.Data()));

				//z momentum sum//
				if (TMath::Abs(p1[i].Py() + p2[j].Py()) < pySumWidth*3 )
				if (TMath::Abs(p1[i].Px() + p2[j].Px()) < pxSumWidth*3 )
					hi.fill(hiOff+16,"PzSumCond_Width3",(p1[i].Pz() + p2[j].Pz()),Form("Pz_{%s} + Pz_{%s} [a.u]",p1.GetName(),p2.GetName()),400,-100,100, Form("%s/MomSums",Hname.Data()));
			}

			///////////////////Momentum sum condition//////////////////////

			//Momenta & Energy & Intensity//
			if (TMath::Abs(p1[i].Px() + p2[j].Px()) < pxSumWidth)
			if (TMath::Abs(p1[i].Py() + p2[j].Py()) < pySumWidth)
			if (TMath::Abs(p1[i].Pz() + p2[j].Pz()) < pzSumWidth)
			{
				//PIPICO//
				hi.fill(hiOff+9,"PIPICOCondXYZ",p1[i].TofCor(),p2[j].TofCor(),Form("tof_{%s} [ns]",p1.GetName()),Form("tof_{%s} [ns]",p2.GetName()),300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());

				//Intensity//
				//for (size_t iI=0; iI<intensity.size();++iI)
				//hi.fill(hiOff+10+iI,"Intensity",intensity[iI],Form("intensity_{%s%s} [arb.units]",p1.GetName(),p2.GetName()),300,0,1000,Form("%s/Intensity",Hname.Data()));
				if (intensity.size()) hi.fill(hiOff+10,Form("intensity%s%s",p1.GetName(),p2.GetName()),intensity[0],Form("intensity_{%s%s} [arb.units]",p1.GetName(),p2.GetName()),50,0,2500);

				//Momentum first Ion//
				int IDX = hiOff+100+intensity.size();
				hi.fill(IDX+0,Form("%sPxPy",p1.GetName()),p1[i].Px(),p1[i].Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p1[i].Pz()) < 30)
				{
					hi.fill(IDX+1,Form("%sPxPySlice",p1.GetName()),p1[i].Px(),p1[i].Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					hi.fill(IDX+2,Form("%sPhiPxPySlice",p1.GetName()),TMath::ATan2(p1[i].Py(),p1[i].Px())*TMath::RadToDeg(),TMath::Sqrt(p1[i].Py()*p1[i].Py() + p1[i].Px()*p1[i].Px()),"#phi [deg]","#sqrt{px^{2} + py^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}
           
				hi.fill(IDX+3,Form("%sPxPz",p1.GetName()),p1[i].Pz(),p1[i].Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p1[i].Py()) < 30)
				{
					hi.fill(IDX+4,Form("%sPxPzSlice",p1.GetName()),p1[i].Pz(),p1[i].Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					hi.fill(IDX+5,Form("%sPhiPxPzSlice",p1.GetName()),TMath::ATan2(p1[i].Px(),p1[i].Pz())*TMath::RadToDeg(),TMath::Sqrt(p1[i].Pz()*p1[i].Pz() + p1[i].Px()*p1[i].Px()),"#phi [deg]","#sqrt{px^{2} + pz^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}
           
				hi.fill(IDX+6,Form("%sPyPz",p1.GetName()),p1[i].Pz(),p1[i].Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p1[i].Px()) < 30)
				{
					hi.fill(IDX+7,Form("%sPyPzSlice",p1.GetName()),p1[i].Pz(),p1[i].Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					hi.fill(IDX+8,Form("%sPhiPyPzSlice",p1.GetName()),TMath::ATan2(p1[i].Py(),p1[i].Pz())*TMath::RadToDeg(),TMath::Sqrt(p1[i].Py()*p1[i].Py() + p1[i].Pz()*p1[i].Pz()),"#phi [deg]","#sqrt{pz^{2} + py^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}

				hi.fill(IDX+9,Form("%sTotalMomentum",p1.GetName()),p1[i].P(),"p [a.u.]",300,0,MomLim,Form("%s/Momenta",Hname.Data()));

				//Energy First Ion//
				hi.fill(IDX+10,Form("%sEnergy",p1.GetName()),p1[i].E(),"Energy [eV]",300,0,60,Form("%s/Energy",Hname.Data()));

				//Raw
				hi.fill(IDX+11,Form("%sTOF",p1.GetName()),p1[i].TofCor(),"Tof [ns]",300,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,Hname.Data());
				hi.fill(IDX+12,Form("%sDetCor",p1.GetName()),p1[i].XCor(),p1[i].YCor(),"x [mm]","y [mm]",300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,300,0-p1.GetCondRad()*1.3,0+p1.GetCondRad()*1.3,Hname.Data());
				hi.fill(IDX+13,Form("%sXPosVsTof",p1.GetName()),p1[i].TofCor(),p1[i].XCorRotScl(),"tof [ns]","x [mm]",500,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p1.GetXcor()-p1.GetCondRad()*1.3,p1.GetXcor()+p1.GetCondRad()*1.3,Hname.Data());
				hi.fill(IDX+14,Form("%sYPosVsTof",p1.GetName()),p1[i].TofCor(),p1[i].YCorRotScl(),"tof [ns]","y [mm]",500,p1.GetCondTofFr()-p1.GetT0()-p1.GetCondTofRange()*0.3,p1.GetCondTofTo()-p1.GetT0()+p1.GetCondTofRange()*0.3,300,p1.GetYcor()-p1.GetCondRad()*1.3,p1.GetYcor()+p1.GetCondRad()*1.3,Hname.Data());

  				//Angular Distribution//added by motomura
				hi.fill(IDX+21,Form("%sThetaX",p1.GetName()),p1[i].ThetaX(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+22,Form("%sThetaY",p1.GetName()),p1[i].ThetaY(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+23,Form("%sThetaZ",p1.GetName()),p1[i].ThetaZ(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+24,Form("%sThetaYvsEnergy",p1.GetName()),p1[i].ThetaY(),p1[i].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
 				hi.fill(IDX+25,Form("%sThetaXvsEnergy",p1.GetName()),p1[i].ThetaX(),p1[i].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+26,Form("%sThetaZvsEnergy",p1.GetName()),p1[i].ThetaZ(),p1[i].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+27,Form("%sPhiXYvsEnergy",p1.GetName()),p1[i].PhiXY(),p1[i].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
 				hi.fill(IDX+28,Form("%sPhiYZvsEnergy",p1.GetName()),p1[i].PhiYZ(),p1[i].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+29,Form("%sPhiZXvsEnergy",p1.GetName()),p1[i].PhiZX(),p1[i].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));

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
				IDX = (p1!=p2)?IDX+100:IDX;
				hi.fill(IDX+0,Form("%sPxPy",p2.GetName()),p2[j].Px(),p2[j].Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p2[j].Pz()) < 30)
				{
					hi.fill(IDX+1,Form("%sPxPySlice",p2.GetName()),p2[j].Px(),p2[j].Py(),"px [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					hi.fill(IDX+2,Form("%sPhiPxPySlice",p2.GetName()),TMath::ATan2(p2[j].Py(),p2[j].Px())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Py()*p2[j].Py() + p2[j].Px()*p2[j].Px()),"#phi [deg]","#sqrt{px^{2} + py^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}
           
				hi.fill(IDX+3,Form("%sPxPz",p2.GetName()),p2[j].Pz(),p2[j].Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p2[j].Py()) < 30)
				{
					hi.fill(IDX+4,Form("%sPxPzSlice",p2.GetName()),p2[j].Pz(),p2[j].Px(),"pz [a.u.]","px [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					hi.fill(IDX+5,Form("%sPhiPxPzSlice",p2.GetName()),TMath::ATan2(p2[j].Px(),p2[j].Pz())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Pz()*p2[j].Pz() + p2[j].Px()*p2[j].Px()),"#phi [deg]","#sqrt{px^{2} + pz^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}
           
				hi.fill(IDX+6,Form("%sPyPz",p2.GetName()),p2[j].Pz(),p2[j].Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
				if (TMath::Abs(p2[j].Px()) < 30)
				{
					hi.fill(IDX+7,Form("%sPyPzSlice",p2.GetName()),p2[j].Pz(),p2[j].Py(),"pz [a.u.]","py [a.u.]",300,-MomLim,MomLim,300,-MomLim,MomLim,Form("%s/Momenta",Hname.Data()));
					hi.fill(IDX+8,Form("%sPhiPyPzSlice",p2.GetName()),TMath::ATan2(p2[j].Py(),p2[j].Pz())*TMath::RadToDeg(),TMath::Sqrt(p2[j].Py()*p2[j].Py() + p2[j].Pz()*p2[j].Pz()),"#phi [deg]","#sqrt{pz^{2} + py^{2}} [a.u.]",360,-180,180,300,0,MomLim,Form("%s/Momenta",Hname.Data()));
				}

				hi.fill(IDX+9,Form("%sTotalMomentum",p2.GetName()),p2[j].P(),"p [a.u.]",300,0,MomLim,Form("%s/Momenta",Hname.Data()));

				//Energy second Ion//
				hi.fill(IDX+10,Form("%sEnergy",p2.GetName()),p2[j].E(),"Energy [eV]",300,0,50,Form("%s/Energy",Hname.Data()));

				//Raw
				hi.fill(IDX+11,Form("%sTOF",p2.GetName()),p2[j].TofCor(),"Tof [ns]",300,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,Hname.Data());
				hi.fill(IDX+12,Form("%sDetCor",p2.GetName()),p2[j].XCor(),p2[j].YCor(),"x [mm]","y [mm]",300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,300,0-p2.GetCondRad()*1.3,0+p2.GetCondRad()*1.3,Hname.Data());
				hi.fill(IDX+13,Form("%sXPosVsTof",p2.GetName()),p2[j].TofCor(),p2[j].XCorRotScl(),"tof [ns]","x [mm]",500,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetXcor()-p2.GetCondRad()*1.3,p2.GetXcor()+p2.GetCondRad()*1.3,Hname.Data());
				hi.fill(IDX+14,Form("%sYPosVsTof",p2.GetName()),p2[j].TofCor(),p2[j].YCorRotScl(),"tof [ns]","y [mm]",500,p2.GetCondTofFr()-p2.GetT0()-p2.GetCondTofRange()*0.3,p2.GetCondTofTo()-p2.GetT0()+p2.GetCondTofRange()*0.3,300,p2.GetYcor()-p2.GetCondRad()*1.3,p2.GetYcor()+p2.GetCondRad()*1.3,Hname.Data());

				//Angular Distribution//added by motomura
				hi.fill(IDX+21,Form("%sThetaX",p2.GetName()),p2[j].ThetaX(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+22,Form("%sThetaY",p2.GetName()),p2[j].ThetaY(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+23,Form("%sThetaZ",p2.GetName()),p2[j].ThetaZ(),"#theta [deg]",36,0,180,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+24,Form("%sThetaYvsEnergy",p2.GetName()),p2[j].ThetaY(),p2[j].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
 				hi.fill(IDX+25,Form("%sThetaXvsEnergy",p2.GetName()),p2[j].ThetaX(),p2[j].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+26,Form("%sThetaZvsEnergy",p2.GetName()),p2[j].ThetaZ(),p2[j].E(),"#theta [deg]","Energy [eV]",360,0,180,200,0,40,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+27,Form("%sPhiXYvsEnergy",p2.GetName()),p2[j].PhiXY(),p2[j].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
 				hi.fill(IDX+28,Form("%sPhiYZvsEnergy",p2.GetName()),p2[j].PhiYZ(),p2[j].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));
				hi.fill(IDX+29,Form("%sPhiZXvsEnergy",p2.GetName()),p2[j].PhiZX(),p2[j].E(),"#phi [deg]","Energy [eV]",360,-180,180,200,0,40,Form("%s/Angular",Hname.Data()));

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
				hi.fill(IDX+50,Form("KER%s",Hname.Data()),p1[i].E()+p2[j].E(),"KER [eV]",200,0,100,Form("%s/KER",Hname.Data()));
				hi.fill(IDX+51,Form("ThetaY1vsKER%s",Hname.Data()),p1[i].ThetaY(),p1[i].E()+p2[j].E(),"#theta [deg]","KER [eV]",360,0,180,200,0,100,Form("%s/KER",Hname.Data()));
				hi.fill(IDX+52,Form("ThetaY2vsKER%s",Hname.Data()),p2[j].ThetaY(),p1[i].E()+p2[j].E(),"#theta [deg]","KER [eV]",360,0,180,200,0,100,Form("%s/KER",Hname.Data()));
				
				//if (((TMath::Abs(p1[i].PhiZX()) > 35)&&(TMath::Abs(p1[i].PhiZX()) < 150))
				//	&&((TMath::Abs(p2[j].PhiZX()) > 35)&&(TMath::Abs(p2[j].PhiZX()) < 150)))
				//{
				//	hi.fill(IDX+60,Form("KER%s_selectPhiXZ",Hname.Data()),p1[i].E()+p2[j].E(),"KER [eV]",200,0,100,Form("%s/KER_masked",Hname.Data()));
				//}

			}
		}
	}
}
//----------------------------------PIPICO ALL-----------------------------------------------------------------//
void fillPIPICO(const MyParticle &p,MyHistos &hi)
{	
	//const double pxSumWidth = 5;//10,9,7,5
	//const double pySumWidth = 5;//10,8,5,4
	//const double pzSumWidth = 5;//5,6,4,2

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




//----------------------------------extraCondition----------------------------------------------------------------//
bool PosCondition(const MyDetektorHit &dh)
{
	//Except residual gas
	//if ((dh.X() < -22) || (dh.X() > 21) || (dh.Y() < -6) || (dh.Y() > 9.5)) return true;
	if ((dh.X() < -23) || (dh.X() > 22) || (dh.Y() < -7) || (dh.Y() > 10.5)) return true;
	return false;
}
bool TofPosCondition(const MyDetektorHit &dh)
{
	////Except Mass 15
	//if (((dh.Tof() > 3922)&&(dh.Tof() < 3957)) && ((dh.Y() > -7.) && (dh.Y() < 10.)) && ((dh.X() > -23) && (dh.X() < 22))) return false;
	////Except Mass 17
	//if (((dh.Tof() > 4121)&&(dh.Tof() < 4161)) && ((dh.Y() > -7.) && (dh.Y() < 10.)) && ((dh.X() > -23) && (dh.X() < 22))) return false;

	//Mass15 N2_2009_Jan
	if (((dh.Tof() > 5600)&&(dh.Tof() < 5660)) && ((dh.Y() > -7.) && (dh.Y() < 8.)) && ((dh.X() > -20) && (dh.X() < 20))) return false;
	//if (((dh.Tof() > 5460)&&(dh.Tof() < 5473)) && ((dh.Y() > -4.) && (dh.Y() < 2.)) && ((dh.X() > -5) && (dh.X() < 4))) return false;

	return true;
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


void fillHistosAfterAnalyzis(const std::vector<MyParticle> &particles, MyHistos &hi,size_t nRegion)
{
}

//extract the intensity from the Channel
double smoothedVal(const short * Data, size_t idx)
{
	return (Data[idx-3]+Data[idx-2]+Data[idx-1]+Data[idx]+Data[idx+1]+Data[idx+2]+Data[idx+3])/7.;
}

//-------------------------------extract the intensity from the Channel-------------------------------------------------//

double Integral(const MyOriginalChannel &oc, const long TRfrom, const long TRto, bool absolute)
{
	//get some infos from the channel and initalize the integral//
	double IntegralTR		= 0;
	//const short baseli		= static_cast<short>(baseline / oc.GetVertGain());
	const double vertGain	= oc.GetVertGain();



	//go through all pulses in this channel//
	for (int puls=0; puls<oc.GetNbrPulses();++puls)
	{
		//get puls and infos from puls and pointer to the array//
		const MyPuls &p				= oc.GetPuls(puls);
		int actPosInEvent			= p.GetIndexToFirstPointOfOriginalWaveform();
		const short *Data			= static_cast<const short*>(oc.GetDataPointerForPuls(p));
		const size_t pLength		= p.GetLength();
		//find baseline//
		//const double baseli			= smoothedVal(Data,pLength-4);
		const double baseli			= oc.GetBaseline();

		//go through the puls and check wether the actual point is still in the range//
		//if that is the case add this point to the integral//
		for (size_t i=0; i<pLength ; ++i, ++actPosInEvent)
		{
			if ((TRfrom < actPosInEvent)&&(actPosInEvent < TRto))
			{
				if (absolute)
				{
					IntegralTR += (TMath::Abs(Data[i]-baseli));
				}
				else
				{
					IntegralTR += (Data[i]-baseli);
				}
			}
		}
	}

	//return the integral//
	return IntegralTR*vertGain;
}

double Average(const MyOriginalChannel &oc, const long TRfrom, const long TRto, bool absolute)
{
	//return the Average//
	return Integral(oc,TRfrom,TRto,absolute)/(TRto-TRfrom-1.);
}
