#include "AnalyzeFunktions.h"
#include "AnalyzeFuncFill.h"
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
#include <TStyle.h> 

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
	//return MyMomentaCalculator::pz(ph.TofCor(),p.GetMass_au(),p.GetCharge_au(),p.GetSpectrometer());
	return MyMomentaCalculator::pz_poly(ph.TofCor(),p.GetSpectrometer());
}
double calcMass(const MyParticle &p, const MyParticleHit &ph)
{
	return MyMomentaCalculator::mass(ph.TofCor(),p.GetSpectrometer());
}
double calcTof(const MyParticle &p, const MyParticle &pIon)
{
	return MyMomentaCalculator::tof(p.GetMass_au()*MyUnitsConv::au2amu()/p.GetCharge_au(),pIon.GetSpectrometer());
}

//___________________________________________________________________________________________________________________________________________________________
void DefineParticlesAndRootFile(MyParticleContainer &particles, MyHistos &hi, const TString &whichParticles)
{
	//setup the particle Infos, this will also load the other settings from the ini files//
	// kind of particle 0:atom 1:molecule 2:electron 3:Ion 4:Ion(simple) 

	particles.Add("Ion",1,1,0);//------------particle 0 --- Do not comment out!


	if (whichParticles=="") return;
	else if(whichParticles=="CH3I") 
	{//---SACLA CH3I molecule
		AddCH3I(particles);
	}
	else if(whichParticles=="IUracil") 
	{//---SACLA I-Uracil
		AddIUracil(particles);
	}
	else
	{
		std::cout << "can not find particles!!" << std::endl;
		return;
	}

	//---SACLA Ar atom
	//AddArgon(particles);
	//---SACLA Xe atom
	//AddXenon(particles);
	//AddXenon132(particles);
	//---SACLA CH3I molecule
	//AddCH3I(particles);
	//AddIUracil(particles);
	//---Test N2 molecule
	//AddNitrogen(particles);
}
//_______SACLA 2012A______________________________________________________________________________
void TofCorrection(MyDetektorHit &dh, const double alpha, const double k2, const double k4, const double xc, const double t0)
{
	const double t = dh.Time()-t0;
	const double x = dh.X()-xc;
	dh.SetTof( t - alpha*t*(k2*x*x+k4*x*x*x*x) );
}
//__________Show mass &ToF spactrum with Particle name_________________________________________
void MyAnalyzer::ShowResult()
{
	canv = new TCanvas("Result","Result",100,100,1000+4,600+28);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	canv->Divide(1,2,0.002,0.002);

	TH1D* massHisto = dynamic_cast<TH1D*>(gFile->GetDirectory("/Ion")->FindObject("Mass"));
	canv->cd(1);
	massHisto->Draw("L");
	txtMass.assign(fParticles.GetNbrOfParticles(),TText());
	for (int i = 1; i < fParticles.GetNbrOfParticles(); i++)
	{
		txtMass[i].SetTextSize(0.04);
		txtMass[i].SetTextAlign(21);
		txtMass[i].SetTextColor(kBlack);
		double massPerQ = fParticles.GetParticle(i).GetMass_au()*MyUnitsConv::au2amu()/fParticles.GetParticle(i).GetCharge_au();
		double height = massHisto->GetBinContent((massHisto->FindBin(massPerQ)));
		txtMass[i].DrawText(massPerQ, height*1.8, fParticles.GetParticle(i).GetName());
	}

	TH1D* tofHisto = dynamic_cast<TH1D*>(gFile->GetDirectory("/Ion")->FindObject("Tof"));
	canv->cd(2);
	tofHisto->Draw("L");
	txtTof.assign(fParticles.GetNbrOfParticles(),TText());
	boxTof.assign(fParticles.GetNbrOfParticles(),TBox());
	for (int i = 1; i < fParticles.GetNbrOfParticles(); i++)
	{
		txtTof[i].SetTextSize(0.04);
		txtTof[i].SetTextAlign(21);
		txtTof[i].SetTextColor(kBlack);
		double t0 = fParticles.GetParticle(0).GetT0();
		double tofPos = calcTof(fParticles.GetParticle(i),fParticles.GetParticle(0)) + t0;
		double height = tofHisto->GetBinContent((tofHisto->FindBin(tofPos)));
		txtTof[i].DrawText(tofPos, height*1.8, fParticles.GetParticle(i).GetName());
		boxTof[i].SetFillColor(kRed);
		boxTof[i].SetFillStyle(3001);
		boxTof[i].DrawBox(fParticles.GetParticle(i).GetCondTofFr(), 0,  fParticles.GetParticle(i).GetCondTofTo(), height *1.6);
	}
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
	if ((intFileName != "")&&(existIntensityData))
	{
		std::map<unsigned int, double>::iterator itTagInt;//BM1
		std::map<unsigned int, double>::iterator itTagInt2;//PD
		itTagInt = tagIntensity.find(TagNumber);
		itTagInt2 = tagIntensity2.find(TagNumber);

		if ((itTagInt == tagIntensity.end())||(itTagInt2 == tagIntensity2.end()))
		{
			std::cout <<"\r"<<TagNumber<< " is not found!!";
			missedTagCount++;
			return;
		}

		fIntensities.push_back(itTagInt2->second * 10000);//[0]
		fIntensities.push_back(itTagInt->second * 10e+9);//[1]

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

		//-----------------------------------
		if (fIntensities.size()) fHi.fill(startIdx,"IntensityBM1",fIntensities[1],"[arb. unit]",1000,0,1000);
		startIdx++;
		if (fIntensities.size()) fHi.fill(startIdx,"IntensityPD",fIntensities[0],"[arb. unit]",1000,0,1000);
		startIdx++;

		//---skip when FEL is stopped
		if (fIntensities.size()){
			if (fIntensities[0]<1) 
			{
				//std::cout << "skip this event" << std::endl;
				return;
			}
		}
		if (fIntensities.size()) fHi.fill(startIdx,"IntensitySelectPD",fIntensities[0],"[arb. unit]",1000,0,500);
		startIdx++;

		if (existIntPartition)
		{
			//---for getting power dep
			if (fIntensities.size() && (intPartition.size()>1)) fHi.fill(startIdx,"IntensityPDForPowDep",fIntensities[0],"[arb. unit]",intPartition.size()-1,&intPartition.front(),"PowerDependence");
			startIdx++;
			int whichRegion = -1;
			for (size_t k=0; k<intPartition.size()-1; k++)
			{
				if ((intPartition[k]<fIntensities[0])&&(fIntensities[0]<intPartition[k+1]))
				{
					whichRegion = k;
					fHi.fill(startIdx+k,Form("Intensity%02d",k),fIntensities[0],"[arb. unit]",1000,0,500,"PowerDependence");
				}
			}
			startIdx += intPartition.size();
		}
	}

	fHi.fill(startIdx,"NumberOfHits",rd.GetNbrOfHits(),"Number of Hits",100,0,100);
	startIdx++;


	//---get the raw mcp events called times//
	const MySignalAnalyzedChannel &sac = fSAE.GetChannel(7-1);
	for (size_t i=0; i<sac.GetNbrPeaks();++i)
	{
		const MyPeak &p = sac.GetPeak(i);
		fHi.fill(startIdx,"Times_MCP",p.GetTime(),"tof [ns]",5000,0,fOE.GetNbrSamples()*fOE.GetSampleInterval()*1e9);
	}
	startIdx++;

	if(fSAE.GetChannel(7-1).GetNbrPeaks())
	{
		const MyPeak &p_t0 = fSAE.GetChannel(7-1).GetPeak(0);
		fHi.fill(startIdx,"Times_MCP_first",p_t0.GetTime(),"tof [ns]",10000,0,2000);
	}
	startIdx++;

	//-----------------------Set jet speed
	//for (size_t i=0;i<fParticles.GetNbrOfParticles();++i)
	//{
	//	MyParticle &p = fParticles.GetParticle(i);
	//	p.SetXVelocity(XVelocity);
	//	p.SetYVelocity(YVelocity);
	//}

	int secondStartIdx=startIdx+20;
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
		//fHi.fill(startIdx+6,"XPosVsTofFine",dh.Tof(),dh.X(),"tof [ns]","x [mm]",static_cast<int>(maxTof/5),0,maxTof,300,-maxPos,maxPos);
		//fHi.fill(startIdx+7,"YPosVsTofFine",dh.Tof(),dh.Y(),"tof [ns]","y [mm]",static_cast<int>(maxTof/5),0,maxTof,300,-maxPos,maxPos);

		
		//-----Particle(0)---Ion---
		//get the particle from the vector//
		MyParticle &p = fParticles.GetParticle(0);
		//if this hit fits both conditions then add the Hit to this Particle and fill the histo for this hit//
		if (p.CheckTofAndPos(dh))
		//select hit by reconstruction method//
		if (dh.RekMeth() < rekmeth)//added by motomura
		{
			const MyParticleHit &ph = p.AddHit(dh);

			fHi.fill(startIdx+8,"Mass",ph.Mass(),"Mass/q",10000,0,200,"Ion");
			fHi.fill(startIdx+9,"TofCor",ph.TofCor(),"tof [ns]",20000,p.GetCondTofFr()-p.GetT0()-p.GetCondTofRange()*0.1,p.GetCondTofTo()-p.GetT0()+p.GetCondTofRange()*0.1,"Ion");
			fHi.fill(startIdx+10,"Det",ph.X(),ph.Y(),"x [mm]","y [mm]",300,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,300,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,"Ion");
			fHi.fill(startIdx+11,"Tof",ph.Tof(),"tof [ns]",20000,0,maxTof,"Ion");
			fHi.fill(startIdx+12,"XPosVsTof",ph.Tof(),ph.X(),"tof [ns]","x [mm]",5000,p.GetCondTofFr()-p.GetCondTofRange()*0.1,p.GetCondTofTo()+p.GetCondTofRange()*0.1,300,p.GetCondRadX()-p.GetCondRad()*1.1,p.GetCondRadX()+p.GetCondRad()*1.1,"Ion");
			fHi.fill(startIdx+13,"YPosVsTof",ph.Tof(),ph.Y(),"tof [ns]","y [mm]",5000,p.GetCondTofFr()-p.GetCondTofRange()*0.1,p.GetCondTofTo()+p.GetCondTofRange()*0.1,300,p.GetCondRadX()-p.GetCondRad()*1.1,p.GetCondRadX()+p.GetCondRad()*1.1,"Ion");
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
			if (p.GetKindParticle() < 0) continue;
			//if this hit fits both conditions then add the Hit to this Particle and fill the histo for this hit//
			if (p.CheckTofAndPos(dh))
			//select hit by reconstruction method//
			if (dh.RekMeth() < rekmeth)
			{
				const MyParticleHit &ph = p.AddHit(dh);
				fillParticleHistograms(p,ph,fIntensities,fHi,secondStartIdx,intPartition);
			}
			//we reserve 100 histograms for one particle//
			secondStartIdx +=100;
			//std::cout <<j<<" "<< secondStartIdx<<std::endl;
		}
	}//------------------------------------------------------------------------------------------------------//
	startIdx += (fParticles.GetNbrOfParticles()*100);

	//fill histograms of hit count
	for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)//particle index j=0 is "ion" 
	{
		fHi.fill(startIdx+j,"NumberOfHits",fParticles.GetParticle(j).GetNbrOfParticleHits(),"Number of Hits",100,0,100,Form("%s",fParticles.GetParticle(j).GetName()));

		for (size_t k=0;k<fParticles.GetParticle(j).GetNbrOfParticleHits();++k)
			fHi.fill(startIdx+fParticles.GetNbrOfParticles(),"NumberOfParticleHits",j,"Particle Number",fParticles.GetNbrOfParticles()+1,0,fParticles.GetNbrOfParticles()+1);
	}

	startIdx += (fParticles.GetNbrOfParticles())+1;

	//now you have found the particles//
	//we can look for coincidences//
	//get coincidence by calcurating aligned momentum-sum
	if (MoleculeAnalysis == 1)
	{
		//clear Coincidence counter
		for (size_t i=0;i<fParticles.GetNbrOfParticles();++i)
			for (size_t j=0;j<fParticles.GetNbrOfParticles();++j)
				molecule[i][j].CoincidenceCount = 0;
		//loop from particle No.1 (except Ion)
		for (size_t i=1;i<fParticles.GetNbrOfParticles();++i)
		{
			const MyParticle &ip = fParticles.GetParticle(i);
			for (size_t j=i;j<fParticles.GetNbrOfParticles();++j)
			{
				const MyParticle &jp = fParticles.GetParticle(j);
				//check if ip and jp is molecule -> particles.Add()
				if ((ip.GetKindParticle() == 1)&&(jp.GetKindParticle() == 1))
				{
					//
					if (
						((ip.GetCoinGroup() != jp.GetCoinGroup()))
						||((ip.GetCoinGroup()==100)&&(jp.GetCoinGroup()==100))
						)
					{
						fillMoleculeHistogram(ip,jp,fIntensities,fHi,startIdx, molecule[i][j], intPartition);
						startIdx += 200;
					}
				}
			}
		}
		//--Fill Number of coincidence
		for (size_t i=1;i<fParticles.GetNbrOfParticles();++i)//from particle No.1
		{
			for (size_t j=i;j<fParticles.GetNbrOfParticles();++j)
			{
				for (size_t k=0; k<molecule[i][j].CoincidenceCount; ++k)
				{
					//Coincidence charge distribution
					fHi.fill(startIdx,"NumberOfCoincidence",i,j,"Particle number","Particle number",
						fParticles.GetNbrOfParticles(),0,fParticles.GetNbrOfParticles(),
						fParticles.GetNbrOfParticles(),0,fParticles.GetNbrOfParticles());
					fHi.fill(startIdx+1,"CoincidentChargeState",fParticles.GetParticle(i).GetCharge_au(),fParticles.GetParticle(j).GetCharge_au(),
						"Carbon Chage","Iodine Chage",
						fParticles.GetNbrOfParticles(),0,fParticles.GetNbrOfParticles(),
						fParticles.GetNbrOfParticles(),0,fParticles.GetNbrOfParticles());
					//sumed chage state
					fHi.fill(startIdx+2, "SumOfChageState", fParticles.GetParticle(i).GetCharge_au()+fParticles.GetParticle(j).GetCharge_au() , 
						"Charge State",	fParticles.GetNbrOfParticles()*2, 0, fParticles.GetNbrOfParticles()*2);
				}
			}
		}
		startIdx += 3;
	}


	//---gate by certain particle hits
	//reserve ID for fillMoleculeHistogram2
	const int nbrOfHistosInfillMol2 = 15;
	secondStartIdx = startIdx + (fParticles.GetNbrOfParticles()+fParticles.GetNbrOfParticles()*fParticles.GetNbrOfParticles())*nbrOfHistosInfillMol2;
	if (MoleculeAnalysis == 2)
	{
		//loop from particle 1 because excepting Ion
		for (size_t i=1;i<fParticles.GetNbrOfParticles();++i)
		{
			const MyParticle &ip = fParticles.GetParticle(i);
			//check whether ip have counts and gatting paticle (I+, I++, ...)
			if ((ip.GetNbrOfParticleHits())&&(ip.GetCoinGroup()==1))
			{
				//fill Ion spectra 
				fillSpectra(ip,fParticles.GetParticle(0),fHi, secondStartIdx+ (i*10));
				//loop to fill the target particles
				for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)
				{
					const MyParticle &jp = fParticles.GetParticle(j);
					//Target particle (C+.N+..O+..CO+..,H+)
					if ((jp.GetKindParticle() == 1)&&(jp.GetCoinGroup()==0))
					{
						fillMoleculeHistogram2(ip,jp,fIntensities,fHi,startIdx+ (i+j*fParticles.GetNbrOfParticles())*nbrOfHistosInfillMol2);
					}
				}
			}
		}
	}

	//Skip already used ID
	startIdx += (fParticles.GetNbrOfParticles()+fParticles.GetNbrOfParticles()*fParticles.GetNbrOfParticles())*nbrOfHistosInfillMol2;
	startIdx += fParticles.GetNbrOfParticles()*10;

	//---Post-analysis---//
	//if (molecule[5][12].CoincidenceCount > 0) 
	//for (size_t i=0; i<rd.GetNbrOfHits();++i)
	//{
	//	MyDetektorHit &dh = rd.GetHit(i);
	//	//the tof is just the timing of the mcp signal//
	//	dh.SetTof(dh.Time());
	//	secondStartIdx=startIdx;
	//	for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)
	//	{
	//		//get the particle from the vector//
	//		MyParticle &p = fParticles.GetParticle(j);
	//		if (p.GetKindParticle() == -1)
	//		if (p.CheckTofAndPos(dh))
	//		//select hit by reconstruction method//
	//		if (dh.RekMeth() < rekmeth)
	//		{
	//			const MyParticleHit &ph = p.AddHit(dh);
	//			fillParticleHistograms(p,ph,fIntensities,fHi,secondStartIdx);
	//		}
	//		secondStartIdx +=100;
	//	}
	//}
	//if (MoleculeAnalysis) fillPIPICO(fParticles.GetParticle(0),fHi);

	//std::cout << startIdx << std::endl;
}
//----------------------------------extraCondition----------------------------------------------------------------//
bool PosCondition(const MyDetektorHit &dh)
{
	//Except residual gas
	//if ((dh.X() < -22) || (dh.X() > 21) || (dh.Y() < -6) || (dh.Y() > 9.5)) return true;
	//if ((dh.X() < -23) || (dh.X() > 22) || (dh.Y() < -7) || (dh.Y() > 10.5)) return true;
	return false;
}
bool TofPosCondition(const MyDetektorHit &dh)
{
	//Mass15 N2_2009_Jan
	//if (((dh.Tof() > 5600)&&(dh.Tof() < 5660)) && ((dh.Y() > -7.) && (dh.Y() < 8.)) && ((dh.X() > -20) && (dh.X() < 20))) return false;
	//if (((dh.Tof() > 5460)&&(dh.Tof() < 5473)) && ((dh.Y() > -4.) && (dh.Y() < 2.)) && ((dh.X() > -5) && (dh.X() < 4))) return false;

	return true;
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
