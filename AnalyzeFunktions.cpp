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
#include <float.h>
#include <TTree.h>
#include <TFile.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TGClient.h>
#include <TStyle.h> 
#include <TDirectory.h>

//#define XVelocity 0.0003704 //(mm/ns)
//#define YVelocity -0.00037775
#define XVelocity 0. //(mm/ns)
#define YVelocity 0.

//_____________________________functions______________________________________________________________________________________________________________________________
//momentum calculation
double calcPx(const MyParticle &p, const MyParticleHit &ph)
{
	//return MyMomentaCalculator::px(ph.XCorRotScl(),ph.YCorRotScl(),ph.TofCor(),p.GetMass_au(),p.GetCharge_au(),p.GetSpectrometer());
	//---vmi model
	return MyMomentaCalculator::pr_VMI(p.GetMass_au(), ph.XCorRotScl(), p.GetSpectrometer());

}
double calcPy(const MyParticle &p, const MyParticleHit &ph)
{
	//return MyMomentaCalculator::py(ph.XCorRotScl(),ph.YCorRotScl(),ph.TofCor(),p.GetMass_au(),p.GetCharge_au(),p.GetSpectrometer());
	//---vmi model
	return MyMomentaCalculator::pr_VMI(p.GetMass_au(), ph.YCorRotScl(), p.GetSpectrometer());
}
double calcPz(const MyParticle &p, const MyParticleHit &ph)
{
	//---coltrims
	//return MyMomentaCalculator::pz(ph.TofCor(),p.GetMass_au(),p.GetCharge_au(),p.GetSpectrometer());
	//---simple poly
	//return MyMomentaCalculator::pz_poly(ph.TofCor(),p.GetSpectrometer());
	//---Polynomial with R
	return MyMomentaCalculator::pz_polyRT(ph.TofCor(), ph.R(), p.GetSpectrometer());
	//---Polynomial with R another version
	//return MyMomentaCalculator::pz_polyRT_Another(ph.TofCor(), ph.R(), p.GetSpectrometer());
}
double calcPr(const MyParticle &p, const MyParticleHit &ph)
{
	//return 0.0;
	//---Polynomial with R and T
	return MyMomentaCalculator::pr_polyRT(ph.TofCor(), ph.R(), p.GetSpectrometer());
	//---vmi model
	//return MyMomentaCalculator::pr_VMI(p.GetMass_au(), ph.R(), p.GetSpectrometer());
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
	else if(whichParticles=="Argon") 
	{
		//---SACLA Ar atom
		AddArgon(particles);
	}
	else if(whichParticles=="XeCluster") 
	{
		//---SACLA Xe cluster
		AddXeCluster(particles);
	}
	else if(whichParticles=="ArCluster") 
	{
		//---SACLA Ar cluster
		AddArCluster(particles);
	}
	else if(whichParticles=="XeIsotope") 
	{
		//---SACLA Xe Isotope
		AddXeIsotope(particles);
	}
	else if(whichParticles=="KrCluster") 
	{
		//---SACLA Kr cluster
		AddKrCluster(particles);
	}
	else if(whichParticles=="KrArCluster") 
	{
		//---SACLA KrAr cluster
		AddKrArCluster(particles);
	}
	else if(whichParticles=="KrIsotope") 
	{
		//---SACLA Kr Isotope
		AddKrIsotope(particles);
	}
	else
	{
		std::cout << "can not find particles!!" << std::endl;
		return;
	}

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

//________Analysis main loop___________________________________________________________________________________________________________________________________________________
void MyAnalyzer::Analyze(MyWaveform &wf)
{
	MyDetektor &rd = fSE.GetDetektor(0);
	int startIdx=10;
	//for Trend Histogram
	const int skipCounter = static_cast<int>(fEntryIterator/trendStep);
	
	//clear the containers//
	fIntensities.clear();
	fParticles.ClearParticles();
	//Get Tag number on this event
	const unsigned int TagNumber = fOE.GetEventID();
	if (optShutMode == 1)
	{
		if(static_cast<int>(DB.GetStatusAndData(TagNumber,2).second+0.1) != 0) 
		{
			fHi.SetMainDir("OptLaserOn");
			startIdx=10;
		}
		if(static_cast<int>(DB.GetStatusAndData(TagNumber,2).second+0.1) == 0) 
		{
			fHi.SetMainDir("OptLaserOff");
			startIdx=5010;
		}
	}
	//select tag number
	//if (TagNumber < 489913416) return;

	//Position information
	//std::map<unsigned int, double>::iterator itPosX;//X
	//std::map<unsigned int, double>::iterator itPosY;//Y
	//itPosX = beamPosX.find(TagNumber);
	//itPosY = beamPosY.find(TagNumber);
	//if ((itPosX == beamPosX.end())||(itPosY == beamPosY.end()))
	//{
	//	std::cout <<"\r"<<TagNumber<< " is not found!!";
	//	missedTagCount++;
	//	return;
	//}
	//const double beamPositionX = itPosX->second;
	//const double beamPositionY = itPosY->second;

	//Get intensity from loaded data map
	if ((zeroDTxtFileName != "")&&(method0D_Data!=0))
	{
		if (method0D_Data==1)
		{
			std::map<unsigned int, double>::iterator itTagDelay;//Delay
			std::map<unsigned int, double>::iterator itTagInt;//PD
			itTagDelay = tagDelay.find(TagNumber);
			itTagInt = tagIntensity.find(TagNumber);

			if ((itTagDelay == tagDelay.end())||(itTagInt == tagIntensity.end()))
			{
				//std::cout <<"\r"<<TagNumber<< " is not found!!";
				missedTagCount++;
				return;
			}

			fIntensities.push_back((factorPMDOffset - itTagDelay->second) / factorPMD);//[0] Delay	BM1:24486*1000000
			fIntensities.push_back(itTagInt->second * factorPD);//[1] PD:
		}

		if (method0D_Data==2)
		{
			
			if(DB.GetStatusAndData(TagNumber,0).first==0)
			{
				//std::cout <<"\r"<<TagNumber<< " is not found!!";
				missedTagCount++;
				return;
			}

			if( _isnan(DB.GetStatusAndData(TagNumber,0).second))
			{
				std::cout <<"Delay of "<<TagNumber<< " is NaN " << std::endl;
			}

			if( _isnan(DB.GetStatusAndData(TagNumber,1).second))
			{
				std::cout <<"Intensity of "<< TagNumber << " is NaN " << std::endl;
			}

			fIntensities.push_back((factorPMDOffset - DB.GetStatusAndData(TagNumber,0).second) / factorPMD);//[0] Delay	BM1:24486*1000000
			fIntensities.push_back(DB.GetStatusAndData(TagNumber,1).second * factorPD);//[1] PD:
		}

		if (method0D_Data==3)
		{
			if(DB.GetStatusAndData(TagNumber,0).first==0)
			{
				//std::cout <<"\r"<<TagNumber<< " is not found!!";
				missedTagCount++;
				return;
			}

			// Delay
			if( _isnan(DB.GetStatusAndData(TagNumber,0).second))
			{
				std::cout <<"Delay of "<<TagNumber<< " is NaN " << std::endl;
				missedTagCount++;
				return;
			}
			
			if( _isnan(DB.GetStatusAndData(TagNumber,1).second))
			{
				std::cout <<"Timing_Valid  of "<<TagNumber<< " is NaN " << std::endl;
				fIntensities.push_back((factorPMDOffset - DB.GetStatusAndData(TagNumber,0).second) / factorPMD);
			}
			
			if(!DB.GetStatusAndData(TagNumber,1).second)
			{
				fIntensities.push_back((factorPMDOffset - DB.GetStatusAndData(TagNumber,0).second) / factorPMD);//[0] Delay	BM1:24486*1000000
			}

			else
			{
				if( _isnan(DB.GetStatusAndData(TagNumber,2).second))
				{
					std::cout <<"Timing_Jitter  of "<<TagNumber<< " is NaN " << std::endl;
				}
				if( _isnan(DB.GetStatusAndData(TagNumber,3).second))
				{
					std::cout <<"xfel_bl_3_st_1_motor_73/position  of "<<TagNumber<< " is NaN " << std::endl;
				}
				else
				{
					fIntensities.push_back((factorPMDOffset - DB.GetStatusAndData(TagNumber,0).second) / factorPMD + DB.GetStatusAndData(TagNumber,2).second);
				}
			}
		}

		// FEL intensity
		if( _isnan(DB.GetStatusAndData(TagNumber,4).second))
		{
			std::cout <<"Intensity of "<< TagNumber << " is NaN " << std::endl;
			missedTagCount++;
			return;
		}
		fIntensities.push_back(DB.GetStatusAndData(TagNumber,4).second * factorPD);//[1] PD:
	}

		////------------------------------------------SACLA 2012A
		////PhotoDiode intensity
		//if ( itTagInt->first < 156088080 )
		//	//before changing Al plate
		//	//fIntensities.push_back((tagIntensity.find(TagNumber)->second)*1000*36.15);
		//	fIntensities.push_back((tagIntensity.find(TagNumber)->second)/1.926E-5);//5.5keV PD->microJ(weak setting)
		//else
		//	fIntensities.push_back((tagIntensity.find(TagNumber)->second)/6.963E-4);//5.5keV PD->microJ
		//	//fIntensities.push_back((tagIntensity.find(TagNumber)->second)/1.17E-3);//5.0keV PD->microJ

		////FEL is attenuated, set tag region
		//const unsigned int tagFrom = 167982762;	//att0p4(5.5keV):159612312		Al25micro(5.0keV):167982762
		//const unsigned int tagTo = 170567544;	//att0p4(5.5keV):161807544		Al25micro(5.0keV):170567544
		//const double attFactor = 0.381;//5.5keV:0.381  5.5keV:0.52(test)  5.0keV:0.271
		//if (( tagFrom <= itTagInt->first ) && ( itTagInt->first <= tagTo ))
		//	fIntensities.push_back(itTagInt->second * attFactor);
		//else 
		//	fIntensities.push_back(itTagInt->second);
		////----------------------------------------------------

		//for (size_t i=0;i<fIntensities.size();++i)
		//{
		//	fHi.fill(startIdx++,Form("Int%d",i),fIntensities[i],Form("Int%d",i),1000,0,1000,"Intensity");
		//	for (size_t j=i+1; j<fIntensities.size();++j)
		//	{
		//		fHi.fill(startIdx++,Form("IntDep(%d)(%d)",i,j),fIntensities[i], fIntensities[j],Form("Int %d",i),Form("Int %d",j),1000,0,1000,1000,0,1000,"Intensity");
		//	}
		//}
		//if (fIntensities.size() > 1)
		//	if ( fIntensities[0]*fIntensities[1]>0)
		//		fHi.fill(startIdx+1,"IntensityDivide0by1",fIntensities[0]/fIntensities[1],"IntensityDivide0by1",1000,0,5,"Intensity");
		//startIdx+=2;

		////-----------------------------------
		fHi.fill(startIdx,"IntensityFEL",fIntensities[1],"[arb. unit]",1000,0,1000);
		startIdx++;
		fHi.fill(startIdx,"Delay",fIntensities[0],"[arb. unit]",1000,-50,50);
		startIdx++;
		////-----Trend plot BM1 & PD intensity
		//fHi.fill(startIdx,"TrendIntensityFEL",skipCounter,Form("[shots/%d]",trendStep),1000,0,1000,"Trend",fIntensities[1]/skipCounter);
		//startIdx++;
		//fHi.fill(startIdx,"TrendDelay",skipCounter,Form("[shots/%d]",trendStep),1000,0,1000,"Trend",fIntensities[0]/skipCounter);
		//startIdx++;
		//fHi.fill(startIdx,"TrendNumberOfHits",skipCounter,Form("[shots/%d]",trendStep),1000,0,1000,"Trend",rd.GetNbrOfHits()/skipCounter);
		//startIdx++;
		//fHi.fill(startIdx,"TrendBeamPosX",skipCounter,Form("[shots/%d]",trendStep),1000,0,1000,"Trend",beamPositionX);
		//startIdx++;
		//fHi.fill(startIdx,"TrendBeamPosY",skipCounter,Form("[shots/%d]",trendStep),1000,0,1000,"Trend",beamPositionY);
		//startIdx++;
		
		//--skip this shot event if FEL is below 5 (FEL is stopped)
		if (fIntensities[1]< 5) return;

		if (selectIntensity)
		{
			//---skip this shot event if FEL is below the lower limit
			if (fIntensities[1]<intensityLowerLimit) return;
			//---skip this shot event if FEL is over the upper limit
			if (fIntensities[1]>intensityUpperLimit) return;
		}
		if (fIntensities.size()) fHi.fill(startIdx,"IntensitySelectPD",fIntensities[1],"[arb. unit]",1000,0,1000);
		startIdx++;

		//if (existIntPartition)
		//{
		//	//---for getting power dep
		//	if (fIntensities.size() && (intPartition.size()>1)) fHi.fill(startIdx,"IntensityPDForPowDep",fIntensities[0],"[arb. unit]",intPartition.size()-1,&intPartition.front(),"PowerDependence");
		//	startIdx++;
		//	int whichRegion = -1;
		//	for (size_t k=0; k<intPartition.size()-1; k++)
		//	{
		//		if ((intPartition[k]<fIntensities[0])&&(fIntensities[0]<intPartition[k+1]))
		//		{
		//			whichRegion = k;
		//			fHi.fill(startIdx+k,Form("Intensity%02d",k),fIntensities[0],"[arb. unit]",1000,0,500,"PowerDependence");
		//		}
		//	}
		//	startIdx += intPartition.size();
		//}

	fHi.fill(startIdx,"NumberOfHits",rd.GetNbrOfHits(),"Number of Hits",100,0,100);
	startIdx++;

	//Analyze MCP intensity

	//fillMCPToFHistograms(fOE, fHi, mcpTofRegion);

	//MCP intensity
	//Ar1p
	//Ar2p
	//Ar3p
	//	double McpIntensityAr1p = Average(fOE.GetChannel(7-1),    7000,7560,true);//Voltage D2:4000
	//	double McpIntensityAr2p = Average(fOE.GetChannel(7-1),    5250,5680,true);//Voltage D2:4000
	//	double McpIntensityAr3p = Average(fOE.GetChannel(7-1),    4400,4740,true);//Voltage D2:4000 // Ar is not yet

	//Xe1p
	//Xe2p
	//Xe3p
	//const double McpIntensityXe1p = Average(fOE.GetChannel(7-1),7000,7600,true);//Voltage D2:4000
	//const double McpIntensityXe2p = Average(fOE.GetChannel(7-1),5250,5680,true);//Voltage D2:4000
	//const double McpIntensityXe3p = Average(fOE.GetChannel(7-1),4400,4740,true);//Voltage D2:4000

	//double McpIntensityXe1p = Average(fOE.GetChannel(7-1),5200,5800,true);//Voltage D2:8000
	//double McpIntensityXe2p = Average(fOE.GetChannel(7-1),4000,4280,true);//Voltage D2:8000
	//double McpIntensityXe3p = Average(fOE.GetChannel(7-1),3400,3550,true);//Voltage D2:8000

	//const double McpIntensityXe1p = Average(fOE.GetChannel(7-1),14127,17150,true);//Voltage D2:800
	//const double McpIntensityXe1pH = Average(fOE.GetChannel(7-1),14000,16200,true);//Voltage D2:800<----High energy region
	//const double McpIntensityXe1pL = Average(fOE.GetChannel(7-1),16329,16554,true);//Voltage D2:800<----Low energy region
	//const double McpIntensityXe2p = Average(fOE.GetChannel(7-1),10254,12300,true);//Voltage D2:800
	//const double McpIntensityXe3p = Average(fOE.GetChannel(7-1),16329,16554,true);//Voltage D2:800

	//Argon D2:800V setting
	//const double McpIntensityXe1p = Average(fOE.GetChannel(7-1),8600,9600,true);//Voltage D2:800
	//const double McpIntensityXe2p = Average(fOE.GetChannel(7-1),6450,6900,true);//Voltage D2:800
	//const double McpIntensityXe3p = Average(fOE.GetChannel(7-1),5400,5600,true);//Voltage D2:800
	//----------------------------------
	//Delay dependence SACLA2014A
	//----------------------------------

	fHi.fill(startIdx+0,"DelayVsShots",fIntensities[0],"Delay [ps]", delayBins, delayFrom, delayTo);
	fHi.fill(startIdx+1,"DelayVsXFELintensity",fIntensities[0],fIntensities[1],"Delay [ps]","XFEL intensity [arb .unit]" ,delayBins, delayFrom, delayTo, 1000 ,0 , 1000);
	//fHi.fill(startIdx+1,"DelayVsXFELintensity",fIntensities[0],fIntensities[1],"Delay [ps]","XFEL intensity [arb .unit]" ,32 , -5, 11, 30, 0, 600);
	//fHi.fill(startIdx+1,"DelayDependenceXe1p2D",fIntensities[0], McpIntensityXe1pH,"Delay [ps]","MCPintensity",200,-50,50,500,0,500);
	//fHi.fill(startIdx+2,"DelayDependenceXe2p2D",fIntensities[0], McpIntensityXe2pL,"Delay [ps]","MCPintensity",200,-50,50,500,0,500);
	//fHi.fill(startIdx+3,"DelayDependenceXe3p2D",fIntensities[0], McpIntensityXe3p,"Delay [ps]","MCPintensity",200,-50,50,500,0,1000);
	//fHi.fill(startIdx+4,"DelayDependenceXe1pH",fIntensities[0],"Delay [ps]", delayBins, delayFrom, delayTo,"",McpIntensityXe1pH);
	//fHi.fill(startIdx+5,"DelayDependenceXe1pL",fIntensities[0],"Delay [ps]", delayBins, delayFrom, delayTo,"",McpIntensityXe1pL);
	//fHi.fill(startIdx+6,"DelayDependenceXe2p",fIntensities[0],"Delay [ps]", delayBins, delayFrom, delayTo,"",McpIntensityXe2p);
	//fHi.fill(startIdx+7,"DelayDependenceXe3p",fIntensities[0],"Delay [ps]", delayBins, delayFrom, delayTo,"",McpIntensityXe3p);
	startIdx += 2;

	//fHi.fill(startIdx+0,"DelayVsShots_UltraWide",fIntensities[0],"Delay [ps]",1000,-50,550);
	//fHi.fill(startIdx+1,"DelayDependenceXe1p2D_UltraWide",fIntensities[0], McpIntensityXe1p,"Delay [ps]","MCPintensity",500,-50,550,500,0,500);
	//fHi.fill(startIdx+2,"DelayDependenceXe2p2D_UltraWide",fIntensities[0], McpIntensityXe2p,"Delay [ps]","MCPintensity",500,-50,550,500,0,500);
	//fHi.fill(startIdx+3,"DelayDependenceXe3p2D_UltraWide",fIntensities[0], McpIntensityXe3p,"Delay [ps]","MCPintensity",500,-50,550,500,0,1000);
	//fHi.fill(startIdx+4,"DelayDependenceXe1p_UltraWide",fIntensities[0],"Delay [ps]",1000,-50,550,"",McpIntensityXe1p);
	//fHi.fill(startIdx+5,"DelayDependenceXe2p_UltraWide",fIntensities[0],"Delay [ps]",1000,-50,550,"",McpIntensityXe2p);
	//fHi.fill(startIdx+6,"DelayDependenceXe3p_UltraWide",fIntensities[0],"Delay [ps]",1000,-50,550,"",McpIntensityXe3p);

	//for (size_t i = 0; i < fOE.GetNbrSamples(); i++)
	//	fHi.fill(startIdx+7,"DelayVsMCPSignal",i,fIntensities[0],"tof [ns]","Delay [ps]",3000,0,fOE.GetNbrSamples(), delayBins, delayFrom, delayTo,"",wf.GetWaveform()[i]);
	//startIdx++;

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
		fHi.fill(startIdx+2,"TofAll",dh.Tof(),"tof [ns]",30000,0,maxTof);
		fHi.fill(startIdx+3,"XPosVsTofAll",dh.Tof(),dh.X(),"tof [ns]","x [mm]",5000,0,maxTof,300,-maxPos,maxPos);
		fHi.fill(startIdx+4,"YPosVsTofAll",dh.Tof(),dh.Y(),"tof [ns]","y [mm]",5000,0,maxTof,300,-maxPos,maxPos);
		fHi.fill(startIdx+5,"TOF",dh.Tof(),"tof [ns]",4000,0,6000);
		fHi.fill(startIdx+6,"XTOF",dh.Tof(),dh.X(),"tof [ns]","x [mm]",4000,0,6000,300,-50,50);
		fHi.fill(startIdx+7,"YTOF",dh.Tof(),dh.Y(),"tof [ns]","y [mm]",4000,0,6000,300,-50,50);

		//-----Particle(0)---Ion---
		//get the particle from the vector//
		MyParticle &p = fParticles.GetParticle(0);
		//if this hit fits both conditions then add the Hit to this Particle and fill the histo for this hit//
		if (p.CheckTofAndPos(dh))
			//select hit by reconstruction method//
				if (dh.RekMeth() < rekmeth)//added by motomura
				{
					const MyParticleHit &ph = p.AddHit(dh);

					fHi.fill(startIdx+8,"Mass",ph.Mass(),"Mass/q",10000,0,500,"Ion");
					//fHi.fill(startIdx+9,"TofCor",ph.TofCor(),"tof [ns]",20000,p.GetCondTofFr()-p.GetT0()-p.GetCondTofRange()*0.1,p.GetCondTofTo()-p.GetT0()+p.GetCondTofRange()*0.1,"Ion");
					fHi.fill(startIdx+9,"TofCor",ph.TofCor(),"tof [ns]",10000,0,maxTof,"Ion");
					fHi.fill(startIdx+10,"Det",ph.X(),ph.Y(),"x [mm]","y [mm]",300,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,300,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,"Ion");
					fHi.fill(startIdx+11,"Tof",ph.Tof(),"tof [ns]",20000,0,maxTof,"Ion");
					fHi.fill(startIdx+12,"XPosVsTof",ph.Tof(),ph.X(),"tof [ns]","x [mm]",5000,p.GetCondTofFr()-p.GetCondTofRange()*0.1,p.GetCondTofTo()+p.GetCondTofRange()*0.1,300,p.GetCondRadX()-p.GetCondRad()*1.1,p.GetCondRadX()+p.GetCondRad()*1.1,"Ion");
					fHi.fill(startIdx+13,"YPosVsTof",ph.Tof(),ph.Y(),"tof [ns]","y [mm]",5000,p.GetCondTofFr()-p.GetCondTofRange()*0.1,p.GetCondTofTo()+p.GetCondTofRange()*0.1,300,p.GetCondRadX()-p.GetCondRad()*1.1,p.GetCondRadX()+p.GetCondRad()*1.1,"Ion");
					fHi.fill(startIdx+14,"DetCorScale",ph.XCorRotScl(),ph.YCorRotScl(),"x [mm]","y [mm]",300,p.GetCondRadX()-p.GetCondRad()*1.3,p.GetCondRadX()+p.GetCondRad()*1.3,300,p.GetCondRadY()-p.GetCondRad()*1.3,p.GetCondRadY()+p.GetCondRad()*1.3,"Ion");
					fHi.fill(startIdx+15,"DelayVsTOF",dh.Tof(),fIntensities[0],"tof [ns]","Delay [ps]",1000,0,maxTof, delayBins, delayFrom, delayTo,"Ion");
					fHi.fill(startIdx+16,"DelayVsTOFCor",ph.TofCor(),fIntensities[0],"tof [ns]","Delay [ps]",1000,0,maxTof, delayBins, delayFrom, delayTo,"Ion");
					//fHi.fill(startIdx+17,"FELIntensityVsTOFCor",ph.TofCor(),fIntensities[1],"tof [ns]","FEL intensity [arb. unit]",1000,0,maxTof, 1000, 0, 1000,"Ion");
					//fHi.fill(startIdx+18,"DelayVsXFELintensityVsTOFCor",ph.TofCor(),fIntensities[0],fIntensities[1],"tof [ns]","Delay [ps]","XFEL intensity [arb. unit]",1000,0,maxTof, delayBins/5, -5, delayTo, 60, 0, 600, "Ion");
					//fHi.fill(startIdx+18,"DelayVsXFELintensityVsTOFCor",ph.TofCor(),fIntensities[0],fIntensities[1],"tof [ns]","Delay [ps]","XFEL intensity [arb. unit]",1000,0,maxTof, 32, -5, 11, 30, 0, 600, "Ion");
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
								//Check the angle of phiZX. if ph is out of condition, it delete.
								if (p.CheckPhiZX(ph))
								{
									fillParticleHistograms(p,ph,fIntensities,fHi,secondStartIdx,intPartition,delayBins,delayFrom,delayTo,limitTheataZ);
								}
							}
							//we reserve 50 histograms for one particle//
							secondStartIdx +=50;
							//std::cout <<j<<" "<< secondStartIdx<<std::endl;
				}
	}//------------------------------------------------------------------------------------------------------//
	startIdx += (fParticles.GetNbrOfParticles()*50+20);
	//fill histograms of hit count
	int iodineCount = 0;
	for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)//particle index j=0 is "ion"  
	{
		fHi.fill(startIdx+j,"NumberOfHits",fParticles.GetParticle(j).GetNbrOfParticleHits(),"Number of Hits",100,0,100,Form("%s",fParticles.GetParticle(j).GetName()));
		fHi.fill(startIdx+fParticles.GetNbrOfParticles()+j,Form("TrendParticleHits%02d", j),skipCounter,Form("[shots/%d]",trendStep),1000,0,1000,"Trend",fParticles.GetParticle(j).GetNbrOfParticleHits());
		//Iodine rate
		if ((fParticles.GetParticle(j).GetMass_au()*MyUnitsConv::au2amu() > 120)&&(fParticles.GetParticle(j).GetCharge_au() < 10))
		{
			iodineCount += fParticles.GetParticle(j).GetNbrOfParticleHits();
		}
		for (size_t k=0;k<fParticles.GetParticle(j).GetNbrOfParticleHits();++k)
		{
			fHi.fill(startIdx+fParticles.GetNbrOfParticles(),"NumberOfParticleHits",j,"Particle Number",fParticles.GetNbrOfParticles()+1,0,fParticles.GetNbrOfParticles()+1);
			fHi.fill(startIdx+2*fParticles.GetNbrOfParticles(),"DelayVsNumberOfParticleHits",j,fIntensities[0],"Particle Number","delay [ps]",fParticles.GetNbrOfParticles()+1,0,fParticles.GetNbrOfParticles()+1,delayBins,delayFrom,delayTo);
		}
	}
	startIdx += (2*fParticles.GetNbrOfParticles()+1);
	fHi.fill(startIdx,"IodineHits",iodineCount,"Number of Hits",100,0,100);
	startIdx++;
	//now you have found the particles//
	//we can look for coincidences//
	//get coincidence by calcurating aligned momentum-sum
	if (MoleculeAnalysis == 1)
	{
		//clear Coincidence counter
		for (size_t i=0;i<fParticles.GetNbrOfParticles();++i)
			for (size_t j=0;j<fParticles.GetNbrOfParticles();++j)
			{
				molecule[i][j].CoincidenceCount = 0;
				molecule[i][j].CoinHitNbrC.clear();
				molecule[i][j].CoinHitNbrI.clear();
			}
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
							startIdx += 120;
							//for proton
							//extraCondition is 0:do not run proton routine, 1:run proton routine, 2:select only one hit event 
							if (extraCondition > 0)
							{
								if (molecule[i][j].CoincidenceCount > 0) 
								{
									int otherHits = 0;
									if (extraCondition == 2)
									{
										for (int k = 2; k < fParticles.GetNbrOfParticles(); k++)
										{
											if ( (k!=i)&&(k!=j) ) 
												otherHits += fParticles.GetParticle(k).GetNbrOfParticleHits();
										} 
									}
									if (otherHits == 0)
									{
										const MyParticle &hp = fParticles.GetParticle(1);
										//if (hp.GetNbrOfParticleHits() < 4)
										fillHydrogenHistogram(hp,ip,jp,fHi,startIdx,molecule[i][j]);
									}
								}
								startIdx += 65;
							}
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
							6,1,7,
							16,1,17);
						if (fParticles.GetParticle(i).GetCoinGroup()==0 && fParticles.GetParticle(j).GetCoinGroup()==1)
						{
							fHi.fill(startIdx+5+i, Form("CoincidentCarbon%dp",static_cast<int>(fParticles.GetParticle(i).GetCharge_au()+0.1)),
								fParticles.GetParticle(j).GetCharge_au(),
								"Iodine charge state", 20, 0, 20);
						}
						//sumed chage state
						fHi.fill(startIdx+2, "SumOfChargeState", fParticles.GetParticle(i).GetCharge_au()+fParticles.GetParticle(j).GetCharge_au() , 
							"Charge State",	fParticles.GetNbrOfParticles()*2, 0, fParticles.GetNbrOfParticles()*2);
					}
				}
			}
			startIdx += (5 + fParticles.GetNbrOfParticles());
	}
	//---gate by certain particle hits
	//reserve ID for fillMoleculeHistogram2
	const int nbrOfHistosInfillMol2 = 27;//
	secondStartIdx = startIdx + (fParticles.GetNbrOfParticles()+fParticles.GetNbrOfParticles()*fParticles.GetNbrOfParticles())*nbrOfHistosInfillMol2;
	if (MoleculeAnalysis == 2)
	{
		//loop from particle 1 for excepting Ion
		for (size_t i=1;i<fParticles.GetNbrOfParticles();++i)
		{
			const MyParticle &ip = fParticles.GetParticle(i);
			//check whether ip have counts and gatting paticle (I+, I++, ...)
			if ((ip.GetNbrOfParticleHits())&&(ip.GetCoinGroup()==1 || ip.GetCoinGroup()==2 ))
			{
				int hitCounter = 0;
				for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)
				{
					if ((j != i) && (fParticles.GetParticle(j).GetCoinGroup() == 1))
						hitCounter += fParticles.GetParticle(j).GetNbrOfParticleHits();
				}
				if (!hitCounter)
				{
					//fill Ion spectra 
					fillSpectra(ip,fParticles.GetParticle(0),fHi, secondStartIdx+ (i*10));
					//loop to fill the target particles
					for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)
					{
						const MyParticle &jp = fParticles.GetParticle(j);
						//Target particle (C+.N+..O+..CO+..,H+)
						if ((jp.GetKindParticle() == 1)&&(jp.GetCoinGroup()==0 || jp.GetCoinGroup()==2))
						{
							int tmpIdx = startIdx + (i+j*fParticles.GetNbrOfParticles())*nbrOfHistosInfillMol2;
							fillMoleculeHistogram2(ip,jp,fIntensities,fHi,tmpIdx);
							fillAngleHistogram(ip, jp, fHi, tmpIdx + 30);	
						}
					}

				}			
			}
		}
		//Skip already used ID
		startIdx += (fParticles.GetNbrOfParticles()+fParticles.GetNbrOfParticles()*fParticles.GetNbrOfParticles())*nbrOfHistosInfillMol2;
		startIdx += fParticles.GetNbrOfParticles()*10;

		//_______________________________//
		//___3-body angular correration__//
		//_______________________________//
		//1st loop from particle 1 for excepting Ion
		for (size_t i=1;i<fParticles.GetNbrOfParticles();++i)
		{
			//2nd loop
			for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)
			{
				//3rd loop
				for (size_t k=1;k<fParticles.GetNbrOfParticles();++k)
				{
					const MyParticle &ip = fParticles.GetParticle(i);
					const MyParticle &jp = fParticles.GetParticle(j);
					const MyParticle &kp = fParticles.GetParticle(k);
					std::string threeParticlesName = ip.GetName();
					threeParticlesName += jp.GetName();
					threeParticlesName += kp.GetName();
					std::vector<std::string>::iterator it = std::find(threeBodyComb.begin(),threeBodyComb.end(), threeParticlesName);
					if (it == threeBodyComb.end()) continue;
					int index = std::distance(threeBodyComb.begin(),it);
					{
						int tmpIdx = startIdx + 2*index;
						fill3BodyHistogram(ip,jp,kp, fHi, tmpIdx);
					}
				}
			}
		}
	}

	//std::cout << startIdx << std::endl;
	//---Post-analysis---//
	////std::cout << molecule[1][1].CoincidenceCount <<":";
	//if (molecule[4][14].CoincidenceCount > 0)  fillHydrogenHistogram(fParticles.GetParticle(1),fParticles.GetParticle(4),fParticles.GetParticle(14),fHi,startIdx);
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
	//			//fillParticleHistograms(p,ph,fIntensities,fHi,secondStartIdx,intPartition);
	//			fillMoleculeHistogram(fParticles.GetParticle(2),fParticles.GetParticle(5),fIntensities,fHi,secondStartIdx, molecule[2][5], intPartition);
	//		}
	//		secondStartIdx +=100;
	//	}
	//}
	//if (MoleculeAnalysis == 1) fillPIPICO(fParticles.GetParticle(0),fHi);


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

///XYZ vector calculation
double calcInnerProduct(const MyParticleHit &ph1,const MyParticleHit &ph2)
{
	return ph1.Px()*ph2.Px()+ph1.Py()*ph2.Py()+ph1.Pz()*ph2.Pz();
}
double calcFormedAngle(const MyParticleHit &ph1,const MyParticleHit &ph2)
{
	const double angle = TMath::ACos(calcInnerProduct(ph1,ph2)/(ph1.P()*ph2.P()));

	return angle * TMath::RadToDeg();
}

//XY stuff
double calcInnerProductXY(const MyParticleHit &ph1,const MyParticleHit &ph2)
{
	return ph1.Px()*ph2.Px()+ph1.Py()*ph2.Py();
}
double calcMagXY(const MyParticleHit &ph)
{
	return TMath::Sqrt(ph.Px()*ph.Px()+ph.Py()*ph.Py());
}
double calcFormedAngleXY(const MyParticleHit &ph1,const MyParticleHit &ph2)
{
	const double angle = TMath::ACos(calcInnerProductXY(ph1,ph2)/(calcMagXY(ph1)*calcMagXY(ph2)));

	return angle * TMath::RadToDeg();
}
//Divide 2D histogram by 1D histogram
void DivideHisto2Dby1D(TH2D *h2d, TH1D *h1d)
{
	if ((h2d->GetNbinsY() != h1d->GetNbinsX()))
	{
		std::cout <<"Histos must have the same binning"<<std::endl;
		return;
	}

	for (int i=1;i<=h2d->GetNbinsX();++i)
	{
		for (int j=1;j<=h2d->GetNbinsY();++j)
		{
			double temp;
			temp = 0.;
			if (fabs(h1d->GetBinContent(j))>1.e-50) temp = h2d->GetBinContent(i,j) / h1d->GetBinContent(j);
			h2d->SetBinContent(i,j,temp);
		}
	}

}