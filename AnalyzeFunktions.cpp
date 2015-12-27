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
#include <cmath>
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
#include <TRandom3.h>

//#define XVelocity 0.0003704 //(mm/ns)
//#define YVelocity -0.00037775
#define XVelocity 0. //(mm/ns)
#define YVelocity 0.


//#define DELAY_RANDOM
#ifdef DELAY_RANDOM
TRandom3 RandomGene;
#endif
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
	else if (whichParticles == "Xenon")
	{
		//---SACLA Ar atom
		AddXenon(particles);
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
	else if (whichParticles == "CH2I2-II")
	{
		//---CH2I2 for SACLA2015B and Lab experiment
		AddCH2I2_II(particles);
	}
	else if (whichParticles == "CH2I2-CI")
	{
		//---CH2I2 for SACLA2015B and Lab experiment
		AddCH2I2_CI(particles);
	}
	else if (whichParticles == "CH2I2-CII")
	{
		//---CH2I2 for SACLA2015B and Lab experiment
		AddCH2I2_CII(particles);
	}
	else if (whichParticles == "CH2BrI")
	{
		//---CH2BrI for SACLA2015B and Lab experiment
		AddCH2BrI(particles);
	}
	else if (whichParticles == "CH2ClI")
	{
		//---CH2ClI for SACLA2015B and Lab experiment
		AddCH2ClI(particles);
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
//________Analysis main loop___________________________________________________________________________________________________________________________________________________
void MyAnalyzer::Analyze(MyWaveform &wf)
{
	MyDetektor &rd = fSE.GetDetektor(0);
	int startIdx=10;
	//for Trend Histogram
	const int skipCounter = static_cast<int>(fEntryIterator/trendStep);
	// clear the containers//
	fIntensities.clear();
	fDelays.clear();
	fFlag.clear();
	fParticles.ClearParticles();
	//Get Tag number on this event
	const unsigned int TagNumber = fOE.GetEventID();
	//Check tag is even (IR + FEL)
	if (TagNumber % 2 != 0) return;
	// Get XFEL intensity, Delay
	if (delayScan == 0) // MySQL databse is not used
	{
		fIntensities.push_back(0.); //[0] Upper PD
		fIntensities.push_back(0.); //[1] Lower PD
		fIntensities.push_back(0.);	//[2] Upper PD + Lower PD
		fDelays.push_back(0.);	//[0] EH Delay	
		fDelays.push_back(0.);	//[1] Jitter 
		fDelays.push_back(0.);	//[2] Cor. Delay
	}
	if (delayScan == 1) // without Jitter from Timing moniter
	{
		// Chack 0D databse
		if (DB0d.GetStatusAndData(TagNumber, 0).first == 0)
		{
			if (missedTagCount % 100 == 0) std::cout << "\r" << TagNumber << " is not found!!";
			missedTagCount++;
			return;
		}
		// Intensity
		if ((_isnan(DB0d.GetStatusAndData(TagNumber, 0).second)))
		{
			std::cout << "Intensity of " << TagNumber << " is NaN " << std::endl;
			missedTagCount++;
			return;
		}
		fIntensities.push_back((DB0d.GetStatusAndData(TagNumber, 0).second) * factorBM1 * 1000000); //[0] BM1
		// Delay
		if (_isnan(DB0d.GetStatusAndData(TagNumber, 1).second))
		{
			std::cout << "Delay of " << TagNumber << " is NaN " << std::endl;
			missedTagCount++;
			return;
		}
		// Optical shutter is open or close
		fFlag.push_back(static_cast<bool>(DB0d.GetStatusAndData(TagNumber, 2).second)); //[0] Optical shutter
		if (fFlag[0] != optShutterOpen)
		{
			//std::cout << "Optical shutter of " << TagNumber << " is close " << std::endl;
			//missedTagCount++;
			return;
		}
		fDelays.push_back((factorPMDOffset - DB0d.GetStatusAndData(TagNumber, 1).second) / factorPMD);	//[0] EH Delay	
		fDelays.push_back(0.);																			//[1] Jitter 
		//fDelays.push_back(fDelays[0]);																	//[2] Cor. Delay
		fDelays.push_back(0.);
	}
	if (delayScan == 2) // with Jitter from Timing moniter, This mode should be chacked.
	{
		// Check 0D database
		if (DBTM.GetStatusAndData(TagNumber, 0).first == 0)
		{
			//if (missedTagCount % 100 == 0) std::cout << "\r" << TagNumber << " is not found!!";
			//std::cout << TagNumber << " is not found!!" << std::endl;
			missedTagCount++;
			return;
		}
		// Intensity
		if ((_isnan(DB0d.GetStatusAndData(TagNumber, 0).second)))
		{
			//std::cout << "Intensity of " << TagNumber << " is NaN " << std::endl;
			missedTagCount++;
			return;
		}
		fIntensities.push_back((DB0d.GetStatusAndData(TagNumber, 0).second) * factorBM1 * 1000000); //[0] BM1
		//
		// Delay
		if (_isnan(DB0d.GetStatusAndData(TagNumber, 1).second))
		{
			//std::cout << "Delay of " << TagNumber << " is NaN " << std::endl;
			missedTagCount++;
			return;
		}
		//
		if ((_isnan(DB0d.GetStatusAndData(TagNumber, 3).second)) || (_isnan(DBTM.GetStatusAndData(TagNumber, 0).second)) || (_isnan(DBTM.GetStatusAndData(TagNumber, 1).second)))
		{
			//std::cout << "Jitter of " << TagNumber << " is NaN " << std::endl;
			missedTagCount++;
			return;
		}
		if (!DBTM.GetStatusAndData(TagNumber, 0).second)
		{
			//std::cout << "Jitter of " << TagNumber << " is not analysed " << std::endl;
			missedTagCount++;
			return;
		}
		fDelays.push_back((DB0d.GetStatusAndData(TagNumber, 1).second - factorPMDOffset) / factorPMD * 1000); //[0] EH Delay [fs]
		//
		fFlag.push_back(static_cast<bool>(DB0d.GetStatusAndData(TagNumber, 2).second)); //[0] Optical shutter
		//std::cout << "Optical shutter valid status (" << TagNumber << "): "<< fFlag[0] << std::endl;
		if (fFlag[0] != optShutterOpen)
		{
			std::cout << "Optical shutter of " << TagNumber << " is close " << std::endl;
			missedTagCount++;
			return;
		}
		fDelays.push_back(factorTM*(DBTM.GetStatusAndData(TagNumber, 1).second - factorTMOffset)); //[1] Jitter 
		// Testing random jitter
#ifdef DELAY_RANDOM
		fDelays[1] = RandomGene.Gaus(-33.45, 304.8);
#endif
		//delay with motor position and jitter
		fDelays.push_back(fDelays[0] + fDelays[1]); // Cor. Delay + (factorTMPMDOffset - DB.GetStatusAndData(TagNumber, 5).second) / factorTMPMD);
		//overlaping deley
		fDelays.push_back(std::fmod(fDelays[0] + fDelays[1] + 10000, 284));

	}
	//
	////-----------------------------------
	//// No selected histgrams //// index = 10
	// XFEL intenisty
	fHi.fill(startIdx + 1, "DelayVsJitter", fDelays[0], fDelays[1], "Delay [fs]", "Jitter [fs]", delayBins, delayFrom, delayTo, delayBins, delayFrom, delayTo); // Delay vs XFEL intenisty
	fHi.fill(startIdx + 2, "IntensityXFEL", fIntensities[0], "[arb. unit]", 1000, 0, 1000); // XFEL intensity BM1
	startIdx += 3; // index = 13
	//
	// Delay
	fHi.fill(startIdx + 0, "Delay_MotorPosition", fDelays[0], "[fs]", delayBins, delayFrom, delayTo); // Delay (inclued jitter)
	fHi.fill(startIdx + 1, "Jitter", fDelays[1], "[fs]", delayBins*100, delayFrom, delayTo); // Jitter
	fHi.fill(startIdx + 2, "Delay", fDelays[2], "[fs]", delayBins, delayFrom, delayTo);  // Delay (no jitter)	
	startIdx += 3; // index = 16
	int maxTrend = 100000;
	// Trend plot XFEL intensity, Delay, number of hit.
	fHi.fill(startIdx + 2, "TrendIntensityXFEL", skipCounter, Form("[shots/%d]", trendStep), maxTrend, 0, maxTrend, "Trend", fIntensities[0]); // XFEL intensity BM1
	fHi.fill(startIdx + 3, "TrendDelay(MotorPosition)", skipCounter, Form("[shots/%d]", trendStep), maxTrend, 0, maxTrend, "Trend", fDelays[0]);  // Delay (inclued jitter)
	fHi.fill(startIdx + 4, "TrendJitter", skipCounter, Form("[shots/%d]", trendStep), maxTrend, 0, maxTrend, "Trend", fDelays[1]); // Jitter
	fHi.fill(startIdx + 5, "TrendDelay", skipCounter, Form("[shots/%d]", trendStep), maxTrend, 0, maxTrend, "Trend", fDelays[2]);  // Delay (no jitter)	
	fHi.fill(startIdx + 6, "TrendNumberOfHits", skipCounter, Form("[shots/%d]", trendStep), maxTrend, 0, maxTrend, "Trend", rd.GetNbrOfHits()); // number of hits
	startIdx += 7; // index = 23
	//---------	
	if (selectIntensity)
	{
	//---skip this shot event if FEL is below the lower limit
		if (fIntensities[0]<intensityLowerLimit) return;
		//---skip this shot event if FEL is over the upper limit
		if (fIntensities[0]>intensityUpperLimit) return;
	}
	if (fIntensities.size())
	{
		fHi.fill(startIdx, "IntensitySelected", fIntensities[0], "[arb. unit]", 1000, 0, 1000); // selected XFEL intensity
	}
	startIdx++;
	//---------		
	//---skip this event if Delay (Delay1 + Delay2) is below the lower limit or over the upper limit
	if (selectDelay)
	{
		if (fDelays[2]<delayLowerLimit) return;
		if (fDelays[2]>delayUpperLimit) return;
	}
	if (fDelays.size())
	{
		fHi.fill(startIdx, "DelaySelected (no correction)", fDelays[0], "[fs]", delayBins, delayFrom, delayTo); // selected Delays (no correction)
	}
	startIdx++;
	//---------		
	//--skip this shot event if Jitter is below the lower limit or over the upper limit
	if (selectJitter)
	{
		if (fDelays[1]<jitterLowerLimit) return;
		if (fDelays[1]>jitterUpperLimit) return;
	}
	if (fDelays.size())
	{
		fHi.fill(startIdx, "JitterSelected", fDelays[1], "[fs]", delayBins, delayFrom, delayTo); // selected Jitter
	}
	startIdx++;
	//---------
	if (fDelays.size())
	{
		fHi.fill(startIdx, "DelaySelected", fDelays[2], "[fs]", delayBins, delayFrom, delayTo); // selected Delays
	}
	startIdx++;
	//---------
	fHi.fill(startIdx,"NumberOfHits",rd.GetNbrOfHits(),"Number of Hits",100,0,100);
	startIdx++; // index = 24
	//
	//if (fSAE.GetChannel(7 - 1).GetNbrPeaks() == 0) return;
	//
	//// Delay vs Shots, XFELintensity, NumberOfHits
	fHi.fill(startIdx + 0, "DelayVsShots", fDelays[2], "Delay [fs]", delayBins, delayFrom, delayTo); // Delay vs Shots
	fHi.fill(startIdx + 1, "DelayVsXFELintensity", fIntensities[0], fDelays[2], "Delay [fs]", "XFEL intensity [arb .unit]", 1000, 0, 1000, delayBins, delayFrom, delayTo); // Delay vs XFEL intenisty
	fHi.fill(startIdx + 2, "DelayVsNumberOfHits", rd.GetNbrOfHits(), fDelays[2], "Delay [fs]", "Number of Hits", 100, 0, 100, delayBins, delayFrom, delayTo); // Delay vs Number of Hits
	startIdx += 3; // index = 27
	//
	int secondStartIdx = startIdx + 20;
	//go through all resorted detektorhits//
	for (size_t i=0; i<rd.GetNbrOfHits();++i)
	{
		MyDetektorHit &dh = rd.GetHit(i);
		//the tof is just the timing of the mcp signal//
		dh.SetTof(dh.Time());
		//detektor & tof for all Hits//
		const double maxPos	= (rd.GetRunTime()+30)*rd.GetSfU();
		MyParticle &LastParticle = fParticles.GetParticle(fParticles.GetNbrOfParticles()-1);
		const double maxTof	= fOE.GetNbrSamples()*fOE.GetSampleInterval()*1e9;
		//
		fHi.fill(startIdx + 0, "Det", dh.X(), dh.Y(), "x [mm]", "y [mm]", 300, -maxPos, maxPos, 300, -maxPos, maxPos);
		fHi.fill(startIdx + 1, "Tof", dh.Tof(), "tof [ns]", 30000, 0, maxTof);
		fHi.fill(startIdx + 2, "XPosVsTof", dh.Tof(), dh.X(), "tof [ns]", "x [mm]", 30000, 0, maxTof, 300, -maxPos, maxPos);
		fHi.fill(startIdx + 3, "YPosVsTof", dh.Tof(), dh.Y(), "tof [ns]", "y [mm]", 30000, 0, maxTof, 300, -maxPos, maxPos);
		fHi.fill(startIdx + 4, "DelayVsTof", dh.Tof(), fDelays[2], "tof [ns]", "Delay [fs]", 5000, 0, maxTof, delayBins, delayFrom, delayTo);
		fHi.fill(startIdx + 5, "XFELIntensityVsTOF", dh.Tof(), fIntensities[2], "tof [ns]", "XFEL intensity [arb. unit]", 5000, 0, maxTof, 1000, 0, 1000, "Ion");
		//
		//-----Particle(0)---Ion---
		//get the particle from the vector//
		MyParticle &p = fParticles.GetParticle(0);
		//if this hit fits both conditions then add the Hit to this Particle and fill the histo for this hit//
		if (p.CheckTofAndPos(dh))
			//select hit by reconstruction method//
				if (dh.RekMeth() < rekmeth)//added by motomura
				{
					const MyParticleHit &ph = p.AddHit(dh);
					fHi.fill(startIdx + 6, "Mass", ph.Mass(), "Mass/q", 10000, 0, 500, "Ion");
					fHi.fill(startIdx + 7, "Det", ph.X(), ph.Y(), "x [mm]", "y [mm]", 300, p.GetCondRadX() - p.GetCondRad()*1.3, p.GetCondRadX() + p.GetCondRad()*1.3, 300, p.GetCondRadY() - p.GetCondRad()*1.3, p.GetCondRadY() + p.GetCondRad()*1.3, "Ion");
					fHi.fill(startIdx + 8, "Tof", ph.Tof(), "tof [ns]", 30000, 0, maxTof, "Ion");
					fHi.fill(startIdx + 9, "XPosVsTof", ph.Tof(), ph.X(), "tof [ns]", "x [mm]", 10000, p.GetCondTofFr() - p.GetCondTofRange()*0.1, p.GetCondTofTo() + p.GetCondTofRange()*0.1, 300, p.GetCondRadX() - p.GetCondRad()*1.1, p.GetCondRadX() + p.GetCondRad()*1.1, "Ion");
					fHi.fill(startIdx + 10, "YPosVsTof", ph.Tof(), ph.Y(), "tof [ns]", "y [mm]", 5000, p.GetCondTofFr() - p.GetCondTofRange()*0.1, p.GetCondTofTo() + p.GetCondTofRange()*0.1, 300, p.GetCondRadX() - p.GetCondRad()*1.1, p.GetCondRadX() + p.GetCondRad()*1.1, "Ion");
					fHi.fill(startIdx + 11, "DelayVsTOF", dh.Tof(), fDelays[2], "tof [ns]", "Delay [ps]", 1000, -p.GetT0(), maxTof - p.GetT0(), delayBins, delayFrom, delayTo, "Ion");
					//
					fHi.fill(startIdx + 12, "DetCorScale", ph.XCorRotScl(), ph.YCorRotScl(), "x [mm]", "y [mm]", 300, p.GetCondRadX() - p.GetCondRad()*1.3, p.GetCondRadX() + p.GetCondRad()*1.3, 300, p.GetCondRadY() - p.GetCondRad()*1.3, p.GetCondRadY() + p.GetCondRad()*1.3, "Ion");
					fHi.fill(startIdx + 13, "TofCor", ph.TofCor(), "tof [ns]", 30000, -p.GetT0(), maxTof - p.GetT0(), "Ion");
					fHi.fill(startIdx + 14, "XPosVsTofCor", ph.TofCor(), ph.X(), "tof [ns]", "x [mm]", 5000, p.GetCondTofFr() - p.GetCondTofRange()*0.1 - p.GetT0(), p.GetCondTofTo() + p.GetCondTofRange()*0.1 - p.GetT0(), 300, p.GetCondRadX() - p.GetCondRad()*1.1, p.GetCondRadX() + p.GetCondRad()*1.1, "Ion");
					fHi.fill(startIdx + 15, "YPosVsTofCor", ph.TofCor(), ph.Y(), "tof [ns]", "y [mm]", 5000, p.GetCondTofFr() - p.GetCondTofRange()*0.1 - p.GetT0(), p.GetCondTofTo() + p.GetCondTofRange()*0.1 - p.GetT0(), 300, p.GetCondRadX() - p.GetCondRad()*1.1, p.GetCondRadX() + p.GetCondRad()*1.1, "Ion");
					fHi.fill(startIdx + 16, "DelayVsTOFCor", ph.TofCor(), fDelays[2], "tof [ns]", "Delay [fs]", 1000, -p.GetT0(), maxTof - p.GetT0(), delayBins, delayFrom, delayTo, "Ion");
				}
		//
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
								//if ((dh.RekMeth() != 8) && (dh.RekMeth() != 11) && (dh.RekMeth() != 12))
							{
								const MyParticleHit &ph = p.AddHit(dh);
								//Check the angle of phiZX. if ph is out of condition, it delete.
								if (p.CheckPhiZX(ph))
								{
									fillParticleHistograms(p, ph, fIntensities, fDelays, fHi, secondStartIdx, intPartition, delayBins, delayFrom, delayTo, selectThetaZ, thetaZLowerLimit, thetaZUpperLimit);
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
	double iTotalChage = 0.0;
	for (size_t j=0;j<fParticles.GetNbrOfParticles();++j)//particle index j=0 is "ion"  
	{
		fHi.fill(startIdx+j,"NumberOfHits",fParticles.GetParticle(j).GetNbrOfParticleHits(),"Number of Hits",100,0,100,Form("%s",fParticles.GetParticle(j).GetName()));
		if (j > 0) fHi.fill(startIdx+fParticles.GetNbrOfParticles()+j,Form("TrendParticleHits%02d", j),skipCounter,Form("[shots/%d]",trendStep),1000,0,1000,"Trend",fParticles.GetParticle(j).GetNbrOfParticleHits());
		//Iodine rate
		if ((fParticles.GetParticle(j).GetMass_au()*MyUnitsConv::au2amu() > 120)&&(fParticles.GetParticle(j).GetCharge_au() < 10))
		{
			iodineCount += fParticles.GetParticle(j).GetNbrOfParticleHits();
			for (size_t k = 0; k < fParticles.GetParticle(j).GetNbrOfParticleHits(); ++k)
			{
				iTotalChage += fParticles.GetParticle(j).GetCharge_au();
			}

		}
		for (size_t k=0;k<fParticles.GetParticle(j).GetNbrOfParticleHits();++k)
		{
			fHi.fill(startIdx+fParticles.GetNbrOfParticles(),"NumberOfParticleHits",j,"Particle Number",fParticles.GetNbrOfParticles()+1,0,fParticles.GetNbrOfParticles()+1);
			fHi.fill(startIdx + 2 * fParticles.GetNbrOfParticles(), "DelayVsNumberOfParticleHits", j, fDelays[2], "Particle Number", "delay [fs]", fParticles.GetNbrOfParticles() + 1, 0, fParticles.GetNbrOfParticles() + 1, delayBins, delayFrom, delayTo);
		}
	}
	startIdx += (2*fParticles.GetNbrOfParticles()+2);
	fHi.fill(startIdx,"IodineHits",iodineCount,"Number of Hits",100,0,100);
	fHi.fill(startIdx + 1, "DelayVsIodineMeanChage", iTotalChage / iodineCount, fDelays[2], "Mean of Charge", "Delay [fs]",  100, 0, 10, delayBins, delayFrom, delayTo);
	startIdx +=2;
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
							fillMoleculeHistogramCH2I2(ip, jp, fIntensities, fHi, startIdx, molecule[i][j], fDelays, delayBins, delayFrom, delayTo);
							startIdx += 115;
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


	//C-I-I 3-body coincidence
	if (MoleculeAnalysis == 3)
	{
		//loop from particle No.1 (except Ion)
		for (size_t i = 1; i < fParticles.GetNbrOfParticles(); ++i)
		{
			//i should be only Carbon (kind of particle 12)
			const MyParticle &ip = fParticles.GetParticle(i);
			if (ip.GetKindParticle() != 12) continue;
			for (size_t j = 1; j < fParticles.GetNbrOfParticles(); ++j)
			{
				//j should be only Iodine (kind of particle 13)
				const MyParticle &jp = fParticles.GetParticle(j);
				if (jp.GetKindParticle() != 13) continue;
				//check if ip and jp is molecule
				for (size_t k = j; k < fParticles.GetNbrOfParticles(); ++k)
				{
					const MyParticle &kp = fParticles.GetParticle(k);
					if (kp.GetKindParticle() != 13) continue;
					{
						//std::cout << ip.GetName() << jp.GetName() << kp.GetName() << std::endl;
						fillMoleculeHistogramCH2I2_3body(ip, jp, kp, fIntensities, fHi, startIdx, /*molecule[i][j],*/ fDelays, delayBins, delayFrom, delayTo);
						startIdx += 115;
					}
				}
			}
		}
	}
	//---gate by certain particle hits
	//reserve ID for fillMoleculeHistogram2
	//const int nbrOfHistosInfillMol2 = 27;//
	//secondStartIdx = startIdx + (fParticles.GetNbrOfParticles()+fParticles.GetNbrOfParticles()*fParticles.GetNbrOfParticles())*nbrOfHistosInfillMol2;
	//if (MoleculeAnalysis == 2)
	//{
	//	//loop from particle 1 for excepting Ion
	//	for (size_t i=1;i<fParticles.GetNbrOfParticles();++i)
	//	{
	//		const MyParticle &ip = fParticles.GetParticle(i);
	//		//check whether ip have counts and gatting paticle (I+, I++, ...)
	//		if ((ip.GetNbrOfParticleHits())&&(ip.GetCoinGroup()==1 || ip.GetCoinGroup()==2 ))
	//		{
	//			int hitCounter = 0;
	//			for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)
	//			{
	//				if ((j != i) && (fParticles.GetParticle(j).GetCoinGroup() == 1))
	//					hitCounter += fParticles.GetParticle(j).GetNbrOfParticleHits();
	//			}
	//			if (!hitCounter)
	//			{
	//				//fill Ion spectra 
	//				fillSpectra(ip,fParticles.GetParticle(0),fHi, secondStartIdx+ (i*10));
	//				//loop to fill the target particles
	//				for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)
	//				{
	//					const MyParticle &jp = fParticles.GetParticle(j);
	//					//Target particle (C+.N+..O+..CO+..,H+)
	//					if ((jp.GetKindParticle() == 1)&&(jp.GetCoinGroup()==0 || jp.GetCoinGroup()==2))
	//					{
	//						int tmpIdx = startIdx + (i+j*fParticles.GetNbrOfParticles())*nbrOfHistosInfillMol2;
	//						fillMoleculeHistogram2(ip,jp,fIntensities,fHi,tmpIdx);
	//						fillAngleHistogram(ip, jp, fHi, tmpIdx + 30);	
	//					}
	//				}

	//			}			
	//		}
	//	}
	//	//Skip already used ID
	//	startIdx += (fParticles.GetNbrOfParticles()+fParticles.GetNbrOfParticles()*fParticles.GetNbrOfParticles())*nbrOfHistosInfillMol2;
	//	startIdx += fParticles.GetNbrOfParticles()*10;

	//	//_______________________________//
	//	//___3-body angular correration__//
	//	//_______________________________//
	//	//1st loop from particle 1 for excepting Ion
	//	for (size_t i=1;i<fParticles.GetNbrOfParticles();++i)
	//	{
	//		//2nd loop
	//		for (size_t j=1;j<fParticles.GetNbrOfParticles();++j)
	//		{
	//			//3rd loop
	//			for (size_t k=1;k<fParticles.GetNbrOfParticles();++k)
	//			{
	//				const MyParticle &ip = fParticles.GetParticle(i);
	//				const MyParticle &jp = fParticles.GetParticle(j);
	//				const MyParticle &kp = fParticles.GetParticle(k);
	//				std::string threeParticlesName = ip.GetName();
	//				threeParticlesName += jp.GetName();
	//				threeParticlesName += kp.GetName();
	//				std::vector<std::string>::iterator it = std::find(threeBodyComb.begin(),threeBodyComb.end(), threeParticlesName);
	//				if (it == threeBodyComb.end()) continue;
	//				int index = std::distance(threeBodyComb.begin(),it);
	//				{
	//					int tmpIdx = startIdx + 2*index;
	//					fill3BodyHistogram(ip,jp,kp, fHi, tmpIdx);
	//				}
	//			}
	//		}
	//	}
	//}

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
	if (MoleculeAnalysis == 1) fillPIPICO(fParticles.GetParticle(0),fHi);


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