#ifndef _MyAnalyzer_h_
#define _MyAnalyzer_h_

#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <TChain.h>
#include <TSystem.h>
#include <TCanvas.h>
#include "TText.h"
#include "TBox.h"


#include "../MyParticle/MyParticleContainer.h"
#include "../FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalEvent.h"
#include "../FilesFromLma2Root/MyEvent/MySortedEvent/MySortedEvent.h"
#include "../FilesFromLma2Root/MyEvent/MySignalAnalyzedEvent/MySignalAnalyzedEvent.h"
#include "../FilesFromLma2Root/MyRootManager/MyHistos.h"
#include "../FilesFromLma2Root/MySettings/MySettings.h"

#include "../MyWaveform.h"
//#include "../SQLiteProcessor/DataBase0d.h"
//_____Information for calc coincidence_____
struct Molecule
{
public:
	Molecule():
		momSumWindowX(0),
		momSumWindowY(0),
		momSumWindowZ(0),
		momSumFactor(0),
		CoincidenceCount(0),
		angleCondition(0),
		momSumFactorLow(0),
		momSumFactorUp(2)
	{};
public:
	double					momSumWindowX;
	double					momSumWindowY;
	double					momSumWindowZ;
	double					momSumFactor;
	double					angleCondition;
	double					momSumFactorLow;
	double					momSumFactorUp;

	size_t					CoincidenceCount;
	std::vector<size_t>		CoinHitNbrC;
	std::vector<size_t>		CoinHitNbrI;
};

class MCPToFRegion
{
public:
	MCPToFRegion():particleName(""),tofFrom(0.0),tofTo(0.0){};
	MCPToFRegion(TString pname, double tfr, double tto)
		:particleName(pname),tofFrom(tfr),tofTo(tto) {};
public:
	TString particleName;
	double tofFrom;
	double tofTo;
};
//___Analyze class___
class MyAnalyzer
{
public:
	//MyAnalyzer();
	MyAnalyzer(MySettings &set);
	~MyAnalyzer();

public:
	void					 Init();
	void					 Init(MySettings &set);
	void					 Run();
	void					 Analyze(MyWaveform&);
	void					 FileOpen();
	void					SetParameter(MySettings &set);
	void					OpenIntensityData();
	void					OpenIntPartition();
	void					OpenMomInfoData();
	void					OpenBeamPositionData();
	void					Open3BodyCombination();
	void					OpenMCPToFRegion();

	void					ShowResult();

	//void test(){std::cout<<"Vals have changed"<<std::endl;	for (size_t i=0; i<ParticleInfos.size();ParticleInfos[i++]->Save());}

	//Set from Setting.txt
	void		SetFileName(const TString& in)		{fileName = in;}

	TString		GetFileName()const			{return fileName;}
	int			GetRekMeth()const			{return rekmeth;}
	int			GetMolecule()const			{return MoleculeAnalysis;}
	int			GetCondition()const			{return extraCondition;}
	TString		GetIntFileName()const		{return intFileName;}

private:
	//the particles//
	MyParticleContainer		 fParticles;
	std::vector<double>		 fIntensities;

	//don't bother with whats below this//
	//the events from the trees//
	MyOriginalEvent			 fOE;
	MyOriginalEvent			*fOEp;
	MySignalAnalyzedEvent	 fSAE;
	MySignalAnalyzedEvent	*fSAEp;
	MySortedEvent			 fSE;
	MySortedEvent			*fSEp;

	//the trees bzw. the chain containing the trees//
	TChain					 fOChain;
	TChain					 fSAChain;
	TChain					 fSChain;
	Long64_t				 fNEntries;
	Long64_t				 fEntryIterator;

	//run control//
	bool					 running;
	TProcessEventTimer		 processTimer;
	TTimer					 runTimer;

	//the histograms//
	MyHistos				 fHi;
	//covariance calcuration stuff//
	MyWaveform				fWf;
	//Canvas for showResult()
	TCanvas					*canv;

	//---Parametrs for analysis (defined by setting.txt)
	//ROOT file name
	TString fileName;
	//File name of intensity data
	TString intFileName;
	//Path name of 0D data
	TString path0D_DataBase;
	//File name of Momentum sum imformation
	TString MomSumInfoName;

	TString whichParticles;
	//reconstruction method (resort parameter)
	int rekmeth;
	//conversion factor for intensity
	double factorBM1;//(to uJ)
	double factorPD; //(to uJ/um^2)
	//step size for trend histogram
	int trendStep;
	//Condition to secelt FEL intensity
	double intensityLowerLimit;
	double intensityUpperLimit;
	//Counter for missed intensity data
	size_t missedTagCount;
	//C-I momentums limits
	double momFactorLowerLimit;
	double momFactorUpperLimit;
	double angleCondition;
	//PM to Delay (2014A)
	double factorPMD;
	double factorPMDOffset;

	int tagFrom;
	int tagTo;

	int delayBins;
	double delayFrom;
	double delayTo;


	//---Analysis frags
	int MoleculeAnalysis;
	int extraCondition;
	bool existIntensityData;
	int method0D_Data;
	bool existIntPartition;
	bool checkingResult;
	bool afterAnalysis;
	bool selectIntensity;


	//BM1 data
	std::map<unsigned int, double>		tagIntensity;
	//PD data
	std::map<unsigned int, double>		tagIntensity2;
	//Intensity partition for making variable bin histogram
	std::vector<double>					intPartition;
	//various molecule data for coincidence (Momentum information, coincidence count, ...)
	std::vector< std::vector<Molecule> > molecule;
	//Position data
	std::map<unsigned int, double>		beamPosX;
	std::map<unsigned int, double>		beamPosY;
	//3-body combination data
	std::vector<std::string>			threeBodyComb;
	//ToF region of MCP intensity 
	std::vector<MCPToFRegion>			mcpTofRegion;

	//---Indicate some information on mass and tof spectrum
	std::vector<TText>					txtMass;
	std::vector<TText>					txtTof;
	std::vector<TBox>					boxTof;

	//DataBase0d db;
};


#endif