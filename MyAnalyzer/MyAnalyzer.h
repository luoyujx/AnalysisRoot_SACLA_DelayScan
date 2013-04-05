#ifndef _MyAnalyzer_h_
#define _MyAnalyzer_h_

#include <iostream>
#include <vector>
#include <map>
#include <TChain.h>
#include <TSystem.h>
#include <TCanvas.h>
#include "TText.h"

#include "../MyParticle/MyParticleContainer.h"
#include "../FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalEvent.h"
#include "../FilesFromLma2Root/MyEvent/MySortedEvent/MySortedEvent.h"
#include "../FilesFromLma2Root/MyEvent/MySignalAnalyzedEvent/MySignalAnalyzedEvent.h"
#include "../FilesFromLma2Root/MyRootManager/MyHistos.h"
#include "../FilesFromLma2Root/MySettings/MySettings.h"

#include "../MyWaveform.h"

//_____Informatin for calc coincidence_____
struct Molecule
{
public:
	Molecule():
		momSumWindowX(0),
		momSumWindowY(0),
		momSumWindowZ(0),
		momSumFactor(0),
		CoincidenceCount(0)
	{};
public:
	double					momSumWindowX;
	double					momSumWindowY;
	double					momSumWindowZ;
	double					momSumFactor;
	size_t					CoincidenceCount;
};


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
	void					 Analyze();
	void					 FileOpen();
	void					SetParameter(MySettings &set);
	void					OpenIntensityData();
	void					OpenIntPartition();
	void					OpenMomInfoData();
	void					ShowResult();

	//void test(){std::cout<<"Vals have changed"<<std::endl;	for (size_t i=0; i<ParticleInfos.size();ParticleInfos[i++]->Save());}

	//---Setting from comand switch---//
	void		SetFileName(const TString& in)		{fileName = in;}
	void		SetRekMeth(int in)					{rekmeth = in;}
	void		SetMolecule(int in)					{MoleculeAnalysis = in;}
	void		SetCondition(int in)				{extraCondition = in;}
	void		SetIntFileName(const TString& in)	{intFileName = in;}

	TString		GetFileName()const			{return fileName;}
	int			GetRekMeth()const			{return rekmeth;}
	int			GetMolecule()const			{return MoleculeAnalysis;}
	int			GetCondition()const			{return extraCondition;}
	TString		GetIntFileName()const			{return intFileName;}

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

	TCanvas					*canv;
	//some parametrs for analysis
	TString fileName;
	TString intFileName;
	int rekmeth;
	size_t missedTagCount;
	//Analysis frags
	int MoleculeAnalysis;
	int extraCondition;
	bool existIntensityData;
	bool existIntPartition;

	std::map<unsigned int, double>		tagIntensity;
	std::map<unsigned int, double>		tagIntensity2;
	std::vector<double>					intPartition;
	std::vector< std::vector<Molecule> >			molecule;
	std::vector<TText*>					txt;
};


#endif