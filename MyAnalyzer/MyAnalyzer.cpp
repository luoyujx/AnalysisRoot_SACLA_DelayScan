#include <iostream>
#include <iomanip>
#include <fstream>
#include <utility>
#include <TTree.h>
#include <TFile.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TGClient.h>



#include "MyAnalyzer.h"
#include "../MyGui/MyGui.h"
#include "../FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalChannel/MyOriginalChannel.h"
#include "../FilesFromLma2Root/MyEvent/MySortedEvent/MyDetektor/MyDetektor.h"
#include "../FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalEventInfo.h"
#include "../FilesFromLma2Root/MyEvent/MySortedEvent/MySortedEventInfo.h"
#include "../FilesFromLma2Root/MyEvent/MySignalAnalyzedEvent/MySignalAnalyzedEventInfo.h"

#include "../AnalyzeFunktions.h"
#include "../AnalyzeFuncFill.h"

//#include "../MyCovariance.h"
//#include "../MyWaveform.h"

using namespace std;

//##################################################################################
//_____________________________The class Members______________________________________________________________________________________________________________________________
MyAnalyzer::MyAnalyzer(MySettings &set):
	fOChain("OriginalEvent"),
	fSAChain("SignalAnalyzedEvent"),
	fSChain("SortedEvent"),
	fHi(false,15000),//Max 10000?
	running(false),
	processTimer(100),
	molecule(0, std::vector<Molecule>(0)),
	missedTagCount(0)
{
	//--------------do not modify this part (unless you know what you're doing START-----------------------------//
	//setup chains//
	int nOrigFiles = fOChain.Add("OrigEvent*.root");
	int nSigAFiles = fSAChain.Add("SigAnaEvent*.root");
	int nSortFiles = fSChain.Add("SortEvent*.root");

	std::cout << "NbrOrigFiles("<<nOrigFiles<<"), NbrSigAnaFiles("<<nSigAFiles<<"),  NbrSortFiles("<<nSortFiles<<")"<<std::endl;
	//make the signal analyzed chain a friend of the sort chain//
	fOChain.AddFriend(&fSAChain);
	fOChain.AddFriend(&fSChain);

	//read the entries of the sort chain (will also load the trees//
	fNEntries = fOChain.GetEntries();
	std::cout << "there are a total of "<<fNEntries<<" Entries in the Chain(Orig)"<<std::endl;
	std::cout << "there are a total of "<<fSAChain.GetEntries()<<" Entries in the Chain(SigAna)"<<std::endl;
	std::cout << "there are a total of "<<fSChain.GetEntries()<<" Entries in the Chain(Sorted)"<<std::endl;

	//set the Branches of the trees in the chain//
	fOEp	= &fOE;
	fSAEp	= &fSAE;
	fSEp	= &fSE;
	fOChain.SetBranchAddress("OriginalEvent", &fOEp);
	fOChain.SetBranchAddress("SignalAnalyzedEvent", &fSAEp);
	fOChain.SetBranchAddress("SortedEvent", &fSEp);

	////get the infos from the trees in the chains//
	MyOriginalEventInfo			*OEInfo		= dynamic_cast<MyOriginalEventInfo*>(fOChain.GetTree()->GetUserInfo()->At(0));
	MySignalAnalyzedEventInfo	*SAEInfo	= dynamic_cast<MySignalAnalyzedEventInfo*>(fSAChain.GetTree()->GetUserInfo()->At(0));
	MySortedEventInfo			*SEInfo		= dynamic_cast<MySortedEventInfo*>(fSChain.GetTree()->GetUserInfo()->At(0));
	//fill the events with infos//
	fOE.ReadFromEventInfo(*OEInfo);
	fSAE.ReadFromEventInfo(*SAEInfo);
	fSE.ReadFromEventInfo(*SEInfo);

	const TString whichPar = set.GetString("WhichParticles","");
	DefineParticlesAndRootFile(fParticles,fHi,whichPar);
	const bool UseGUI=static_cast<int>(set.GetValue("UseGUI", true)+0.1);
	if (UseGUI)
	{
		//create gui//
		//run the gui and give it the particle infos of the particles that need to be calibrated//
		//and then setup the connection//
		MyGui *gui = new MyGui(gClient->GetRoot(),fParticles.GetParticleInfos());
		gui->Connect("ValsChanged()","MyAnalyzer",this,"Init()");
	}
		//start run//
		runTimer.Connect("Timeout()","MyAnalyzer",this,"Run()");
		runTimer.Start(1000);
}
MyAnalyzer::~MyAnalyzer()
{
	if (canv) delete canv;
}

//___________________________________________________________________________________________________________________________________________________________
void MyAnalyzer::FileOpen()
{
	if (fileName == "")
		fHi.OpenRootFile("analysis.root");
	else 
		fHi.OpenRootFile(fileName);
}
	
//___________________________________________________________________________________________________________________________________________________________
void MyAnalyzer::Init()
{
	//save the current settings to the ini files//
	fParticles.SaveParticleInfos();

	//reset the iterator//
	fEntryIterator=0;

	//clear all histograms//
	fHi.ResetAll();

	//read the particle info to the actual particle//
	fParticles.Init();

	//Init Covariance stuff//
	//fWf.Init(fOE,fHi);

	if (MoleculeAnalysis == 1)OpenMomInfoData();
}
void MyAnalyzer::Init(MySettings &set)
{
	//set parameters from setting.ext
	SetParameter(set);

	//call Initialize
	Init();
}

//________________________This should not be modified___________________________________________________________________________________________________________________________________
void MyAnalyzer::Run()
{
	//std::cout <<"enter Run"<<std::endl;
	//turn the timer off//
	runTimer.TurnOff();
	//create a flag that shows wether we have analyzed things//
	bool WasRunningBefore=false;
	bool realyBreak=false;
	//run while we are still analysing the entries from the tree//
	size_t firstTAG = 0;
	while(fEntryIterator < fNEntries)
	{
		if (fEntryIterator % 1000 == 0)  std::cout << "\r" << "Entry Number :"<< std::setw(7) << std::setfill(' ') << fEntryIterator;

		//Clear the events//
		fOE.Clear();
		fSAE.Clear();
		fSE.Clear();
		//set the flag//
		WasRunningBefore=true;
		//read the entry from the chain//
		fOChain.GetEntry(fEntryIterator);
		fSAChain.GetEntry(fEntryIterator);
		fSChain.GetEntry(fEntryIterator);

		//check EventID//
		if (fEntryIterator == 0) firstTAG = fOE.GetEventID();
		if ((fOE.GetEventID()!=fSAE.GetEventID())||(fOE.GetEventID()!=fSE.GetEventID())) 
		{
			std::cout << std::endl << "******Error!!! EventID Mismatch !!!! ******"<<std::endl;
			realyBreak=true;
			break;
		}
		//analyze the event//
		Analyze();
		//calc covariance map//
		//fWf.ExtractWaveform(fOE,fHi,7-1);
		//increase the counter//
		fEntryIterator++;
		//if(fEntryIterator > 1000) {std::cout << "user requested break"<<std::endl;realyBreak=true;break;}
		//the timer will only process events when it has timed out//
		if (processTimer.ProcessEvents()) {std::cout << "user requested break"<<std::endl;realyBreak=true;break;}
	}
	if (WasRunningBefore)
	{
		if (afterAnalysis) fillHistosAfterAnalyzis(fParticles.GetParticles(),fHi,intPartition.size()-1);
		//fWf.FillHist(fHi);
		std::cout << "<- Done, now saving Histograms!!!!"<<std::endl;
		std::cout << "First TAG: "<<firstTAG << " Last Tag: " << fOE.GetEventID()<< std::endl;

		if (missedTagCount) std::cout << "Can not find "<< missedTagCount << " intensity data." << std::endl;
		fHi.FlushRootFile();
		if (checkingResult) ShowResult();
	}
	
	//restart run at this time//
	if (!realyBreak) runTimer.Start(3000);
	//std::cout << "leaving run"<<std::endl;
}

//_____________
void MyAnalyzer::SetParameter(MySettings &set)
{
	//set parameters
	intFileName = set.GetString("IntensityFile","Intensity.txt");
	MomSumInfoName = set.GetString("MomSumInfoFile","MomentumInfo.txt");
	rekmeth = static_cast<int>(set.GetValue("ReconstructionMethod", 20)+0.1);
	MoleculeAnalysis = static_cast<int>(set.GetValue("Molecule", 0)+0.1);
	extraCondition = static_cast<int>(set.GetValue("ExtraCondition", false)+0.1);
	existIntensityData=static_cast<int>(set.GetValue("IntensityData", false)+0.1);
	existIntPartition=static_cast<int>(set.GetValue("IntensityPartition", false)+0.1);
	factorBM1=set.GetValue("ConversionFactorBM1", 10e+9);
	factorPD=set.GetValue("ConversionFactorPD", 10000);
	selectIntensity=static_cast<int>(set.GetValue("SelectIntensity", false)+0.1);
	intensityLowerLimit=set.GetValue("IntensityLowerLimit", 0.0);
	intensityUpperLimit=set.GetValue("IntensityUpperLimit", 100000);
	checkingResult=static_cast<int>(set.GetValue("CheckResult", false)+0.1);
	afterAnalysis=static_cast<int>(set.GetValue("AfterAnalysis", false)+0.1);
	trendStep=static_cast<int>(set.GetValue("TrendStep", 100)+0.1);
}


//
//_____Read Intensity DATA
void MyAnalyzer::OpenIntensityData()
{
	if ((intFileName == "")||(!existIntensityData)) return;
	
	std::ifstream ifs(intFileName,std::ios::in);
	if (ifs.fail()){
		std::cout<<"Can not open "<<intFileName<<std::endl;
		return;
	}

	unsigned int uintBuf = 0;
	double doubleBuf1;
	double doubleBuf2;
	char tmp[256];
	while (!ifs.eof())
	{
		//read the data Tag and Intensity (uint/double)
		ifs >> uintBuf >> doubleBuf1 >> doubleBuf2;
		//go to nextline
		ifs.getline(tmp,256);
		if (uintBuf % 6 != 0)
			std::cout<< "wrong Tag number!! "<< uintBuf;
		if (!ifs.fail())
		{
			//add to map (tagIntensity)
			tagIntensity.insert(pair<unsigned int, double>(uintBuf,doubleBuf1));
			tagIntensity2.insert(pair<unsigned int, double>(uintBuf,doubleBuf2));
		}
	}

	std::cout << "Intensity data: "<< tagIntensity.size() << " records have been loaded." << std::endl;
	std::map<unsigned int, double>::iterator itbegin = tagIntensity.begin();
	std::map<unsigned int, double>::iterator itend = tagIntensity.end();
	itend--;
	std::cout << "Tag number is from " << itbegin->first << " to " << itend->first << ". total records should be " << (itend->first-itbegin->first)/6 +1 << std::endl;
}	
//_____Read Intensity region DATA
void MyAnalyzer::OpenIntPartition()
{
	if (!existIntPartition) return;

	std::ifstream ifs("IntensityPartition.txt",std::ios::in);
	if (ifs.fail())
	{
		std::cout<<"Can not open "<<"IntensityPartition.txt"<<std::endl;
		return;
	}
	double doubleBuf;
	char tmp[256];
	while (!ifs.eof())
	{
		//read the data Intensity Region (double)
		ifs >> doubleBuf;
		//go to nextline
		ifs.getline(tmp,256);
		if (!ifs.fail())
		{
			//add to vector
			intPartition.push_back(doubleBuf);
		}
	}
	std::cout << intPartition.size()-1 << " partitions have been set." << std::endl;
}
//_____Read Molecule momentumsum DATA
void MyAnalyzer::OpenMomInfoData()
{
	if (MoleculeAnalysis!=1) return;
	//initialize 2D vector
	molecule.resize(fParticles.GetNbrOfParticles());
	for (int i=0; i<fParticles.GetNbrOfParticles(); ++i)
		molecule[i].resize(fParticles.GetNbrOfParticles());

	std::ifstream ifs(MomSumInfoName,std::ios::in);
	if (ifs.fail())
	{
		std::cout<<"Can not open MomentumInfo.txt. Use (Make) default value."<<std::endl;
		std::ofstream ofs("MomentumInfo.txt",std::ios::out);
		for (int i=1; i<fParticles.GetNbrOfParticles(); ++i)
		{
			for (int j=i+1; j<fParticles.GetNbrOfParticles(); ++j)
			{
				if ((fParticles.GetParticle(i).GetKindParticle() == 1)&&(fParticles.GetParticle(j).GetKindParticle() == 1))
				{
					if (
						((fParticles.GetParticle(i).GetCoinGroup() != fParticles.GetParticle(j).GetCoinGroup()))
						||((fParticles.GetParticle(i).GetCoinGroup()==100)&&(fParticles.GetParticle(j).GetCoinGroup()==100))
						)

					{
						molecule[i][j].momSumWindowX = 100;
						molecule[i][j].momSumWindowY = 100;
						molecule[i][j].momSumWindowZ = 100;
						molecule[i][j].momSumFactor = 1;
						//std::cout<<molecule[i].size() << ":" << molecule[i][j].momSumWindowX<<std::endl;
						string molName(fParticles.GetParticle(i).GetName());
						molName += fParticles.GetParticle(j).GetName();
						ofs << molName; 
						ofs << "\t" << molecule[i][j].momSumWindowX;
						ofs << "\t" << molecule[i][j].momSumWindowY;
						ofs << "\t" << molecule[i][j].momSumWindowZ;
						ofs << "\t" << molecule[i][j].momSumFactor;
						ofs << "\t" << i;
						ofs << "\t" << j;
						ofs << std::endl;
						std::cout << molName << std::endl;
					}
				}
			}
		}
		ofs.close();
		return;
	}
	//-----can open MomentumInfo
	double doubleBuf[6];
	string strBuf;
	string tmp;
	Molecule molBuf;
	map<string,Molecule> bufMap;
	//---Load Momentum sum data from "MomentumInfo.txt"
	while (!ifs.eof())
	{
		//read the data (string, double)
		ifs >> strBuf;
		if (ifs.eof()) break;
		for (int i=0; i<4; ++i)
			ifs >> doubleBuf[i];
		if (ifs.fail())
		{
			std::cout << "MomentumSumInfo: Data read error!" <<std::endl;
			break;
		}
		//--set to buffer structure
		molBuf.momSumWindowX = doubleBuf[0];
		molBuf.momSumWindowY = doubleBuf[1];
		molBuf.momSumWindowZ = doubleBuf[2];
		molBuf.momSumFactor = doubleBuf[3];
		//add to Map
		bufMap.insert(pair<string,Molecule>(strBuf, molBuf));
		//go to nextline
		std::getline(ifs, tmp);
	}
	std::cout << "MomentumSumInfo: "<< bufMap.size() << " records" <<std::endl;
	for (int i=1; i<fParticles.GetNbrOfParticles(); ++i)
	{
		for (int j=i+1; j<fParticles.GetNbrOfParticles(); ++j)
		{
			if ((fParticles.GetParticle(i).GetKindParticle() == 1)&&(fParticles.GetParticle(j).GetKindParticle() == 1))
			{
				if (
					((fParticles.GetParticle(i).GetCoinGroup() != fParticles.GetParticle(j).GetCoinGroup()))
					||((fParticles.GetParticle(i).GetCoinGroup()==100)&&(fParticles.GetParticle(j).GetCoinGroup()==100))
					)
				{
					string molName(fParticles.GetParticle(i).GetName());
					molName += fParticles.GetParticle(j).GetName();
					map<string,Molecule>::iterator it = bufMap.find(molName);
					if (it != bufMap.end())
					{
						molecule[i][j].momSumWindowX = it->second.momSumWindowX;
						molecule[i][j].momSumWindowY = it->second.momSumWindowY;
						molecule[i][j].momSumWindowZ = it->second.momSumWindowZ;
						molecule[i][j].momSumFactor = it->second.momSumFactor;
					}
					else
					{
						std::cout << "MomentumSumInfo: Can not find " << molName << " data!!"<< std::endl;
						std::cout << "MomentumSumInfo: Set window size at 0."<<std::endl;

						molecule[i][j].momSumWindowX = 0;
						molecule[i][j].momSumWindowY = 0;
						molecule[i][j].momSumWindowZ = 0;
						molecule[i][j].momSumFactor = 1;
					}
				}
			}
		}
	}
}

void MyAnalyzer::OpenBeamPositionData()
{
	//if ((intFileName == "")||(!existIntensityData)) return;
	const TString posFileName("BeamPosition.txt");
	std::ifstream ifs(posFileName,std::ios::in);
	if (ifs.fail()){
		std::cout<<"Can not open "<<posFileName<<std::endl;
		return;
	}

	unsigned int uintBuf = 0;
	double doubleBuf1;
	double doubleBuf2;
	char tmp[256];
	while (!ifs.eof())
	{
		//read the data Tag and Intensity (uint/double)
		ifs >> uintBuf >> doubleBuf1 >> doubleBuf2;
		//go to nextline
		ifs.getline(tmp,256);
		if (uintBuf % 6 != 0)
			std::cout<< "wrong Tag number!! "<< uintBuf;
		if (!ifs.fail())
		{
			//add to map (tagIntensity)
			beamPosX.insert(pair<unsigned int, double>(uintBuf,doubleBuf1));
			beamPosY.insert(pair<unsigned int, double>(uintBuf,doubleBuf2));
		}
	}

	std::cout << "Position data: "<< beamPosX.size() << " records have been loaded." << std::endl;
}	