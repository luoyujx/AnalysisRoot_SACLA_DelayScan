#include <iostream>
#include <iomanip>
#include <fstream>
#include <TTree.h>
#include <TFile.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TGClient.h>
#include <utility>

#include "MyAnalyzer.h"
#include "../MyGui/MyGui.h"
#include "../FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalChannel/MyOriginalChannel.h"
#include "../FilesFromLma2Root/MyEvent/MySortedEvent/MyDetektor/MyDetektor.h"
#include "../FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalEventInfo.h"
#include "../FilesFromLma2Root/MyEvent/MySortedEvent/MySortedEventInfo.h"
#include "../FilesFromLma2Root/MyEvent/MySignalAnalyzedEvent/MySignalAnalyzedEventInfo.h"

#include "../AnalyzeFunktions.h"

//#include "../MyCovariance.h"
//#include "../MyWaveform.h"

using namespace std;

//##################################################################################
//_____________________________The class Members______________________________________________________________________________________________________________________________
MyAnalyzer::MyAnalyzer(int UseGUI):
	fOChain("OriginalEvent"),
	fSAChain("SignalAnalyzedEvent"),
	fSChain("SortedEvent"),
	fHi(false),
	running(false),
	processTimer(100)
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

	DefineParticlesAndRootFile(fParticles,fHi);

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
		runTimer.Start(100);
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
	fWf.Init(fOE,fHi);
}
//Read Intensity DATA
void MyAnalyzer::OpenIntensityData()
{
	if (intFileName == "") return;
	
	std::ifstream ifs(intFileName,std::ios::in);
	if (ifs.fail()){
		std::cout<<"Can not open "<<intFileName<<std::endl;
		return;
	}

	unsigned int uintBuf;
	double doubleBuf1;
	double doubleBuf2;
	char tmp[256];
	while (!ifs.eof())
	{
		//read the data Tag and Intensity (uint/double)
		ifs >> uintBuf >> doubleBuf1 >> doubleBuf2;
		//go to nextline
		ifs.getline(tmp,256);
		if (uintBuf % 6 != 0) std::cout<< "wrong Tag number!! "<< uintBuf;
		if (!ifs.fail())
		{
			//add to map (tagIntensity)
			tagIntensity.insert(pair<unsigned int, double>(uintBuf,doubleBuf1));
			tagIntensity2.insert(pair<unsigned int, double>(uintBuf,doubleBuf2));
		}
	}
	std::cout << tagIntensity.size() << " records have been loaded." << std::endl;
	std::map<unsigned int, double>::iterator itbegin = tagIntensity.begin();
	std::map<unsigned int, double>::iterator itend = tagIntensity.end();
	itend--;
	std::cout << "Tag number is from " << itbegin->first << " to " << itend->first << ". total records should be " << (itend->first-itbegin->first)/6 +1 << std::endl;
}
//Read Intensity region DATA
void MyAnalyzer::OpenIntRegionData()
{
	std::ifstream ifs("IntensityRegion.txt",std::ios::in);
	if (ifs.fail()){
		std::cout<<"Can not open "<<"IntensityRegion.txt"<<std::endl;
		return;
	}
	double doubleBuf;
	char tmp[256];
	while (!ifs.eof())
	{
		//read the data Intensity Region (double)
		ifs >> doubleBuf;
		//go to nextline
		if (!ifs.fail())
		{
			//add to vector
			intRegion.push_back(doubleBuf);
		}
	}
	std::cout << intRegion.size()-1 << " region have been loaded." << std::endl;
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

	while(fEntryIterator < fNEntries)
	{
		std::cout << "\r" << "Entry Number :"<< std::setw(7) << std::setfill(' ') << fEntryIterator;

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

		std::cout << "   *** " << "EventID :"<< static_cast<unsigned int>(fOE.GetEventID()) <<" : "<< static_cast<unsigned int>(fSAE.GetEventID()) <<" : "<< static_cast<unsigned int>(fSE.GetEventID());
		//check EventID//
		if ((fOE.GetEventID()!=fSAE.GetEventID())||(fOE.GetEventID()!=fSE.GetEventID())) 
		{
			std::cout << std::endl << "******Error!!! EventID Mismatch !!!! ******"<<std::endl;
			realyBreak=true;
			break;
		}
		
		//analyze the event//
		Analyze();
		//calc covariance map//
		fWf.ExtractWaveform(fOE,fHi,7-1);
		//increase the counter//
		fEntryIterator++;
		//if(fEntryIterator > 1000) {std::cout << "user requested break"<<std::endl;realyBreak=true;break;}
		//the timer will only process events when it has timed out//
		if (processTimer.ProcessEvents()) {std::cout << "user requested break"<<std::endl;realyBreak=true;break;}
	}
	if (WasRunningBefore)
	{
		fillHistosAfterAnalyzis(fParticles.GetParticles(),fHi,intRegion.size()-1);
		fWf.FillHist(fHi);
		std::cout << "<- Done, now saving Histograms!!!!"<<std::endl;
		fHi.FlushRootFile();
	}
	
	//restart run at this time//
	if (!realyBreak) runTimer.Start(1000);
	//std::cout << "leaving run"<<std::endl;
}
