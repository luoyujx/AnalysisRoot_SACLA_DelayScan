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
	processTimer(100),
	molecule(0, std::vector<Molecule>(0))
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
//_____Read Intensity DATA
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
//_____Read Intensity region DATA
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
		ifs.getline(tmp,256);
		if (!ifs.fail())
		{
			//add to vector
			intRegion.push_back(doubleBuf);
		}
	}
	std::cout << intRegion.size()-1 << " region have been loaded." << std::endl;
}
//_____Read Molecule momentumsum DATA
void MyAnalyzer::OpenMoleculeData()
{
	//initialize 2D vector
	molecule.resize(fParticles.GetNbrOfParticles());
	for (int i=0; i<fParticles.GetNbrOfParticles(); ++i)
		molecule[i].resize(fParticles.GetNbrOfParticles());

	std::ifstream ifs("MomentumInfo.txt",std::ios::in);
	if (ifs.fail())
	{
		std::cout<<"Can not open MomentumInfo.txt. Use default value."<<std::endl;
		std::ofstream ofs("MomentumInfo.txt",std::ios::out);
		for (int i=1; i<fParticles.GetNbrOfParticles(); ++i)
			for (int j=i+1; j<fParticles.GetNbrOfParticles(); ++j)
				if ((fParticles.GetParticle(i).GetKindParticle() == 1)&&(fParticles.GetParticle(j).GetKindParticle() == 1))
					if (
						((fParticles.GetParticle(i).GetCoinGroup() != fParticles.GetParticle(j).GetCoinGroup()))
						||((fParticles.GetParticle(i).GetCoinGroup()==100)&&(fParticles.GetParticle(j).GetCoinGroup()==100))
						)

					{
						molecule[i][j].momSumWindowX = 50;
						molecule[i][j].momSumWindowY = 50;
						molecule[i][j].momSumWindowZ = 50;
						molecule[i][j].momSumFactorX = 1;
						molecule[i][j].momSumFactorY = 1;
						molecule[i][j].momSumFactorZ = 1;
						//std::cout<<molecule[i].size() << ":" << molecule[i][j].momSumWindowX<<std::endl;
						string molName(fParticles.GetParticle(i).GetName());
						molName += fParticles.GetParticle(j).GetName();
						ofs << molName; 
						ofs << "\t" << molecule[i][j].momSumWindowX;
						ofs << "\t" << molecule[i][j].momSumWindowY;
						ofs << "\t" << molecule[i][j].momSumWindowZ;
						ofs << "\t" << molecule[i][j].momSumFactorX;
						ofs << "\t" << molecule[i][j].momSumFactorY;
						ofs << "\t" << molecule[i][j].momSumFactorZ;
						ofs << std::endl;
						std::cout << molName << std::endl;
					}
					ofs.close();
					return;
	}

	double doubleBuf[6];
	string strBuf;
	char tmp[256];
	Molecule molBuf;
	map<string,Molecule> bufMap;
	while (!ifs.eof())
	{
		//read the data (double)
		ifs >> strBuf;
		for (int i=0; i<6; ++i)
			ifs >> doubleBuf[i];
		//go to nextline
		ifs.getline(tmp,256);
		if (!ifs.fail())
		{
			molBuf.momSumWindowX = doubleBuf[0];
			molBuf.momSumWindowY = doubleBuf[1];
			molBuf.momSumWindowZ = doubleBuf[2];
			molBuf.momSumFactorX = doubleBuf[3];
			molBuf.momSumFactorY = doubleBuf[4];
			molBuf.momSumFactorZ = doubleBuf[5];
			//add to Map
			bufMap.insert(pair<string,Molecule>(strBuf, molBuf));
		}
		else
		{
			molBuf.momSumWindowX = 50;
			molBuf.momSumWindowY = 50;
			molBuf.momSumWindowZ = 50;
			molBuf.momSumFactorX = 1;
			molBuf.momSumFactorY = 1;
			molBuf.momSumFactorZ = 1;
			bufMap.insert(pair<string,Molecule>(strBuf, molBuf));
		}
	}

	for (int i=1; i<fParticles.GetNbrOfParticles(); ++i)
		for (int j=i+1; j<fParticles.GetNbrOfParticles(); ++j)
			if ((fParticles.GetParticle(i).GetKindParticle() == 1)&&(fParticles.GetParticle(j).GetKindParticle() == 1))
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
						molecule[i][j].momSumFactorX = it->second.momSumFactorX;
						molecule[i][j].momSumFactorY = it->second.momSumFactorY;
						molecule[i][j].momSumFactorZ = it->second.momSumFactorZ;
					}
					else
					{
						molecule[i][j].momSumWindowX = 50;
						molecule[i][j].momSumWindowY = 50;
						molecule[i][j].momSumWindowZ = 50;
						molecule[i][j].momSumFactorX = 1;
						molecule[i][j].momSumFactorY = 1;
						molecule[i][j].momSumFactorZ = 1;
					}
				}
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

		//std::cout << "   *** " << "EventID :"<< static_cast<unsigned int>(fOE.GetEventID()) <<" : "<< static_cast<unsigned int>(fSAE.GetEventID()) <<" : "<< static_cast<unsigned int>(fSE.GetEventID());
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
