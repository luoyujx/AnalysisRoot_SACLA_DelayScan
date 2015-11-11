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
#include "../MyMomentaCalculator/MyMomentaCalculator.h"

//#include "../MyCovariance.h"
//#include "../MyWaveform.h"

using namespace std;

//Global object 
DataBase0d DB0d;
DataBase0d DBTM;

//##################################################################################
//_____________________________The class Members______________________________________________________________________________________________________________________________
MyAnalyzer::MyAnalyzer(MySettings &set):
	fOChain("OriginalEvent"),
	fSAChain("SignalAnalyzedEvent"),
	fSChain("SortedEvent"),
	fHi(false,10000),//Max 10000?
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

	whichParticles = set.GetString("WhichParticles","");
	DefineParticlesAndRootFile(fParticles,fHi,whichParticles);
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

	//Init Raw waveform stuff//
	fWf.Init(fOE,fHi);

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
		if (fEntryIterator % 1000 == 0) std::cout << "\r" << "Entry Number :"<< std::setw(7) << std::setfill(' ') << fEntryIterator;
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
		//Analyze Raw waveform//
		fWf.ExtractWaveform(fOE,fHi,7-1);
		//analyze the event//
		Analyze(fWf);
		//increase the counter//
		fEntryIterator++;
		//if(fEntryIterator > 1000) {std::cout << "user requested break"<<std::endl;realyBreak=true;break;}
		//the timer will only process events when it has timed out//
		if (processTimer.ProcessEvents()) {std::cout << "user requested break"<<std::endl;realyBreak=true;break;}
	}
	if (WasRunningBefore)
	{
		if (afterAnalysis) fillHistosAfterAnalyzis(fParticles.GetParticles(),fHi,intPartition.size()-1,delayBins,delayFrom,delayTo);
		fWf.FillHist(fHi);
		std::cout << "<- Done, now saving Histograms!!!!"<<std::endl;
		std::cout << "First TAG: "<<firstTAG << " Last Tag: " << fOE.GetEventID()<< std::endl;
		if (missedTagCount) std::cout << "Can not find "<< missedTagCount << " intensity data." << std::endl;
		fHi.FlushRootFile();
		//if (checkingResult) ShowResult();
	}
	//restart run at this time//
	if (!realyBreak) runTimer.Start(3000);
	//std::cout << "leaving run"<<std::endl;
}
//_____________
void MyAnalyzer::SetParameter(MySettings &set)
{
	//set parameters
	MomSumInfoName = set.GetString("MomSumInfoFile","MomentumInfo.txt");
	//Should be 18 or 21 (resort method)
	rekmeth = static_cast<int>(set.GetValue("ReconstructionMethod", 20)+0.1);
	//The way for extracting coincidence(1:Momentum sum 2 : gated by particle)
	MoleculeAnalysis = static_cast<int>(set.GetValue("Molecule", 0)+0.1);
	//Investigate proton(0:none 1 : all 2 : only)
	extraCondition = static_cast<int>(set.GetValue("ExtraCondition", false)+0.1);
	//Information for the gate of momentum sums(only need if Molecule = 1)
	MomSumInfoName = set.GetString("MomSumInfoFile", "MomentumInfo.txt");
	//Coincidence condition
	angleCondition = set.GetValue("AngleCondition", 0.0);
	momFactorLowerLimit = set.GetValue("MomFactorLowerLimit", 0.0);
	momFactorUpperLimit = set.GetValue("MomFactorUpperLimit", 2);
	// delay scan mode (0: no 0d data, 1: use MySQL, 2:with Timing monitor)
	delayScan = static_cast<int>(set.GetValue("delayScan", false)+0.1);
	if (delayScan == 0)
	{
		std::cout << "Fixed delay. 0D databese is not used." << std::endl;
		std::cout << "First Tag: " << tagFrom << " Last Tag: " << tagTo << std::endl;
	}
	// for delay scan
	else
	{
		// MySQL information
		hostMySQL = set.GetString("hostMySQL", "");
		userMySQL = set.GetString("userMySQL", "");
		passMySQL = set.GetString("passMySQL", "");
		nameMySQL = set.GetString("nameMySQL", "");
	}
	// Table & Field names for 0d Database 
	tableBL = set.GetString("tableBL", "");
	BM1FN = set.GetString("BM1FieldName", ""); //0
	delayFN = set.GetString("delayFieldName", ""); //1
	// Table & Field names for Timing moniter
	if (delayScan == 2)
	{
		std::cout << "Delay scan with jitter extracted from Timing moniter" << std::endl;
		tableTM = set.GetString("tableTM", ""); 
		flagTMFN = set.GetString("timingValidName", ""); //0
		jitterFN = set.GetString("jitterFieldName", ""); //1			
		delayTMFN = set.GetString("timingMoniterDelayName", ""); //2
	}
	else
	{
		std::cout << "Delay scan without jitter extracted from Timing moniter" << std::endl;
	}
	//Tag from, Tag to
	tagFrom = static_cast<int>(set.GetValue("TagFrom", 0)+0.1);
	tagTo = static_cast<int>(set.GetValue("TagTo", 0)+0.1);
	//Intensity informaion
	factorBM1 = set.GetValue("ConversionFactorBM1", 10000);
	selectIntensity = static_cast<int>(set.GetValue("SelectIntensity", false) + 0.1);
	intensityLowerLimit = set.GetValue("IntensityLowerLimit", 0.0);
	intensityUpperLimit = set.GetValue("IntensityUpperLimit", 100000);
	//Delay informatyion
	factorPMD = set.GetValue("ConversionPMtoDelay", 150);
	factorPMDOffset = set.GetValue("PMOffset", 0);
	factorTM = set.GetValue("ConversionPIXtoJitter", 3.8);
	factorTMOffset = set.GetValue("TimingMoniterOffset", 1200);
	delayBins = static_cast<int>(set.GetValue("DelayBins", 80));
	delayFrom = set.GetValue("DelayFrom", -4000);
	delayTo = set.GetValue("DelayTo", 4000);
	afterAnalysis = static_cast<int>(set.GetValue("AfterAnalysis", false)+0.1);
	trendStep = static_cast<int>(set.GetValue("TrendStep", 100)+0.1);
	limitOfThetaZ = static_cast<int>(set.GetValue("limitOfThetaZ", 180) + 0.1);
	//
	lowerDelay = static_cast<int>(set.GetValue("lowerDelay", -1000) + 0.1);
	upperDelay = static_cast<int>(set.GetValue("upperDelay", 1000) + 0.1);
}
//__________Show mass &ToF spactrum with Particle name_________________________________________
void MyAnalyzer::ShowResult()
{
	canv = new TCanvas("Result","Result",100,100,1000+4,600+28);
	//gStyle->SetOptStat(0);
	//gStyle->SetOptFit(0);
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
		//txtMass[i].DrawText(massPerQ, height*1.8, fParticles.GetParticle(i).GetName());
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
		//txtTof[i].DrawText(tofPos, height*1.8, fParticles.GetParticle(i).GetName());
		boxTof[i].SetFillColor(kRed);
		boxTof[i].SetFillStyle(3001);
		boxTof[i].DrawBox(fParticles.GetParticle(i).GetCondTofFr(), 0,  fParticles.GetParticle(i).GetCondTofTo(), height *1.6);
	}
}

//
//_____Read Intensity DATA
void MyAnalyzer::OpenIntensityData()
{
	std::cout << "Method of 0D Data: " << delayScan << std::endl;
	if (delayScan == 0)
	{
		//std::cout << "Fixed delay. 0D databese is not used." << std::endl;
		return;
	}
	else if (delayScan == 1) // without Jitter from Timing Moniter
	{
		std::cout << "Delay scan without jitter extracted from Timing moniter" << std::endl;
		DB0d.Connect(hostMySQL, userMySQL, passMySQL, nameMySQL);
		vector<string> fields0d;
		fields0d.push_back(BM1FN);		//0
		fields0d.push_back(delayFN);	//1
		DB0d.LoadDataM(tagFrom, tagTo, fields0d, tableBL);
		DB0d.CloseMySQL();
		//DB0d.ShowTable();
	}
	else if (delayScan == 2) // with Jitter form Timing Moniter
	{
		std::cout << "Delay scan with jitter extracted from Timing moniter" << std::endl;
		DB0d.Connect(hostMySQL, userMySQL, passMySQL, nameMySQL);
		vector<string> fields0d;
		fields0d.push_back(BM1FN);		//0
		fields0d.push_back(delayFN);	//1
		fields0d.push_back(delayTMFN);	//2
		DB0d.LoadDataM(tagFrom, tagTo, fields0d, tableBL);
		DB0d.CloseMySQL();

		DBTM.Connect(hostMySQL, userMySQL, passMySQL, nameMySQL);
		vector<string> fieldsTM;
		fieldsTM.push_back(flagTMFN);	//0
		fieldsTM.push_back(jitterFN);	//1
		DBTM.LoadDataM(tagFrom, tagTo, fieldsTM, tableTM);
		DBTM.CloseMySQL();

		std::cout << "OK" << std::endl;
		//DBTM.ShowTable();
	}
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
						molecule[i][j].angleCondition = angleCondition;
						molecule[i][j].momSumFactorLow = momFactorLowerLimit;
						molecule[i][j].momSumFactorUp = momFactorUpperLimit;
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
	//if ((zeroDTxtFileName == "")||(method0D_Data==0)) return;
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
		if (uintBuf % 2 != 0)
			std::cout<< "wrong Tag number!! "<< uintBuf;
		if (!ifs.fail())
		{
			//add to map (tagDelay)
			beamPosX.insert(pair<unsigned int, double>(uintBuf,doubleBuf1));
			beamPosY.insert(pair<unsigned int, double>(uintBuf,doubleBuf2));
		}
	}

	std::cout << "Position data: "<< beamPosX.size() << " records have been loaded." << std::endl;
}
//_____Open 3-body combination data
void MyAnalyzer::Open3BodyCombination()
{
	//if ((zeroDTxtFileName == "")||(method0D_Data==0)) return;
	const TString posFileName("3bodyCombination.txt");
	std::ifstream ifs(posFileName,std::ios::in);
	if (ifs.fail())
	{
		std::cout<<"Can not open "<<posFileName<<std::endl;
		return;
	}

	std::string strBuff;
	char tmp[256];
	while (!ifs.eof())
	{
		//read the data
		ifs >> strBuff;
		//go to nextline
		ifs.getline(tmp,256);
		if (!ifs.fail())
		{
			//add to vector
			threeBodyComb.push_back(strBuff);
		}
	}

	std::cout << "3-body combination data: "<< threeBodyComb.size() << " records have been loaded." << std::endl;
}

void MyAnalyzer::OpenMCPToFRegion()
{
	//if ((zeroDTxtFileName == "")||(method0D_Data == 0)) return;
	const TString posFileName("MCPToFRegion.txt");
	std::ifstream ifs(posFileName,std::ios::in);
	if (ifs.fail())
	{
		std::cout<<"Can not open "<<posFileName<<std::endl;
		return;
	}

	TString strBuff;
	double doubleBuff1;
	double doubleBuff2;
	char tmp[256];
	while (!ifs.eof())
	{
		//read the data
		ifs >> strBuff >> doubleBuff1 >> doubleBuff2;
		//go to nextline
		ifs.getline(tmp,256);
		if (!ifs.fail())
		{
			//add to vector
			//MCPToFRegion mtr (strBuff,doubleBuff1, doubleBuff2);
			//mtr.particleName = strBuff;
			//mtr.tofFrom = doubleBuff1;
			//mtr.tofTo = doubleBuff2;
			mcpTofRegion.push_back(MCPToFRegion(strBuff,doubleBuff1, doubleBuff2));
		}
	}

	std::cout << "ToF (MCP) region data: "<< mcpTofRegion.size() << " records have been loaded." << std::endl;
}