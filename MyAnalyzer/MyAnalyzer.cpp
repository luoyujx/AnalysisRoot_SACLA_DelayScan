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
DataBase0d DB;

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
		//Analyze Raw waveform//
		fWf.ExtractWaveform(fOE,fHi,7-1);
		//analyze the event//
		if (optShutMode == 0)
		{
			Analyze(fWf);
		}
		if (optShutMode == 1)
		{
			if(DB.GetStatusAndData(fOE.GetEventID(),2).second == 1) Analyze(fWf); //memo UV-on
			if(DB.GetStatusAndData(fOE.GetEventID(),2).second == 0) Analyze(fWf); //memo UV-off
		}
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
	zeroDTxtFileName = set.GetString("0D_DataTxtFile","Scan.txt");
	MomSumInfoName = set.GetString("MomSumInfoFile","MomentumInfo.txt");
	rekmeth = static_cast<int>(set.GetValue("ReconstructionMethod", 20)+0.1);
	MoleculeAnalysis = static_cast<int>(set.GetValue("Molecule", 0)+0.1);
	extraCondition = static_cast<int>(set.GetValue("ExtraCondition", false)+0.1);
	method0D_Data = static_cast<int>(set.GetValue("0D_Data", false)+0.1);
	pathSQLite = set.GetString("pathSQLite","");
	hostMySQL = set.GetString("hostMySQL","");
	userMySQL = set.GetString("userMySQL","");
	passMySQL = set.GetString("passMySQL","");
	nameMySQL = set.GetString("nameMySQL","");
	runMode = static_cast<int>(set.GetValue("runMode", false)+0.1);
	runNum = static_cast<int>(set.GetValue("runNumber", false)+0.1);
	tagFrom = static_cast<int>(set.GetValue("TagFrom", 0)+0.1);
	tagTo = static_cast<int>(set.GetValue("TagTo", 0)+0.1);
	intfield = set.GetString("intensityFieldName","");
	existIntPartition=static_cast<int>(set.GetValue("IntensityPartition", false)+0.1);
	factorBM1=set.GetValue("ConversionFactorBM1", 10e+9);
	factorPD=set.GetValue("ConversionFactorPD", 10000);
	selectIntensity=static_cast<int>(set.GetValue("SelectIntensity", false)+0.1);
	intensityLowerLimit=set.GetValue("IntensityLowerLimit", 0.0);
	intensityUpperLimit=set.GetValue("IntensityUpperLimit", 100000);
	checkingResult=static_cast<int>(set.GetValue("CheckResult", false)+0.1);
	afterAnalysis=static_cast<int>(set.GetValue("AfterAnalysis", false)+0.1);
	trendStep=static_cast<int>(set.GetValue("TrendStep", 100)+0.1);
	momFactorLowerLimit=set.GetValue("MomFactorLowerLimit", 0.0);
	momFactorUpperLimit=set.GetValue("MomFactorUpperLimit", 2);
	angleCondition=set.GetValue("AngleCondition", 0.0);
	optShutMode = static_cast<int>(set.GetValue("OpticalLaser_OnOff", false)+0.1);
	delayfield = set.GetString("delayFieldName","");;
	factorPMDOffset=set.GetValue("ConversionPMOffset", -1750);
	factorPMD=set.GetValue("ConversionPMtoDelay", 150);
	delayBins=static_cast<int>(set.GetValue("DelayBins", 300));
	delayFrom=set.GetValue("DelayFrom", -10);
	delayTo=set.GetValue("DelayTo", 20);

	limitTheataZ = set.GetValue("LimitOfTheataZ", 90);
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
	if ((zeroDTxtFileName == "")||(method0D_Data==0)) return;
	{
		if (method0D_Data==1)
		{
			std::ifstream ifs(zeroDTxtFileName,std::ios::in);
			if (ifs.fail())
			{
				std::cout<<"Can not open "<<zeroDTxtFileName<<std::endl;
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
					tagDelay.insert(pair<unsigned int, double>(uintBuf,doubleBuf1));
					tagIntensity.insert(pair<unsigned int, double>(uintBuf,doubleBuf2));
				}
			}

			std::cout << "Intensity data: "<< tagDelay.size() << " records have been loaded." << std::endl;
			std::map<unsigned int, double>::iterator itbegin = tagDelay.begin();
			std::map<unsigned int, double>::iterator itend = tagDelay.end();
			itend--;
			std::cout << "Tag number is from " << itbegin->first << " to " << itend->first << ". total records should be " << (itend->first-itbegin->first)/6 +1 << std::endl;
		}

		if (method0D_Data==2)
		{
			DB.Open(pathSQLite);
			vector<string> fields;
			fields.push_back(delayfield);
			fields.push_back(intfield);
			DB.LoadDataL(tagFrom, tagTo, fields);
			//DB.ShowTable();
		}

		if (method0D_Data==3)
		{
			DB.Connect(hostMySQL, userMySQL, passMySQL, nameMySQL); //memo: 最終的にはpath0D_DataBaseMへ書き換え(kuma)
			vector<string> fields;
			fields.push_back(delayfield);
			fields.push_back(intfield);
			if (runMode == 1) 
			{
				std::cout << "Run mode" << std::endl;
				DB.LoadDataM(runNum, fields);
			}
			else
			{	
				std::cout << "Tag mode" << std::endl;
				DB.LoadDataM(tagFrom, tagTo, fields);
			}
			//DB.ShowTable();
		}
	}
}

//void MyAnalyzer::Open0D_Data()
//{
//	if ((zeroDTxtFileName == "")||(method0D_Data==0)) return;
//	
//	std::ifstream ifs(zeroDTxtFileName,std::ios::in);
//	if (ifs.fail()){
//		std::cout<<"Can not open "<<zeroDTxtFileName<<std::endl;
//		return;
//	}
//
//	unsigned int uintBuf = 0;
//	double doubleBuf1;
//	double doubleBuf2;
//	char tmp[256];
//	while (!ifs.eof())
//	{
//		//read the data Tag and Intensity (uint/double)
//		ifs >> uintBuf >> doubleBuf1 >> doubleBuf2;
//		//go to nextline
//		ifs.getline(tmp,256);
//		if (uintBuf % 2 != 0)
//			std::cout<< "wrong Tag number!! "<< uintBuf;
//		if (!ifs.fail())
//		{
//			//add to map (tagDelay)
//			tagDelay.insert(pair<unsigned int, double>(uintBuf,doubleBuf1));
//			tagIntensity.insert(pair<unsigned int, double>(uintBuf,doubleBuf2));
//		}
//	}
//
//	std::cout << "Intensity data: "<< tagDelay.size() << " records have been loaded." << std::endl;
//	std::map<unsigned int, double>::iterator itbegin = tagDelay.begin();
//	std::map<unsigned int, double>::iterator itend = tagDelay.end();
//	itend--;
//	std::cout << "Tag number is from " << itbegin->first << " to " << itend->first << ". total records should be " << (itend->first-itbegin->first)/6 +1 << std::endl;
//}

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