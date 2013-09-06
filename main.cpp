#include <TRint.h>
#include <stdio.h>
#include <stdlib.h>
#include <TString.h>
#include <fstream>

#include "./MyAnalyzer/MyAnalyzer.h"
#include "./FilesFromLma2Root/MySettings/MySettings.h"

bool LoadSettings(TString &filename, MySettings &set);

int main(int argc, char *argv[])
{
	TRint theApp("App", &argc, argv);
		
	MySettings set(true);
	TString filename("setting.txt");
	LoadSettings(filename,set);

	MyAnalyzer fAn(set);
	fAn.SetFileName(set.GetString("OutputROOTFile","Analysis.root"));

	fAn.FileOpen();
	fAn.Init(set);

	fAn.OpenIntensityData();
	fAn.OpenIntPartition();
	fAn.Open3BodyCombination();
	//test
	//fAn.OpenBeamPositionData();

	std::cout<<"Root File Name : "<<fAn.GetFileName()<<std::endl;
	std::cout<<"Reconstruction Method : "<<fAn.GetRekMeth()<<std::endl;
	std::cout<<"Analyze molecule : "<<fAn.GetMolecule()<<std::endl;
	std::cout<<"Intensity data file name : "<<fAn.GetIntFileName()<<std::endl;

	theApp.Run();
	return 0;
}

bool LoadSettings(TString &filename, MySettings &set)
{
	TString LeftHand;
	TString RightHand;
	std::ifstream ifs(filename,std::ios::in);
	if (ifs.fail())
	{
		std::cout<<"Can not open "<<filename<<std::endl;
		std::cout << "Analyze default settings." << std::endl;
		return false;
	}
	while(!ifs.eof())
	{
		char tmp[256];
		ifs.getline(tmp,256);
		TString param = tmp;
		//std::cout << filename.Data()<<std::endl;
		//--when q return true--//
		if (!param.CompareTo("q"))
			return true;
		
		//--if there is a '#--' its a commented line that should be printed--//
		else if ( param.Contains("#") && param.Contains("--") )
		{
			std::cout <<std::endl<<param.Data()<<std::endl;
		}
		//--if there is a '#' its a commented line--//
		else if (param.Contains("#"))
		{
			//do nothing//
		}
		//--when there is a "=" then check wether a Parameter is changed--//
		else if (param.Contains("="))
		{
			//--separate input into the lefthand and righthand side of the "="--//
			LeftHand = param(0,param.Index("="));
			RightHand = param(param.Index("=")+1,(param.Length()-param.Index("=")));
			if (RightHand.IsFloat())
				set.Change(LeftHand.Data(),RightHand.Atof());
			else 
				set.Change(LeftHand.Data(),RightHand.Data());
		}
		//--else it must have been a filename, so return--//
		else
		{
			return false;
		}
	}
	return true;
}
