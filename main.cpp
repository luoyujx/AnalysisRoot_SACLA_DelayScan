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

	//------for Setting-------//
	MyAnalyzer fAn(static_cast<int>(set.GetValue("UseGUI", true)+0.1));
	fAn.SetRekMeth(static_cast<int>(set.GetValue("ReconstructionMethod", 20)+0.1));
	fAn.SetMolecule(static_cast<int>(set.GetValue("Molecule", false)+0.1));
	fAn.SetCondition(static_cast<int>(set.GetValue("ExtraCondition", false)+0.1));
	fAn.SetFileName(set.GetString("OutputROOTFile","Analysis.root"));
	fAn.SetIntFileName(set.GetString("IntensityFile","Intensity.txt"));

	std::cout<<"Root File Name : "<<fAn.GetFileName()<<std::endl;
	std::cout<<"Reconstruction Method : "<<fAn.GetRekMeth()<<std::endl;
	std::cout<<"Analyze molecule : "<<fAn.GetMolecule()<<std::endl;
	std::cout<<"Extra condition : "<<fAn.GetCondition()<<std::endl;
	std::cout<<"Intensity data file name : "<<fAn.GetIntFileName()<<std::endl;
	//------------------------//

	fAn.FileOpen();
	fAn.Init(set);

	fAn.OpenIntensityData();
	fAn.OpenIntRegionData();
	fAn.OpenMoleculeData();

	theApp.Run();
	return 0;
}

bool LoadSettings(TString &filename, MySettings &set)
{
	TString LeftHand;
	TString RightHand;
	std::ifstream ifs(filename,std::ios::in);
	if (ifs.fail()){
		std::cout<<"Can not open "<<filename<<std::endl;
		std::cout << "Analyze default settings." << std::endl;
		return false;
	}
	while(!ifs.eof())
	{
		char tmp[128];
		ifs >> tmp;
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
