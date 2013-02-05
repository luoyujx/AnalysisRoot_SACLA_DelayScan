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

	int UseGUI = 1;
	if ( argc > 1 )
	for (size_t i = 1; i <= argc; i++) 
	{
		TString progArg(argv[i]);
		if (progArg(0,3) == "/g:")
		{
			TString ComandArg(progArg(3,progArg.Length()));
			UseGUI = ComandArg.Atoi();
		}
	}	
		
	MyAnalyzer fAn(UseGUI);
	MySettings set(false);

	//------for Setting-------//
	fAn.SetRekMeth(20);
	fAn.SetMolecule(false);
	fAn.SetCondition(false);

	if ( argc > 1 )
	for (size_t i = 1; i <= argc; i++) 
	{
		TString progArg(argv[i]);
		//std::cout<<progArg(0,3)<<std::endl;

		if (progArg(0,3) == "/r:")
		{
			TString ComandArg(progArg(3,progArg.Length()));
			fAn.SetRekMeth(ComandArg.Atoi());
		}
		else if (progArg(0,3) == "/f:")
		{
			TString ComandArg(progArg(3,progArg.Length()));
			fAn.SetFileName(ComandArg);
		}
		else if (progArg(0,3) == "/m:")
		{
			TString ComandArg(progArg(3,progArg.Length()));
			fAn.SetMolecule(ComandArg.Atoi());
		}

		else if (progArg(0,3) == "/e:")
		{
			TString ComandArg(progArg(3,progArg.Length()));
			fAn.SetCondition(ComandArg.Atoi());
		}
		else if (progArg(0,3) == "/i:")
		{
			TString ComandArg(progArg(3,progArg.Length()));
			fAn.SetIntFileName(ComandArg);
		}

	}

	std::cout<<"Root File Name : "<<fAn.GetFileName()<<std::endl;
	std::cout<<"Reconstruction Method : "<<fAn.GetRekMeth()<<std::endl;
	std::cout<<"Analyze molecule : "<<fAn.GetMolecule()<<std::endl;
	std::cout<<"Extra condition : "<<fAn.GetCondition()<<std::endl;
	std::cout<<"Intensity data file name : "<<fAn.GetIntFileName()<<std::endl;
	std::cout<<"Use GUI : "<<UseGUI<<std::endl;
	//------------------------//
	TString filename("setting.txt");
	LoadSettings(filename,set);
	fAn.FileOpen();
	fAn.Init(set);

	fAn.OpenIntensityData();
//	fAn.OpenIntRegionData();
	fAn.OpenMoleculeData();

	theApp.Run();
	theApp.Terminate(0);
	return 0;
}

bool LoadSettings(TString &filename, MySettings &set)
{
	TString LeftHand;
	TString RightHand;
	std::ifstream ifs(filename,std::ios::in);
	if (ifs.fail()){
		std::cout<<"Can not open "<<filename<<std::endl;
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
