#include <TRint.h>
#include <stdio.h>
#include <stdlib.h>
#include <TString.h>

#include "./MyAnalyzer/MyAnalyzer.h"

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

	fAn.FileOpen();

	fAn.Init();

	fAn.OpenIntensityData();
//	fAn.OpenIntRegionData();
	fAn.OpenMoleculeData();

	theApp.Run();
	theApp.Terminate(0);
	return 0;
}