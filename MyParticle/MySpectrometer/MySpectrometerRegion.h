#ifndef __MySpectrometerRegion_H_
#define __MySpectrometerRegion_H_

#include <TObject.h>
//-------------Class describing a Spectrometer Region---------------------
class MySpectrometerRegion
{
public:
	MySpectrometerRegion():fL(0),fF(0)	{}
	MySpectrometerRegion(double length_mm, double efield_Vpcm):
		fL(length_mm),
		fF(efield_Vpcm)	{}

	double  EField_Vpcm()const			{return fF;}
	double &EField_Vpcm()				{return fF;}
	double  Length_mm()const			{return fL;}
	double &Length_mm()					{return fL;}

private:
	double  fL;					//the length of the spectrometer region
	double  fF;					//the electric field inside this region

	ClassDef(MySpectrometerRegion,1)		//one part of the spectrometer
};

#endif