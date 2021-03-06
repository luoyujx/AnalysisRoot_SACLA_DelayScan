#ifndef __MySpectrometer_H_
#define __MySpectrometer_H_

#include <vector>

#include "MySpectrometerRegion.h"

typedef std::vector<MySpectrometerRegion> SpecRegions;
class MySpectrometer	//-------------Class describing a Spectrometer-Part ---------------------
{
public:
	MySpectrometer():fMFieldOn(false),fCw(false),fCycFreq(0)								{}

	double				 CyclotronPeriod_ns()const											{return fCycFreq;}
	double				&CyclotronPeriod_ns()												{return fCycFreq;}
	bool				 MagneticFieldIsOn()const											{return fMFieldOn;}
	bool				&MagneticFieldIsOn()												{return fMFieldOn;}
	bool				 RotationClockWise()const											{return fCw;}
	bool				&RotationClockWise()												{return fCw;}

	void				 AddSpectrometerRegion(double Length_mm, double EField_Vpcm)		{fSr.push_back(MySpectrometerRegion(Length_mm,EField_Vpcm));}
	const SpecRegions	&GetSpectrometerRegions() const										{return fSr;}
	SpecRegions			&GetSpectrometerRegions() 											{return fSr;}
	void				 Clear()															{fSr.clear();}

private:
	SpecRegions			 fSr;																//the spectrometer has several regions
	bool				 fMFieldOn;															//flag showing wether a magnetic field was turned on
	double				 fCycFreq;															//the cyclotron frequency of electrons, defines the strength of magnetic field
	bool				 fCw;																//flag showing the direction of the magnetic field (by saying wether electron turn clock or counterclockwise

	ClassDef(MySpectrometer,1)																//the spektrometer for a particle
};
#endif