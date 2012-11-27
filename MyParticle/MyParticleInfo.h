#ifndef __MyParticleInfo_H__
#define __MyParticleInfo_H__

#include <TString.h>
#include <TObject.h>
#include <TMath.h>

#include "./MySpectrometer/MySpectrometer.h"

//----------------------------------the particle info class-------------------------------------------------------------------
class MyParticleInfo : public TObject
{
public:
	MyParticleInfo()														{}

public:
	void							 Load();
	void							 Save();
	void							 SetProperties(const char * name, double charge_au, double mass_amu, int kindParticle = 0);

public:
	const char						*GetName()const						{return fName.Data();}
	const MySpectrometer			&GetSpectrometer()const				{return fSp;}
	MySpectrometer					&GetSpectrometer()					{return fSp;}

	double							 GetXcor()const						{return fXcor;}
	double							&GetXcor()							{return fXcor;}
	double							 GetYcor()const						{return fYcor;}
	double							&GetYcor()							{return fYcor;}
	double							 GetSfx()const						{return fSfx;}
	double							&GetSfx()							{return fSfx;}
	double							 GetSfy()const						{return fSfy;}
	double							&GetSfy()							{return fSfy;}
	double							 GetAngle()const					{return fAngle;}
	double							&GetAngle()							{return fAngle;}
	double							 GetT0()const						{return fT0;}
	double							&GetT0()							{return fT0;}

	double							 GetCondTofFr()const				{return fCondTofFr;}
	double							&GetCondTofFr()						{return fCondTofFr;}
	double							 GetCondTofTo()const				{return fCondTofTo;}
	double							&GetCondTofTo()						{return fCondTofTo;}
	double							 GetCondRad()const					{return fCondRad;}
	double							&GetCondRad()						{return fCondRad;}
	double							 GetCondRadX()const					{return fCondRadX;}
	double							&GetCondRadX()						{return fCondRadX;}
	double							 GetCondRadY()const					{return fCondRadY;}
	double							&GetCondRadY()						{return fCondRadY;}
	double							 GetCondWidthX()const				{return fCondWidthX;}
	double							&GetCondWidthX()					{return fCondWidthX;}
	double							 GetCondWidthY()const				{return fCondWidthY;}
	double							&GetCondWidthY()					{return fCondWidthY;}
	bool							 GetPosFlag()const					{return fPosFlag;}
	bool							&GetPosFlag()						{return fPosFlag;}

	double							 GetCharge_au()const				{return fCharge_au;}
	double							 GetMass_amu()const					{return fMass_amu;}
	int								 GetKindParticle()const				{return fKindParticle;}

private:
	//these informations will be in the info class//
	double							 fCondTofFr;						//the lower edge of the time of flight that this particle will be in
	double							 fCondTofTo;						//the upper edge of the time of flight that this particle will be in
	double							 fCondRad;							//the Radius of the circle on the detektor that this particle will be in
	double							 fCondRadX;							//the x-center of the circle on the detektor that this particle will be in
	double							 fCondRadY;							//the y-center of the circle on the detektor that this particle will be in
	double							 fCondWidthX;						//the width x when choosing a quad condition
	double							 fCondWidthY;						//the width y when choosing a quad condition
	bool							 fPosFlag;							//flag to tell wether you want a radius or quad condition

	double							 fAngle;							//the angle to turn the detector in RAD; will be setted in deg but converted to rad while setting it
	double							 fXcor;								//to move the center of the distribution to center of raw detektor
	double							 fYcor;								//to move the center of the distribution to center of raw detektor
	double							 fSfx;								//to scale the x such that it is realy mm
	double							 fSfy;								//to scale the y such that it is realy mm
	double							 fT0;								//to correct the time of flight

	double							 fMass_amu;							//the Mass of this Particle in a.u.
	double							 fCharge_au;						//the Charge of this Particle in a.u.
	TString							 fName;								//how is this particle called

	int								 fKindParticle;						//kind of perticle 0:atom 1:molecule 2:electron 3:Ion 4:Ion(simple) ;motomura

	MySpectrometer					 fSp;								//the Spectrometer Properties this Particle flies through

	ClassDef(MyParticleInfo,1)											//the infos about a particle
};

#endif