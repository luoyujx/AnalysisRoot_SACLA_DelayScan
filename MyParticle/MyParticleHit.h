#ifndef __MyParticleHit_H__
#define __MyParticleHit_H__

#include "../FilesFromLma2Root/MyEvent/MySortedEvent/MyDetektor/MyDetektorHit.h"

class MyHistos;
class MyParticle;

class MyParticleHit : public MyDetektorHit
{
public:
	MyParticleHit()					{}
	MyParticleHit(const MyDetektorHit&, const MyParticle&);


public:
	double		XCor()const			{return fXCor;}
	double		YCor()const			{return fYCor;}
	double		XCorRot()const		{return fXCorRot;}
	double		YCorRot()const		{return fYCorRot;}
	double		XCorRotScl()const	{return fXCorRotScl;}
	double		YCorRotScl()const	{return fYCorRotScl;}
	double		TofCor()const		{return fTofC;}
	double		Mass()const			{return fMassCalc;}//added by motomura


	double		Px()const			{return fPx;}
	double		Py()const			{return fPy;}
	double		Pz()const			{return fPz;}
	double		P()const			{return fP;}
	double		E()const			{return fE;}
	double		ThetaX()const		{return fThetaX;}//added by motomura
	double		ThetaY()const		{return fThetaY;}//added by motomura
	double		ThetaZ()const		{return fThetaZ;}//added by motomura
	double		PhiXY()const		{return fPhiXY;}//added by motomura
	double		PhiYZ()const		{return fPhiYZ;}//added by motomura
	double		PhiZX()const		{return fPhiZX;}//added by motomura

private:
	double		fPx;				//the momentum in x
	double		fPy;				//the momentum in y
	double		fPz;				//the momentum in z
	double		fP;					//the total momentum
	double		fE;					//the energy
	double		fXCor;				//corrected x position on the detektor
	double		fYCor;				//corrected y position on the detektor
	double		fXCorRot;			//rotated and corrected x position on the detektor
	double		fYCorRot;			//rotated and corrected y position on the detektor
	double		fXCorRotScl;		//rotated, corrected and scaled x position on the detektor
	double		fYCorRotScl;		//rotated, corrected and scaled y position on the detektor
	double		fTofC;				//the corrected time of flight
	double		fMassCalc;			//added by motomura
	double		fThetaX;			//added by motomura
	double		fThetaY;			//added by motomura
	double		fThetaZ;			//added by motomura
	double		fPhiXY;				//added by motomura				
	double		fPhiYZ;				//added by motomura				
	double		fPhiZX;				//added by motomura				


	ClassDef(MyParticleHit,1)		//a hit of a particle on the detektor
};

#endif