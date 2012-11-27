#include <TMath.h>

#include "MyParticleHit.h"
#include "MyParticle.h"

double calcPx(const MyParticle&, const MyParticleHit&);
double calcPy(const MyParticle&, const MyParticleHit&);
double calcPz(const MyParticle&, const MyParticleHit&);
double calcMass(const MyParticle&, const MyParticleHit&);

MyParticleHit::MyParticleHit(const MyDetektorHit &dh, const MyParticle &p):
	MyDetektorHit(dh)
{
	//calculate the values for this particlehit from the detektorhit//
	fTofC		= fTof - p.GetT0();
	fXCor		= fX_mm - p.GetXcor() - (p.GetXVelocity() * fTofC);
	fYCor		= fY_mm - p.GetYcor() - (p.GetYVelocity() * fTofC);
	fXCorRot	= ( TMath::Cos(p.GetAngle())*fXCor + TMath::Sin(p.GetAngle())*fYCor );
	fYCorRot	= (-TMath::Sin(p.GetAngle())*fXCor + TMath::Cos(p.GetAngle())*fYCor );
	fXCorRotScl	= fXCorRot * p.GetSfx();
	fYCorRotScl	= fYCorRot * p.GetSfy();
	fMassCalc	= calcMass(p,*this);//added by motomura


	//calculate the momenta and energy of this Particle//
	fPx = calcPx(p,*this);
	fPy = calcPy(p,*this);
	fPz = calcPz(p,*this);

	fP  = TMath::Sqrt(fPx*fPx + fPy*fPy + fPz*fPz);
	fE  = (13.6* fP*fP / p.GetMass_au());

	fThetaX = acos(fPx / fP )/TMath::Pi()*180.;//added by motomura
	fThetaY = acos(fPy / fP )/TMath::Pi()*180.;//added by motomura
	fThetaZ = acos(fPz / fP )/TMath::Pi()*180.;//added by motomura

	fPhiXY = atan2(fPy, fPx)/TMath::Pi()*180.;//added by motomura				
	fPhiYZ = atan2(fPz, fPy)/TMath::Pi()*180.;//added by motomura				
	fPhiZX = atan2(fPx, fPz)/TMath::Pi()*180.;//added by motomura				

} 
