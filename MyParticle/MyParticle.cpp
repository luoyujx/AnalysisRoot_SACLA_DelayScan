#include "MyParticle.h"
#include "MyParticleInfo.h"

#include "../FilesFromLma2Root/MyRootManager/MyHistos.h"
#include "../MyMomentaCalculator/MyMomentaCalculator.h"

//___________________________________________________________________________________________________________________________________________________________
void MyParticle::ReadFromInfo(const MyParticleInfo &pi)
{
	fCondTofFr	= pi.GetCondTofFr();
	fCondTofTo	= pi.GetCondTofTo();	
	fCondRad	= pi.GetPosFlag() ? TMath::Max(pi.GetCondWidthX(),pi.GetCondWidthY()) : pi.GetCondRad();
	fCondRadX	= pi.GetCondRadX();
	fCondRadY	= pi.GetCondRadY();
	fCondWidthX	= pi.GetCondWidthX();
	fCondWidthY	= pi.GetCondWidthY();
	fPosFlag	= pi.GetPosFlag();
	fAngle		= pi.GetAngle()*TMath::DegToRad();
	fXcor		= pi.GetXcor();
	fYcor		= pi.GetYcor();
	fSfx		= pi.GetSfx();
	fSfy		= pi.GetSfy();
	fT0			= pi.GetT0();
	fMass_au	= pi.GetMass_amu()*MyUnitsConv::amu2au();
	fCharge_au	= pi.GetCharge_au();
	fName		= pi.GetName();
	fSp			= pi.GetSpectrometer();
	fEnergyFrom = pi.GetEnergyFrom();
	fEnergyTo	= pi.GetEnergyTo();
	fPhiZXFrom	= pi.GetPhiZXFrom();
	fPhiZXTo	= pi.GetPhiZXTo();

	fKindParticle = pi.GetKindParticle();
	fCoinGroup = pi.GetCoinGroup();

}

//___________________________________________________________________________________________________________________________________________________________
const MyParticleHit &MyParticle::AddHit(const MyDetektorHit &dh)
{
	//add a hit to this particle
	fPh.push_back(MyParticleHit(dh,*this));
	return fPh.back();
}


bool MyParticle::CheckPhiZX(const MyParticleHit &ph)
{
	if ((fPhiZXFrom == 0.)&&(fPhiZXTo == 0.)) return true;
	if (!((ph.PhiZX() > fPhiZXFrom && ph.PhiZX() < fPhiZXTo) || (ph.PhiZX() > -fPhiZXTo && ph.PhiZX() < -fPhiZXFrom)))
	{
		fPh.pop_back();
		return false;
	}
	return true;
}
bool MyParticle::CheckPhiZX(const MyParticleHit &ph, double from, double to)
{
	if (!((ph.PhiZX() > from && ph.PhiZX() < to) || (ph.PhiZX() > -to && ph.PhiZX() < -from)))
	{
		fPh.pop_back();
		return false;
	}
	return true;
}
void MyParticle::EraseHit(size_t idx)
{
	//delete hit
	fPh.erase(fPh.begin() + idx);
}