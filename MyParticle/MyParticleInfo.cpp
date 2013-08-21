#include <TMath.h>
#include <TEnv.h>

#include "MyParticleInfo.h"

//___________________________________________________________________________________________________________________________________________________________
void MyParticleInfo::SetProperties(const char * name, double charge_au, double mass_amu,  int kindParticle, int coinGroup)
{
	//set the infos provided by the user//
	fName		= name;
	fCharge_au	= charge_au;
	fMass_amu	= mass_amu;
	fKindParticle = kindParticle;
	fCoinGroup = coinGroup;
	//load the other settings from the ini file//
	Load();
}
//___________________________________________________________________________________________________________________________________________________________
void MyParticleInfo::Load()
{
	TEnv e(Form("%s.ini",fName.Data()));
	//correction Values//
	fXcor		= e.GetValue("Corrections.XCorrection_mm",0.);
	fYcor		= e.GetValue("Corrections.YCorrection_mm",0.);
	fSfx		= e.GetValue("Corrections.XScaleFactor",1.);
	fSfy		= e.GetValue("Corrections.YScaleFactor",1.);
	fAngle		= e.GetValue("Corrections.Angle_deg",0.);
	fT0			= e.GetValue("Corrections.T0_ns",0.);
	//condition Values//
	fCondTofFr	= e.GetValue("Conditions.TofFrom_ns",0.);
	fCondTofTo	= e.GetValue("Conditions.TofTo_ns",0.);
	fCondRadX	= e.GetValue("Conditions.CenterOfRadiusX_mm",0.);
	fCondRadY	= e.GetValue("Conditions.CenterOfRadiusY_mm",0.);
	fCondRad	= e.GetValue("Conditions.Radius_mm",0.);
	fCondWidthX	= e.GetValue("Conditions.WidthOfX_mm",0.);
	fCondWidthY	= e.GetValue("Conditions.WidthOfY_mm",0.);
	fPosFlag	= e.GetValue("Conditions.UseQuadForPosCondition",false);
	fEnergyFrom = e.GetValue("Conditions.EnegyFrom",0.);
	fEnergyTo	= e.GetValue("Conditions.EnegyTo",10000);

	//Spectrometer//
	fSp.Clear();
	fSp.CyclotronPeriod_ns()	= e.GetValue("Spectrometer.CyclotronFrequency_ns",0.);
	fSp.MagneticFieldIsOn()		= e.GetValue("Spectrometer.MagneticFieldIsOn",false);
	fSp.RotationClockWise()		= e.GetValue("Spectrometer.RotationIsClockwise",true);
	int nSpecRegions			= e.GetValue("Spectrometer.NbrOfSpectrometerRegions",1);
	for (int i=0;i<nSpecRegions;++i)
	{
		double l = e.GetValue(Form("Spectrometer.Region%d.Length_mm",i+1),10.);
		double f = e.GetValue(Form("Spectrometer.Region%d.EField_VPercm",i+1),10.);
		fSp.AddSpectrometerRegion(l,f);
	}

}

//___________________________________________________________________________________________________________________________________________________________
void MyParticleInfo::Save()
{
	TEnv e(Form("%s.ini",fName.Data()));
	//correction Values//
	e.SetValue("Corrections.XCorrection_mm",fXcor);
	e.SetValue("Corrections.YCorrection_mm",fYcor);
	e.SetValue("Corrections.XScaleFactor",fSfx);
	e.SetValue("Corrections.YScaleFactor",fSfy);
	e.SetValue("Corrections.Angle_deg",fAngle);
	e.SetValue("Corrections.T0_ns",fT0);
	//condition Values//
	e.SetValue("Conditions.TofFrom_ns",fCondTofFr);
	e.SetValue("Conditions.TofTo_ns",fCondTofTo);
	e.SetValue("Conditions.CenterOfRadiusX_mm",fCondRadX);
	e.SetValue("Conditions.CenterOfRadiusY_mm",fCondRadY);
	e.SetValue("Conditions.Radius_mm",fCondRad);
	e.SetValue("Conditions.WidthOfX_mm",fCondWidthX);
	e.SetValue("Conditions.WidthOfY_mm",fCondWidthY);
	e.SetValue("Conditions.UseQuadForPosCondition",fPosFlag);
	//Spectrometer//
	e.SetValue("Spectrometer.CyclotronFrequency_ns",fSp.CyclotronPeriod_ns());
	e.SetValue("Spectrometer.MagneticFieldIsOn",fSp.MagneticFieldIsOn());
	e.SetValue("Spectrometer.RotationIsClockwise",fSp.RotationClockWise());
	//SpectrometerRegions//
	SpecRegions sr = fSp.GetSpectrometerRegions();
	e.SetValue("Spectrometer.NbrOfSpectrometerRegions",static_cast<Int_t>(sr.size()));
	for (size_t i=0;i<sr.size();++i)
	{
		e.SetValue(Form("Spectrometer.Region%d.Length_mm",i+1),sr[i].Length_mm());
		e.SetValue(Form("Spectrometer.Region%d.EField_VPercm",i+1),sr[i].EField_Vpcm());
	}

	//save the settings to file//
	e.WriteFile(Form("%s.ini",fName.Data()));
}
