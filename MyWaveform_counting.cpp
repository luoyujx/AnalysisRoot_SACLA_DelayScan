#include <algorithm>
#include "TMath.h"

#include "MyWaveform.h"

//#include "../lma2Root/MyEvent/MyOriginalEvent/MyOriginalEvent.h"
#include "../lma2Root/MyRootManager/MyHistos.h"
#include "../lma2Root/MyEvent/MyOriginalEvent/MyOriginalChannel/MyOriginalChannel.h"
#include "../lma2Root/MyEvent/MySortedEvent/MyDetektor/MyDetektor.h"
#include "../lma2Root/MyEvent/MyOriginalEvent/MyOriginalEventInfo.h"
#include "../lma2Root/MyEvent/MySortedEvent/MySortedEventInfo.h"
#include "../lma2Root/MyEvent/MySignalAnalyzedEvent/MySignalAnalyzedEventInfo.h"

MyWaveform::MyWaveform():IDOffset(1000), eventCounter(0),aveIntensity(0.),varIntensity(0.)
{
}
MyWaveform::~MyWaveform()
{
}
void MyWaveform::Init(const MyOriginalEvent &oe, MyHistos &rm)
{
	length= oe.GetNbrSamples();
	arraySize = 500;
	rebin = 5;
	TimeRangeFrom = 2500;

	waveform.resize(length,0.);
	aveWaveform.resize(length,0.);
	rebinWaveform.resize(arraySize,0.);
	aveRebinWaveform.resize(arraySize,0.);

	// make histgram ---ToF, Intensity,,,,
	rm.create1d(IDOffset+1,"Tof","ToF",length,0,length);
	rm.create1d(IDOffset+2,"Tof_rebin","ToF",arraySize,TimeRangeFrom,TimeRangeFrom+arraySize*rebin);
	rm.create1d(IDOffset+10,"FEL_intensity","pulse energy (microJ)",200,0,50);

}
void MyWaveform::Clear()
{
	std::fill(waveform.begin(),waveform.end(),0.);
	std::fill(rebinWaveform.begin(),rebinWaveform.end(),0.);
	std::fill(aveWaveform.begin(),aveWaveform.end(),0.);
	std::fill(aveRebinWaveform.begin(),aveRebinWaveform.end(),0.);
}

void MyWaveform::ExtractWaveform(const MyOriginalEvent &fOE, const MySortedEvent &fSE, MyHistos &fHi)
{
	const MyDetektor &rd = fSE.GetDetektor(0);


	//--Filling the waveform array--//	
	std::fill(waveform.begin(),waveform.end(),0.);


	for (size_t i=0; i<rd.GetNbrOfHits();++i)
	{
		const MyDetektorHit &dh = rd.GetHit(i);

		//the tof is just the timing of the mcp signal//
		//dh.SetTof(dh.Time());
		waveform[static_cast<long>(dh.Time())]++;
		
		const double maxTof	= fOE.GetNbrSamples()*fOE.GetSampleInterval()*1e9;

		fHi.fill(IDOffset+100,"TofAll",dh.Time(),"tof [ns]",10000,0,maxTof);
		//if (dh.RekMeth() < 15) fHi.fill(IDOffset+101,"Tof_Rekmeth15",dh.Tof(),"tof [ns]",10000,0,maxTof);
		//if (dh.RekMeth() < 18) fHi.fill(IDOffset+102,"Tof_Rekmeth18",dh.Tof(),"tof [ns]",10000,0,maxTof);
		//if (dh.RekMeth() < 20) fHi.fill(IDOffset+103,"Tof_Rekmeth20",dh.Tof(),"tof [ns]",10000,0,maxTof);

	}
	std::fill(rebinWaveform.begin(),rebinWaveform.end(),0.);
	for (size_t i=0; i<arraySize; ++i)
		for (size_t j=0; j<rebin; ++j)
			rebinWaveform[i] += waveform[TimeRangeFrom+i*rebin+j];//fill rebined waveform

	eventCounter++;
	const double scale = 1./eventCounter;

	std::transform(waveform.begin(),waveform.end(),
		aveWaveform.begin(),aveWaveform.begin(),
		Average(scale));

	std::transform(rebinWaveform.begin(),rebinWaveform.end(),
		aveRebinWaveform.begin(),aveRebinWaveform.begin(),
		Average(scale));
}

void MyWaveform::ExtractIntensity(const MyOriginalEvent &oe, MyHistos &fHi, int chan)
{
	intensity = CalcGMDIntensity(oe.GetChannel(chan),950,14000,false)
	-(CalcGMDIntensity(oe.GetChannel(chan),200,700,false)
	+CalcGMDIntensity(oe.GetChannel(chan),200,700,false))/2.;
	intensity *= 0.05717;

	fHi.fill1d(IDOffset+10,intensity);

	const double scale = 1./eventCounter;
	const double preAveIntensity = aveIntensity;

	Average ave(scale); 
	aveIntensity = ave(intensity,aveIntensity);

	varIntensity = ((varIntensity*(eventCounter-1)) + (intensity-preAveIntensity)* (intensity-aveIntensity))*scale;
}

void MyWaveform::FillHist(MyHistos &fHi)
{
	for (size_t i=0; i<aveWaveform.size(); ++i)
		fHi.plot1d(IDOffset+1, i+1,aveWaveform[i]);

	for (size_t i=0; i<aveRebinWaveform.size(); ++i)
		fHi.plot1d(IDOffset+2, i+1,aveRebinWaveform[i]);

	//std::cout<<std::endl;
	//std::cout<<"average of intensity = "<<aveIntensity<<std::endl;
	//std::cout<<"Standard Deviation of intensity = "<<TMath::Sqrt(varIntensity)<<std::endl;

}




//------------------------------------------------------------------------------------------------------------------------------------
//---Integral waveform from TRfrom to TRto. 
double MyWaveform::IntegralWaveform(const MyOriginalChannel &oc, const long TRfrom, const long TRto, bool absolute)
{
	//get some infos from the channel and initalize the integral//
	double IntegralTR		= 0;
	const double vertGain	= oc.GetVertGain();
	const short baseli		= oc.GetBaseline();

	//go through all pulses in this channel//
	for (int puls=0; puls<oc.GetNbrPulses();++puls)
	{
		//get puls and infos from puls and pointer to the array//
		const MyPuls &p				= oc.GetPuls(puls);
		int actPosInEvent			= p.GetIndexToFirstPointOfOriginalWaveform();
		const short *Data			= static_cast<const short*>(oc.GetDataPointerForPuls(p));
		const size_t pLength		= p.GetLength();

		//go through the puls and check wether the actual point is still in the range//
		//if that is the case add this point to the integral//
		for (size_t i=0; i<pLength ; ++i, ++actPosInEvent)
		{
			if ((TRfrom <= actPosInEvent)&&(actPosInEvent <= TRto))
			{
				if (absolute)
				{
					IntegralTR += (TMath::Abs(Data[i]-baseli));
				}
				else
				{
					IntegralTR += (Data[i]-baseli);
				}
			}
		}
	}

	//return the integral//
	return IntegralTR*vertGain;
}

