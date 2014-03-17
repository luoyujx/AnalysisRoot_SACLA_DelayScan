#include <algorithm>
#include "TMath.h"

#include "MyWaveform.h"

//#include "FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalEvent.h"
#include "FilesFromLma2Root/MyRootManager/MyHistos.h"
#include "FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalChannel/MyOriginalChannel.h"
#include "FilesFromLma2Root/MyEvent/MySortedEvent/MyDetektor/MyDetektor.h"
#include "FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalEventInfo.h"
#include "FilesFromLma2Root/MyEvent/MySortedEvent/MySortedEventInfo.h"
#include "FilesFromLma2Root/MyEvent/MySignalAnalyzedEvent/MySignalAnalyzedEventInfo.h"

MyWaveform::MyWaveform():IDOffset(9900), eventCounter(0),aveIntensity(0.),varIntensity(0.)
{
}
MyWaveform::~MyWaveform()
{
}
void MyWaveform::Init(const MyOriginalEvent &oe, MyHistos &rm)
{
	length= oe.GetNbrSamples();
//	arraySize = 500;
	//rebin = 5;
	//TimeRangeFrom = 2500;

	//--allocate array and clear---//
	waveform.resize(length,0.);
	aveWaveform.resize(length,0.);
	//rebinWaveform.resize(arraySize,0.);
	//aveRebinWaveform.resize(arraySize,0.);

	// make histgram ---ToF, Intensity,,,,
	rm.create1d(IDOffset+1,"WaveTrace","SampleInt",length,0,length);
	//rm.create1d(IDOffset+2,"Tof_rebin","ToF",arraySize,TimeRangeFrom,TimeRangeFrom+arraySize*rebin);

	//---check signal distribution---//
	//const MyOriginalChannel &oc = oe.GetChannel(7-1);
	//const double vertGain = oc.GetVertGain();
	//const short offset = oc.GetBaseline();
	//const short scale = oc.GetFullScale();
	//const double range = scale*rebin;

	//rm.create1d(IDOffset+100,"CheckSignalHight96","signal hight (mV*rebin)",200,-range,range,"test");
	//rm.create1d(IDOffset+101,"CheckSignalHight20","signal hight (mV*rebin)",200,-range,range,"test");
	//rm.create1d(IDOffset+102,"CheckSignalHight85","signal hight (mV*rebin)",200,-range,range,"test");
	//rm.create1d(IDOffset+103,"CheckSignalHight48","signal hight (mV*rebin)",200,-range,range,"test");
	//rm.create1d(IDOffset+110,"CheckSignalOver","ToF",length,0,length,"test");

	//std::cout<<scale<<std::endl;
	//std::cout<<offset<<std::endl;
	//std::cout<<vertGain<<std::endl;
	//std::cout<<range<<std::endl;

}
void MyWaveform::Clear()
{
	std::fill(waveform.begin(),waveform.end(),0.);
	//std::fill(rebinWaveform.begin(),rebinWaveform.end(),0.);
	std::fill(aveWaveform.begin(),aveWaveform.end(),0.);
	//std::fill(aveRebinWaveform.begin(),aveRebinWaveform.end(),0.);
}

void MyWaveform::ExtractWaveform(const MyOriginalEvent &oe, MyHistos &rm, int chan)
{
	const short* data;
	const MyOriginalChannel &oc = oe.GetChannel(chan);
	const double vertGain = oc.GetVertGain();
	const short offset = oc.GetBaseline();

	//--Filling the waveform array--//	
	std::fill(waveform.begin(),waveform.end(),0.);
	for (size_t i=0;i<oc.GetNbrPulses();++i)
	{
		const MyPuls &p = oc.GetPuls(i);
		data			= static_cast<const short*>(oc.GetDataPointerForPuls(p));
		const long pulsLength = p.GetLength();
		const long pulsStartPoint = p.GetIndexToFirstPointOfOriginalWaveform();
		for (size_t j=0; j<pulsLength;++j)
		{
			const double data_mV = (data[j]-offset)*vertGain;
			waveform[(j+pulsStartPoint)] = data_mV;
			//if (data_mV > 1400) rm.fill1d(IDOffset+110,(j+pulsStartPoint));//FAMP max output voltage is 1.5V
		}
	}
	////-------test-------------------
	//waveform_t BLData(length);
	//BaseLineCorr(&waveform[0], &BLData[0],0,length,50,2);
	//for (size_t i=0; i<length; i++) waveform[i]=waveform[i]-BLData[i];
	////-------test-------------------

	//---rebin waveform---//
	//std::fill(rebinWaveform.begin(),rebinWaveform.end(),0.);
	//for (size_t i=0; i<arraySize; ++i)
	//	for (size_t j=0; j<rebin; ++j)
	//		rebinWaveform[i] += waveform[TimeRangeFrom+i*rebin+j];//fill rebined waveform

	eventCounter++;

	//---average waveform---//
	const double scale = 1./eventCounter;
	std::transform(waveform.begin(),waveform.end(),
		aveWaveform.begin(),aveWaveform.begin(),
		Average(scale));
	//std::transform(rebinWaveform.begin(),rebinWaveform.end(),
	//	aveRebinWaveform.begin(),aveRebinWaveform.begin(),
	//	Average(scale));
}


void MyWaveform::FillHist(MyHistos &rm)
{
	//-fill averaged raw waveform---//
	for (size_t i=0; i<aveWaveform.size(); ++i)
		rm.plot1d(IDOffset+1, i+1,aveWaveform[i]);

	//-fill averaged rebin waveform---//
	//for (size_t i=0; i<aveRebinWaveform.size(); ++i)
	//	rm.plot1d(IDOffset+2, i+1,aveRebinWaveform[i]);

}



