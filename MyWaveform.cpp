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

void BaseLineCorr(double* data, double* BLdata, double offset, const long ndata, const int aveWindow, const int nAve);

MyWaveform::MyWaveform():IDOffset(9900), eventCounter(0)
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
	waveform_t BLData(length);
	BaseLineCorr(&waveform[0], &BLData[0],0,length,50, 1);
	for (size_t i=0; i<length; i++) waveform[i]=waveform[i]-BLData[i];
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


double MyWaveform::GetIntegral(const long TRfrom, const long TRto, bool absolute)
{
	double signalSum = 0.0;
	for (size_t i = TRfrom; i < TRto; i++)
	{
		if (absolute) signalSum += TMath::Abs(waveform[i]);
		else signalSum += waveform[i];
	}
	return signalSum;
}

double MyWaveform::GetAverage(const long TRfrom, const long TRto, bool absolute)
{
	//return the Average//
	return GetIntegral(TRfrom, TRto, absolute) / (TRto - TRfrom);
}

//################BaseLine Correction (test)#############################
void BaseLineCorr(double* data, double* BLdata, double offset, const long ndata, const int aveWindow, const int nAve)//positive
{
	double temp;
	int n;
	double* tempData = new double[ndata];

	for (int i = 0; i<ndata; i++)
	{
		temp = 0.;
		for (int j = i - aveWindow; j <= i + aveWindow; ++j)
		{
			if ((j<0) || (j >= ndata)) temp += offset;
			else temp += data[j];
		}
		BLdata[i] = temp / (2 * aveWindow + 1);

		for (int k = 0; k<nAve; ++k)
		{
			temp = 0.;
			n = 0;
			for (int j = i - aveWindow; j <= i + aveWindow; ++j)
			{
				if (((j >= 0) && (j<ndata)) && (BLdata[i]>data[j]))
				{
					temp += data[j];
					n++;
				}
			}
			if (n>0) BLdata[i] = temp / n;
			else break;
		}
	}

	for (int i = 0; i<ndata; i++)
	{
		temp = 0.;
		for (int j = i - aveWindow; j <= i + aveWindow; ++j)
		{
			if ((j<0) || (j >= ndata)) temp += offset;
			else temp += BLdata[j];
		}
		tempData[i] = temp / (2 * aveWindow + 1);
	}

	for (int i = 0; i<ndata; i++)
		BLdata[i] = tempData[i];

	delete[] tempData;
}

