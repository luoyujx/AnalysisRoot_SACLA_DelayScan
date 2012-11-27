#include <iostream>
#include <fstream>
#include <TH1.h>
#include <algorithm>
#include <conio.h>

#include "MyCovariance.h"
#include "MyWaveform.h"

#include "FilesFromLma2Root/MyRootManager/MyHistos.h"
#include "FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalEvent.h"
#include "FilesFromLma2Root/MyEvent/MyOriginalEvent/MyPuls/MyPuls.h"
#include "FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalChannel/MyOriginalChannel.h"

MyCovariance::MyCovariance():
IDOffset(2000)
{}
MyCovariance::~MyCovariance()
{}

void MyCovariance::Init(const MyOriginalEvent &oe, const MyWaveform &wf)
{
	length= oe.GetNbrSamples();
	arraySize = wf.GetArraySize();
	rebin = wf.GetRebinSize();
	TimeRangeFrom = wf.GetTimeRangeFrom();

	//allocate and clear array
	CorrectionVec.resize(arraySize,0.);
	CovarianceMap.resize(arraySize*arraySize,0.);
	CorrectionMap.resize(arraySize*arraySize,0.);
}
void MyCovariance::Clear()
{
	std::fill(CorrectionVec.begin(),CorrectionVec.end(),0.);
	std::fill(CovarianceMap.begin(),CovarianceMap.end(),0.);
	std::fill(CorrectionMap.begin(),CorrectionMap.end(),0.);
}
void MyCovariance::MakeCovMap(MyHistos &rm)
{
	rm.create2d(IDOffset+1,"CovarianceMap","ToF","ToF",arraySize,TimeRangeFrom,TimeRangeFrom+arraySize*rebin,arraySize,TimeRangeFrom,TimeRangeFrom+arraySize*rebin,"Covariance");
	rm.create2d(IDOffset+3,"IntensityCorrMap","ToF","ToF",arraySize,TimeRangeFrom,TimeRangeFrom+arraySize*rebin,arraySize,TimeRangeFrom,TimeRangeFrom+arraySize*rebin,"Covariance");
	rm.create2d(IDOffset+4,"PartialCovarianceMap","ToF","ToF",arraySize,TimeRangeFrom,TimeRangeFrom+arraySize*rebin,arraySize,TimeRangeFrom,TimeRangeFrom+arraySize*rebin,"Covariance");
	rm.create1d(IDOffset+5,"IntensityCorrVec","ToF",arraySize,TimeRangeFrom,TimeRangeFrom+arraySize*rebin,"Covariance");
}
//------------------------------------------------------------------------------------------------------------------------------------
//---Filling Covariance Map
//---this function is only called on last shot event
void MyCovariance::FillCovMap(MyHistos &fHi)
{
	for (size_t i=0; i<arraySize; i++)
	{
		for (size_t j=0; j<arraySize; j++)
		{
				fHi.plot2d(IDOffset+1,i+1,j+1,(CovarianceMap[i*arraySize+j]));
		}
	}
	//for (size_t i=0; i<arraySize; i++)
	//	for (size_t j=i+1; j<arraySize; j++)
	//	{
	//			rm.plot2d(IDOffset+1,j+1,i+1,(CovarianceMap[i*arraySize+j]));//plot
	//	}

	//std::cout<<std::endl<<" Total shots: "<<eventCounter<<std::endl;
}
//------------------------------------------------------------------------------------------------------------------------------------
//---calculate (<I*X'> - <I>*<X'>) * (<I*X> - <I>*<X>) / (<I^2> - <I>^2)
//---Filling Intensity Correction Map
//---this function is only called on last shot event
void MyCovariance::FillCorrMap(const MyWaveform &wf, MyHistos &fHi)
{
	for (size_t i=0; i<arraySize; ++i)
		fHi.plot1d(IDOffset+5,i+1,CorrectionVec[i]);

	for (size_t i=0; i<arraySize; i++)
		for (size_t j=0; j<arraySize; j++)
			CorrectionMap[i*arraySize+j] = CorrectionVec[i]*CorrectionVec[j]/wf.GetVarIntensity();
	
	for (size_t i=0; i<arraySize; i++)
		for (size_t j=0; j<arraySize; j++)
		{
				fHi.plot2d(IDOffset+3,i+1,j+1,(CorrectionMap[i*arraySize+j]));
				fHi.plot2d(IDOffset+4,i+1,j+1,(CovarianceMap[i*arraySize+j] - CorrectionMap[i*arraySize+j]));
		}

	//for (size_t i=0; i<arraySize; i++)
	//	for (size_t j=i+1; j<arraySize; j++)
	//	{
	//			rm.plot2d(IDOffset+3,j+1,i+1,(CorrectionMap[i*arraySize+j]));
	//			rm.plot2d(IDOffset+4,j+1,i+1,(CovarianceMap[i*arraySize+j] - CorrectionMap[i*arraySize+j]));
	//	}

}
//------------------------------------------------------------------------------------------------------------------------------------
//---calculate covariance map
//---this function is called on each shot event
void MyCovariance::CalcCovMap(const MyWaveform &wf)
{
	const waveform_t &rebinWaveform = wf.GetRebinWaveform();
	const waveform_t &aveRebinWaveform = wf.GetAveRebinWaveform();

	waveform_t preAveRebinWaveform(wf.GetAveRebinWaveform().size());

	double scale(1./wf.GetEventCount());

	std::transform(rebinWaveform.begin(),rebinWaveform.end(),
		aveRebinWaveform.begin(),preAveRebinWaveform.begin(),
		PreAverage(scale));

	for (unsigned int i=0; i<arraySize; ++i)
        for (unsigned int j=0; j<arraySize; ++j)
        {
            CovarianceMap[i*arraySize+j] = ((CovarianceMap[i*arraySize+j]*(wf.GetEventCount()-1)) 
				+ (rebinWaveform[i]-preAveRebinWaveform[i])
				* (rebinWaveform[j]-aveRebinWaveform[j]))*scale;
        }


}
//------------------------------------------------------------------------------------------------------------------------------------
//---calculate Intensity Correction array
//---this function is called on each shot event
void MyCovariance::CalcCorrMap(const MyWaveform &wf)
{
	double scale(1./wf.GetEventCount());
	const waveform_t &rebinWaveform = wf.GetRebinWaveform();
	const waveform_t &aveRebinWaveform = wf.GetAveRebinWaveform();

	PreAverage preave(scale); 
	const double PreAveIntensity = preave(wf.GetIntensity(),wf.GetAveIntensity());
	
    for (size_t i=0; i<arraySize; ++i)
        {
			CorrectionVec[i] = (CorrectionVec[i]*(wf.GetEventCount()-1) 
				+ (rebinWaveform[i]-aveRebinWaveform[i])
				* (wf.GetIntensity()-PreAveIntensity))*scale;
        }
}
