#ifndef __MyWaveform_H__
#define __MyWaveform_H__

#include<vector>
#include <functional>

#include "FilesFromLma2Root/MyEvent/MyOriginalEvent/MyOriginalEvent.h"
#include "FilesFromLma2Root/MyEvent/MySortedEvent/MySortedEvent.h"
#include "FilesFromLma2Root/MyEvent/MySignalAnalyzedEvent/MySignalAnalyzedEvent.h"
#include "FilesFromLma2Root/MyRootManager/MyHistos.h"

class MyHistos;
class MyOriginalEvent;
class MyOriginalEventInfo;
class MyOriginalChannel;
class MyOriginalEvent;
//class MySortedEvent;

class MyWaveform
{
public:
	MyWaveform();
	~MyWaveform();

	typedef std::vector<double> waveform_t;

	void Init(const MyOriginalEvent& ,MyHistos&);
	void Clear();
	void ExtractWaveform(const MyOriginalEvent &oe, MyHistos &rm, int chan);
	void FillHist(MyHistos &fHi);

public:
	const waveform_t &GetWaveform()const				{return waveform;}
	const waveform_t &GetAveWaveform()const				{return aveWaveform;}
	//const waveform_t &GetRebinWaveform()const			{return rebinWaveform;}
	//const waveform_t &GetAveRebinWaveform()const		{return aveRebinWaveform;}

	int GetArraySize()const						{return arraySize;}
	int GetRebinSize()const						{return rebin;}
	int GetTimeRangeFrom()const					{return TimeRangeFrom;}
	long GetEventCount()const					{return eventCounter;}

	double GetIntegral(const long TRfrom, const long TRto, bool absolute);
	double GetAverage(const long TRfrom, const long TRto, bool absolute);

private:

	waveform_t waveform;
	waveform_t aveWaveform;
	waveform_t rebinWaveform;
	waveform_t aveRebinWaveform;

	const int IDOffset;
	long length;
	int arraySize;
	int rebin;
	long eventCounter;
	int TimeRangeFrom;
};


//--------------------------------------------------------------------------------------------
//-------------------------- binary function -------------------------------------------------
//--------------------------------------------------------------------------------------------
/** binary function for averaging.*/

class Average : std::binary_function<double,double,double>
{
public:
  explicit Average(double alpha)
	:_alpha(alpha)
  {}

  /** operator.
   *
   * the operator calculates the average using the function
   * \f$Y_N = Y_{N-1} + \alpha(y-Y_{N-1})\f$
   * where when \f$\alpha\f$ is equal to N it is a cumulative moving average,
   * otherwise it will be a exponential moving average.
   */
  double operator()(double currentValue, double Average_Nm1)
  {
	return Average_Nm1 + _alpha*(currentValue - Average_Nm1);
  }

protected:
  double _alpha;
};

//***Squre Average
class SqAverage : std::binary_function<double,double,double>
{
public:
  explicit SqAverage(double alpha)
	:_alpha(alpha)
  {}

  double operator()(double currentValue, double Average_Nm1)
  {
	return Average_Nm1 + _alpha*(currentValue*currentValue - Average_Nm1);
  }

protected:
  double _alpha;
};

//*** calc previous value
class PreAverage : std::binary_function<double,double,double>
{
public:
  explicit PreAverage(double alpha)
	:_alpha(alpha)
  {}

  float operator()(double currentValue, double Average_Nm1)
  {
	return Average_Nm1 - _alpha*(currentValue - Average_Nm1);
  }

protected:
  double _alpha;
};
#endif