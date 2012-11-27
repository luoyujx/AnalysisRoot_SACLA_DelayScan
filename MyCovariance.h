#ifndef __MyCovariance_H__
#define __MyCovariance_H__

#include<vector>

class MyHistos;
class MyOriginalEvent;
class MyOriginalEventInfo;
class MySignalAnalyzedEvent;
class MySignalAnalyzedEventInfo;
class MyOriginalChannel;
class MySettings;
class MyWaveform;

class MyCovariance
{
public:
	MyCovariance();
	~MyCovariance();

	typedef std::vector<double> waveform_t;
	typedef std::vector<waveform_t> waveforms_t;
	typedef std::vector<std::vector<double> > matrix_t;
	typedef std::vector<double> intensities_t;
	typedef std::vector<double> vector_t;

public:
	void Init(const MyOriginalEvent&, const MyWaveform&);
	void Clear();
	void MakeCovMap(MyHistos &);
	void CalcCovMap(const MyWaveform &wf);
	void CalcCorrMap(const MyWaveform &wf);
	void FillCovMap(MyHistos&);
	void FillCorrMap(const MyWaveform &wf, MyHistos &rm);

private:
	const int IDOffset;
	long length;
	int arraySize;
	int rebin;
	int TimeRangeFrom;
	vector_t CorrectionVec;
	vector_t CovarianceMap;
	vector_t CorrectionMap;

};
#endif