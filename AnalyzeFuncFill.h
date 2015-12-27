#ifndef __AnalyzeFunFuncFill_h_
#define __AnalyzeFunFuncFill_h_

#include "./MyAnalyzer/MyAnalyzer.h"


//void fillParticleHistograms(const MyParticle &p, const MyParticleHit &ph, std::vector<double>& intensity, MyHistos &hi, int hiOff, std::vector<double>& intPart);
void fillParticleHistograms(const MyParticle &p, const MyParticleHit &ph, std::vector<double>& intensity, std::vector<double>& delay, MyHistos &hi, int hiOff, std::vector<double>& intPart, int& delayBins, double& delayFrom, double& delayTo, bool& selectThetaZ, double& thetaZLowerLimit, double& thetaZUpperLimit);
void fillParticleConditionsPos(const MyOriginalEvent &oe, const MyDetektor &det, const MyParticle &p, const MyDetektorHit &dh, std::vector<double>& intensity, MyHistos &hi, int hiOff);
void fillParticleConditionsTof(const MyOriginalEvent &oe, const MyDetektor &det,const MyParticle &p, const MyDetektorHit &dh, std::vector<double>& intensity, MyHistos &hi, int hiOff);
//void fillHistosAfterAnalyzis(const std::vector<MyParticle> &particles, MyHistos &hi,size_t);
void fillHistosAfterAnalyzis(const std::vector<MyParticle> &particles, MyHistos &hi,size_t, int, double, double);
void fillMoleculeHistogram(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff, Molecule &mol, std::vector<double>& intPart);
void fillMoleculeHistogram2(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff);
void fillSpectra(const MyParticle &p1, const MyParticle &p2, MyHistos &hi, int hiOff);
void fillHydrogenHistogram(const MyParticle &p1, const MyParticle &p2,const MyParticle &p3, MyHistos &hi, int hiOff, Molecule &mol);
void fillAngleHistogram(const MyParticle &p1, const MyParticle &p2, MyHistos &hi, int hiOff);
void fill3BodyHistogram(const MyParticle &p1, const MyParticle &p2,const MyParticle &p3, MyHistos &hi, int hiOff);
void fillMCPToFHistograms(const MyOriginalEvent &oe, MyHistos &hi, std::vector<MCPToFRegion> &mcpTofRegion);
void fillMoleculeHistogramCH2I2(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff, Molecule &mol, std::vector<double>& delay, int& delayBins, double& delayFrom, double& delayTo);
void fillMoleculeHistogramCH2I2_3body(const MyParticle &p1, const MyParticle &p2, const MyParticle &p3, std::vector<double>& intensity, MyHistos &hi, int hiOff,/* Molecule &mol,*/ std::vector<double>& delay, int& delayBins, double& delayFrom, double& delayTo);
//void fillAnalogHistogram(const MyOriginalEvent &oe, MyHistos &hi, const MyParticle &p, double delay, int delayBins, double delayFrom, double delayTo);

#endif