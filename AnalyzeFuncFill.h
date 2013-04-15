#ifndef __AnalyzeFunFuncFill_h_
#define __AnalyzeFunFuncFill_h_

#include "./MyAnalyzer/MyAnalyzer.h"


void fillParticleHistograms(const MyParticle &p, const MyParticleHit &ph, std::vector<double>& intensity, MyHistos &hi, int hiOff, std::vector<double>& intPart);
void fillParticleConditionsPos(const MyOriginalEvent &oe, const MyDetektor &det, const MyParticle &p, const MyDetektorHit &dh, std::vector<double>& intensity, MyHistos &hi, int hiOff);
void fillParticleConditionsTof(const MyOriginalEvent &oe, const MyDetektor &det,const MyParticle &p, const MyDetektorHit &dh, std::vector<double>& intensity, MyHistos &hi, int hiOff);
void fillHistosAfterAnalyzis(const std::vector<MyParticle> &particles, MyHistos &hi,size_t);
void fillMoleculeHistogram(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff, Molecule &mol, std::vector<double>& intPart);
void fillMoleculeHistogram2(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff);
void fillSpectra(const MyParticle &p1, const MyParticle &p2, MyHistos &hi, int hiOff);

#endif