#ifndef __AnalyzeFunktions_h_
#define __AnalyzeFunktions_h_

#include <vector>
#include <time.h>
#include "./MyAnalyzer/MyAnalyzer.h"

class MyParticleContainer;
class MyHistos;
class MyParticle;
class MyParticleHit;
class MyDetektor;
class MyDetektorHit;
class MyOriginalEvent;
class MyOriginalChannel;

void DefineParticlesAndRootFile(MyParticleContainer &particles, MyHistos &hi);

void fillParticleHistograms(const MyParticle &p, const MyParticleHit &ph, std::vector<double>& intensity, MyHistos &hi, int hiOff);
void fillParticleConditionsPos(const MyOriginalEvent &oe, const MyDetektor &det, const MyParticle &p, const MyDetektorHit &dh, std::vector<double>& intensity, MyHistos &hi, int hiOff);
void fillParticleConditionsTof(const MyOriginalEvent &oe, const MyDetektor &det,const MyParticle &p, const MyDetektorHit &dh, std::vector<double>& intensity, MyHistos &hi, int hiOff);
void fillHistosAfterAnalyzis(const std::vector<MyParticle> &particles, MyHistos &hi,size_t);
void fillMoleculeHistogram(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff, Molecule &mol, std::vector<double>& intPart);
void fillMoleculeHistogram2(const MyParticle &p1, const MyParticle &p2, std::vector<double>& intensity, MyHistos &hi, int hiOff);
void fillSpectra(const MyParticle &p1, const MyParticle &p2, MyHistos &hi, int hiOff);

void fillPIPICO(const MyParticle &p,MyHistos &hi);

double Integral(const MyOriginalChannel &oc, const long TRfrom, const long TRto, bool absolute);
double Average(const MyOriginalChannel &oc, const long TRfrom, const long TRto, bool absolute);

double calcPx(const MyParticle &p, const MyParticleHit &ph);
double calcPy(const MyParticle &p, const MyParticleHit &ph);
double calcPy(const MyParticle &p, const MyParticleHit &ph);

bool PosCondition(const MyDetektorHit &dh);
bool TofPosCondition(const MyDetektorHit &dh);

#endif