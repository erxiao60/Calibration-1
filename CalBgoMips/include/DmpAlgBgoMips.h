#ifndef DmpAlgBgoMips_H
#define DmpAlgBgoMips_H

#include "DmpVAlg.h"
//#include "DmpEvtMCPrimaryParticle.h"
#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"
#include "DmpEvtBgoMips.h"
#include <map>
class DmpEvtBgoRaw;
class DmpEvtHeader;
class DmpEvtBgoMips;

#include "TH1F.h"
class DmpAlgBgoMips : public DmpVAlg{
/* 
 *  DmpAlgBgoMips
 *
 */
public:
  DmpAlgBgoMips();
  ~DmpAlgBgoMips();
  bool Reset(); 
  bool Initialize();
  bool ProcessThisEvent();    // only for algorithm
  bool RawTrack();    // only for algorithm
//  Double_t langaufun(Double_t *x,Double_t *par);
  bool Finalize();
public:
//DmpEvtMCPrimaryParticle *fEvtMCPrimaryParticle;
DmpEvtHeader *fEvtHeader;
DmpEvtBgoRaw *fBgoRaw;
DmpEvtBgoMips *fBgoMips;
std::map<short,TH1F*>fMipsHist;
std::map<short,std::map<short,TH1F*> >fSideMap;
short MaxBar[14];//The maximum bar of each layer
short MaxAdc[14][3];
//bool goodtrack;

};

#endif
