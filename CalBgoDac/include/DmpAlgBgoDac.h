/*
 *  $Id: DmpAlgBgoPed.cc, 2014-09-25 14:00:22 DAMPE $
 *  Author(s):
 *  yong (zhzhy@mail.ustc.edu.cn) 25/09/2014
*/
#ifndef DmpAlgBgoDac_H
#define DmpAlgBgoDac_H

#include "TH1F.h"
#include "TH2F.h"
#include "DmpVAlg.h"
#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"
#include "DmpEvtBgoDac.h"
#include <map>
class DmpEvtHeader;
class DmpEvtBgoRaw;
class DmpEvtBgoDac;

class DmpAlgBgoDac : public DmpVAlg{
/*
 *  DmpAlgBgoDac
 *
 */
public:
  DmpAlgBgoDac();
  ~DmpAlgBgoDac();

  //void Set(const std::string &type,const std::string &value);
  // if you need to set some options for your algorithm at run time. Overload Set()
  bool Initialize();
  bool ProcessThisEvent();    // only for algorithm
  bool Finalize();
  bool PrepareCalVol();
  void DacHist2F_fill();
  Float_t GetVoltage(const short &ipoint,const short &gid);
  short GetFeeID(const short &gid);


private:
  DmpEvtHeader          *fEvtHeader;
  DmpEvtBgoRaw          *fBgoRaw;
  DmpEvtBgoDac          *fBgoDac;
  std::map<short,TH1F*>  fDacHist1F;
  std::map<short,TH2F*>  fDacHist2F;
  Float_t Cal_Vol[16][15];
  short DacCalPot;
  short TriggerNub;

  // key is global dynode ID

};

#endif
