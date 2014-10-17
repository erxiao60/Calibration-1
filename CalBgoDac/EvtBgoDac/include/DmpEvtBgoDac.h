/*   $Id: DmpEvtBgoDac.h, 2014-09-25 17:12:51+08:00 DAMPE $
 *--------------------------------------------------------
 *  Author(s):
 *
 *--------------------------------------------------------
*/

#ifndef DmpEvtBgoDac_H
#define DmpEvtBgoDac_H

#include "TObject.h"

class DmpEvtBgoDac : public TObject{
/*
 *  DmpEvtBgoDac
 *
 */
public:
  DmpEvtBgoDac();
  ~DmpEvtBgoDac();

  void Reset();

public:
  
  std::string   UsedFileName;      // file name of raw data
  int   StartTime;     // the time of the first event used to cal. DAC
  int   StopTime;      // 
  std::vector<short>    GlobalDynodeID;
  std::vector<double>    QdcCoe;//Quadratic Coefficients
  std::vector<double>    Slope;
  std::vector<double>    Cst;// the Constant
  
  ClassDef(DmpEvtBgoDac,1)

};

#endif
