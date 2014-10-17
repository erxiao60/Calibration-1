/*   $Id: DmpEvtBgoMips.h, 2014-10-05 15:30:26+08:00 DAMPE $
 *--------------------------------------------------------
 *  Author(s):
 *
 *--------------------------------------------------------
*/

#ifndef DmpEvtBgoMips_H
#define DmpEvtBgoMips_H

#include "TObject.h"

class DmpEvtBgoMips : public TObject{
/*
 *  DmpEvtBgoMips
 *
 */
public:
  DmpEvtBgoMips();
  ~DmpEvtBgoMips();

  void Reset();

public:
std::string   UsedFileName;      // file name of raw data
int   StartTime;     // the time of the first event used to cal. ped
int   StopTime;      // 
std::vector<short> GlobalBarID;
std::vector<short> BgoSide;//0,1 side and combined 2; 

std::vector<double> MPV;
std::vector<double> Gsigma;
std::vector<double> Lwidth;

std::vector<double> MPV_Error;
std::vector<double> Gsigma_Error;
std::vector<double> Lwidth_Error;

std::vector<double> ChiS;

  /*
   *  Define your data members at here
   *
   */
  
  ClassDef(DmpEvtBgoMips,1)

};

#endif
