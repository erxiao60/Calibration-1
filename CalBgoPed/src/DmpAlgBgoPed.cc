/*
 *  $Id: DmpAlgBgoPed.cc, 2014-09-01 14:52:22 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 19/07/2014
*/

#include <stdio.h>

#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"

#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"
#include "DmpEvtBgoPed.h"
#include "DmpAlgBgoPed.h"
#include "DmpDataBuffer.h"
#include "DmpParameterBgo.h"
#include "DmpBgoBase.h"
#include "DmpCore.h"
//-------------------------------------------------------------------
DmpAlgBgoPed::DmpAlgBgoPed()
 :DmpVAlg("Cal/Bgo/Ped"),
  fEvtHeader(0),
  fBgoRaw(0),
  fBgoPed(0)
{
	fTgap=0;
	fTlast=0;
}

//-------------------------------------------------------------------
DmpAlgBgoPed::~DmpAlgBgoPed(){
}

//-------------------------------------------------------------------
bool DmpAlgBgoPed::Initialize(){
 // gRootIOSvc->Set("Output/FileName",gRootIOSvc->GetInputStem()+"_ped.root");
  // read input data
  fEvtHeader = dynamic_cast<DmpEvtHeader*>(gDataBuffer->ReadObject("Event/Dacor/EventHeader"));
  if(fEvtHeader==0){
    fEvtHeader= new DmpEvtHeader();
    gDataBuffer->LinkRootFile("Event/Rdc/EventHeader",fEvtHeader);
  }
  fBgoRaw = dynamic_cast<DmpEvtBgoRaw*>(gDataBuffer->ReadObject("Event/Dacor/Bgo"));
  if(fBgoRaw==0){
    fBgoRaw= new DmpEvtBgoRaw();
    gDataBuffer->LinkRootFile("Event/Rdc/Bgo",fBgoRaw);
  }
  // create output data holder
  fBgoPed = new DmpEvtBgoPed();
  gDataBuffer->RegisterObject("Calibration/Bgo/Pedestal",fBgoPed,"DmpEvtBgoPed");

  fBgoPed->UsedFileName = gRootIOSvc->GetInputFileName();
  gRootIOSvc->PrepareEvent(0);
  fBgoPed->StartTime = fEvtHeader->fSecond;
  fTlast=fEvtHeader->fSecond;
  // create Hist map
  short layerNo = DmpParameterBgo::kPlaneNo*2;
  short barNo = DmpParameterBgo::kBarNo+DmpParameterBgo::kRefBarNo;
  for(short l=0;l<layerNo;++l){
    for(short b=0;b<barNo;++b){
      for(short s=0;s<DmpParameterBgo::kSideNo;++s){
        for(short d=0;d<DmpParameterBgo::kDyNo;++d){
          char name[50];
          short gid_dy = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,d*3+2);
          snprintf(name,50,"BgoPed_%05d-L%02d_B%02d_Dy%02d",gid_dy,l,b,s*10+d*3+2);
          fPedHist.insert(std::make_pair(gid_dy,new TH1F(name,name,1000,-500,500)));
        }
      }
    }
  }
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoPed::ProcessThisEvent(){
  short nSignal = fBgoRaw->fGlobalDynodeID.size();
  short gid = 0;
  double adc = 0.;
  for(short i=0;i<nSignal;++i){
    gid=fBgoRaw->fGlobalDynodeID[i];
    adc=fBgoRaw->fADC[i]; 
      fPedHist[gid]->Fill(adc);
   }
  //check time gap to cut last package
  fTgap =fEvtHeader->fSecond-fTlast;
  fTlast=fEvtHeader->fSecond;
  if(fTgap>60){
    // long entries=fEvtHeader 
    long entries= gCore->GetCurrentEventID();
    if(entries<500){
      std::cout<<"Cut last package & Reset histograms!";
      for(std::map<short,TH1F*>::iterator aHist=fPedHist.begin();aHist!=fPedHist.end();++aHist){
        aHist->second->Reset();
      }
      fBgoPed->StartTime = fEvtHeader->fSecond;
    }
  } 

  return true;
} 

//-------------------------------------------------------------------
bool DmpAlgBgoPed::Finalize(){
  TF1 *gausFit = new TF1("GausFit","gaus",-500,500);
  std::string histFileName = "./Pedestal/Histograms/"+gRootIOSvc->GetInputStem()+"_ped_Hist.root";
  TFile *histFile = new TFile(histFileName.c_str(),"RECREATE");
  fBgoPed->StopTime = fEvtHeader->fSecond;
  for(std::map<short,TH1F*>::iterator aHist=fPedHist.begin();aHist!=fPedHist.end();++aHist){
      fBgoPed->GlobalDynodeID.push_back(aHist->first);
// *
// *  TODO: fit and save output data 
// *
      float mean = aHist->second->GetMean(), sigma = aHist->second->GetRMS();
      for(short i = 0;i<3;++i){
        gausFit->SetRange(mean-2*sigma,mean+2*sigma);
        aHist->second->Fit(gausFit,"RQ");
        mean = gausFit->GetParameter(1);
        sigma = gausFit->GetParameter(2);
      }
      fBgoPed->Mean.push_back(mean);
      fBgoPed->Sigma.push_back(sigma);
      if((mean > 160 || mean<-160) && sigma >30){
         DmpLogWarning<<"GID = "<<aHist->first<<"\tmean = "<<mean<<"\tsigma = "<<sigma<<DmpLogEndl;
      }
      aHist->second->Write();
      delete aHist->second;
  }
  delete histFile;
  return true;
}


