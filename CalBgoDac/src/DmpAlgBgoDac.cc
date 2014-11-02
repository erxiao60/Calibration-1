/*
 *  $Id: DmpAlgBgoPed.cc, 2014-09-25 14:00:22 DAMPE $
 *  Author(s):
 *  yong (zhzhy@mail.ustc.edu.cn) 25/09/2014
*/
#include <stdio.h>

#include "TF1.h"
#include "TFile.h"

#include "DmpAlgBgoDac.h"
#include "DmpDataBuffer.h"
#include "DmpParameterBgo.h"
#include "DmpBgoBase.h"
#include "DmpCore.h"
//-------------------------------------------------------------------
DmpAlgBgoDac::DmpAlgBgoDac()
 :DmpVAlg("BgoDac"),
	fEvtHeader(0),
	fBgoRaw(0),
	fBgoDac(0)
{
	DacCalPot=0;
	TriggerNub=0;

}

//-------------------------------------------------------------------
DmpAlgBgoDac::~DmpAlgBgoDac(){

}

//-------------------------------------------------------------------
bool DmpAlgBgoDac::Initialize(){
//  gRootIOSvc->Set("Output/FileName","./"+gRootIOSvc->GetInputFileName()+"_Dac.root");
  // read input data
  fEvtHeader= new DmpEvtHeader();
  fBgoRaw=new DmpEvtBgoRaw();
  gDataBuffer->LinkRootFile("Event/Rdc/EventHeader",fEvtHeader);
  gDataBuffer->LinkRootFile("Event/Rdc/Bgo",fBgoRaw);
  fEvtHeader = dynamic_cast<DmpEvtHeader*>(gDataBuffer->ReadObject("Event/Rdc/EventHeader"));
  fBgoRaw = dynamic_cast<DmpEvtBgoRaw*>(gDataBuffer->ReadObject("Event/Rdc/Bgo"));
  // create output data holder
  fBgoDac = new DmpEvtBgoDac();
  gDataBuffer->RegisterObject("Calibration/Bgo/Dac",fBgoDac,"DmpEvtBgoDac"); 
  fBgoDac->UsedFileName = gRootIOSvc->GetInputFileName();
  gRootIOSvc->PrepareEvent(0);
  fBgoDac->StartTime = fEvtHeader->GetSecond();
  //fTlast=fEvtHeader->GetSecond();
  // create Hist map
  short layerNo = DmpParameterBgo::kPlaneNo*2;
  short barNo = DmpParameterBgo::kBarNo+DmpParameterBgo::kRefBarNo;
  for(short l=0;l<layerNo;++l){
    for(short b=0;b<barNo;++b){
      for(short s=0;s<DmpParameterBgo::kSideNo;++s){
        for(short d=0;d<DmpParameterBgo::kDyNo;++d){
          char title1F[50];
	  char title2F[50];
	  char name2F[80];
          short gid_dy = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,d*3+2);
          snprintf(title1F,50,"BgoDac_%05d_TH1F",gid_dy);
          snprintf(title2F,50,"BgoDac_L%02d_B%02d_Dy%02d",l,b,s*10+d*3+2);
          snprintf(name2F,80,"BgoDac_L%02d_B%02d_Dy%02d;Input(mV);ADC channels",l,b,s*10+d*3+2);
          fDacHist1F.insert(std::make_pair(gid_dy,new TH1F(title1F,title1F,1700,-1000,16000)));
          fDacHist2F.insert(std::make_pair(gid_dy,new TH2F(title2F,name2F,250,0,2500,1700,-1000,16000)));
        }
      }
    }
  }
  //prepare  DAC Input Voltage from DACcode
  if(!PrepareCalVol()){
  return false;
  }
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoDac::PrepareCalVol(){
 // Float_t Cal_DACcode[15]={160,320,480,640,800,960,1120,1280,1440,1600,1760,1920,2080,2240,2400};
  Float_t Cal_DACcode[15]={160,320,640,800,1120,1280,1600,1760,2080,2240,2400};
  //Float_t Cal_Vol[16][15];
  for(int i=0;i<15;i++){
    Cal_Vol[0][i]=0.76575*Cal_DACcode[i]+5.7013;
    Cal_Vol[1][i]=0.76288*Cal_DACcode[i]+3.28139;
    Cal_Vol[2][i]=0.76615*Cal_DACcode[i]+1.8961;
    Cal_Vol[3][i]=0.75676*Cal_DACcode[i]-0.74459;
    Cal_Vol[4][i]=0.75278*Cal_DACcode[i]-2.12987;
    Cal_Vol[5][i]=0.76327*Cal_DACcode[i]+1.44156;
    Cal_Vol[6][i]=0.76291*Cal_DACcode[i]+0.38528;
    Cal_Vol[7][i]=0.75692*Cal_DACcode[i]+2.47619;
    Cal_Vol[8][i]=0.76002*Cal_DACcode[i]+1.31169;
    Cal_Vol[9][i]=0.75968*Cal_DACcode[i]+1.41991;
    Cal_Vol[10][i]=0.75951*Cal_DACcode[i]-0.26407;
    Cal_Vol[11][i]=0.75845*Cal_DACcode[i]+3.08658;
    Cal_Vol[12][i]=0.76415*Cal_DACcode[i]+1.93074;
    Cal_Vol[13][i]=0.76422*Cal_DACcode[i]-0.58874;
    Cal_Vol[14][i]=0.76073*Cal_DACcode[i]+1.59307;
    Cal_Vol[15][i]=0.75579*Cal_DACcode[i]+3.91775;
  }
  return true;
}
//-------------------------------------------------------------------
bool DmpAlgBgoDac::ProcessThisEvent(){
    //check run mode
  if(fBgoRaw->GetRunMode()!=2)
    return true;
     //Fill
  short nSignal = fBgoRaw->fGlobalDynodeID.size();
  short gid = 0;
  double adc =0.;
  for(short i=0;i<nSignal;++i){
    gid=fBgoRaw->fGlobalDynodeID[i];
    adc=fBgoRaw->fADC[i]; 
      fDacHist1F[gid]->Fill(adc);
  }
  TriggerNub++;
  //check DAC cal Event Number
  if(TriggerNub==512||gCore->GetCurrentEventID()==gCore->GetMaxEventNumber()-1){
    DacHist2F_fill();
  //DacHist1F_Reset();//conformed in DacHist2F_fill();
    TriggerNub=0;
    DacCalPot++;
   }
   return true;
} 

//-------------------------------------------------------------------
void DmpAlgBgoDac::DacHist2F_fill(){
  Float_t DacMean=0;
  short gid=0;
  Float_t InputVoltage=0;
  for(std::map<short,TH1F*>::iterator aHist=fDacHist1F.begin();aHist!=fDacHist1F.end();++aHist){
  DacMean=(Float_t)aHist->second->GetMean();
  gid=aHist->first;
  InputVoltage=GetVoltage(DacCalPot,gid);
  fDacHist2F[gid]->Fill(InputVoltage,DacMean);
  aHist->second->Reset();
  }
}
//-------------------------------------------------------------------
Float_t DmpAlgBgoDac::GetVoltage(const short &ipoint,const short &gid){
        if(ipoint>=15){
		std::cout<<"~~~~~~~~~~~~~~~~~~~"<<std::endl;
		std::cout<<"Dac Cal data event overflow,please chech the trigger Number of one Cal point!";
	exit(0);
	}
	Int_t feeid=(Int_t)GetFeeID(gid);
	return Cal_Vol[feeid][ipoint];
}
short DmpAlgBgoDac::GetFeeID(const short &gid){
short ilayer=DmpBgoBase::GetLayerID(gid);
short iside=DmpBgoBase::GetSideID(gid);
if((ilayer==0||ilayer==2)&&iside==0)
  return 0;
else if((ilayer==4||ilayer==6)&&iside==0)
  return 1;
else if(ilayer==8&&iside==0)
  return 2;
else if((ilayer==10||ilayer==12)&&iside==0)
  return 3;
if((ilayer==0||ilayer==2)&&iside==1)
  return 4;
else if((ilayer==4||ilayer==6)&&iside==1)
  return 5;
else if(ilayer==8&&iside==1)
  return 6;
else if((ilayer==10||ilayer==12)&&iside==1)
  return 7;
if((ilayer==1||ilayer==3)&&iside==0)
  return 8;
else if((ilayer==5||ilayer==7)&&iside==0)
  return 9;
else if(ilayer==9&&iside==0)
  return 10;
else if((ilayer==11||ilayer==13)&&iside==0)
  return 11;
if((ilayer==1||ilayer==3)&&iside==1)
  return 12;
else if((ilayer==5||ilayer==7)&&iside==1)
  return 13;
else if(ilayer==9&&iside==1)
  return 14;
else if((ilayer==11||ilayer==13)&&iside==1)
  return 15;
  else {std::cout<<"gid Error(Layer,Side):"<<ilayer<<","<<iside<<std::endl;
  exit(0);}
}
//-------------------------------------------------------------------
bool DmpAlgBgoDac::Finalize(){
  TF1 *Dacfit=new TF1("DacFit","[0]*x*x+[1]*x+[2]",80,1400);
  std::string histFileName = "./DAC/Histograms/"+gRootIOSvc->GetInputStem()+"_dac_Hist.root";
  TFile *histFile = new TFile(histFileName.c_str(),"RECREATE");
  fBgoDac->StopTime = fEvtHeader->GetSecond();
  Double_t DacPar[3];
  for(std::map<short,TH2F*>::iterator aHist=fDacHist2F.begin();aHist!=fDacHist2F.end();++aHist){
      fBgoDac->GlobalDynodeID.push_back(aHist->first);
      aHist->second->SetMarkerStyle(25);
      aHist->second->Fit(Dacfit,"RQ");
      Dacfit->GetParameters(DacPar);
      fBgoDac->QdcCoe.push_back(DacPar[0]);
      fBgoDac->Slope.push_back(DacPar[1]);
      fBgoDac->Cst.push_back(DacPar[2]);

      if(DacPar[1] > 12 || DacPar[1]<8){
         DmpLogWarning<<"Layer="<<DmpBgoBase::GetLayerID(aHist->first)<<" Bar="<<DmpBgoBase::GetBarID(aHist->first)<<" Side="<<DmpBgoBase::GetSideID(aHist->first)<<"\tSlope= "<<DacPar[1]<<DmpLogEndl;
       }
      aHist->second->Write();
      delete aHist->second;
   }
  delete histFile;
  return true;
}

