/* for MIPs calibration, fitted by langaus*/
/* yong zhzhy@mail.ustc.edu.cn 04-10-2014 */

#include "DmpAlgBgoMips.h"
#include "DmpRootIOSvc.h"
#include "DmpDataBuffer.h"
#include "DmpBgoBase.h"
#include "DmpParameterBgo.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

//-------------------------------------------------------------------
//define laugaufun
Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}
//-------------------------------------------------------------------
DmpAlgBgoMips::DmpAlgBgoMips()
 :DmpVAlg("CalBgoMips"),
  fEvtHeader(0),
  fBgoRaw(0),
  fBgoMips(0)
{ 
Reset();
}

bool DmpAlgBgoMips::Reset(){
	memset(MaxBar,0,sizeof(MaxBar));
	for(int i=0;i<14;++i){
	  memset(MaxAdc[i],0,sizeof(MaxAdc[i]));
	}
	return true;
}
//-------------------------------------------------------------------
DmpAlgBgoMips::~DmpAlgBgoMips(){
}
//-------------------------------------------------------------------
bool DmpAlgBgoMips::Initialize(){
//  gRootIOSvc->Set("OutData/FileName","./"+gRootIOSvc->GetInputFileName()+"_mip.root");
  // read input data
// gDataBuffer->LinkRootFile("Event/Rdc/EventHeader",fEvtHeader);
// gDataBuffer->LinkRootFile("Event/Rdc/Bgo",fBgoRaw);
  fEvtHeader = dynamic_cast<DmpEvtHeader*>(gDataBuffer->ReadObject("Event/Cutped/EventHeader"));
  fBgoRaw = dynamic_cast<DmpEvtBgoRaw*>(gDataBuffer->ReadObject("Event/Cutped/Bgo"));
    
  // create output data holder
  fBgoMips = new DmpEvtBgoMips();
  gDataBuffer->RegisterObject("Calibration/Bgo/Mips",fBgoMips,"DmpEvtBgoMips");
  fBgoMips->UsedFileName = gRootIOSvc->GetInputFileName();
//  std::cout<<"0000000000000000000000000"<<std::endl;
  gRootIOSvc->PrepareEvent(0);
//  std::cout<<"0000000000000000000000000"<<std::endl;
  fBgoMips->StartTime = fEvtHeader->fSecond;
//  std::cout<<"0000000000000000000000000"<<std::endl;
  // create Hist map
  short layerNo = DmpParameterBgo::kPlaneNo*2;
  short barNo = DmpParameterBgo::kBarNo;
  for(short s=0;s<3;++s){//Side 0,1 and 2(combined)
    for(short l=0;l<layerNo;++l){
      for(short b=0;b<barNo;++b){
        char name[50];
        short gid=DmpBgoBase::ConstructGlobalBarID(l,b); 
          sprintf(name,"BgoMips_L%02d_B%02d_S%02d",l,b,s); 
    	  fMipsHist.insert(std::make_pair(gid,new TH1F(name,name,300,15,5100)));
       }
     }
	  std::cout<<std::endl;
     fSideMap.insert(std::make_pair(s,fMipsHist));
     fMipsHist.clear();
   }
 return true;
} 

//-------------------------------------------------------------------
bool DmpAlgBgoMips::ProcessThisEvent(){
  //find a raw track
  bool goodtrack=RawTrack();
  //fll
  if(goodtrack){
   // std::cout<<"~~~~~~~~~~~~~~~~~~~"<<std::endl;
    for(short is=0;is<3;++is){
      for(short il=0;il<14;++il){
        short gid=DmpBgoBase::ConstructGlobalBarID(il,MaxBar[il]);
        if(MaxAdc[il][is]>0)	
        fSideMap[is][gid]->Fill(MaxAdc[il][is]);
    //    std::cout<<MaxAdc[il][is]<<"\t";
      }
    }
  }
  Reset();
  return true;
}
bool DmpAlgBgoMips::RawTrack(){
//  Lmax=(Int_t*)malloc(20);
  double AdcBuffer[14][22][2];
  short gid=0,l=0,b=0,s=0,d=0;
  short nHits=0;
  double adc =0.;
  short nSignal = fBgoRaw->fGlobalDynodeID.size();
//  std::cout<<"nSignal:"<<nSignal<<std::endl;
  for(short i=0; i<nSignal;i++){
    gid=fBgoRaw->fGlobalDynodeID[i];
    adc=fBgoRaw->fADC[i]; 
    l=DmpBgoBase::GetLayerID(gid);
    b=DmpBgoBase::GetBarID(gid);
    s=DmpBgoBase::GetSideID(gid);
    d=DmpBgoBase::GetDynodeID(gid);
    if(b>=22){continue;}//spare channels
     if(s==0&&d==8){//only check side 0;
	nHits++;
        if(adc>MaxAdc[l][0]){
	MaxAdc[l][0]=adc;
	MaxBar[l]=b;//0--21
       }
      }
    else{AdcBuffer[l][b][s]=adc;}
     } //single event Loop end
   for(short j=0;j<14;j++){
  MaxAdc[j][1]=AdcBuffer[j][MaxBar[j]][1];
  MaxAdc[j][2]=TMath::Sqrt(MaxAdc[j][0]*MaxAdc[j][1]);
  }
  //MIPs Tracker check
  if(nHits<10||nHits>45){
  return false;
  }
  if((MaxAdc[0][0]<100&&MaxAdc[1][0]<100)||(MaxAdc[12][0]<100&&MaxAdc[13][0]<100)) {
  return false;
  }
  for(int i=0;i<14;i++){
  if(MaxAdc[i][0]>14000){return false;}
  }
  return true;
}   

//-------------------------------------------------------------------
 bool  DmpAlgBgoMips::Finalize(){
 // TF1 *gausFit = new TF1("GausFit","gaus",-500,500);
  std::string histFileName ="./MIPs/Histograms/"+gRootIOSvc->GetInputStem()+"_mip_Hist.root";
  fBgoMips->StopTime = fEvtHeader->fSecond;
  Double_t mean=0.,par[4]={0,0,0,0},chis=0;
  Long64_t nentries;
  Double_t *par_error;
  TF1 *myMIPs=new TF1("myMIPs",langaufun,0.,5800.,4);
  myMIPs->SetParNames("Width","MP","Area","GSigma");
  TFile *histFile = new TFile(histFileName.c_str(),"RECREATE");
  for(std::map<short,std::map<short,TH1F*> >::iterator aSide=fSideMap.begin();aSide!=fSideMap.end();++aSide){
    for(std::map<short,TH1F*>::iterator aHist=aSide->second.begin();aHist!=aSide->second.end();++aHist){
      mean=aHist->second->GetMean();
      nentries=aHist->second->GetEntries();
      myMIPs->SetRange(mean*0.42,mean*2.0);
	  //myMIPs->SetRange(par[1]-2.8*par[2],par[1]+6.50*par[2]);
	   //myMIPs->SetParameters(par);
	   // Once again, here are the Landau * Gaussian parameters:
	   //   par[0]=Width (scale) parameter of Landau density
	   //   par[1]=Most Probable (MP, location) parameter of Landau density
	   //   par[2]=Total area (integral -inf to inf, normalization constant)
	   //   par[3]=Width (sigma) of convoluted Gaussian function
      myMIPs->SetParLimits(3,0,0.8*mean);
      myMIPs->SetParLimits(1,mean*0.6,mean*1.1);
  //    myMIPs->SetParLimits(2,nentries*0.6,nentries*1.4);
      myMIPs->SetParameter(1,mean*0.8);
      myMIPs->SetParameter(3,mean*0.3);
    //  myMIPs->SetParameter(2,nentries);
      myMIPs->SetParameter(0,mean*80);
      aHist->second->Fit(myMIPs,"RQ");
      //fill MIPs event class
      fBgoMips->GlobalBarID.push_back(aHist->first);
      fBgoMips->BgoSide.push_back(aSide->first);
    //using maximum instead of MPV 
      myMIPs->GetParameters(par);
      double peak=myMIPs->GetMaximumX(0.8*par[1],1.5*par[1]);
      par[1]=peak;
      myMIPs->SetParameters(par);
      myMIPs->SetRange(par[1]*0.5,par[1]*3);
      

      par_error=myMIPs->GetParErrors();
      chis=myMIPs->GetChisquare();
      fBgoMips->MPV.push_back(par[1]);
      fBgoMips->MPV_Error.push_back(par_error[1]);
      fBgoMips->Gsigma.push_back(par[3]);
      fBgoMips->Gsigma_Error.push_back(par_error[3]);
      fBgoMips->Lwidth.push_back(par[0]);
      fBgoMips->Lwidth_Error.push_back(par_error[0]);
      fBgoMips->ChiS.push_back(chis);
      
      aHist->second->Write();
      delete aHist->second;
     } 
   }
  histFile->Close();
  return true;
} 





