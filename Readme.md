Bgo ECal Calibration of DMPSW:    DmpBgoCal_v1.0.0
>--------------------------------------------------

>--------------------------------------------------
Packages:                                   Authors:
RacEQM3Det(Rdc), CalBgoPed;                 Chi Wang
CalBgoDac, CorBgoDacor, CalBgoHits;         Zhiyong Zhang  Yifeng Wei
CorBgoCutped, CalBgoMips;                   Zhiyong Zhang, Chi Wang
CalBgoDyCoe, CalBgoAttLen(under testing);   Jianing Dong
JobOpt                                      Chi Wang Zhiyong Zhang



Usages:
      ./xxx_job inputpath/inputfilename;

      for example:
      1) ./Raw2ROOT_job ./binary/A2data.frd
        Results are saved in ./RawROOT/A2data.root   
      2) ./CalBgoMips_job ./RawROOT/A2data.root
      Parameters are saved in ./MIPs/A2data_BgoMips.root,
      and histograms are saved in ./MIPs/Histograms/A2data_BgoMips_Hist.root


notices:
      The Parameters used in the packages are in the fixed path.
      for example:
      MIPs Parameters shoud be stored as ./MIPs/MipsPar.root, when they are called.


      
                 
