#include "DYMvaInCMSSW/GetDYMVA/interface/GetDYMVA.h"

using namespace std;

GetDYMVA::GetDYMVA():isInit(false){
  GetDYMVA(-1);
}

GetDYMVA::GetDYMVA(int version):isInit(false),version_(version){
  
  if (version_==0) {
    init("BDTB", 
	 getenv("CMSSW_BASE")+string("/src/DYMvaInCMSSW/GetDYMVA/data/TMVA_0j_BDTB.weights.xml"), 
	 getenv("CMSSW_BASE")+string("/src/DYMvaInCMSSW/GetDYMVA/data/TMVA_1j_BDTB.weights.xml") );
  } else if (version_==1) {
    init("BDTG", 
	 getenv("CMSSW_BASE")+string("/src/DYMvaInCMSSW/GetDYMVA/data/TMVA_BDTG_0j_MCtrain.weights.xml"),  
	 getenv("CMSSW_BASE")+string("/src/DYMvaInCMSSW/GetDYMVA/data/TMVA_BDTG_1j_MCtrain.weights.xml") );
  } else {
    cerr << "version not supported: use 0 for data train and 1 for MC train" << endl;
  }

}

GetDYMVA::~GetDYMVA(){
  delete theReader_0j;
  delete theReader_1j;
}

void GetDYMVA::init(string methodName, string weightsfile_0j, string weightsfile_1j){

  if (isInit) return;

  theReader_0j = new TMVA::Reader();
  theReader_1j = new TMVA::Reader();

  if (version_==0) {
    theReader_0j->AddVariable( "met"          , &met_           );
    theReader_0j->AddVariable( "trackMet"     , &trackMet_      );
    theReader_0j->AddVariable( "jet1pt"       , &ptjet1_        );
    theReader_0j->AddVariable( "metSig"       , &metSig_        );
    theReader_0j->AddVariable( "dPhiDiLepJet1", &dPhiDiLepJet1_ );
    theReader_0j->AddVariable( "dPhiJet1MET"  , &dPhiJet1MET_   );
    theReader_0j->AddVariable( "mt"           , &mt_            );
    
    theReader_1j->AddVariable( "met"          , &met_           );
    theReader_1j->AddVariable( "trackMet"     , &trackMet_      );
    theReader_1j->AddVariable( "jet1pt"       , &ptjet1_        );
    theReader_1j->AddVariable( "metSig"       , &metSig_        );
    theReader_1j->AddVariable( "dPhiDiLepJet1", &dPhiDiLepJet1_ );
    theReader_1j->AddVariable( "dPhiJet1MET"  , &dPhiJet1MET_   );
    theReader_1j->AddVariable( "mt"           , &mt_            );
  } else if (version_==1) {
    theReader_0j->AddVariable( "pmet"         , &pmet_          );
    theReader_0j->AddVariable( "pTrackMet"    , &pTrackMet_     );
    theReader_0j->AddVariable( "nvtx"         , &nvtx_          );
    theReader_0j->AddVariable( "dilpt"        , &dilpt_         );
    theReader_0j->AddVariable( "jet1pt"       , &ptjet1_        );
    theReader_0j->AddVariable( "metSig"       , &metSig_        );
    theReader_0j->AddVariable( "dPhiDiLepJet1", &dPhiDiLepJet1_ );
    theReader_0j->AddVariable( "dPhiDiLepMET" , &dPhiDiLepMET_  );
    theReader_0j->AddVariable( "dPhiJet1MET"  , &dPhiJet1MET_   );
    theReader_0j->AddVariable( "recoil"       , &recoil_        );
    theReader_0j->AddVariable( "mt"           , &mt_            );
    
    theReader_1j->AddVariable( "pmet"         , &pmet_          );
    theReader_1j->AddVariable( "pTrackMet"    , &pTrackMet_     );
    theReader_1j->AddVariable( "nvtx"         , &nvtx_          );
    theReader_1j->AddVariable( "dilpt"        , &dilpt_         );
    theReader_1j->AddVariable( "jet1pt"       , &ptjet1_        );
    theReader_1j->AddVariable( "metSig"       , &metSig_        );
    theReader_1j->AddVariable( "dPhiDiLepJet1", &dPhiDiLepJet1_ );
    theReader_1j->AddVariable( "dPhiDiLepMET" , &dPhiDiLepMET_  );
    theReader_1j->AddVariable( "dPhiJet1MET"  , &dPhiJet1MET_   );
    theReader_1j->AddVariable( "recoil"       , &recoil_        );
    theReader_1j->AddVariable( "mt"           , &mt_            );
  } else {
    isInit=false;
    return;
  }

  theReader_0j->BookMVA(methodName,weightsfile_0j);
  theReader_1j->BookMVA(methodName,weightsfile_1j);

  methodname_ = methodName;
  isInit=true;

}


double GetDYMVA::getValue(int njets, double met, double trackMet, double ptjet1, double metSig, 
			  double dPhiDiLepJet1, double dPhiJet1MET, double mt){

  if (version_!=0){
    cerr << "MVA version mismatch!" << endl;
    return -999.;
  }
  if (!isInit) {
    cerr << "MVA reader not initialized!" << endl;
    return -999.;
  }

  if (ptjet1<15.){
    ptjet1 = 0.;
    dPhiDiLepJet1 = -0.1;
    dPhiJet1MET = -0.1;
  }

  met_           = met; 
  trackMet_      = trackMet; 
  ptjet1_        = ptjet1; 
  metSig_        = metSig; 
  dPhiDiLepJet1_ = dPhiDiLepJet1; 
  dPhiJet1MET_   = dPhiJet1MET; 
  mt_            = mt;

  if (njets==0) return theReader_0j->EvaluateMVA( methodname_ );
  if (njets==1) return theReader_1j->EvaluateMVA( methodname_ );
  return -999;

}

double GetDYMVA::getValue(int njets, double pmet, double pTrackMet, int nvtx, double dilpt, double ptjet1, double metSig, 
			  double dPhiDiLepJet1, double dPhiDiLepMET, double dPhiJet1MET, double recoil, double mt){

  if (version_!=1){
    cerr << "MVA version mismatch!" << endl;
    return -999.;
  }
  if (!isInit) {
    cerr << "MVA reader not initialized!" << endl;
    return -999.;
  }

  if (ptjet1<15.){
    ptjet1 = 15.;
    dPhiDiLepJet1 = -0.1;
    dPhiJet1MET = -0.1;
  }

  pmet_          = pmet;
  pTrackMet_     = pTrackMet;
  nvtx_          = nvtx;
  dilpt_         = dilpt;     
  ptjet1_        = ptjet1; 
  metSig_        = metSig; 
  dPhiDiLepJet1_ = dPhiDiLepJet1; 
  dPhiDiLepMET_  = dPhiDiLepMET;
  dPhiJet1MET_   = dPhiJet1MET; 
  recoil_        = recoil;
  mt_            = mt;

  if (njets==0) return theReader_0j->EvaluateMVA( methodname_ );
  if (njets==1) return theReader_1j->EvaluateMVA( methodname_ );
  return -999;

}
