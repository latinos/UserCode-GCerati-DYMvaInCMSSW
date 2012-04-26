#ifndef GetDYMVA_h
#define GetDYMVA_h

#include "TMVA/Reader.h"

class GetDYMVA {

 public:
  GetDYMVA();
  GetDYMVA(int version);
  ~GetDYMVA();

  void init(std::string methodName, std::string weightsfile_0j, std::string weightsfile_1j);
  double getValue(int njets, double met, double trackMet, double ptjet1, double metSig, double dPhiDiLepJet1, double dPhiJet1MET, double mt);
  double getValue(int njets, double pmet, double pTrackMet, int nvtx, double dilpt, double ptjet1, double metSig, 
		  double dPhiDiLepJet1, double dPhiDiLepMET, double dPhiJet1MET, double recoil, double mt);


 private:

  TMVA::Reader* theReader_0j;
  TMVA::Reader* theReader_1j;

  //common variables
  float ptjet1_; 
  float metSig_;//this is in common but the meaning is different for v0 and v1!≈ß
  float dPhiDiLepJet1_; 
  float dPhiJet1MET_; 
  float mt_;
  //variables for data train only
  float met_; 
  float trackMet_; 
  //variables for MC train only
  float pmet_;
  float pTrackMet_;
  float nvtx_;//this is an int but the reader does not like it
  float dilpt_;     
  float dPhiDiLepMET_;
  float recoil_;

  std::string methodname_;

  bool isInit;
  int version_;

};

#endif
