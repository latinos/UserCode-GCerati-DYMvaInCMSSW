// -*- C++ -*-
//
// Package:    TestGetDYMVA
// Class:      TestGetDYMVA
// 
/**\class TestGetDYMVA TestGetDYMVA.cc DYMvaInCMSSW/TestGetDYMVA/src/TestGetDYMVA.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giuseppe Cerati,28 S-012,+41227678302,
//         Created:  Wed Apr 25 18:06:59 CEST 2012
// $Id$
//
//

#include <memory>

#include <TMath.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DYMvaInCMSSW/GetDYMVA/interface/GetDYMVA.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class TestGetDYMVA : public edm::EDAnalyzer {
   public:
      explicit TestGetDYMVA(const edm::ParameterSet&);
      ~TestGetDYMVA();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      GetDYMVA* getDYMVA_v0;
      GetDYMVA* getDYMVA_v1;
};

TestGetDYMVA::TestGetDYMVA(const edm::ParameterSet& iConfig){
  getDYMVA_v0 = new GetDYMVA(0);
  getDYMVA_v1 = new GetDYMVA(1);
}

TestGetDYMVA::~TestGetDYMVA(){
  delete getDYMVA_v0;
  delete getDYMVA_v1;
}

void TestGetDYMVA::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;

  //this is just an example: consider only events with 2 muons
  Handle<View<Muon> > muon_h;
  iEvent.getByLabel( "muons" , muon_h );
  if (muon_h->size()!=2) return;
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > dilep = muon_h->at(0).p4()+muon_h->at(1).p4();
  
  Handle<View<PFMET> > met_h;
  iEvent.getByLabel("pfMet", met_h);
  float met    = ( met_h->front() ).et();
  float metPhi = ( met_h->front() ).phi();
  float metSig_v0 = ( met_h->front() ).significance();
  float metSig_v1 = ( met_h->front() ).mEtSig();
  
  //compute track met
  edm::Handle<reco::PFCandidateCollection> pfCand_h;
  iEvent.getByLabel("particleFlow", pfCand_h);
  reco::PFCandidateCollection pfCandCollection = *(pfCand_h.product());
  //Get the PV collection
  Handle<VertexCollection> pvCollection;
  iEvent.getByLabel("offlinePrimaryVertices", pvCollection);
  VertexCollection::const_iterator vertex = pvCollection->begin();
  double pX(0), pY(0);
  pX -= muon_h->at(0).px();
  pY -= muon_h->at(0).py();
  pX -= muon_h->at(1).px();
  pY -= muon_h->at(1).py();
  reco::PFCandidateCollection::const_iterator pf;
  for (pf = pfCandCollection.begin(); pf != pfCandCollection.end(); ++pf) {
    //reject neutrals
    if (pf->charge()==0) continue;
    //avoid overlap with leptons
    float dR0 = reco::deltaR(pf->eta(), pf->phi(), muon_h->at(0).p4().eta(), muon_h->at(0).p4().phi());
    if ( fabs(dR0)<0.1 ) continue;
    float dR1 = reco::deltaR(pf->eta(), pf->phi(), muon_h->at(1).p4().eta(), muon_h->at(1).p4().phi());
    if ( fabs(dR1)<0.1 ) continue;
    //dz cut
    const reco::TrackRef pfTrack  = pf->trackRef();
    if (pfTrack.isNonnull()==false) continue;
    if( fabs(pfTrack->dz(vertex->position()))>0.1 ) continue;
    pX -= pf->px();
    pY -= pf->py();
  }
  double trackMet    = sqrt(pX * pX + pY * pY);
  double trackMetPhi = atan2(pY, pX);

  Handle<View<PFJet> > pfJetsHandle;
  iEvent.getByLabel("ak5PFJets", pfJetsHandle);
  const JetCorrector* correctorL1FastL2L3 = JetCorrector::getJetCorrector ( "ak5PFL1FastL2L3" , iSetup );
  float jet1pt = 0.;
  float jet1phi = 0.;
  int njets=0;
  for(View<PFJet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {
    if (deltaR(*pfjet_it, muon_h->at(0))<0.3) continue;
    if (deltaR(*pfjet_it, muon_h->at(1))<0.3) continue;
    float L1FastL2L3JetScale = correctorL1FastL2L3->correction( *pfjet_it,  iEvent, iSetup );
    float jpt=pfjet_it->p4().pt()*L1FastL2L3JetScale;
    if ( fabs(pfjet_it->p4().eta())>5.0) continue;
    if ( jpt>30.) njets++;
    if ( jpt < jet1pt ) continue;
    jet1pt = jpt;
    jet1phi = pfjet_it->p4().Phi();     
  }
  
  double mt = 2*sqrt(dilep.pt()*met)*fabs(sin(deltaPhi(dilep.phi(),metPhi)/2));   
  double dPhiDiLepJet1 = fabs(deltaPhi(dilep.phi(),jet1phi));
  double dPhiJet1MET = fabs(deltaPhi(jet1phi,metPhi));
  if (jet1pt<15.){
    jet1pt = 0.;
    dPhiDiLepJet1 = -0.1;
    dPhiJet1MET = -0.1;
  }

  double pmet = met;
  float minDPhiMet = TMath::Min(fabs(deltaPhi(muon_h->at(0).p4().phi(),metPhi)), fabs(deltaPhi(muon_h->at(1).p4().phi(),metPhi)));
  if (minDPhiMet < TMath::Pi()/2) pmet = met*TMath::Sin(minDPhiMet);

  double pTrackMet = trackMet;
  float minDPhiTkMet = TMath::Min(fabs(deltaPhi(muon_h->at(0).p4().phi(),trackMetPhi)), fabs(deltaPhi(muon_h->at(1).p4().phi(),trackMetPhi)));
  if (minDPhiTkMet < TMath::Pi()/2) pTrackMet = trackMet*TMath::Sin(minDPhiTkMet);

  int nvtx = 0;
  for ( unsigned int ivtx = 0; ivtx < pvCollection->size(); ++ivtx ){
    if (pvCollection->at(ivtx).isFake()) continue;
    if (pvCollection->at(ivtx).ndof() <= 4.) continue;
    if (pvCollection->at(ivtx).position().Rho() > 2.0) continue;
    if (fabs(pvCollection->at(ivtx).position().Z()) > 24.0) continue;
    nvtx++;
  }

  double dilpt        = (muon_h->at(0).p4()+muon_h->at(1).p4()).pt();
  double dPhiDiLepMET = fabs(deltaPhi(dilep.phi(),metPhi));;

  float px_rec = met*cos(metPhi) + (muon_h->at(0).p4()+muon_h->at(1).p4()).px();       
  float py_rec = met*sin(metPhi) + (muon_h->at(0).p4()+muon_h->at(1).p4()).py();
  double recoil = sqrt(px_rec*px_rec+py_rec*py_rec);

  double mvaValue_v0 = getDYMVA_v0->getValue(njets, met, trackMet, jet1pt, metSig_v0,dPhiDiLepJet1, dPhiJet1MET, mt);
  double mvaValue_v1 = getDYMVA_v1->getValue(njets, pmet, pTrackMet, nvtx, dilpt, jet1pt, metSig_v1,dPhiDiLepJet1, dPhiDiLepMET, dPhiJet1MET, recoil, mt);
  
  cout << iEvent.id().event() << " " << njets << " " 
       << met  << " " << trackMet  << " " 
       << pmet  << " " << pTrackMet << " " 
       << nvtx << " " << dilpt << " " 
       << jet1pt  << " " << metSig_v0  << " " << metSig_v1 << " " 
       << dPhiDiLepJet1 << " " << dPhiDiLepMET  << " " << dPhiJet1MET  << " " 
       << recoil << " " << mt  << " "
       << mvaValue_v0 << " " << mvaValue_v1
       << endl;
   
}


void TestGetDYMVA::beginJob(){}

void TestGetDYMVA::endJob() {}

void TestGetDYMVA::beginRun(edm::Run const&, edm::EventSetup const&){}

void TestGetDYMVA::endRun(edm::Run const&, edm::EventSetup const&){}

void TestGetDYMVA::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}

void TestGetDYMVA::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}

void TestGetDYMVA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TestGetDYMVA);
