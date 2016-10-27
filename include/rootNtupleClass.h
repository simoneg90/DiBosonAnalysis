//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 17 14:43:23 2016 by ROOT version 6.06/01
// from TChain tree/
//////////////////////////////////////////////////////////

#ifndef rootNtupleClass_h
#define rootNtupleClass_h

//// Lines added by make_rootNtupleClass.sh - BEGIN 
#include <vector> 
using namespace std; 
//// Lines added by make_rootNtupleClass.sh - END 

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class rootNtupleClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       evt;
   Int_t           isData;
   Float_t         xsec;
   Float_t         puWeight;
   Float_t         nTrueInt;
   Float_t         genWeight;
   Float_t         rho;
   Int_t           nVert;
   Int_t           nJet25;
   Int_t           nBJetLoose25;
   Int_t           nBJetMedium25;
   Int_t           nBJetTight25;
   Int_t           nJet30;
   Int_t           nJet30a;
   Int_t           nBJetLoose30;
   Int_t           nBJetMedium30;
   Int_t           nBJetTight30;
   Int_t           nJet40;
   Int_t           nJet40a;
   Int_t           nBJetLoose40;
   Int_t           nBJetMedium40;
   Int_t           nBJetTight40;
   Int_t           nLepGood20;
   Int_t           nLepGood15;
   Int_t           nLepGood10;
   Float_t         htJet25;
   Float_t         mhtJet25;
   Float_t         htJet40j;
   Float_t         htJet40ja;
   Float_t         htJet40;
   Float_t         htJet40a;
   Float_t         mhtJet40;
   Float_t         mhtJet40a;
   Float_t         apcjetmetmin;
   Float_t         mZ1;
   Float_t         mZ1SFSS;
   Float_t         minMllSFOS;
   Float_t         maxMllSFOS;
   Float_t         minMllAFOS;
   Float_t         maxMllAFOS;
   Float_t         minMllAFSS;
   Float_t         maxMllAFSS;
   Float_t         minMllAFAS;
   Float_t         maxMllAFAS;
   Float_t         m2l;
   Int_t           nBJet25;
   Int_t           nFatJet100;
   Int_t           nMuons10;
   Int_t           nElectrons10;
   Int_t           nTaus20;
   Int_t           nGammas20;
   Float_t         rhoCN;
   Int_t           hbheFilterNew50ns;
   Int_t           hbheFilterNew25ns;
   Int_t           hbheFilterIso;
   Float_t         met_trkPt;
   Float_t         met_trkPhi;
   Int_t           HLT_BIT_HLT_PFMET170_NoiseCleaned_v;
   Int_t           HLT_BIT_HLT_PFMET170_v;
   Int_t           HLT_BIT_HLT_PFMET170_HBHECleaned_v;
   Int_t           HLT_BIT_HLT_PFMET170_JetIdCleaned_v;
   Int_t           HLT_Met170;
   Int_t           HLT_BIT_HLT_Photon155_v;
   Int_t           HLT_BIT_HLT_Photon165_HE10_v;
   Int_t           HLT_BIT_HLT_Photon175_v;
   Int_t           HLT_BIT_HLT_PFJet40_v;
   Int_t           HLT_BIT_HLT_PFJet60_v;
   Int_t           HLT_BIT_HLT_PFJet80_v;
   Int_t           HLT_BIT_HLT_PFJet140_v;
   Int_t           HLT_BIT_HLT_PFJet200_v;
   Int_t           HLT_BIT_HLT_PFJet260_v;
   Int_t           HLT_BIT_HLT_PFJet320_v;
   Int_t           HLT_BIT_HLT_PFJet400_v;
   Int_t           HLT_BIT_HLT_PFJet450_v;
   Int_t           HLT_BIT_HLT_PFJet500_v;
   Int_t           HLT_SinglePho;
   Int_t           HLT_BIT_HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v;
   Int_t           HLT_BIT_HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v;
   Int_t           HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;
   Int_t           HLT_MonoJetMetNoMuMHT120;
   Int_t           HLT_BIT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v;
   Int_t           HLT_BIT_HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v;
   Int_t           HLT_BIT_HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v;
   Int_t           HLT_MonoJetMetNoMuMHT90;
   Int_t           HLT_BIT_HLT_PFHT800_v;
   Int_t           HLT_BIT_HLT_ECALHT800_v;
   Int_t           HLT_Jet_HT;
   Int_t           HLT_BIT_HLT_PFMET300_NoiseCleaned_v;
   Int_t           HLT_BIT_HLT_PFMET300_v;
   Int_t           HLT_BIT_HLT_PFMET300_JetIdCleaned_v;
   Int_t           HLT_Met300;
   Int_t           HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
   Int_t           HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
   Int_t           HLT_DoubleEl;
   Int_t           HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v;
   Int_t           HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v;
   Int_t           HLT_DoubleMu;
   Int_t           HLT_BIT_HLT_IsoMu24_eta2p1_v;
   Int_t           HLT_BIT_HLT_IsoTkMu24_eta2p1_v;
   Int_t           HLT_BIT_HLT_IsoMu18_v;
   Int_t           HLT_BIT_HLT_IsoMu20_v;
   Int_t           HLT_BIT_HLT_IsoTkMu20_v;
   Int_t           HLT_BIT_HLT_IsoMu27_v;
   Int_t           HLT_BIT_HLT_IsoTkMu27_v;
   Int_t           HLT_SingleMu;
   Int_t           HLT_BIT_HLT_Ele23_WPLoose_Gsf_v;
   Int_t           HLT_BIT_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v;
   Int_t           HLT_BIT_HLT_Ele27_WPLoose_Gsf_v;
   Int_t           HLT_BIT_HLT_Ele27_eta2p1_WPLoose_Gsf_v;
   Int_t           HLT_BIT_HLT_Ele32_eta2p1_WPLoose_Gsf_v;
   Int_t           HLT_BIT_HLT_Ele27_WP85_Gsf_v;
   Int_t           HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_v;
   Int_t           HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_v;
   Int_t           HLT_SingleEl;
   Int_t           Flag_hcalLaserEventFilter;
   Int_t           Flag_trkPOG_logErrorTooManyClusters;
   Int_t           Flag_EcalDeadCellTriggerPrimitiveFilter;
   Int_t           Flag_trkPOGFilters;
   Int_t           Flag_trackingFailureFilter;
   Int_t           Flag_CSCTightHaloFilter;
   Int_t           Flag_HBHENoiseFilter;
   Int_t           Flag_HBHENoiseIsoFilter;
   Int_t           Flag_goodVertices;
   Int_t           Flag_trkPOG_manystripclus53X;
   Int_t           Flag_METFilters;
   Int_t           Flag_ecalLaserCorrFilter;
   Int_t           Flag_trkPOG_toomanystripclus53X;
   Int_t           Flag_CSCTightHalo2015Filter;
   Int_t           Flag_eeBadScFilter;
   Int_t           Flag_globalTightHalo2016Filter;
   Float_t         met_pt;
   Float_t         met_eta;
   Float_t         met_phi;
   Float_t         met_mass;
   Float_t         met_sumEt;
   Float_t         met_rawPt;
   Float_t         met_rawPhi;
   Float_t         met_rawSumEt;
   Float_t         met_genPt;
   Float_t         met_genPhi;
   Float_t         met_genEta;
   Float_t         metNoMu_pt;
   Float_t         metNoMu_eta;
   Float_t         metNoMu_phi;
   Float_t         metNoMu_mass;
   Float_t         metNoMu_sumEt;
   Float_t         metNoMu_rawPt;
   Float_t         metNoMu_rawPhi;
   Float_t         metNoMu_rawSumEt;
   Float_t         metNoMu_genPt;
   Float_t         metNoMu_genPhi;
   Float_t         metNoMu_genEta;
   Float_t         metPuppi_pt;
   Float_t         metPuppi_eta;
   Float_t         metPuppi_phi;
   Float_t         metPuppi_mass;
   Float_t         metPuppi_sumEt;
   Float_t         metPuppi_rawPt;
   Float_t         metPuppi_rawPhi;
   Float_t         metPuppi_rawSumEt;
   Float_t         metPuppi_genPt;
   Float_t         metPuppi_genPhi;
   Float_t         metPuppi_genEta;
   Float_t         metNoPU_pt;
   Float_t         metNoPU_eta;
   Float_t         metNoPU_phi;
   Float_t         metNoPU_mass;
   Int_t           nFatJet;
   Int_t           FatJet_id[6];   //[nFatJet]
   Int_t           FatJet_puId[6];   //[nFatJet]
   Float_t         FatJet_btagCSV[6];   //[nFatJet]
   Float_t         FatJet_btagCMVA[6];   //[nFatJet]
   Float_t         FatJet_rawPt[6];   //[nFatJet]
   Float_t         FatJet_mcPt[6];   //[nFatJet]
   Int_t           FatJet_mcFlavour[6];   //[nFatJet]
   Int_t           FatJet_partonFlavour[6];   //[nFatJet]
   Int_t           FatJet_hadronFlavour[6];   //[nFatJet]
   Int_t           FatJet_mcMatchId[6];   //[nFatJet]
   Float_t         FatJet_corr_JECUp[6];   //[nFatJet]
   Float_t         FatJet_corr_JECDown[6];   //[nFatJet]
   Float_t         FatJet_corr[6];   //[nFatJet]
   Float_t         FatJet_corr_JERUp[6];   //[nFatJet]
   Float_t         FatJet_corr_JERDown[6];   //[nFatJet]
   Float_t         FatJet_corr_JER[6];   //[nFatJet]
   Float_t         FatJet_pt[6];   //[nFatJet]
   Float_t         FatJet_eta[6];   //[nFatJet]
   Float_t         FatJet_phi[6];   //[nFatJet]
   Float_t         FatJet_mass[6];   //[nFatJet]
   Float_t         FatJet_softDropMass[6];   //[nFatJet]
   Float_t         FatJet_tau1[6];   //[nFatJet]
   Float_t         FatJet_tau2[6];   //[nFatJet]
   Float_t         FatJet_tau3[6];   //[nFatJet]
   Float_t         FatJet_topMass[6];   //[nFatJet]
   Float_t         FatJet_minMass[6];   //[nFatJet]
   Float_t         FatJet_nSubJets[6];   //[nFatJet]
   Float_t         FatJet_puMva[6];   //[nFatJet]
   Float_t         FatJet_prunedMass[6];   //[nFatJet]
   Float_t         FatJet_chHEF[6];   //[nFatJet]
   Float_t         FatJet_neHEF[6];   //[nFatJet]
   Float_t         FatJet_chEmEF[6];   //[nFatJet]
   Float_t         FatJet_neEmEF[6];   //[nFatJet]
   Float_t         FatJet_chMult[6];   //[nFatJet]
   Float_t         FatJet_neMult[6];   //[nFatJet]
   Float_t         FatJet_puppiPt[6];   //[nFatJet]
   Float_t         FatJet_puppiEta[6];   //[nFatJet]
   Float_t         FatJet_puppiPhi[6];   //[nFatJet]
   Float_t         FatJet_puppiMass[6];   //[nFatJet]
   Float_t         FatJet_puppiTau1[6];   //[nFatJet]
   Float_t         FatJet_puppiTau2[6];   //[nFatJet]
   Float_t         FatJet_puppiTau3[6];   //[nFatJet]
   Int_t           nisoTrack;
   Int_t           isoTrack_pdgId[4];   //[nisoTrack]
   Float_t         isoTrack_pt[4];   //[nisoTrack]
   Float_t         isoTrack_eta[4];   //[nisoTrack]
   Float_t         isoTrack_phi[4];   //[nisoTrack]
   Float_t         isoTrack_mass[4];   //[nisoTrack]
   Int_t           isoTrack_charge[4];   //[nisoTrack]
   Float_t         isoTrack_dz[4];   //[nisoTrack]
   Float_t         isoTrack_absIso[4];   //[nisoTrack]
   Float_t         isoTrack_relIsoAn04[4];   //[nisoTrack]
   Int_t           isoTrack_mcMatchId[4];   //[nisoTrack]
   Int_t           ncustomPuppiAK8;
   Float_t         customPuppiAK8_pt[5];   //[ncustomPuppiAK8]
   Float_t         customPuppiAK8_eta[5];   //[ncustomPuppiAK8]
   Float_t         customPuppiAK8_phi[5];   //[ncustomPuppiAK8]
   Float_t         customPuppiAK8_mass[5];   //[ncustomPuppiAK8]
   Float_t         customPuppiAK8_massCorrected[5];   //[ncustomPuppiAK8]
   Int_t           nGenPart;
   Int_t           GenPart_motherId[26];   //[nGenPart]
   Int_t           GenPart_grandmotherId[26];   //[nGenPart]
   Int_t           GenPart_sourceId[26];   //[nGenPart]
   Float_t         GenPart_charge[26];   //[nGenPart]
   Int_t           GenPart_status[26];   //[nGenPart]
   Int_t           GenPart_isPromptHard[26];   //[nGenPart]
   Int_t           GenPart_pdgId[26];   //[nGenPart]
   Float_t         GenPart_pt[26];   //[nGenPart]
   Float_t         GenPart_eta[26];   //[nGenPart]
   Float_t         GenPart_phi[26];   //[nGenPart]
   Float_t         GenPart_mass[26];   //[nGenPart]
   Int_t           GenPart_motherIndex[26];   //[nGenPart]
   Int_t           nLepGood;
   Int_t           LepGood_charge[3];   //[nLepGood]
   Int_t           LepGood_tightId[3];   //[nLepGood]
   Int_t           LepGood_eleCutIdCSA14_25ns_v1[3];   //[nLepGood]
   Int_t           LepGood_eleCutIdCSA14_50ns_v1[3];   //[nLepGood]
   Int_t           LepGood_eleCutIdSpring15_25ns_v1[3];   //[nLepGood]
   Float_t         LepGood_dxy[3];   //[nLepGood]
   Float_t         LepGood_dz[3];   //[nLepGood]
   Float_t         LepGood_edxy[3];   //[nLepGood]
   Float_t         LepGood_edz[3];   //[nLepGood]
   Float_t         LepGood_ip3d[3];   //[nLepGood]
   Float_t         LepGood_sip3d[3];   //[nLepGood]
   Int_t           LepGood_convVeto[3];   //[nLepGood]
   Int_t           LepGood_lostHits[3];   //[nLepGood]
   Float_t         LepGood_relIso03[3];   //[nLepGood]
   Float_t         LepGood_relIso04[3];   //[nLepGood]
   Float_t         LepGood_miniRelIso[3];   //[nLepGood]
   Float_t         LepGood_relIsoAn04[3];   //[nLepGood]
   Int_t           LepGood_tightCharge[3];   //[nLepGood]
   Int_t           LepGood_mcMatchId[3];   //[nLepGood]
   Int_t           LepGood_mcMatchAny[3];   //[nLepGood]
   Int_t           LepGood_mcMatchTau[3];   //[nLepGood]
   Float_t         LepGood_mcPt[3];   //[nLepGood]
   Int_t           LepGood_mediumMuonId[3];   //[nLepGood]
   Int_t           LepGood_ICHEPmediumMuonId[3];   //[nLepGood]
   Int_t           LepGood_pdgId[3];   //[nLepGood]
   Float_t         LepGood_pt[3];   //[nLepGood]
   Float_t         LepGood_eta[3];   //[nLepGood]
   Float_t         LepGood_phi[3];   //[nLepGood]
   Float_t         LepGood_mass[3];   //[nLepGood]
   Int_t           LepGood_looseIdOnly[3];   //[nLepGood]
   Float_t         LepGood_chargedHadRelIso03[3];   //[nLepGood]
   Float_t         LepGood_chargedHadRelIso04[3];   //[nLepGood]
   Int_t           LepGood_softMuonId[3];   //[nLepGood]
   Int_t           LepGood_pfMuonId[3];   //[nLepGood]
   Int_t           LepGood_eleCutId2012_full5x5[3];   //[nLepGood]
   Int_t           LepGood_trackerLayers[3];   //[nLepGood]
   Int_t           LepGood_pixelLayers[3];   //[nLepGood]
   Int_t           LepGood_trackerHits[3];   //[nLepGood]
   Int_t           LepGood_lostOuterHits[3];   //[nLepGood]
   Float_t         LepGood_innerTrackValidHitFraction[3];   //[nLepGood]
   Float_t         LepGood_innerTrackChi2[3];   //[nLepGood]
   Float_t         LepGood_nStations[3];   //[nLepGood]
   Float_t         LepGood_caloCompatibility[3];   //[nLepGood]
   Float_t         LepGood_globalTrackChi2[3];   //[nLepGood]
   Float_t         LepGood_trkKink[3];   //[nLepGood]
   Float_t         LepGood_segmentCompatibility[3];   //[nLepGood]
   Float_t         LepGood_chi2LocalPosition[3];   //[nLepGood]
   Float_t         LepGood_chi2LocalMomentum[3];   //[nLepGood]
   Float_t         LepGood_glbTrackProbability[3];   //[nLepGood]
   Int_t           LepGood_TMOneStationTightMuonId[3];   //[nLepGood]
   Int_t           LepGood_trackHighPurityMuon[3];   //[nLepGood]
   Int_t           LepGood_isGlobalMuon[3];   //[nLepGood]
   Int_t           LepGood_isTrackerMuon[3];   //[nLepGood]
   Float_t         LepGood_sigmaIEtaIEta[3];   //[nLepGood]
   Float_t         LepGood_dEtaScTrkIn[3];   //[nLepGood]
   Float_t         LepGood_dPhiScTrkIn[3];   //[nLepGood]
   Float_t         LepGood_hadronicOverEm[3];   //[nLepGood]
   Float_t         LepGood_eInvMinusPInv[3];   //[nLepGood]
   Float_t         LepGood_eInvMinusPInv_tkMom[3];   //[nLepGood]
   Float_t         LepGood_etaSc[3];   //[nLepGood]
   Int_t           LepGood_mcMatchPdgId[3];   //[nLepGood]
   Float_t         LepGood_r9[3];   //[nLepGood]
   Float_t         LepGood_e5x5[3];   //[nLepGood]
   Float_t         LepGood_sigmaIetaIeta[3];   //[nLepGood]
   Float_t         LepGood_sigmaIphiIphi[3];   //[nLepGood]
   Float_t         LepGood_hcalOverEcal[3];   //[nLepGood]
   Float_t         LepGood_full5x5_e5x5[3];   //[nLepGood]
   Float_t         LepGood_full5x5_r9[3];   //[nLepGood]
   Float_t         LepGood_full5x5_sigmaIetaIeta[3];   //[nLepGood]
   Float_t         LepGood_full5x5_sigmaIphiIphi[3];   //[nLepGood]
   Float_t         LepGood_full5x5_hcalOverEcal[3];   //[nLepGood]
   Float_t         LepGood_correctedEcalEnergy[3];   //[nLepGood]
   Float_t         LepGood_eSuperClusterOverP[3];   //[nLepGood]
   Float_t         LepGood_ecalEnergy[3];   //[nLepGood]
   Float_t         LepGood_superCluster_rawEnergy[3];   //[nLepGood]
   Float_t         LepGood_superCluster_preshowerEnergy[3];   //[nLepGood]
   Float_t         LepGood_superCluster_correctedEnergy[3];   //[nLepGood]
   Float_t         LepGood_superCluster_energy[3];   //[nLepGood]
   Float_t         LepGood_superCluster_clustersSize[3];   //[nLepGood]
   Float_t         LepGood_superCluster_seed_energy[3];   //[nLepGood]
   Float_t         LepGood_e2x5Max[3];   //[nLepGood]
   Float_t         LepGood_e1x5[3];   //[nLepGood]
   Float_t         LepGood_isolTrkPt[3];   //[nLepGood]
   Float_t         LepGood_isolEmHadDepth1[3];   //[nLepGood]
   Float_t         LepGood_eleClusterDEta[3];   //[nLepGood]
   Int_t           LepGood_isEcalDriven[3];   //[nLepGood]
   Float_t         LepGood_muonDB[3];   //[nLepGood]
   Float_t         LepGood_pixelHits[3];   //[nLepGood]
   Float_t         LepGood_muTrackIso[3];   //[nLepGood]
   Float_t         LepGood_muon_dz[3];   //[nLepGood]
   Float_t         LepGood_nChamberHits[3];   //[nLepGood]
   Float_t         LepGood_muonPtRatio[3];   //[nLepGood]
   Int_t           nTauGood;
   Int_t           TauGood_charge[3];   //[nTauGood]
   Int_t           TauGood_decayMode[3];   //[nTauGood]
   Int_t           TauGood_idDecayMode[3];   //[nTauGood]
   Int_t           TauGood_idDecayModeNewDMs[3];   //[nTauGood]
   Float_t         TauGood_dxy[3];   //[nTauGood]
   Float_t         TauGood_dz[3];   //[nTauGood]
   Int_t           TauGood_idMVA[3];   //[nTauGood]
   Int_t           TauGood_idMVANewDM[3];   //[nTauGood]
   Int_t           TauGood_idCI3hit[3];   //[nTauGood]
   Int_t           TauGood_idAntiMu[3];   //[nTauGood]
   Int_t           TauGood_idAntiE[3];   //[nTauGood]
   Float_t         TauGood_isoCI3hit[3];   //[nTauGood]
   Int_t           TauGood_mcMatchId[3];   //[nTauGood]
   Int_t           TauGood_pdgId[3];   //[nTauGood]
   Float_t         TauGood_pt[3];   //[nTauGood]
   Float_t         TauGood_eta[3];   //[nTauGood]
   Float_t         TauGood_phi[3];   //[nTauGood]
   Float_t         TauGood_mass[3];   //[nTauGood]
   Int_t           TauGood_idMVAOldDMRun2[3];   //[nTauGood]
   Int_t           TauGood_idMVAOldDMRun2dR03[3];   //[nTauGood]
   Int_t           ngenJet;
   Int_t           genJet_pdgId[10];   //[ngenJet]
   Float_t         genJet_pt[10];   //[ngenJet]
   Float_t         genJet_eta[10];   //[ngenJet]
   Float_t         genJet_phi[10];   //[ngenJet]
   Float_t         genJet_mass[10];   //[ngenJet]
   Float_t         genJet_charge[10];   //[ngenJet]
   Int_t           genJet_status[10];   //[ngenJet]
   Int_t           genJet_isPromptHard[10];   //[ngenJet]
   Int_t           nSoftDropSubJet;
   Float_t         SoftDropSubJet_pt[9];   //[nSoftDropSubJet]
   Float_t         SoftDropSubJet_eta[9];   //[nSoftDropSubJet]
   Float_t         SoftDropSubJet_phi[9];   //[nSoftDropSubJet]
   Float_t         SoftDropSubJet_mass[9];   //[nSoftDropSubJet]
   Int_t           nPuppiSubJet;
   Float_t         PuppiSubJet_pt[9];   //[nPuppiSubJet]
   Float_t         PuppiSubJet_eta[9];   //[nPuppiSubJet]
   Float_t         PuppiSubJet_phi[9];   //[nPuppiSubJet]
   Float_t         PuppiSubJet_mass[9];   //[nPuppiSubJet]
   Int_t           nJet;
   Float_t         Jet_CorrFactor_L1[16];   //[nJet]
   Float_t         Jet_CorrFactor_L1L2[16];   //[nJet]
   Float_t         Jet_CorrFactor_L1L2L3[16];   //[nJet]
   Float_t         Jet_CorrFactor_L1L2L3Res[16];   //[nJet]
   Int_t           Jet_mcMatchFlav[16];   //[nJet]
   Float_t         Jet_charge[16];   //[nJet]
   Float_t         Jet_ctagCsvL[16];   //[nJet]
   Float_t         Jet_ctagCsvB[16];   //[nJet]
   Float_t         Jet_area[16];   //[nJet]
   Float_t         Jet_qgl[16];   //[nJet]
   Float_t         Jet_ptd[16];   //[nJet]
   Float_t         Jet_axis2[16];   //[nJet]
   Int_t           Jet_mult[16];   //[nJet]
   Int_t           Jet_partonId[16];   //[nJet]
   Int_t           Jet_partonMotherId[16];   //[nJet]
   Float_t         Jet_nLeptons[16];   //[nJet]
   Int_t           Jet_id[16];   //[nJet]
   Int_t           Jet_puId[16];   //[nJet]
   Float_t         Jet_btagCSV[16];   //[nJet]
   Float_t         Jet_btagCMVA[16];   //[nJet]
   Float_t         Jet_rawPt[16];   //[nJet]
   Float_t         Jet_mcPt[16];   //[nJet]
   Int_t           Jet_mcFlavour[16];   //[nJet]
   Int_t           Jet_partonFlavour[16];   //[nJet]
   Int_t           Jet_hadronFlavour[16];   //[nJet]
   Int_t           Jet_mcMatchId[16];   //[nJet]
   Float_t         Jet_corr_JECUp[16];   //[nJet]
   Float_t         Jet_corr_JECDown[16];   //[nJet]
   Float_t         Jet_corr[16];   //[nJet]
   Float_t         Jet_corr_JERUp[16];   //[nJet]
   Float_t         Jet_corr_JERDown[16];   //[nJet]
   Float_t         Jet_corr_JER[16];   //[nJet]
   Float_t         Jet_pt[16];   //[nJet]
   Float_t         Jet_eta[16];   //[nJet]
   Float_t         Jet_phi[16];   //[nJet]
   Float_t         Jet_mass[16];   //[nJet]
   Float_t         Jet_prunedMass[16];   //[nJet]
   Int_t           Jet_mcNumPartons[16];   //[nJet]
   Int_t           Jet_mcNumLeptons[16];   //[nJet]
   Int_t           Jet_mcNumTaus[16];   //[nJet]
   Float_t         Jet_mcAnyPartonMass[16];   //[nJet]
   Int_t           Jet_nSubJets[16];   //[nJet]
   Int_t           Jet_nSubJets25[16];   //[nJet]
   Int_t           Jet_nSubJets30[16];   //[nJet]
   Int_t           Jet_nSubJets40[16];   //[nJet]
   Int_t           Jet_nSubJetsZ01[16];   //[nJet]
   Float_t         Jet_chHEF[16];   //[nJet]
   Float_t         Jet_neHEF[16];   //[nJet]
   Float_t         Jet_phEF[16];   //[nJet]
   Float_t         Jet_eEF[16];   //[nJet]
   Float_t         Jet_muEF[16];   //[nJet]
   Float_t         Jet_HFHEF[16];   //[nJet]
   Float_t         Jet_HFEMEF[16];   //[nJet]
   Int_t           Jet_chHMult[16];   //[nJet]
   Int_t           Jet_neHMult[16];   //[nJet]
   Int_t           Jet_phMult[16];   //[nJet]
   Int_t           Jet_eMult[16];   //[nJet]
   Int_t           Jet_muMult[16];   //[nJet]
   Int_t           Jet_HFHMult[16];   //[nJet]
   Int_t           Jet_HFEMMult[16];   //[nJet]
   Float_t         Jet_puMva[16];   //[nJet]
   Int_t           ngenLep;
   Int_t           genLep_motherId[10];   //[ngenLep]
   Int_t           genLep_grandmotherId[10];   //[ngenLep]
   Int_t           genLep_sourceId[10];   //[ngenLep]
   Float_t         genLep_charge[10];   //[ngenLep]
   Int_t           genLep_status[10];   //[ngenLep]
   Int_t           genLep_isPromptHard[10];   //[ngenLep]
   Int_t           genLep_pdgId[10];   //[ngenLep]
   Float_t         genLep_pt[10];   //[ngenLep]
   Float_t         genLep_eta[10];   //[ngenLep]
   Float_t         genLep_phi[10];   //[ngenLep]
   Float_t         genLep_mass[10];   //[ngenLep]
   Int_t           genLep_motherIndex[10];   //[ngenLep]
   Int_t           nLHEweight;
   Int_t           LHEweight_id[1];   //[nLHEweight]
   Float_t         LHEweight_wgt[1];   //[nLHEweight]
   Int_t           ncustomPuppiSoftDropAK8;
   Float_t         customPuppiSoftDropAK8_pt[5];   //[ncustomPuppiSoftDropAK8]
   Float_t         customPuppiSoftDropAK8_eta[5];   //[ncustomPuppiSoftDropAK8]
   Float_t         customPuppiSoftDropAK8_phi[5];   //[ncustomPuppiSoftDropAK8]
   Float_t         customPuppiSoftDropAK8_mass[5];   //[ncustomPuppiSoftDropAK8]
   Float_t         customPuppiSoftDropAK8_massCorrected[5];   //[ncustomPuppiSoftDropAK8]
   Int_t           nPuppiJet;
   Float_t         PuppiJet_prunedMass[9];   //[nPuppiJet]
   Int_t           PuppiJet_mcNumPartons[9];   //[nPuppiJet]
   Int_t           PuppiJet_mcNumLeptons[9];   //[nPuppiJet]
   Int_t           PuppiJet_mcNumTaus[9];   //[nPuppiJet]
   Float_t         PuppiJet_mcAnyPartonMass[9];   //[nPuppiJet]
   Int_t           PuppiJet_nSubJets[9];   //[nPuppiJet]
   Int_t           PuppiJet_nSubJets25[9];   //[nPuppiJet]
   Int_t           PuppiJet_nSubJets30[9];   //[nPuppiJet]
   Int_t           PuppiJet_nSubJets40[9];   //[nPuppiJet]
   Int_t           PuppiJet_nSubJetsZ01[9];   //[nPuppiJet]
   Float_t         PuppiJet_chHEF[9];   //[nPuppiJet]
   Float_t         PuppiJet_neHEF[9];   //[nPuppiJet]
   Float_t         PuppiJet_phEF[9];   //[nPuppiJet]
   Float_t         PuppiJet_eEF[9];   //[nPuppiJet]
   Float_t         PuppiJet_muEF[9];   //[nPuppiJet]
   Float_t         PuppiJet_HFHEF[9];   //[nPuppiJet]
   Float_t         PuppiJet_HFEMEF[9];   //[nPuppiJet]
   Int_t           PuppiJet_chHMult[9];   //[nPuppiJet]
   Int_t           PuppiJet_neHMult[9];   //[nPuppiJet]
   Int_t           PuppiJet_phMult[9];   //[nPuppiJet]
   Int_t           PuppiJet_eMult[9];   //[nPuppiJet]
   Int_t           PuppiJet_muMult[9];   //[nPuppiJet]
   Int_t           PuppiJet_HFHMult[9];   //[nPuppiJet]
   Int_t           PuppiJet_HFEMMult[9];   //[nPuppiJet]
   Float_t         PuppiJet_puMva[9];   //[nPuppiJet]
   Float_t         PuppiJet_CorrFactor_L1[9];   //[nPuppiJet]
   Float_t         PuppiJet_CorrFactor_L1L2[9];   //[nPuppiJet]
   Float_t         PuppiJet_CorrFactor_L1L2L3[9];   //[nPuppiJet]
   Float_t         PuppiJet_CorrFactor_L1L2L3Res[9];   //[nPuppiJet]
   Int_t           PuppiJet_mcMatchFlav[9];   //[nPuppiJet]
   Float_t         PuppiJet_charge[9];   //[nPuppiJet]
   Float_t         PuppiJet_ctagCsvL[9];   //[nPuppiJet]
   Float_t         PuppiJet_ctagCsvB[9];   //[nPuppiJet]
   Float_t         PuppiJet_area[9];   //[nPuppiJet]
   Float_t         PuppiJet_qgl[9];   //[nPuppiJet]
   Float_t         PuppiJet_ptd[9];   //[nPuppiJet]
   Float_t         PuppiJet_axis2[9];   //[nPuppiJet]
   Int_t           PuppiJet_mult[9];   //[nPuppiJet]
   Int_t           PuppiJet_partonId[9];   //[nPuppiJet]
   Int_t           PuppiJet_partonMotherId[9];   //[nPuppiJet]
   Float_t         PuppiJet_nLeptons[9];   //[nPuppiJet]
   Int_t           PuppiJet_id[9];   //[nPuppiJet]
   Int_t           PuppiJet_puId[9];   //[nPuppiJet]
   Float_t         PuppiJet_btagCSV[9];   //[nPuppiJet]
   Float_t         PuppiJet_btagCMVA[9];   //[nPuppiJet]
   Float_t         PuppiJet_rawPt[9];   //[nPuppiJet]
   Float_t         PuppiJet_mcPt[9];   //[nPuppiJet]
   Int_t           PuppiJet_mcFlavour[9];   //[nPuppiJet]
   Int_t           PuppiJet_partonFlavour[9];   //[nPuppiJet]
   Int_t           PuppiJet_hadronFlavour[9];   //[nPuppiJet]
   Int_t           PuppiJet_mcMatchId[9];   //[nPuppiJet]
   Float_t         PuppiJet_corr_JECUp[9];   //[nPuppiJet]
   Float_t         PuppiJet_corr_JECDown[9];   //[nPuppiJet]
   Float_t         PuppiJet_corr[9];   //[nPuppiJet]
   Float_t         PuppiJet_corr_JERUp[9];   //[nPuppiJet]
   Float_t         PuppiJet_corr_JERDown[9];   //[nPuppiJet]
   Float_t         PuppiJet_corr_JER[9];   //[nPuppiJet]
   Float_t         PuppiJet_pt[9];   //[nPuppiJet]
   Float_t         PuppiJet_eta[9];   //[nPuppiJet]
   Float_t         PuppiJet_phi[9];   //[nPuppiJet]
   Float_t         PuppiJet_mass[9];   //[nPuppiJet]
   Float_t         PuppiJet_massCorrected[9];   //[nPuppiJet]
   Int_t           ngenFatJet;
   Int_t           genFatJet_pdgId[7];   //[ngenFatJet]
   Float_t         genFatJet_pt[7];   //[ngenFatJet]
   Float_t         genFatJet_eta[7];   //[ngenFatJet]
   Float_t         genFatJet_phi[7];   //[ngenFatJet]
   Float_t         genFatJet_mass[7];   //[ngenFatJet]
   Float_t         genFatJet_charge[7];   //[ngenFatJet]
   Int_t           genFatJet_status[7];   //[ngenFatJet]
   Int_t           genFatJet_isPromptHard[7];   //[ngenFatJet]
   Int_t           nJetFwd;
   Float_t         JetFwd_CorrFactor_L1[5];   //[nJetFwd]
   Float_t         JetFwd_CorrFactor_L1L2[5];   //[nJetFwd]
   Float_t         JetFwd_CorrFactor_L1L2L3[5];   //[nJetFwd]
   Float_t         JetFwd_CorrFactor_L1L2L3Res[5];   //[nJetFwd]
   Int_t           JetFwd_mcMatchFlav[5];   //[nJetFwd]
   Float_t         JetFwd_charge[5];   //[nJetFwd]
   Float_t         JetFwd_ctagCsvL[5];   //[nJetFwd]
   Float_t         JetFwd_ctagCsvB[5];   //[nJetFwd]
   Float_t         JetFwd_area[5];   //[nJetFwd]
   Float_t         JetFwd_qgl[5];   //[nJetFwd]
   Float_t         JetFwd_ptd[5];   //[nJetFwd]
   Float_t         JetFwd_axis2[5];   //[nJetFwd]
   Int_t           JetFwd_mult[5];   //[nJetFwd]
   Int_t           JetFwd_partonId[5];   //[nJetFwd]
   Int_t           JetFwd_partonMotherId[5];   //[nJetFwd]
   Float_t         JetFwd_nLeptons[5];   //[nJetFwd]
   Int_t           JetFwd_id[5];   //[nJetFwd]
   Int_t           JetFwd_puId[5];   //[nJetFwd]
   Float_t         JetFwd_btagCSV[5];   //[nJetFwd]
   Float_t         JetFwd_btagCMVA[5];   //[nJetFwd]
   Float_t         JetFwd_rawPt[5];   //[nJetFwd]
   Float_t         JetFwd_mcPt[5];   //[nJetFwd]
   Int_t           JetFwd_mcFlavour[5];   //[nJetFwd]
   Int_t           JetFwd_partonFlavour[5];   //[nJetFwd]
   Int_t           JetFwd_hadronFlavour[5];   //[nJetFwd]
   Int_t           JetFwd_mcMatchId[5];   //[nJetFwd]
   Float_t         JetFwd_corr_JECUp[5];   //[nJetFwd]
   Float_t         JetFwd_corr_JECDown[5];   //[nJetFwd]
   Float_t         JetFwd_corr[5];   //[nJetFwd]
   Float_t         JetFwd_corr_JERUp[5];   //[nJetFwd]
   Float_t         JetFwd_corr_JERDown[5];   //[nJetFwd]
   Float_t         JetFwd_corr_JER[5];   //[nJetFwd]
   Float_t         JetFwd_pt[5];   //[nJetFwd]
   Float_t         JetFwd_eta[5];   //[nJetFwd]
   Float_t         JetFwd_phi[5];   //[nJetFwd]
   Float_t         JetFwd_mass[5];   //[nJetFwd]
   Float_t         JetFwd_prunedMass[5];   //[nJetFwd]
   Int_t           JetFwd_mcNumPartons[5];   //[nJetFwd]
   Int_t           JetFwd_mcNumLeptons[5];   //[nJetFwd]
   Int_t           JetFwd_mcNumTaus[5];   //[nJetFwd]
   Float_t         JetFwd_mcAnyPartonMass[5];   //[nJetFwd]
   Int_t           JetFwd_nSubJets[5];   //[nJetFwd]
   Int_t           JetFwd_nSubJets25[5];   //[nJetFwd]
   Int_t           JetFwd_nSubJets30[5];   //[nJetFwd]
   Int_t           JetFwd_nSubJets40[5];   //[nJetFwd]
   Int_t           JetFwd_nSubJetsZ01[5];   //[nJetFwd]
   Float_t         JetFwd_chHEF[5];   //[nJetFwd]
   Float_t         JetFwd_neHEF[5];   //[nJetFwd]
   Float_t         JetFwd_phEF[5];   //[nJetFwd]
   Float_t         JetFwd_eEF[5];   //[nJetFwd]
   Float_t         JetFwd_muEF[5];   //[nJetFwd]
   Float_t         JetFwd_HFHEF[5];   //[nJetFwd]
   Float_t         JetFwd_HFEMEF[5];   //[nJetFwd]
   Int_t           JetFwd_chHMult[5];   //[nJetFwd]
   Int_t           JetFwd_neHMult[5];   //[nJetFwd]
   Int_t           JetFwd_phMult[5];   //[nJetFwd]
   Int_t           JetFwd_eMult[5];   //[nJetFwd]
   Int_t           JetFwd_muMult[5];   //[nJetFwd]
   Int_t           JetFwd_HFHMult[5];   //[nJetFwd]
   Int_t           JetFwd_HFEMMult[5];   //[nJetFwd]
   Float_t         JetFwd_puMva[5];   //[nJetFwd]
   Int_t           nGammaGood;
   Float_t         GammaGood_etaSc[3];   //[nGammaGood]
   Int_t           GammaGood_idCutBased[3];   //[nGammaGood]
   Float_t         GammaGood_hOverE[3];   //[nGammaGood]
   Float_t         GammaGood_r9[3];   //[nGammaGood]
   Float_t         GammaGood_sigmaIetaIeta[3];   //[nGammaGood]
   Float_t         GammaGood_chHadIso04[3];   //[nGammaGood]
   Float_t         GammaGood_chHadIso[3];   //[nGammaGood]
   Float_t         GammaGood_phIso[3];   //[nGammaGood]
   Float_t         GammaGood_neuHadIso[3];   //[nGammaGood]
   Float_t         GammaGood_relIso[3];   //[nGammaGood]
   Int_t           GammaGood_mcMatchId[3];   //[nGammaGood]
   Float_t         GammaGood_mcPt[3];   //[nGammaGood]
   Int_t           GammaGood_pdgId[3];   //[nGammaGood]
   Float_t         GammaGood_pt[3];   //[nGammaGood]
   Float_t         GammaGood_eta[3];   //[nGammaGood]
   Float_t         GammaGood_phi[3];   //[nGammaGood]
   Float_t         GammaGood_mass[3];   //[nGammaGood]
   Float_t         GammaGood_genIso04[3];   //[nGammaGood]
   Float_t         GammaGood_genIso03[3];   //[nGammaGood]
   Float_t         GammaGood_chHadIsoRC04[3];   //[nGammaGood]
   Float_t         GammaGood_chHadIsoRC[3];   //[nGammaGood]
   Float_t         GammaGood_drMinParton[3];   //[nGammaGood]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_xsec;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nVert;   //!
   TBranch        *b_nJet25;   //!
   TBranch        *b_nBJetLoose25;   //!
   TBranch        *b_nBJetMedium25;   //!
   TBranch        *b_nBJetTight25;   //!
   TBranch        *b_nJet30;   //!
   TBranch        *b_nJet30a;   //!
   TBranch        *b_nBJetLoose30;   //!
   TBranch        *b_nBJetMedium30;   //!
   TBranch        *b_nBJetTight30;   //!
   TBranch        *b_nJet40;   //!
   TBranch        *b_nJet40a;   //!
   TBranch        *b_nBJetLoose40;   //!
   TBranch        *b_nBJetMedium40;   //!
   TBranch        *b_nBJetTight40;   //!
   TBranch        *b_nLepGood20;   //!
   TBranch        *b_nLepGood15;   //!
   TBranch        *b_nLepGood10;   //!
   TBranch        *b_htJet25;   //!
   TBranch        *b_mhtJet25;   //!
   TBranch        *b_htJet40j;   //!
   TBranch        *b_htJet40ja;   //!
   TBranch        *b_htJet40;   //!
   TBranch        *b_htJet40a;   //!
   TBranch        *b_mhtJet40;   //!
   TBranch        *b_mhtJet40a;   //!
   TBranch        *b_apcjetmetmin;   //!
   TBranch        *b_mZ1;   //!
   TBranch        *b_mZ1SFSS;   //!
   TBranch        *b_minMllSFOS;   //!
   TBranch        *b_maxMllSFOS;   //!
   TBranch        *b_minMllAFOS;   //!
   TBranch        *b_maxMllAFOS;   //!
   TBranch        *b_minMllAFSS;   //!
   TBranch        *b_maxMllAFSS;   //!
   TBranch        *b_minMllAFAS;   //!
   TBranch        *b_maxMllAFAS;   //!
   TBranch        *b_m2l;   //!
   TBranch        *b_nBJet25;   //!
   TBranch        *b_nFatJet100;   //!
   TBranch        *b_nMuons10;   //!
   TBranch        *b_nElectrons10;   //!
   TBranch        *b_nTaus20;   //!
   TBranch        *b_nGammas20;   //!
   TBranch        *b_rhoCN;   //!
   TBranch        *b_hbheFilterNew50ns;   //!
   TBranch        *b_hbheFilterNew25ns;   //!
   TBranch        *b_hbheFilterIso;   //!
   TBranch        *b_met_trkPt;   //!
   TBranch        *b_met_trkPhi;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET170_NoiseCleaned_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET170_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET170_HBHECleaned_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET170_JetIdCleaned_v;   //!
   TBranch        *b_HLT_Met170;   //!
   TBranch        *b_HLT_BIT_HLT_Photon155_v;   //!
   TBranch        *b_HLT_BIT_HLT_Photon165_HE10_v;   //!
   TBranch        *b_HLT_BIT_HLT_Photon175_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet40_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet60_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet80_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet140_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet200_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet260_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet320_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet400_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet450_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet500_v;   //!
   TBranch        *b_HLT_SinglePho;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;   //!
   TBranch        *b_HLT_MonoJetMetNoMuMHT120;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v;   //!
   TBranch        *b_HLT_MonoJetMetNoMuMHT90;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT800_v;   //!
   TBranch        *b_HLT_BIT_HLT_ECALHT800_v;   //!
   TBranch        *b_HLT_Jet_HT;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET300_NoiseCleaned_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET300_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET300_JetIdCleaned_v;   //!
   TBranch        *b_HLT_Met300;   //!
   TBranch        *b_HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_HLT_DoubleEl;   //!
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v;   //!
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v;   //!
   TBranch        *b_HLT_DoubleMu;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu24_eta2p1_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoTkMu24_eta2p1_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu18_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu20_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoTkMu20_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu27_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoTkMu27_v;   //!
   TBranch        *b_HLT_SingleMu;   //!
   TBranch        *b_HLT_BIT_HLT_Ele23_WPLoose_Gsf_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele27_WPLoose_Gsf_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele27_eta2p1_WPLoose_Gsf_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele32_eta2p1_WPLoose_Gsf_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele27_WP85_Gsf_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_v;   //!
   TBranch        *b_HLT_SingleEl;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_trackingFailureFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_METFilters;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_mass;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_met_rawPt;   //!
   TBranch        *b_met_rawPhi;   //!
   TBranch        *b_met_rawSumEt;   //!
   TBranch        *b_met_genPt;   //!
   TBranch        *b_met_genPhi;   //!
   TBranch        *b_met_genEta;   //!
   TBranch        *b_metNoMu_pt;   //!
   TBranch        *b_metNoMu_eta;   //!
   TBranch        *b_metNoMu_phi;   //!
   TBranch        *b_metNoMu_mass;   //!
   TBranch        *b_metNoMu_sumEt;   //!
   TBranch        *b_metNoMu_rawPt;   //!
   TBranch        *b_metNoMu_rawPhi;   //!
   TBranch        *b_metNoMu_rawSumEt;   //!
   TBranch        *b_metNoMu_genPt;   //!
   TBranch        *b_metNoMu_genPhi;   //!
   TBranch        *b_metNoMu_genEta;   //!
   TBranch        *b_metPuppi_pt;   //!
   TBranch        *b_metPuppi_eta;   //!
   TBranch        *b_metPuppi_phi;   //!
   TBranch        *b_metPuppi_mass;   //!
   TBranch        *b_metPuppi_sumEt;   //!
   TBranch        *b_metPuppi_rawPt;   //!
   TBranch        *b_metPuppi_rawPhi;   //!
   TBranch        *b_metPuppi_rawSumEt;   //!
   TBranch        *b_metPuppi_genPt;   //!
   TBranch        *b_metPuppi_genPhi;   //!
   TBranch        *b_metPuppi_genEta;   //!
   TBranch        *b_metNoPU_pt;   //!
   TBranch        *b_metNoPU_eta;   //!
   TBranch        *b_metNoPU_phi;   //!
   TBranch        *b_metNoPU_mass;   //!
   TBranch        *b_nFatJet;   //!
   TBranch        *b_FatJet_id;   //!
   TBranch        *b_FatJet_puId;   //!
   TBranch        *b_FatJet_btagCSV;   //!
   TBranch        *b_FatJet_btagCMVA;   //!
   TBranch        *b_FatJet_rawPt;   //!
   TBranch        *b_FatJet_mcPt;   //!
   TBranch        *b_FatJet_mcFlavour;   //!
   TBranch        *b_FatJet_partonFlavour;   //!
   TBranch        *b_FatJet_hadronFlavour;   //!
   TBranch        *b_FatJet_mcMatchId;   //!
   TBranch        *b_FatJet_corr_JECUp;   //!
   TBranch        *b_FatJet_corr_JECDown;   //!
   TBranch        *b_FatJet_corr;   //!
   TBranch        *b_FatJet_corr_JERUp;   //!
   TBranch        *b_FatJet_corr_JERDown;   //!
   TBranch        *b_FatJet_corr_JER;   //!
   TBranch        *b_FatJet_pt;   //!
   TBranch        *b_FatJet_eta;   //!
   TBranch        *b_FatJet_phi;   //!
   TBranch        *b_FatJet_mass;   //!
   TBranch        *b_FatJet_softDropMass;   //!
   TBranch        *b_FatJet_tau1;   //!
   TBranch        *b_FatJet_tau2;   //!
   TBranch        *b_FatJet_tau3;   //!
   TBranch        *b_FatJet_topMass;   //!
   TBranch        *b_FatJet_minMass;   //!
   TBranch        *b_FatJet_nSubJets;   //!
   TBranch        *b_FatJet_puMva;   //!
   TBranch        *b_FatJet_prunedMass;   //!
   TBranch        *b_FatJet_chHEF;   //!
   TBranch        *b_FatJet_neHEF;   //!
   TBranch        *b_FatJet_chEmEF;   //!
   TBranch        *b_FatJet_neEmEF;   //!
   TBranch        *b_FatJet_chMult;   //!
   TBranch        *b_FatJet_neMult;   //!
   TBranch        *b_FatJet_puppiPt;   //!
   TBranch        *b_FatJet_puppiEta;   //!
   TBranch        *b_FatJet_puppiPhi;   //!
   TBranch        *b_FatJet_puppiMass;   //!
   TBranch        *b_FatJet_puppiTau1;   //!
   TBranch        *b_FatJet_puppiTau2;   //!
   TBranch        *b_FatJet_puppiTau3;   //!
   TBranch        *b_nisoTrack;   //!
   TBranch        *b_isoTrack_pdgId;   //!
   TBranch        *b_isoTrack_pt;   //!
   TBranch        *b_isoTrack_eta;   //!
   TBranch        *b_isoTrack_phi;   //!
   TBranch        *b_isoTrack_mass;   //!
   TBranch        *b_isoTrack_charge;   //!
   TBranch        *b_isoTrack_dz;   //!
   TBranch        *b_isoTrack_absIso;   //!
   TBranch        *b_isoTrack_relIsoAn04;   //!
   TBranch        *b_isoTrack_mcMatchId;   //!
   TBranch        *b_ncustomPuppiAK8;   //!
   TBranch        *b_customPuppiAK8_pt;   //!
   TBranch        *b_customPuppiAK8_eta;   //!
   TBranch        *b_customPuppiAK8_phi;   //!
   TBranch        *b_customPuppiAK8_mass;   //!
   TBranch        *b_customPuppiAK8_massCorrected;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_motherId;   //!
   TBranch        *b_GenPart_grandmotherId;   //!
   TBranch        *b_GenPart_sourceId;   //!
   TBranch        *b_GenPart_charge;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_isPromptHard;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_motherIndex;   //!
   TBranch        *b_nLepGood;   //!
   TBranch        *b_LepGood_charge;   //!
   TBranch        *b_LepGood_tightId;   //!
   TBranch        *b_LepGood_eleCutIdCSA14_25ns_v1;   //!
   TBranch        *b_LepGood_eleCutIdCSA14_50ns_v1;   //!
   TBranch        *b_LepGood_eleCutIdSpring15_25ns_v1;   //!
   TBranch        *b_LepGood_dxy;   //!
   TBranch        *b_LepGood_dz;   //!
   TBranch        *b_LepGood_edxy;   //!
   TBranch        *b_LepGood_edz;   //!
   TBranch        *b_LepGood_ip3d;   //!
   TBranch        *b_LepGood_sip3d;   //!
   TBranch        *b_LepGood_convVeto;   //!
   TBranch        *b_LepGood_lostHits;   //!
   TBranch        *b_LepGood_relIso03;   //!
   TBranch        *b_LepGood_relIso04;   //!
   TBranch        *b_LepGood_miniRelIso;   //!
   TBranch        *b_LepGood_relIsoAn04;   //!
   TBranch        *b_LepGood_tightCharge;   //!
   TBranch        *b_LepGood_mcMatchId;   //!
   TBranch        *b_LepGood_mcMatchAny;   //!
   TBranch        *b_LepGood_mcMatchTau;   //!
   TBranch        *b_LepGood_mcPt;   //!
   TBranch        *b_LepGood_mediumMuonId;   //!
   TBranch        *b_LepGood_ICHEPmediumMuonId;   //!
   TBranch        *b_LepGood_pdgId;   //!
   TBranch        *b_LepGood_pt;   //!
   TBranch        *b_LepGood_eta;   //!
   TBranch        *b_LepGood_phi;   //!
   TBranch        *b_LepGood_mass;   //!
   TBranch        *b_LepGood_looseIdOnly;   //!
   TBranch        *b_LepGood_chargedHadRelIso03;   //!
   TBranch        *b_LepGood_chargedHadRelIso04;   //!
   TBranch        *b_LepGood_softMuonId;   //!
   TBranch        *b_LepGood_pfMuonId;   //!
   TBranch        *b_LepGood_eleCutId2012_full5x5;   //!
   TBranch        *b_LepGood_trackerLayers;   //!
   TBranch        *b_LepGood_pixelLayers;   //!
   TBranch        *b_LepGood_trackerHits;   //!
   TBranch        *b_LepGood_lostOuterHits;   //!
   TBranch        *b_LepGood_innerTrackValidHitFraction;   //!
   TBranch        *b_LepGood_innerTrackChi2;   //!
   TBranch        *b_LepGood_nStations;   //!
   TBranch        *b_LepGood_caloCompatibility;   //!
   TBranch        *b_LepGood_globalTrackChi2;   //!
   TBranch        *b_LepGood_trkKink;   //!
   TBranch        *b_LepGood_segmentCompatibility;   //!
   TBranch        *b_LepGood_chi2LocalPosition;   //!
   TBranch        *b_LepGood_chi2LocalMomentum;   //!
   TBranch        *b_LepGood_glbTrackProbability;   //!
   TBranch        *b_LepGood_TMOneStationTightMuonId;   //!
   TBranch        *b_LepGood_trackHighPurityMuon;   //!
   TBranch        *b_LepGood_isGlobalMuon;   //!
   TBranch        *b_LepGood_isTrackerMuon;   //!
   TBranch        *b_LepGood_sigmaIEtaIEta;   //!
   TBranch        *b_LepGood_dEtaScTrkIn;   //!
   TBranch        *b_LepGood_dPhiScTrkIn;   //!
   TBranch        *b_LepGood_hadronicOverEm;   //!
   TBranch        *b_LepGood_eInvMinusPInv;   //!
   TBranch        *b_LepGood_eInvMinusPInv_tkMom;   //!
   TBranch        *b_LepGood_etaSc;   //!
   TBranch        *b_LepGood_mcMatchPdgId;   //!
   TBranch        *b_LepGood_r9;   //!
   TBranch        *b_LepGood_e5x5;   //!
   TBranch        *b_LepGood_sigmaIetaIeta;   //!
   TBranch        *b_LepGood_sigmaIphiIphi;   //!
   TBranch        *b_LepGood_hcalOverEcal;   //!
   TBranch        *b_LepGood_full5x5_e5x5;   //!
   TBranch        *b_LepGood_full5x5_r9;   //!
   TBranch        *b_LepGood_full5x5_sigmaIetaIeta;   //!
   TBranch        *b_LepGood_full5x5_sigmaIphiIphi;   //!
   TBranch        *b_LepGood_full5x5_hcalOverEcal;   //!
   TBranch        *b_LepGood_correctedEcalEnergy;   //!
   TBranch        *b_LepGood_eSuperClusterOverP;   //!
   TBranch        *b_LepGood_ecalEnergy;   //!
   TBranch        *b_LepGood_superCluster_rawEnergy;   //!
   TBranch        *b_LepGood_superCluster_preshowerEnergy;   //!
   TBranch        *b_LepGood_superCluster_correctedEnergy;   //!
   TBranch        *b_LepGood_superCluster_energy;   //!
   TBranch        *b_LepGood_superCluster_clustersSize;   //!
   TBranch        *b_LepGood_superCluster_seed_energy;   //!
   TBranch        *b_LepGood_e2x5Max;   //!
   TBranch        *b_LepGood_e1x5;   //!
   TBranch        *b_LepGood_isolTrkPt;   //!
   TBranch        *b_LepGood_isolEmHadDepth1;   //!
   TBranch        *b_LepGood_eleClusterDEta;   //!
   TBranch        *b_LepGood_isEcalDriven;   //!
   TBranch        *b_LepGood_muonDB;   //!
   TBranch        *b_LepGood_pixelHits;   //!
   TBranch        *b_LepGood_muTrackIso;   //!
   TBranch        *b_LepGood_muon_dz;   //!
   TBranch        *b_LepGood_nChamberHits;   //!
   TBranch        *b_LepGood_muonPtRatio;   //!
   TBranch        *b_nTauGood;   //!
   TBranch        *b_TauGood_charge;   //!
   TBranch        *b_TauGood_decayMode;   //!
   TBranch        *b_TauGood_idDecayMode;   //!
   TBranch        *b_TauGood_idDecayModeNewDMs;   //!
   TBranch        *b_TauGood_dxy;   //!
   TBranch        *b_TauGood_dz;   //!
   TBranch        *b_TauGood_idMVA;   //!
   TBranch        *b_TauGood_idMVANewDM;   //!
   TBranch        *b_TauGood_idCI3hit;   //!
   TBranch        *b_TauGood_idAntiMu;   //!
   TBranch        *b_TauGood_idAntiE;   //!
   TBranch        *b_TauGood_isoCI3hit;   //!
   TBranch        *b_TauGood_mcMatchId;   //!
   TBranch        *b_TauGood_pdgId;   //!
   TBranch        *b_TauGood_pt;   //!
   TBranch        *b_TauGood_eta;   //!
   TBranch        *b_TauGood_phi;   //!
   TBranch        *b_TauGood_mass;   //!
   TBranch        *b_TauGood_idMVAOldDMRun2;   //!
   TBranch        *b_TauGood_idMVAOldDMRun2dR03;   //!
   TBranch        *b_ngenJet;   //!
   TBranch        *b_genJet_pdgId;   //!
   TBranch        *b_genJet_pt;   //!
   TBranch        *b_genJet_eta;   //!
   TBranch        *b_genJet_phi;   //!
   TBranch        *b_genJet_mass;   //!
   TBranch        *b_genJet_charge;   //!
   TBranch        *b_genJet_status;   //!
   TBranch        *b_genJet_isPromptHard;   //!
   TBranch        *b_nSoftDropSubJet;   //!
   TBranch        *b_SoftDropSubJet_pt;   //!
   TBranch        *b_SoftDropSubJet_eta;   //!
   TBranch        *b_SoftDropSubJet_phi;   //!
   TBranch        *b_SoftDropSubJet_mass;   //!
   TBranch        *b_nPuppiSubJet;   //!
   TBranch        *b_PuppiSubJet_pt;   //!
   TBranch        *b_PuppiSubJet_eta;   //!
   TBranch        *b_PuppiSubJet_phi;   //!
   TBranch        *b_PuppiSubJet_mass;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_CorrFactor_L1;   //!
   TBranch        *b_Jet_CorrFactor_L1L2;   //!
   TBranch        *b_Jet_CorrFactor_L1L2L3;   //!
   TBranch        *b_Jet_CorrFactor_L1L2L3Res;   //!
   TBranch        *b_Jet_mcMatchFlav;   //!
   TBranch        *b_Jet_charge;   //!
   TBranch        *b_Jet_ctagCsvL;   //!
   TBranch        *b_Jet_ctagCsvB;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_qgl;   //!
   TBranch        *b_Jet_ptd;   //!
   TBranch        *b_Jet_axis2;   //!
   TBranch        *b_Jet_mult;   //!
   TBranch        *b_Jet_partonId;   //!
   TBranch        *b_Jet_partonMotherId;   //!
   TBranch        *b_Jet_nLeptons;   //!
   TBranch        *b_Jet_id;   //!
   TBranch        *b_Jet_puId;   //!
   TBranch        *b_Jet_btagCSV;   //!
   TBranch        *b_Jet_btagCMVA;   //!
   TBranch        *b_Jet_rawPt;   //!
   TBranch        *b_Jet_mcPt;   //!
   TBranch        *b_Jet_mcFlavour;   //!
   TBranch        *b_Jet_partonFlavour;   //!
   TBranch        *b_Jet_hadronFlavour;   //!
   TBranch        *b_Jet_mcMatchId;   //!
   TBranch        *b_Jet_corr_JECUp;   //!
   TBranch        *b_Jet_corr_JECDown;   //!
   TBranch        *b_Jet_corr;   //!
   TBranch        *b_Jet_corr_JERUp;   //!
   TBranch        *b_Jet_corr_JERDown;   //!
   TBranch        *b_Jet_corr_JER;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_prunedMass;   //!
   TBranch        *b_Jet_mcNumPartons;   //!
   TBranch        *b_Jet_mcNumLeptons;   //!
   TBranch        *b_Jet_mcNumTaus;   //!
   TBranch        *b_Jet_mcAnyPartonMass;   //!
   TBranch        *b_Jet_nSubJets;   //!
   TBranch        *b_Jet_nSubJets25;   //!
   TBranch        *b_Jet_nSubJets30;   //!
   TBranch        *b_Jet_nSubJets40;   //!
   TBranch        *b_Jet_nSubJetsZ01;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_phEF;   //!
   TBranch        *b_Jet_eEF;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_HFHEF;   //!
   TBranch        *b_Jet_HFEMEF;   //!
   TBranch        *b_Jet_chHMult;   //!
   TBranch        *b_Jet_neHMult;   //!
   TBranch        *b_Jet_phMult;   //!
   TBranch        *b_Jet_eMult;   //!
   TBranch        *b_Jet_muMult;   //!
   TBranch        *b_Jet_HFHMult;   //!
   TBranch        *b_Jet_HFEMMult;   //!
   TBranch        *b_Jet_puMva;   //!
   TBranch        *b_ngenLep;   //!
   TBranch        *b_genLep_motherId;   //!
   TBranch        *b_genLep_grandmotherId;   //!
   TBranch        *b_genLep_sourceId;   //!
   TBranch        *b_genLep_charge;   //!
   TBranch        *b_genLep_status;   //!
   TBranch        *b_genLep_isPromptHard;   //!
   TBranch        *b_genLep_pdgId;   //!
   TBranch        *b_genLep_pt;   //!
   TBranch        *b_genLep_eta;   //!
   TBranch        *b_genLep_phi;   //!
   TBranch        *b_genLep_mass;   //!
   TBranch        *b_genLep_motherIndex;   //!
   TBranch        *b_nLHEweight;   //!
   TBranch        *b_LHEweight_id;   //!
   TBranch        *b_LHEweight_wgt;   //!
   TBranch        *b_ncustomPuppiSoftDropAK8;   //!
   TBranch        *b_customPuppiSoftDropAK8_pt;   //!
   TBranch        *b_customPuppiSoftDropAK8_eta;   //!
   TBranch        *b_customPuppiSoftDropAK8_phi;   //!
   TBranch        *b_customPuppiSoftDropAK8_mass;   //!
   TBranch        *b_customPuppiSoftDropAK8_massCorrected;   //!
   TBranch        *b_nPuppiJet;   //!
   TBranch        *b_PuppiJet_prunedMass;   //!
   TBranch        *b_PuppiJet_mcNumPartons;   //!
   TBranch        *b_PuppiJet_mcNumLeptons;   //!
   TBranch        *b_PuppiJet_mcNumTaus;   //!
   TBranch        *b_PuppiJet_mcAnyPartonMass;   //!
   TBranch        *b_PuppiJet_nSubJets;   //!
   TBranch        *b_PuppiJet_nSubJets25;   //!
   TBranch        *b_PuppiJet_nSubJets30;   //!
   TBranch        *b_PuppiJet_nSubJets40;   //!
   TBranch        *b_PuppiJet_nSubJetsZ01;   //!
   TBranch        *b_PuppiJet_chHEF;   //!
   TBranch        *b_PuppiJet_neHEF;   //!
   TBranch        *b_PuppiJet_phEF;   //!
   TBranch        *b_PuppiJet_eEF;   //!
   TBranch        *b_PuppiJet_muEF;   //!
   TBranch        *b_PuppiJet_HFHEF;   //!
   TBranch        *b_PuppiJet_HFEMEF;   //!
   TBranch        *b_PuppiJet_chHMult;   //!
   TBranch        *b_PuppiJet_neHMult;   //!
   TBranch        *b_PuppiJet_phMult;   //!
   TBranch        *b_PuppiJet_eMult;   //!
   TBranch        *b_PuppiJet_muMult;   //!
   TBranch        *b_PuppiJet_HFHMult;   //!
   TBranch        *b_PuppiJet_HFEMMult;   //!
   TBranch        *b_PuppiJet_puMva;   //!
   TBranch        *b_PuppiJet_CorrFactor_L1;   //!
   TBranch        *b_PuppiJet_CorrFactor_L1L2;   //!
   TBranch        *b_PuppiJet_CorrFactor_L1L2L3;   //!
   TBranch        *b_PuppiJet_CorrFactor_L1L2L3Res;   //!
   TBranch        *b_PuppiJet_mcMatchFlav;   //!
   TBranch        *b_PuppiJet_charge;   //!
   TBranch        *b_PuppiJet_ctagCsvL;   //!
   TBranch        *b_PuppiJet_ctagCsvB;   //!
   TBranch        *b_PuppiJet_area;   //!
   TBranch        *b_PuppiJet_qgl;   //!
   TBranch        *b_PuppiJet_ptd;   //!
   TBranch        *b_PuppiJet_axis2;   //!
   TBranch        *b_PuppiJet_mult;   //!
   TBranch        *b_PuppiJet_partonId;   //!
   TBranch        *b_PuppiJet_partonMotherId;   //!
   TBranch        *b_PuppiJet_nLeptons;   //!
   TBranch        *b_PuppiJet_id;   //!
   TBranch        *b_PuppiJet_puId;   //!
   TBranch        *b_PuppiJet_btagCSV;   //!
   TBranch        *b_PuppiJet_btagCMVA;   //!
   TBranch        *b_PuppiJet_rawPt;   //!
   TBranch        *b_PuppiJet_mcPt;   //!
   TBranch        *b_PuppiJet_mcFlavour;   //!
   TBranch        *b_PuppiJet_partonFlavour;   //!
   TBranch        *b_PuppiJet_hadronFlavour;   //!
   TBranch        *b_PuppiJet_mcMatchId;   //!
   TBranch        *b_PuppiJet_corr_JECUp;   //!
   TBranch        *b_PuppiJet_corr_JECDown;   //!
   TBranch        *b_PuppiJet_corr;   //!
   TBranch        *b_PuppiJet_corr_JERUp;   //!
   TBranch        *b_PuppiJet_corr_JERDown;   //!
   TBranch        *b_PuppiJet_corr_JER;   //!
   TBranch        *b_PuppiJet_pt;   //!
   TBranch        *b_PuppiJet_eta;   //!
   TBranch        *b_PuppiJet_phi;   //!
   TBranch        *b_PuppiJet_mass;   //!
   TBranch        *b_PuppiJet_massCorrected;   //!
   TBranch        *b_ngenFatJet;   //!
   TBranch        *b_genFatJet_pdgId;   //!
   TBranch        *b_genFatJet_pt;   //!
   TBranch        *b_genFatJet_eta;   //!
   TBranch        *b_genFatJet_phi;   //!
   TBranch        *b_genFatJet_mass;   //!
   TBranch        *b_genFatJet_charge;   //!
   TBranch        *b_genFatJet_status;   //!
   TBranch        *b_genFatJet_isPromptHard;   //!
   TBranch        *b_nJetFwd;   //!
   TBranch        *b_JetFwd_CorrFactor_L1;   //!
   TBranch        *b_JetFwd_CorrFactor_L1L2;   //!
   TBranch        *b_JetFwd_CorrFactor_L1L2L3;   //!
   TBranch        *b_JetFwd_CorrFactor_L1L2L3Res;   //!
   TBranch        *b_JetFwd_mcMatchFlav;   //!
   TBranch        *b_JetFwd_charge;   //!
   TBranch        *b_JetFwd_ctagCsvL;   //!
   TBranch        *b_JetFwd_ctagCsvB;   //!
   TBranch        *b_JetFwd_area;   //!
   TBranch        *b_JetFwd_qgl;   //!
   TBranch        *b_JetFwd_ptd;   //!
   TBranch        *b_JetFwd_axis2;   //!
   TBranch        *b_JetFwd_mult;   //!
   TBranch        *b_JetFwd_partonId;   //!
   TBranch        *b_JetFwd_partonMotherId;   //!
   TBranch        *b_JetFwd_nLeptons;   //!
   TBranch        *b_JetFwd_id;   //!
   TBranch        *b_JetFwd_puId;   //!
   TBranch        *b_JetFwd_btagCSV;   //!
   TBranch        *b_JetFwd_btagCMVA;   //!
   TBranch        *b_JetFwd_rawPt;   //!
   TBranch        *b_JetFwd_mcPt;   //!
   TBranch        *b_JetFwd_mcFlavour;   //!
   TBranch        *b_JetFwd_partonFlavour;   //!
   TBranch        *b_JetFwd_hadronFlavour;   //!
   TBranch        *b_JetFwd_mcMatchId;   //!
   TBranch        *b_JetFwd_corr_JECUp;   //!
   TBranch        *b_JetFwd_corr_JECDown;   //!
   TBranch        *b_JetFwd_corr;   //!
   TBranch        *b_JetFwd_corr_JERUp;   //!
   TBranch        *b_JetFwd_corr_JERDown;   //!
   TBranch        *b_JetFwd_corr_JER;   //!
   TBranch        *b_JetFwd_pt;   //!
   TBranch        *b_JetFwd_eta;   //!
   TBranch        *b_JetFwd_phi;   //!
   TBranch        *b_JetFwd_mass;   //!
   TBranch        *b_JetFwd_prunedMass;   //!
   TBranch        *b_JetFwd_mcNumPartons;   //!
   TBranch        *b_JetFwd_mcNumLeptons;   //!
   TBranch        *b_JetFwd_mcNumTaus;   //!
   TBranch        *b_JetFwd_mcAnyPartonMass;   //!
   TBranch        *b_JetFwd_nSubJets;   //!
   TBranch        *b_JetFwd_nSubJets25;   //!
   TBranch        *b_JetFwd_nSubJets30;   //!
   TBranch        *b_JetFwd_nSubJets40;   //!
   TBranch        *b_JetFwd_nSubJetsZ01;   //!
   TBranch        *b_JetFwd_chHEF;   //!
   TBranch        *b_JetFwd_neHEF;   //!
   TBranch        *b_JetFwd_phEF;   //!
   TBranch        *b_JetFwd_eEF;   //!
   TBranch        *b_JetFwd_muEF;   //!
   TBranch        *b_JetFwd_HFHEF;   //!
   TBranch        *b_JetFwd_HFEMEF;   //!
   TBranch        *b_JetFwd_chHMult;   //!
   TBranch        *b_JetFwd_neHMult;   //!
   TBranch        *b_JetFwd_phMult;   //!
   TBranch        *b_JetFwd_eMult;   //!
   TBranch        *b_JetFwd_muMult;   //!
   TBranch        *b_JetFwd_HFHMult;   //!
   TBranch        *b_JetFwd_HFEMMult;   //!
   TBranch        *b_JetFwd_puMva;   //!
   TBranch        *b_nGammaGood;   //!
   TBranch        *b_GammaGood_etaSc;   //!
   TBranch        *b_GammaGood_idCutBased;   //!
   TBranch        *b_GammaGood_hOverE;   //!
   TBranch        *b_GammaGood_r9;   //!
   TBranch        *b_GammaGood_sigmaIetaIeta;   //!
   TBranch        *b_GammaGood_chHadIso04;   //!
   TBranch        *b_GammaGood_chHadIso;   //!
   TBranch        *b_GammaGood_phIso;   //!
   TBranch        *b_GammaGood_neuHadIso;   //!
   TBranch        *b_GammaGood_relIso;   //!
   TBranch        *b_GammaGood_mcMatchId;   //!
   TBranch        *b_GammaGood_mcPt;   //!
   TBranch        *b_GammaGood_pdgId;   //!
   TBranch        *b_GammaGood_pt;   //!
   TBranch        *b_GammaGood_eta;   //!
   TBranch        *b_GammaGood_phi;   //!
   TBranch        *b_GammaGood_mass;   //!
   TBranch        *b_GammaGood_genIso04;   //!
   TBranch        *b_GammaGood_genIso03;   //!
   TBranch        *b_GammaGood_chHadIsoRC04;   //!
   TBranch        *b_GammaGood_chHadIsoRC;   //!
   TBranch        *b_GammaGood_drMinParton;   //!

   rootNtupleClass(TTree *tree=0);
   virtual ~rootNtupleClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef rootNtupleClass_cxx
rootNtupleClass::rootNtupleClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("tree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("tree","");
      chain->Add("ZGamma2050.root/tree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

rootNtupleClass::~rootNtupleClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t rootNtupleClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t rootNtupleClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void rootNtupleClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("xsec", &xsec, &b_xsec);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nVert", &nVert, &b_nVert);
   fChain->SetBranchAddress("nJet25", &nJet25, &b_nJet25);
   fChain->SetBranchAddress("nBJetLoose25", &nBJetLoose25, &b_nBJetLoose25);
   fChain->SetBranchAddress("nBJetMedium25", &nBJetMedium25, &b_nBJetMedium25);
   fChain->SetBranchAddress("nBJetTight25", &nBJetTight25, &b_nBJetTight25);
   fChain->SetBranchAddress("nJet30", &nJet30, &b_nJet30);
   fChain->SetBranchAddress("nJet30a", &nJet30a, &b_nJet30a);
   fChain->SetBranchAddress("nBJetLoose30", &nBJetLoose30, &b_nBJetLoose30);
   fChain->SetBranchAddress("nBJetMedium30", &nBJetMedium30, &b_nBJetMedium30);
   fChain->SetBranchAddress("nBJetTight30", &nBJetTight30, &b_nBJetTight30);
   fChain->SetBranchAddress("nJet40", &nJet40, &b_nJet40);
   fChain->SetBranchAddress("nJet40a", &nJet40a, &b_nJet40a);
   fChain->SetBranchAddress("nBJetLoose40", &nBJetLoose40, &b_nBJetLoose40);
   fChain->SetBranchAddress("nBJetMedium40", &nBJetMedium40, &b_nBJetMedium40);
   fChain->SetBranchAddress("nBJetTight40", &nBJetTight40, &b_nBJetTight40);
   fChain->SetBranchAddress("nLepGood20", &nLepGood20, &b_nLepGood20);
   fChain->SetBranchAddress("nLepGood15", &nLepGood15, &b_nLepGood15);
   fChain->SetBranchAddress("nLepGood10", &nLepGood10, &b_nLepGood10);
   fChain->SetBranchAddress("htJet25", &htJet25, &b_htJet25);
   fChain->SetBranchAddress("mhtJet25", &mhtJet25, &b_mhtJet25);
   fChain->SetBranchAddress("htJet40j", &htJet40j, &b_htJet40j);
   fChain->SetBranchAddress("htJet40ja", &htJet40ja, &b_htJet40ja);
   fChain->SetBranchAddress("htJet40", &htJet40, &b_htJet40);
   fChain->SetBranchAddress("htJet40a", &htJet40a, &b_htJet40a);
   fChain->SetBranchAddress("mhtJet40", &mhtJet40, &b_mhtJet40);
   fChain->SetBranchAddress("mhtJet40a", &mhtJet40a, &b_mhtJet40a);
   fChain->SetBranchAddress("apcjetmetmin", &apcjetmetmin, &b_apcjetmetmin);
   fChain->SetBranchAddress("mZ1", &mZ1, &b_mZ1);
   fChain->SetBranchAddress("mZ1SFSS", &mZ1SFSS, &b_mZ1SFSS);
   fChain->SetBranchAddress("minMllSFOS", &minMllSFOS, &b_minMllSFOS);
   fChain->SetBranchAddress("maxMllSFOS", &maxMllSFOS, &b_maxMllSFOS);
   fChain->SetBranchAddress("minMllAFOS", &minMllAFOS, &b_minMllAFOS);
   fChain->SetBranchAddress("maxMllAFOS", &maxMllAFOS, &b_maxMllAFOS);
   fChain->SetBranchAddress("minMllAFSS", &minMllAFSS, &b_minMllAFSS);
   fChain->SetBranchAddress("maxMllAFSS", &maxMllAFSS, &b_maxMllAFSS);
   fChain->SetBranchAddress("minMllAFAS", &minMllAFAS, &b_minMllAFAS);
   fChain->SetBranchAddress("maxMllAFAS", &maxMllAFAS, &b_maxMllAFAS);
   fChain->SetBranchAddress("m2l", &m2l, &b_m2l);
   fChain->SetBranchAddress("nBJet25", &nBJet25, &b_nBJet25);
   fChain->SetBranchAddress("nFatJet100", &nFatJet100, &b_nFatJet100);
   fChain->SetBranchAddress("nMuons10", &nMuons10, &b_nMuons10);
   fChain->SetBranchAddress("nElectrons10", &nElectrons10, &b_nElectrons10);
   fChain->SetBranchAddress("nTaus20", &nTaus20, &b_nTaus20);
   fChain->SetBranchAddress("nGammas20", &nGammas20, &b_nGammas20);
   fChain->SetBranchAddress("rhoCN", &rhoCN, &b_rhoCN);
   fChain->SetBranchAddress("hbheFilterNew50ns", &hbheFilterNew50ns, &b_hbheFilterNew50ns);
   fChain->SetBranchAddress("hbheFilterNew25ns", &hbheFilterNew25ns, &b_hbheFilterNew25ns);
   fChain->SetBranchAddress("hbheFilterIso", &hbheFilterIso, &b_hbheFilterIso);
   fChain->SetBranchAddress("met_trkPt", &met_trkPt, &b_met_trkPt);
   fChain->SetBranchAddress("met_trkPhi", &met_trkPhi, &b_met_trkPhi);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET170_NoiseCleaned_v", &HLT_BIT_HLT_PFMET170_NoiseCleaned_v, &b_HLT_BIT_HLT_PFMET170_NoiseCleaned_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET170_v", &HLT_BIT_HLT_PFMET170_v, &b_HLT_BIT_HLT_PFMET170_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET170_HBHECleaned_v", &HLT_BIT_HLT_PFMET170_HBHECleaned_v, &b_HLT_BIT_HLT_PFMET170_HBHECleaned_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET170_JetIdCleaned_v", &HLT_BIT_HLT_PFMET170_JetIdCleaned_v, &b_HLT_BIT_HLT_PFMET170_JetIdCleaned_v);
   fChain->SetBranchAddress("HLT_Met170", &HLT_Met170, &b_HLT_Met170);
   fChain->SetBranchAddress("HLT_BIT_HLT_Photon155_v", &HLT_BIT_HLT_Photon155_v, &b_HLT_BIT_HLT_Photon155_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Photon165_HE10_v", &HLT_BIT_HLT_Photon165_HE10_v, &b_HLT_BIT_HLT_Photon165_HE10_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Photon175_v", &HLT_BIT_HLT_Photon175_v, &b_HLT_BIT_HLT_Photon175_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet40_v", &HLT_BIT_HLT_PFJet40_v, &b_HLT_BIT_HLT_PFJet40_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet60_v", &HLT_BIT_HLT_PFJet60_v, &b_HLT_BIT_HLT_PFJet60_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet80_v", &HLT_BIT_HLT_PFJet80_v, &b_HLT_BIT_HLT_PFJet80_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet140_v", &HLT_BIT_HLT_PFJet140_v, &b_HLT_BIT_HLT_PFJet140_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet200_v", &HLT_BIT_HLT_PFJet200_v, &b_HLT_BIT_HLT_PFJet200_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet260_v", &HLT_BIT_HLT_PFJet260_v, &b_HLT_BIT_HLT_PFJet260_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet320_v", &HLT_BIT_HLT_PFJet320_v, &b_HLT_BIT_HLT_PFJet320_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet400_v", &HLT_BIT_HLT_PFJet400_v, &b_HLT_BIT_HLT_PFJet400_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet450_v", &HLT_BIT_HLT_PFJet450_v, &b_HLT_BIT_HLT_PFJet450_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet500_v", &HLT_BIT_HLT_PFJet500_v, &b_HLT_BIT_HLT_PFJet500_v);
   fChain->SetBranchAddress("HLT_SinglePho", &HLT_SinglePho, &b_HLT_SinglePho);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v", &HLT_BIT_HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v, &b_HLT_BIT_HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v", &HLT_BIT_HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v, &b_HLT_BIT_HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v", &HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v, &b_HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v);
   fChain->SetBranchAddress("HLT_MonoJetMetNoMuMHT120", &HLT_MonoJetMetNoMuMHT120, &b_HLT_MonoJetMetNoMuMHT120);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v", &HLT_BIT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v, &b_HLT_BIT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v", &HLT_BIT_HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v, &b_HLT_BIT_HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v", &HLT_BIT_HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v, &b_HLT_BIT_HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v);
   fChain->SetBranchAddress("HLT_MonoJetMetNoMuMHT90", &HLT_MonoJetMetNoMuMHT90, &b_HLT_MonoJetMetNoMuMHT90);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT800_v", &HLT_BIT_HLT_PFHT800_v, &b_HLT_BIT_HLT_PFHT800_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_ECALHT800_v", &HLT_BIT_HLT_ECALHT800_v, &b_HLT_BIT_HLT_ECALHT800_v);
   fChain->SetBranchAddress("HLT_Jet_HT", &HLT_Jet_HT, &b_HLT_Jet_HT);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET300_NoiseCleaned_v", &HLT_BIT_HLT_PFMET300_NoiseCleaned_v, &b_HLT_BIT_HLT_PFMET300_NoiseCleaned_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET300_v", &HLT_BIT_HLT_PFMET300_v, &b_HLT_BIT_HLT_PFMET300_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET300_JetIdCleaned_v", &HLT_BIT_HLT_PFMET300_JetIdCleaned_v, &b_HLT_BIT_HLT_PFMET300_JetIdCleaned_v);
   fChain->SetBranchAddress("HLT_Met300", &HLT_Met300, &b_HLT_Met300);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("HLT_DoubleEl", &HLT_DoubleEl, &b_HLT_DoubleEl);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", &HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v);
   fChain->SetBranchAddress("HLT_DoubleMu", &HLT_DoubleMu, &b_HLT_DoubleMu);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu24_eta2p1_v", &HLT_BIT_HLT_IsoMu24_eta2p1_v, &b_HLT_BIT_HLT_IsoMu24_eta2p1_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoTkMu24_eta2p1_v", &HLT_BIT_HLT_IsoTkMu24_eta2p1_v, &b_HLT_BIT_HLT_IsoTkMu24_eta2p1_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu18_v", &HLT_BIT_HLT_IsoMu18_v, &b_HLT_BIT_HLT_IsoMu18_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu20_v", &HLT_BIT_HLT_IsoMu20_v, &b_HLT_BIT_HLT_IsoMu20_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoTkMu20_v", &HLT_BIT_HLT_IsoTkMu20_v, &b_HLT_BIT_HLT_IsoTkMu20_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu27_v", &HLT_BIT_HLT_IsoMu27_v, &b_HLT_BIT_HLT_IsoMu27_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoTkMu27_v", &HLT_BIT_HLT_IsoTkMu27_v, &b_HLT_BIT_HLT_IsoTkMu27_v);
   fChain->SetBranchAddress("HLT_SingleMu", &HLT_SingleMu, &b_HLT_SingleMu);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele23_WPLoose_Gsf_v", &HLT_BIT_HLT_Ele23_WPLoose_Gsf_v, &b_HLT_BIT_HLT_Ele23_WPLoose_Gsf_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v", &HLT_BIT_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_BIT_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele27_WPLoose_Gsf_v", &HLT_BIT_HLT_Ele27_WPLoose_Gsf_v, &b_HLT_BIT_HLT_Ele27_WPLoose_Gsf_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele27_eta2p1_WPLoose_Gsf_v", &HLT_BIT_HLT_Ele27_eta2p1_WPLoose_Gsf_v, &b_HLT_BIT_HLT_Ele27_eta2p1_WPLoose_Gsf_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele32_eta2p1_WPLoose_Gsf_v", &HLT_BIT_HLT_Ele32_eta2p1_WPLoose_Gsf_v, &b_HLT_BIT_HLT_Ele32_eta2p1_WPLoose_Gsf_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele27_WP85_Gsf_v", &HLT_BIT_HLT_Ele27_WP85_Gsf_v, &b_HLT_BIT_HLT_Ele27_WP85_Gsf_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_v", &HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_v, &b_HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_v", &HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_v, &b_HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_v);
   fChain->SetBranchAddress("HLT_SingleEl", &HLT_SingleEl, &b_HLT_SingleEl);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, &b_Flag_trackingFailureFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_mass", &met_mass, &b_met_mass);
   fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   fChain->SetBranchAddress("met_rawPt", &met_rawPt, &b_met_rawPt);
   fChain->SetBranchAddress("met_rawPhi", &met_rawPhi, &b_met_rawPhi);
   fChain->SetBranchAddress("met_rawSumEt", &met_rawSumEt, &b_met_rawSumEt);
   fChain->SetBranchAddress("met_genPt", &met_genPt, &b_met_genPt);
   fChain->SetBranchAddress("met_genPhi", &met_genPhi, &b_met_genPhi);
   fChain->SetBranchAddress("met_genEta", &met_genEta, &b_met_genEta);
   fChain->SetBranchAddress("metNoMu_pt", &metNoMu_pt, &b_metNoMu_pt);
   fChain->SetBranchAddress("metNoMu_eta", &metNoMu_eta, &b_metNoMu_eta);
   fChain->SetBranchAddress("metNoMu_phi", &metNoMu_phi, &b_metNoMu_phi);
   fChain->SetBranchAddress("metNoMu_mass", &metNoMu_mass, &b_metNoMu_mass);
   fChain->SetBranchAddress("metNoMu_sumEt", &metNoMu_sumEt, &b_metNoMu_sumEt);
   fChain->SetBranchAddress("metNoMu_rawPt", &metNoMu_rawPt, &b_metNoMu_rawPt);
   fChain->SetBranchAddress("metNoMu_rawPhi", &metNoMu_rawPhi, &b_metNoMu_rawPhi);
   fChain->SetBranchAddress("metNoMu_rawSumEt", &metNoMu_rawSumEt, &b_metNoMu_rawSumEt);
   fChain->SetBranchAddress("metNoMu_genPt", &metNoMu_genPt, &b_metNoMu_genPt);
   fChain->SetBranchAddress("metNoMu_genPhi", &metNoMu_genPhi, &b_metNoMu_genPhi);
   fChain->SetBranchAddress("metNoMu_genEta", &metNoMu_genEta, &b_metNoMu_genEta);
   fChain->SetBranchAddress("metPuppi_pt", &metPuppi_pt, &b_metPuppi_pt);
   fChain->SetBranchAddress("metPuppi_eta", &metPuppi_eta, &b_metPuppi_eta);
   fChain->SetBranchAddress("metPuppi_phi", &metPuppi_phi, &b_metPuppi_phi);
   fChain->SetBranchAddress("metPuppi_mass", &metPuppi_mass, &b_metPuppi_mass);
   fChain->SetBranchAddress("metPuppi_sumEt", &metPuppi_sumEt, &b_metPuppi_sumEt);
   fChain->SetBranchAddress("metPuppi_rawPt", &metPuppi_rawPt, &b_metPuppi_rawPt);
   fChain->SetBranchAddress("metPuppi_rawPhi", &metPuppi_rawPhi, &b_metPuppi_rawPhi);
   fChain->SetBranchAddress("metPuppi_rawSumEt", &metPuppi_rawSumEt, &b_metPuppi_rawSumEt);
   fChain->SetBranchAddress("metPuppi_genPt", &metPuppi_genPt, &b_metPuppi_genPt);
   fChain->SetBranchAddress("metPuppi_genPhi", &metPuppi_genPhi, &b_metPuppi_genPhi);
   fChain->SetBranchAddress("metPuppi_genEta", &metPuppi_genEta, &b_metPuppi_genEta);
   fChain->SetBranchAddress("metNoPU_pt", &metNoPU_pt, &b_metNoPU_pt);
   fChain->SetBranchAddress("metNoPU_eta", &metNoPU_eta, &b_metNoPU_eta);
   fChain->SetBranchAddress("metNoPU_phi", &metNoPU_phi, &b_metNoPU_phi);
   fChain->SetBranchAddress("metNoPU_mass", &metNoPU_mass, &b_metNoPU_mass);
   fChain->SetBranchAddress("nFatJet", &nFatJet, &b_nFatJet);
   fChain->SetBranchAddress("FatJet_id", FatJet_id, &b_FatJet_id);
   fChain->SetBranchAddress("FatJet_puId", FatJet_puId, &b_FatJet_puId);
   fChain->SetBranchAddress("FatJet_btagCSV", FatJet_btagCSV, &b_FatJet_btagCSV);
   fChain->SetBranchAddress("FatJet_btagCMVA", FatJet_btagCMVA, &b_FatJet_btagCMVA);
   fChain->SetBranchAddress("FatJet_rawPt", FatJet_rawPt, &b_FatJet_rawPt);
   fChain->SetBranchAddress("FatJet_mcPt", FatJet_mcPt, &b_FatJet_mcPt);
   fChain->SetBranchAddress("FatJet_mcFlavour", FatJet_mcFlavour, &b_FatJet_mcFlavour);
   fChain->SetBranchAddress("FatJet_partonFlavour", FatJet_partonFlavour, &b_FatJet_partonFlavour);
   fChain->SetBranchAddress("FatJet_hadronFlavour", FatJet_hadronFlavour, &b_FatJet_hadronFlavour);
   fChain->SetBranchAddress("FatJet_mcMatchId", FatJet_mcMatchId, &b_FatJet_mcMatchId);
   fChain->SetBranchAddress("FatJet_corr_JECUp", FatJet_corr_JECUp, &b_FatJet_corr_JECUp);
   fChain->SetBranchAddress("FatJet_corr_JECDown", FatJet_corr_JECDown, &b_FatJet_corr_JECDown);
   fChain->SetBranchAddress("FatJet_corr", FatJet_corr, &b_FatJet_corr);
   fChain->SetBranchAddress("FatJet_corr_JERUp", FatJet_corr_JERUp, &b_FatJet_corr_JERUp);
   fChain->SetBranchAddress("FatJet_corr_JERDown", FatJet_corr_JERDown, &b_FatJet_corr_JERDown);
   fChain->SetBranchAddress("FatJet_corr_JER", FatJet_corr_JER, &b_FatJet_corr_JER);
   fChain->SetBranchAddress("FatJet_pt", FatJet_pt, &b_FatJet_pt);
   fChain->SetBranchAddress("FatJet_eta", FatJet_eta, &b_FatJet_eta);
   fChain->SetBranchAddress("FatJet_phi", FatJet_phi, &b_FatJet_phi);
   fChain->SetBranchAddress("FatJet_mass", FatJet_mass, &b_FatJet_mass);
   fChain->SetBranchAddress("FatJet_softDropMass", FatJet_softDropMass, &b_FatJet_softDropMass);
   fChain->SetBranchAddress("FatJet_tau1", FatJet_tau1, &b_FatJet_tau1);
   fChain->SetBranchAddress("FatJet_tau2", FatJet_tau2, &b_FatJet_tau2);
   fChain->SetBranchAddress("FatJet_tau3", FatJet_tau3, &b_FatJet_tau3);
   fChain->SetBranchAddress("FatJet_topMass", FatJet_topMass, &b_FatJet_topMass);
   fChain->SetBranchAddress("FatJet_minMass", FatJet_minMass, &b_FatJet_minMass);
   fChain->SetBranchAddress("FatJet_nSubJets", FatJet_nSubJets, &b_FatJet_nSubJets);
   fChain->SetBranchAddress("FatJet_puMva", FatJet_puMva, &b_FatJet_puMva);
   fChain->SetBranchAddress("FatJet_prunedMass", FatJet_prunedMass, &b_FatJet_prunedMass);
   fChain->SetBranchAddress("FatJet_chHEF", FatJet_chHEF, &b_FatJet_chHEF);
   fChain->SetBranchAddress("FatJet_neHEF", FatJet_neHEF, &b_FatJet_neHEF);
   fChain->SetBranchAddress("FatJet_chEmEF", FatJet_chEmEF, &b_FatJet_chEmEF);
   fChain->SetBranchAddress("FatJet_neEmEF", FatJet_neEmEF, &b_FatJet_neEmEF);
   fChain->SetBranchAddress("FatJet_chMult", FatJet_chMult, &b_FatJet_chMult);
   fChain->SetBranchAddress("FatJet_neMult", FatJet_neMult, &b_FatJet_neMult);
   fChain->SetBranchAddress("FatJet_puppiPt", FatJet_puppiPt, &b_FatJet_puppiPt);
   fChain->SetBranchAddress("FatJet_puppiEta", FatJet_puppiEta, &b_FatJet_puppiEta);
   fChain->SetBranchAddress("FatJet_puppiPhi", FatJet_puppiPhi, &b_FatJet_puppiPhi);
   fChain->SetBranchAddress("FatJet_puppiMass", FatJet_puppiMass, &b_FatJet_puppiMass);
   fChain->SetBranchAddress("FatJet_puppiTau1", FatJet_puppiTau1, &b_FatJet_puppiTau1);
   fChain->SetBranchAddress("FatJet_puppiTau2", FatJet_puppiTau2, &b_FatJet_puppiTau2);
   fChain->SetBranchAddress("FatJet_puppiTau3", FatJet_puppiTau3, &b_FatJet_puppiTau3);
   fChain->SetBranchAddress("nisoTrack", &nisoTrack, &b_nisoTrack);
   fChain->SetBranchAddress("isoTrack_pdgId", isoTrack_pdgId, &b_isoTrack_pdgId);
   fChain->SetBranchAddress("isoTrack_pt", isoTrack_pt, &b_isoTrack_pt);
   fChain->SetBranchAddress("isoTrack_eta", isoTrack_eta, &b_isoTrack_eta);
   fChain->SetBranchAddress("isoTrack_phi", isoTrack_phi, &b_isoTrack_phi);
   fChain->SetBranchAddress("isoTrack_mass", isoTrack_mass, &b_isoTrack_mass);
   fChain->SetBranchAddress("isoTrack_charge", isoTrack_charge, &b_isoTrack_charge);
   fChain->SetBranchAddress("isoTrack_dz", isoTrack_dz, &b_isoTrack_dz);
   fChain->SetBranchAddress("isoTrack_absIso", isoTrack_absIso, &b_isoTrack_absIso);
   fChain->SetBranchAddress("isoTrack_relIsoAn04", isoTrack_relIsoAn04, &b_isoTrack_relIsoAn04);
   fChain->SetBranchAddress("isoTrack_mcMatchId", isoTrack_mcMatchId, &b_isoTrack_mcMatchId);
   fChain->SetBranchAddress("ncustomPuppiAK8", &ncustomPuppiAK8, &b_ncustomPuppiAK8);
   fChain->SetBranchAddress("customPuppiAK8_pt", customPuppiAK8_pt, &b_customPuppiAK8_pt);
   fChain->SetBranchAddress("customPuppiAK8_eta", customPuppiAK8_eta, &b_customPuppiAK8_eta);
   fChain->SetBranchAddress("customPuppiAK8_phi", customPuppiAK8_phi, &b_customPuppiAK8_phi);
   fChain->SetBranchAddress("customPuppiAK8_mass", customPuppiAK8_mass, &b_customPuppiAK8_mass);
   fChain->SetBranchAddress("customPuppiAK8_massCorrected", customPuppiAK8_massCorrected, &b_customPuppiAK8_massCorrected);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_motherId", GenPart_motherId, &b_GenPart_motherId);
   fChain->SetBranchAddress("GenPart_grandmotherId", GenPart_grandmotherId, &b_GenPart_grandmotherId);
   fChain->SetBranchAddress("GenPart_sourceId", GenPart_sourceId, &b_GenPart_sourceId);
   fChain->SetBranchAddress("GenPart_charge", GenPart_charge, &b_GenPart_charge);
   fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_isPromptHard", GenPart_isPromptHard, &b_GenPart_isPromptHard);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_motherIndex", GenPart_motherIndex, &b_GenPart_motherIndex);
   fChain->SetBranchAddress("nLepGood", &nLepGood, &b_nLepGood);
   fChain->SetBranchAddress("LepGood_charge", LepGood_charge, &b_LepGood_charge);
   fChain->SetBranchAddress("LepGood_tightId", LepGood_tightId, &b_LepGood_tightId);
   fChain->SetBranchAddress("LepGood_eleCutIdCSA14_25ns_v1", LepGood_eleCutIdCSA14_25ns_v1, &b_LepGood_eleCutIdCSA14_25ns_v1);
   fChain->SetBranchAddress("LepGood_eleCutIdCSA14_50ns_v1", LepGood_eleCutIdCSA14_50ns_v1, &b_LepGood_eleCutIdCSA14_50ns_v1);
   fChain->SetBranchAddress("LepGood_eleCutIdSpring15_25ns_v1", LepGood_eleCutIdSpring15_25ns_v1, &b_LepGood_eleCutIdSpring15_25ns_v1);
   fChain->SetBranchAddress("LepGood_dxy", LepGood_dxy, &b_LepGood_dxy);
   fChain->SetBranchAddress("LepGood_dz", LepGood_dz, &b_LepGood_dz);
   fChain->SetBranchAddress("LepGood_edxy", LepGood_edxy, &b_LepGood_edxy);
   fChain->SetBranchAddress("LepGood_edz", LepGood_edz, &b_LepGood_edz);
   fChain->SetBranchAddress("LepGood_ip3d", LepGood_ip3d, &b_LepGood_ip3d);
   fChain->SetBranchAddress("LepGood_sip3d", LepGood_sip3d, &b_LepGood_sip3d);
   fChain->SetBranchAddress("LepGood_convVeto", LepGood_convVeto, &b_LepGood_convVeto);
   fChain->SetBranchAddress("LepGood_lostHits", LepGood_lostHits, &b_LepGood_lostHits);
   fChain->SetBranchAddress("LepGood_relIso03", LepGood_relIso03, &b_LepGood_relIso03);
   fChain->SetBranchAddress("LepGood_relIso04", LepGood_relIso04, &b_LepGood_relIso04);
   fChain->SetBranchAddress("LepGood_miniRelIso", LepGood_miniRelIso, &b_LepGood_miniRelIso);
   fChain->SetBranchAddress("LepGood_relIsoAn04", LepGood_relIsoAn04, &b_LepGood_relIsoAn04);
   fChain->SetBranchAddress("LepGood_tightCharge", LepGood_tightCharge, &b_LepGood_tightCharge);
   fChain->SetBranchAddress("LepGood_mcMatchId", LepGood_mcMatchId, &b_LepGood_mcMatchId);
   fChain->SetBranchAddress("LepGood_mcMatchAny", LepGood_mcMatchAny, &b_LepGood_mcMatchAny);
   fChain->SetBranchAddress("LepGood_mcMatchTau", LepGood_mcMatchTau, &b_LepGood_mcMatchTau);
   fChain->SetBranchAddress("LepGood_mcPt", LepGood_mcPt, &b_LepGood_mcPt);
   fChain->SetBranchAddress("LepGood_mediumMuonId", LepGood_mediumMuonId, &b_LepGood_mediumMuonId);
   fChain->SetBranchAddress("LepGood_ICHEPmediumMuonId", LepGood_ICHEPmediumMuonId, &b_LepGood_ICHEPmediumMuonId);
   fChain->SetBranchAddress("LepGood_pdgId", LepGood_pdgId, &b_LepGood_pdgId);
   fChain->SetBranchAddress("LepGood_pt", LepGood_pt, &b_LepGood_pt);
   fChain->SetBranchAddress("LepGood_eta", LepGood_eta, &b_LepGood_eta);
   fChain->SetBranchAddress("LepGood_phi", LepGood_phi, &b_LepGood_phi);
   fChain->SetBranchAddress("LepGood_mass", LepGood_mass, &b_LepGood_mass);
   fChain->SetBranchAddress("LepGood_looseIdOnly", LepGood_looseIdOnly, &b_LepGood_looseIdOnly);
   fChain->SetBranchAddress("LepGood_chargedHadRelIso03", LepGood_chargedHadRelIso03, &b_LepGood_chargedHadRelIso03);
   fChain->SetBranchAddress("LepGood_chargedHadRelIso04", LepGood_chargedHadRelIso04, &b_LepGood_chargedHadRelIso04);
   fChain->SetBranchAddress("LepGood_softMuonId", LepGood_softMuonId, &b_LepGood_softMuonId);
   fChain->SetBranchAddress("LepGood_pfMuonId", LepGood_pfMuonId, &b_LepGood_pfMuonId);
   fChain->SetBranchAddress("LepGood_eleCutId2012_full5x5", LepGood_eleCutId2012_full5x5, &b_LepGood_eleCutId2012_full5x5);
   fChain->SetBranchAddress("LepGood_trackerLayers", LepGood_trackerLayers, &b_LepGood_trackerLayers);
   fChain->SetBranchAddress("LepGood_pixelLayers", LepGood_pixelLayers, &b_LepGood_pixelLayers);
   fChain->SetBranchAddress("LepGood_trackerHits", LepGood_trackerHits, &b_LepGood_trackerHits);
   fChain->SetBranchAddress("LepGood_lostOuterHits", LepGood_lostOuterHits, &b_LepGood_lostOuterHits);
   fChain->SetBranchAddress("LepGood_innerTrackValidHitFraction", LepGood_innerTrackValidHitFraction, &b_LepGood_innerTrackValidHitFraction);
   fChain->SetBranchAddress("LepGood_innerTrackChi2", LepGood_innerTrackChi2, &b_LepGood_innerTrackChi2);
   fChain->SetBranchAddress("LepGood_nStations", LepGood_nStations, &b_LepGood_nStations);
   fChain->SetBranchAddress("LepGood_caloCompatibility", LepGood_caloCompatibility, &b_LepGood_caloCompatibility);
   fChain->SetBranchAddress("LepGood_globalTrackChi2", LepGood_globalTrackChi2, &b_LepGood_globalTrackChi2);
   fChain->SetBranchAddress("LepGood_trkKink", LepGood_trkKink, &b_LepGood_trkKink);
   fChain->SetBranchAddress("LepGood_segmentCompatibility", LepGood_segmentCompatibility, &b_LepGood_segmentCompatibility);
   fChain->SetBranchAddress("LepGood_chi2LocalPosition", LepGood_chi2LocalPosition, &b_LepGood_chi2LocalPosition);
   fChain->SetBranchAddress("LepGood_chi2LocalMomentum", LepGood_chi2LocalMomentum, &b_LepGood_chi2LocalMomentum);
   fChain->SetBranchAddress("LepGood_glbTrackProbability", LepGood_glbTrackProbability, &b_LepGood_glbTrackProbability);
   fChain->SetBranchAddress("LepGood_TMOneStationTightMuonId", LepGood_TMOneStationTightMuonId, &b_LepGood_TMOneStationTightMuonId);
   fChain->SetBranchAddress("LepGood_trackHighPurityMuon", LepGood_trackHighPurityMuon, &b_LepGood_trackHighPurityMuon);
   fChain->SetBranchAddress("LepGood_isGlobalMuon", LepGood_isGlobalMuon, &b_LepGood_isGlobalMuon);
   fChain->SetBranchAddress("LepGood_isTrackerMuon", LepGood_isTrackerMuon, &b_LepGood_isTrackerMuon);
   fChain->SetBranchAddress("LepGood_sigmaIEtaIEta", LepGood_sigmaIEtaIEta, &b_LepGood_sigmaIEtaIEta);
   fChain->SetBranchAddress("LepGood_dEtaScTrkIn", LepGood_dEtaScTrkIn, &b_LepGood_dEtaScTrkIn);
   fChain->SetBranchAddress("LepGood_dPhiScTrkIn", LepGood_dPhiScTrkIn, &b_LepGood_dPhiScTrkIn);
   fChain->SetBranchAddress("LepGood_hadronicOverEm", LepGood_hadronicOverEm, &b_LepGood_hadronicOverEm);
   fChain->SetBranchAddress("LepGood_eInvMinusPInv", LepGood_eInvMinusPInv, &b_LepGood_eInvMinusPInv);
   fChain->SetBranchAddress("LepGood_eInvMinusPInv_tkMom", LepGood_eInvMinusPInv_tkMom, &b_LepGood_eInvMinusPInv_tkMom);
   fChain->SetBranchAddress("LepGood_etaSc", LepGood_etaSc, &b_LepGood_etaSc);
   fChain->SetBranchAddress("LepGood_mcMatchPdgId", LepGood_mcMatchPdgId, &b_LepGood_mcMatchPdgId);
   fChain->SetBranchAddress("LepGood_r9", LepGood_r9, &b_LepGood_r9);
   fChain->SetBranchAddress("LepGood_e5x5", LepGood_e5x5, &b_LepGood_e5x5);
   fChain->SetBranchAddress("LepGood_sigmaIetaIeta", LepGood_sigmaIetaIeta, &b_LepGood_sigmaIetaIeta);
   fChain->SetBranchAddress("LepGood_sigmaIphiIphi", LepGood_sigmaIphiIphi, &b_LepGood_sigmaIphiIphi);
   fChain->SetBranchAddress("LepGood_hcalOverEcal", LepGood_hcalOverEcal, &b_LepGood_hcalOverEcal);
   fChain->SetBranchAddress("LepGood_full5x5_e5x5", LepGood_full5x5_e5x5, &b_LepGood_full5x5_e5x5);
   fChain->SetBranchAddress("LepGood_full5x5_r9", LepGood_full5x5_r9, &b_LepGood_full5x5_r9);
   fChain->SetBranchAddress("LepGood_full5x5_sigmaIetaIeta", LepGood_full5x5_sigmaIetaIeta, &b_LepGood_full5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("LepGood_full5x5_sigmaIphiIphi", LepGood_full5x5_sigmaIphiIphi, &b_LepGood_full5x5_sigmaIphiIphi);
   fChain->SetBranchAddress("LepGood_full5x5_hcalOverEcal", LepGood_full5x5_hcalOverEcal, &b_LepGood_full5x5_hcalOverEcal);
   fChain->SetBranchAddress("LepGood_correctedEcalEnergy", LepGood_correctedEcalEnergy, &b_LepGood_correctedEcalEnergy);
   fChain->SetBranchAddress("LepGood_eSuperClusterOverP", LepGood_eSuperClusterOverP, &b_LepGood_eSuperClusterOverP);
   fChain->SetBranchAddress("LepGood_ecalEnergy", LepGood_ecalEnergy, &b_LepGood_ecalEnergy);
   fChain->SetBranchAddress("LepGood_superCluster_rawEnergy", LepGood_superCluster_rawEnergy, &b_LepGood_superCluster_rawEnergy);
   fChain->SetBranchAddress("LepGood_superCluster_preshowerEnergy", LepGood_superCluster_preshowerEnergy, &b_LepGood_superCluster_preshowerEnergy);
   fChain->SetBranchAddress("LepGood_superCluster_correctedEnergy", LepGood_superCluster_correctedEnergy, &b_LepGood_superCluster_correctedEnergy);
   fChain->SetBranchAddress("LepGood_superCluster_energy", LepGood_superCluster_energy, &b_LepGood_superCluster_energy);
   fChain->SetBranchAddress("LepGood_superCluster_clustersSize", LepGood_superCluster_clustersSize, &b_LepGood_superCluster_clustersSize);
   fChain->SetBranchAddress("LepGood_superCluster_seed.energy", LepGood_superCluster_seed_energy, &b_LepGood_superCluster_seed_energy);
   fChain->SetBranchAddress("LepGood_e2x5Max", LepGood_e2x5Max, &b_LepGood_e2x5Max);
   fChain->SetBranchAddress("LepGood_e1x5", LepGood_e1x5, &b_LepGood_e1x5);
   fChain->SetBranchAddress("LepGood_isolTrkPt", LepGood_isolTrkPt, &b_LepGood_isolTrkPt);
   fChain->SetBranchAddress("LepGood_isolEmHadDepth1", LepGood_isolEmHadDepth1, &b_LepGood_isolEmHadDepth1);
   fChain->SetBranchAddress("LepGood_eleClusterDEta", LepGood_eleClusterDEta, &b_LepGood_eleClusterDEta);
   fChain->SetBranchAddress("LepGood_isEcalDriven", LepGood_isEcalDriven, &b_LepGood_isEcalDriven);
   fChain->SetBranchAddress("LepGood_muonDB", LepGood_muonDB, &b_LepGood_muonDB);
   fChain->SetBranchAddress("LepGood_pixelHits", LepGood_pixelHits, &b_LepGood_pixelHits);
   fChain->SetBranchAddress("LepGood_muTrackIso", LepGood_muTrackIso, &b_LepGood_muTrackIso);
   fChain->SetBranchAddress("LepGood_muon_dz", LepGood_muon_dz, &b_LepGood_muon_dz);
   fChain->SetBranchAddress("LepGood_nChamberHits", LepGood_nChamberHits, &b_LepGood_nChamberHits);
   fChain->SetBranchAddress("LepGood_muonPtRatio", LepGood_muonPtRatio, &b_LepGood_muonPtRatio);
   fChain->SetBranchAddress("nTauGood", &nTauGood, &b_nTauGood);
   fChain->SetBranchAddress("TauGood_charge", TauGood_charge, &b_TauGood_charge);
   fChain->SetBranchAddress("TauGood_decayMode", TauGood_decayMode, &b_TauGood_decayMode);
   fChain->SetBranchAddress("TauGood_idDecayMode", TauGood_idDecayMode, &b_TauGood_idDecayMode);
   fChain->SetBranchAddress("TauGood_idDecayModeNewDMs", TauGood_idDecayModeNewDMs, &b_TauGood_idDecayModeNewDMs);
   fChain->SetBranchAddress("TauGood_dxy", TauGood_dxy, &b_TauGood_dxy);
   fChain->SetBranchAddress("TauGood_dz", TauGood_dz, &b_TauGood_dz);
   fChain->SetBranchAddress("TauGood_idMVA", TauGood_idMVA, &b_TauGood_idMVA);
   fChain->SetBranchAddress("TauGood_idMVANewDM", TauGood_idMVANewDM, &b_TauGood_idMVANewDM);
   fChain->SetBranchAddress("TauGood_idCI3hit", TauGood_idCI3hit, &b_TauGood_idCI3hit);
   fChain->SetBranchAddress("TauGood_idAntiMu", TauGood_idAntiMu, &b_TauGood_idAntiMu);
   fChain->SetBranchAddress("TauGood_idAntiE", TauGood_idAntiE, &b_TauGood_idAntiE);
   fChain->SetBranchAddress("TauGood_isoCI3hit", TauGood_isoCI3hit, &b_TauGood_isoCI3hit);
   fChain->SetBranchAddress("TauGood_mcMatchId", TauGood_mcMatchId, &b_TauGood_mcMatchId);
   fChain->SetBranchAddress("TauGood_pdgId", TauGood_pdgId, &b_TauGood_pdgId);
   fChain->SetBranchAddress("TauGood_pt", TauGood_pt, &b_TauGood_pt);
   fChain->SetBranchAddress("TauGood_eta", TauGood_eta, &b_TauGood_eta);
   fChain->SetBranchAddress("TauGood_phi", TauGood_phi, &b_TauGood_phi);
   fChain->SetBranchAddress("TauGood_mass", TauGood_mass, &b_TauGood_mass);
   fChain->SetBranchAddress("TauGood_idMVAOldDMRun2", TauGood_idMVAOldDMRun2, &b_TauGood_idMVAOldDMRun2);
   fChain->SetBranchAddress("TauGood_idMVAOldDMRun2dR03", TauGood_idMVAOldDMRun2dR03, &b_TauGood_idMVAOldDMRun2dR03);
   fChain->SetBranchAddress("ngenJet", &ngenJet, &b_ngenJet);
   fChain->SetBranchAddress("genJet_pdgId", genJet_pdgId, &b_genJet_pdgId);
   fChain->SetBranchAddress("genJet_pt", genJet_pt, &b_genJet_pt);
   fChain->SetBranchAddress("genJet_eta", genJet_eta, &b_genJet_eta);
   fChain->SetBranchAddress("genJet_phi", genJet_phi, &b_genJet_phi);
   fChain->SetBranchAddress("genJet_mass", genJet_mass, &b_genJet_mass);
   fChain->SetBranchAddress("genJet_charge", genJet_charge, &b_genJet_charge);
   fChain->SetBranchAddress("genJet_status", genJet_status, &b_genJet_status);
   fChain->SetBranchAddress("genJet_isPromptHard", genJet_isPromptHard, &b_genJet_isPromptHard);
   fChain->SetBranchAddress("nSoftDropSubJet", &nSoftDropSubJet, &b_nSoftDropSubJet);
   fChain->SetBranchAddress("SoftDropSubJet_pt", SoftDropSubJet_pt, &b_SoftDropSubJet_pt);
   fChain->SetBranchAddress("SoftDropSubJet_eta", SoftDropSubJet_eta, &b_SoftDropSubJet_eta);
   fChain->SetBranchAddress("SoftDropSubJet_phi", SoftDropSubJet_phi, &b_SoftDropSubJet_phi);
   fChain->SetBranchAddress("SoftDropSubJet_mass", SoftDropSubJet_mass, &b_SoftDropSubJet_mass);
   fChain->SetBranchAddress("nPuppiSubJet", &nPuppiSubJet, &b_nPuppiSubJet);
   fChain->SetBranchAddress("PuppiSubJet_pt", PuppiSubJet_pt, &b_PuppiSubJet_pt);
   fChain->SetBranchAddress("PuppiSubJet_eta", PuppiSubJet_eta, &b_PuppiSubJet_eta);
   fChain->SetBranchAddress("PuppiSubJet_phi", PuppiSubJet_phi, &b_PuppiSubJet_phi);
   fChain->SetBranchAddress("PuppiSubJet_mass", PuppiSubJet_mass, &b_PuppiSubJet_mass);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_CorrFactor_L1", Jet_CorrFactor_L1, &b_Jet_CorrFactor_L1);
   fChain->SetBranchAddress("Jet_CorrFactor_L1L2", Jet_CorrFactor_L1L2, &b_Jet_CorrFactor_L1L2);
   fChain->SetBranchAddress("Jet_CorrFactor_L1L2L3", Jet_CorrFactor_L1L2L3, &b_Jet_CorrFactor_L1L2L3);
   fChain->SetBranchAddress("Jet_CorrFactor_L1L2L3Res", Jet_CorrFactor_L1L2L3Res, &b_Jet_CorrFactor_L1L2L3Res);
   fChain->SetBranchAddress("Jet_mcMatchFlav", Jet_mcMatchFlav, &b_Jet_mcMatchFlav);
   fChain->SetBranchAddress("Jet_charge", Jet_charge, &b_Jet_charge);
   fChain->SetBranchAddress("Jet_ctagCsvL", Jet_ctagCsvL, &b_Jet_ctagCsvL);
   fChain->SetBranchAddress("Jet_ctagCsvB", Jet_ctagCsvB, &b_Jet_ctagCsvB);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   fChain->SetBranchAddress("Jet_ptd", Jet_ptd, &b_Jet_ptd);
   fChain->SetBranchAddress("Jet_axis2", Jet_axis2, &b_Jet_axis2);
   fChain->SetBranchAddress("Jet_mult", Jet_mult, &b_Jet_mult);
   fChain->SetBranchAddress("Jet_partonId", Jet_partonId, &b_Jet_partonId);
   fChain->SetBranchAddress("Jet_partonMotherId", Jet_partonMotherId, &b_Jet_partonMotherId);
   fChain->SetBranchAddress("Jet_nLeptons", Jet_nLeptons, &b_Jet_nLeptons);
   fChain->SetBranchAddress("Jet_id", Jet_id, &b_Jet_id);
   fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   fChain->SetBranchAddress("Jet_btagCSV", Jet_btagCSV, &b_Jet_btagCSV);
   fChain->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA, &b_Jet_btagCMVA);
   fChain->SetBranchAddress("Jet_rawPt", Jet_rawPt, &b_Jet_rawPt);
   fChain->SetBranchAddress("Jet_mcPt", Jet_mcPt, &b_Jet_mcPt);
   fChain->SetBranchAddress("Jet_mcFlavour", Jet_mcFlavour, &b_Jet_mcFlavour);
   fChain->SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
   fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   fChain->SetBranchAddress("Jet_mcMatchId", Jet_mcMatchId, &b_Jet_mcMatchId);
   fChain->SetBranchAddress("Jet_corr_JECUp", Jet_corr_JECUp, &b_Jet_corr_JECUp);
   fChain->SetBranchAddress("Jet_corr_JECDown", Jet_corr_JECDown, &b_Jet_corr_JECDown);
   fChain->SetBranchAddress("Jet_corr", Jet_corr, &b_Jet_corr);
   fChain->SetBranchAddress("Jet_corr_JERUp", Jet_corr_JERUp, &b_Jet_corr_JERUp);
   fChain->SetBranchAddress("Jet_corr_JERDown", Jet_corr_JERDown, &b_Jet_corr_JERDown);
   fChain->SetBranchAddress("Jet_corr_JER", Jet_corr_JER, &b_Jet_corr_JER);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_prunedMass", Jet_prunedMass, &b_Jet_prunedMass);
   fChain->SetBranchAddress("Jet_mcNumPartons", Jet_mcNumPartons, &b_Jet_mcNumPartons);
   fChain->SetBranchAddress("Jet_mcNumLeptons", Jet_mcNumLeptons, &b_Jet_mcNumLeptons);
   fChain->SetBranchAddress("Jet_mcNumTaus", Jet_mcNumTaus, &b_Jet_mcNumTaus);
   fChain->SetBranchAddress("Jet_mcAnyPartonMass", Jet_mcAnyPartonMass, &b_Jet_mcAnyPartonMass);
   fChain->SetBranchAddress("Jet_nSubJets", Jet_nSubJets, &b_Jet_nSubJets);
   fChain->SetBranchAddress("Jet_nSubJets25", Jet_nSubJets25, &b_Jet_nSubJets25);
   fChain->SetBranchAddress("Jet_nSubJets30", Jet_nSubJets30, &b_Jet_nSubJets30);
   fChain->SetBranchAddress("Jet_nSubJets40", Jet_nSubJets40, &b_Jet_nSubJets40);
   fChain->SetBranchAddress("Jet_nSubJetsZ01", Jet_nSubJetsZ01, &b_Jet_nSubJetsZ01);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_phEF", Jet_phEF, &b_Jet_phEF);
   fChain->SetBranchAddress("Jet_eEF", Jet_eEF, &b_Jet_eEF);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_HFHEF", Jet_HFHEF, &b_Jet_HFHEF);
   fChain->SetBranchAddress("Jet_HFEMEF", Jet_HFEMEF, &b_Jet_HFEMEF);
   fChain->SetBranchAddress("Jet_chHMult", Jet_chHMult, &b_Jet_chHMult);
   fChain->SetBranchAddress("Jet_neHMult", Jet_neHMult, &b_Jet_neHMult);
   fChain->SetBranchAddress("Jet_phMult", Jet_phMult, &b_Jet_phMult);
   fChain->SetBranchAddress("Jet_eMult", Jet_eMult, &b_Jet_eMult);
   fChain->SetBranchAddress("Jet_muMult", Jet_muMult, &b_Jet_muMult);
   fChain->SetBranchAddress("Jet_HFHMult", Jet_HFHMult, &b_Jet_HFHMult);
   fChain->SetBranchAddress("Jet_HFEMMult", Jet_HFEMMult, &b_Jet_HFEMMult);
   fChain->SetBranchAddress("Jet_puMva", Jet_puMva, &b_Jet_puMva);
   fChain->SetBranchAddress("ngenLep", &ngenLep, &b_ngenLep);
   fChain->SetBranchAddress("genLep_motherId", genLep_motherId, &b_genLep_motherId);
   fChain->SetBranchAddress("genLep_grandmotherId", genLep_grandmotherId, &b_genLep_grandmotherId);
   fChain->SetBranchAddress("genLep_sourceId", genLep_sourceId, &b_genLep_sourceId);
   fChain->SetBranchAddress("genLep_charge", genLep_charge, &b_genLep_charge);
   fChain->SetBranchAddress("genLep_status", genLep_status, &b_genLep_status);
   fChain->SetBranchAddress("genLep_isPromptHard", genLep_isPromptHard, &b_genLep_isPromptHard);
   fChain->SetBranchAddress("genLep_pdgId", genLep_pdgId, &b_genLep_pdgId);
   fChain->SetBranchAddress("genLep_pt", genLep_pt, &b_genLep_pt);
   fChain->SetBranchAddress("genLep_eta", genLep_eta, &b_genLep_eta);
   fChain->SetBranchAddress("genLep_phi", genLep_phi, &b_genLep_phi);
   fChain->SetBranchAddress("genLep_mass", genLep_mass, &b_genLep_mass);
   fChain->SetBranchAddress("genLep_motherIndex", genLep_motherIndex, &b_genLep_motherIndex);
   fChain->SetBranchAddress("nLHEweight", &nLHEweight, &b_nLHEweight);
   fChain->SetBranchAddress("LHEweight_id", &LHEweight_id, &b_LHEweight_id);
   fChain->SetBranchAddress("LHEweight_wgt", &LHEweight_wgt, &b_LHEweight_wgt);
   fChain->SetBranchAddress("ncustomPuppiSoftDropAK8", &ncustomPuppiSoftDropAK8, &b_ncustomPuppiSoftDropAK8);
   fChain->SetBranchAddress("customPuppiSoftDropAK8_pt", customPuppiSoftDropAK8_pt, &b_customPuppiSoftDropAK8_pt);
   fChain->SetBranchAddress("customPuppiSoftDropAK8_eta", customPuppiSoftDropAK8_eta, &b_customPuppiSoftDropAK8_eta);
   fChain->SetBranchAddress("customPuppiSoftDropAK8_phi", customPuppiSoftDropAK8_phi, &b_customPuppiSoftDropAK8_phi);
   fChain->SetBranchAddress("customPuppiSoftDropAK8_mass", customPuppiSoftDropAK8_mass, &b_customPuppiSoftDropAK8_mass);
   fChain->SetBranchAddress("customPuppiSoftDropAK8_massCorrected", customPuppiSoftDropAK8_massCorrected, &b_customPuppiSoftDropAK8_massCorrected);
   fChain->SetBranchAddress("nPuppiJet", &nPuppiJet, &b_nPuppiJet);
   fChain->SetBranchAddress("PuppiJet_prunedMass", PuppiJet_prunedMass, &b_PuppiJet_prunedMass);
   fChain->SetBranchAddress("PuppiJet_mcNumPartons", PuppiJet_mcNumPartons, &b_PuppiJet_mcNumPartons);
   fChain->SetBranchAddress("PuppiJet_mcNumLeptons", PuppiJet_mcNumLeptons, &b_PuppiJet_mcNumLeptons);
   fChain->SetBranchAddress("PuppiJet_mcNumTaus", PuppiJet_mcNumTaus, &b_PuppiJet_mcNumTaus);
   fChain->SetBranchAddress("PuppiJet_mcAnyPartonMass", PuppiJet_mcAnyPartonMass, &b_PuppiJet_mcAnyPartonMass);
   fChain->SetBranchAddress("PuppiJet_nSubJets", PuppiJet_nSubJets, &b_PuppiJet_nSubJets);
   fChain->SetBranchAddress("PuppiJet_nSubJets25", PuppiJet_nSubJets25, &b_PuppiJet_nSubJets25);
   fChain->SetBranchAddress("PuppiJet_nSubJets30", PuppiJet_nSubJets30, &b_PuppiJet_nSubJets30);
   fChain->SetBranchAddress("PuppiJet_nSubJets40", PuppiJet_nSubJets40, &b_PuppiJet_nSubJets40);
   fChain->SetBranchAddress("PuppiJet_nSubJetsZ01", PuppiJet_nSubJetsZ01, &b_PuppiJet_nSubJetsZ01);
   fChain->SetBranchAddress("PuppiJet_chHEF", PuppiJet_chHEF, &b_PuppiJet_chHEF);
   fChain->SetBranchAddress("PuppiJet_neHEF", PuppiJet_neHEF, &b_PuppiJet_neHEF);
   fChain->SetBranchAddress("PuppiJet_phEF", PuppiJet_phEF, &b_PuppiJet_phEF);
   fChain->SetBranchAddress("PuppiJet_eEF", PuppiJet_eEF, &b_PuppiJet_eEF);
   fChain->SetBranchAddress("PuppiJet_muEF", PuppiJet_muEF, &b_PuppiJet_muEF);
   fChain->SetBranchAddress("PuppiJet_HFHEF", PuppiJet_HFHEF, &b_PuppiJet_HFHEF);
   fChain->SetBranchAddress("PuppiJet_HFEMEF", PuppiJet_HFEMEF, &b_PuppiJet_HFEMEF);
   fChain->SetBranchAddress("PuppiJet_chHMult", PuppiJet_chHMult, &b_PuppiJet_chHMult);
   fChain->SetBranchAddress("PuppiJet_neHMult", PuppiJet_neHMult, &b_PuppiJet_neHMult);
   fChain->SetBranchAddress("PuppiJet_phMult", PuppiJet_phMult, &b_PuppiJet_phMult);
   fChain->SetBranchAddress("PuppiJet_eMult", PuppiJet_eMult, &b_PuppiJet_eMult);
   fChain->SetBranchAddress("PuppiJet_muMult", PuppiJet_muMult, &b_PuppiJet_muMult);
   fChain->SetBranchAddress("PuppiJet_HFHMult", PuppiJet_HFHMult, &b_PuppiJet_HFHMult);
   fChain->SetBranchAddress("PuppiJet_HFEMMult", PuppiJet_HFEMMult, &b_PuppiJet_HFEMMult);
   fChain->SetBranchAddress("PuppiJet_puMva", PuppiJet_puMva, &b_PuppiJet_puMva);
   fChain->SetBranchAddress("PuppiJet_CorrFactor_L1", PuppiJet_CorrFactor_L1, &b_PuppiJet_CorrFactor_L1);
   fChain->SetBranchAddress("PuppiJet_CorrFactor_L1L2", PuppiJet_CorrFactor_L1L2, &b_PuppiJet_CorrFactor_L1L2);
   fChain->SetBranchAddress("PuppiJet_CorrFactor_L1L2L3", PuppiJet_CorrFactor_L1L2L3, &b_PuppiJet_CorrFactor_L1L2L3);
   fChain->SetBranchAddress("PuppiJet_CorrFactor_L1L2L3Res", PuppiJet_CorrFactor_L1L2L3Res, &b_PuppiJet_CorrFactor_L1L2L3Res);
   fChain->SetBranchAddress("PuppiJet_mcMatchFlav", PuppiJet_mcMatchFlav, &b_PuppiJet_mcMatchFlav);
   fChain->SetBranchAddress("PuppiJet_charge", PuppiJet_charge, &b_PuppiJet_charge);
   fChain->SetBranchAddress("PuppiJet_ctagCsvL", PuppiJet_ctagCsvL, &b_PuppiJet_ctagCsvL);
   fChain->SetBranchAddress("PuppiJet_ctagCsvB", PuppiJet_ctagCsvB, &b_PuppiJet_ctagCsvB);
   fChain->SetBranchAddress("PuppiJet_area", PuppiJet_area, &b_PuppiJet_area);
   fChain->SetBranchAddress("PuppiJet_qgl", PuppiJet_qgl, &b_PuppiJet_qgl);
   fChain->SetBranchAddress("PuppiJet_ptd", PuppiJet_ptd, &b_PuppiJet_ptd);
   fChain->SetBranchAddress("PuppiJet_axis2", PuppiJet_axis2, &b_PuppiJet_axis2);
   fChain->SetBranchAddress("PuppiJet_mult", PuppiJet_mult, &b_PuppiJet_mult);
   fChain->SetBranchAddress("PuppiJet_partonId", PuppiJet_partonId, &b_PuppiJet_partonId);
   fChain->SetBranchAddress("PuppiJet_partonMotherId", PuppiJet_partonMotherId, &b_PuppiJet_partonMotherId);
   fChain->SetBranchAddress("PuppiJet_nLeptons", PuppiJet_nLeptons, &b_PuppiJet_nLeptons);
   fChain->SetBranchAddress("PuppiJet_id", PuppiJet_id, &b_PuppiJet_id);
   fChain->SetBranchAddress("PuppiJet_puId", PuppiJet_puId, &b_PuppiJet_puId);
   fChain->SetBranchAddress("PuppiJet_btagCSV", PuppiJet_btagCSV, &b_PuppiJet_btagCSV);
   fChain->SetBranchAddress("PuppiJet_btagCMVA", PuppiJet_btagCMVA, &b_PuppiJet_btagCMVA);
   fChain->SetBranchAddress("PuppiJet_rawPt", PuppiJet_rawPt, &b_PuppiJet_rawPt);
   fChain->SetBranchAddress("PuppiJet_mcPt", PuppiJet_mcPt, &b_PuppiJet_mcPt);
   fChain->SetBranchAddress("PuppiJet_mcFlavour", PuppiJet_mcFlavour, &b_PuppiJet_mcFlavour);
   fChain->SetBranchAddress("PuppiJet_partonFlavour", PuppiJet_partonFlavour, &b_PuppiJet_partonFlavour);
   fChain->SetBranchAddress("PuppiJet_hadronFlavour", PuppiJet_hadronFlavour, &b_PuppiJet_hadronFlavour);
   fChain->SetBranchAddress("PuppiJet_mcMatchId", PuppiJet_mcMatchId, &b_PuppiJet_mcMatchId);
   fChain->SetBranchAddress("PuppiJet_corr_JECUp", PuppiJet_corr_JECUp, &b_PuppiJet_corr_JECUp);
   fChain->SetBranchAddress("PuppiJet_corr_JECDown", PuppiJet_corr_JECDown, &b_PuppiJet_corr_JECDown);
   fChain->SetBranchAddress("PuppiJet_corr", PuppiJet_corr, &b_PuppiJet_corr);
   fChain->SetBranchAddress("PuppiJet_corr_JERUp", PuppiJet_corr_JERUp, &b_PuppiJet_corr_JERUp);
   fChain->SetBranchAddress("PuppiJet_corr_JERDown", PuppiJet_corr_JERDown, &b_PuppiJet_corr_JERDown);
   fChain->SetBranchAddress("PuppiJet_corr_JER", PuppiJet_corr_JER, &b_PuppiJet_corr_JER);
   fChain->SetBranchAddress("PuppiJet_pt", PuppiJet_pt, &b_PuppiJet_pt);
   fChain->SetBranchAddress("PuppiJet_eta", PuppiJet_eta, &b_PuppiJet_eta);
   fChain->SetBranchAddress("PuppiJet_phi", PuppiJet_phi, &b_PuppiJet_phi);
   fChain->SetBranchAddress("PuppiJet_mass", PuppiJet_mass, &b_PuppiJet_mass);
   fChain->SetBranchAddress("PuppiJet_massCorrected", PuppiJet_massCorrected, &b_PuppiJet_massCorrected);
   fChain->SetBranchAddress("ngenFatJet", &ngenFatJet, &b_ngenFatJet);
   fChain->SetBranchAddress("genFatJet_pdgId", genFatJet_pdgId, &b_genFatJet_pdgId);
   fChain->SetBranchAddress("genFatJet_pt", genFatJet_pt, &b_genFatJet_pt);
   fChain->SetBranchAddress("genFatJet_eta", genFatJet_eta, &b_genFatJet_eta);
   fChain->SetBranchAddress("genFatJet_phi", genFatJet_phi, &b_genFatJet_phi);
   fChain->SetBranchAddress("genFatJet_mass", genFatJet_mass, &b_genFatJet_mass);
   fChain->SetBranchAddress("genFatJet_charge", genFatJet_charge, &b_genFatJet_charge);
   fChain->SetBranchAddress("genFatJet_status", genFatJet_status, &b_genFatJet_status);
   fChain->SetBranchAddress("genFatJet_isPromptHard", genFatJet_isPromptHard, &b_genFatJet_isPromptHard);
   fChain->SetBranchAddress("nJetFwd", &nJetFwd, &b_nJetFwd);
   fChain->SetBranchAddress("JetFwd_CorrFactor_L1", JetFwd_CorrFactor_L1, &b_JetFwd_CorrFactor_L1);
   fChain->SetBranchAddress("JetFwd_CorrFactor_L1L2", JetFwd_CorrFactor_L1L2, &b_JetFwd_CorrFactor_L1L2);
   fChain->SetBranchAddress("JetFwd_CorrFactor_L1L2L3", JetFwd_CorrFactor_L1L2L3, &b_JetFwd_CorrFactor_L1L2L3);
   fChain->SetBranchAddress("JetFwd_CorrFactor_L1L2L3Res", JetFwd_CorrFactor_L1L2L3Res, &b_JetFwd_CorrFactor_L1L2L3Res);
   fChain->SetBranchAddress("JetFwd_mcMatchFlav", JetFwd_mcMatchFlav, &b_JetFwd_mcMatchFlav);
   fChain->SetBranchAddress("JetFwd_charge", JetFwd_charge, &b_JetFwd_charge);
   fChain->SetBranchAddress("JetFwd_ctagCsvL", JetFwd_ctagCsvL, &b_JetFwd_ctagCsvL);
   fChain->SetBranchAddress("JetFwd_ctagCsvB", JetFwd_ctagCsvB, &b_JetFwd_ctagCsvB);
   fChain->SetBranchAddress("JetFwd_area", JetFwd_area, &b_JetFwd_area);
   fChain->SetBranchAddress("JetFwd_qgl", JetFwd_qgl, &b_JetFwd_qgl);
   fChain->SetBranchAddress("JetFwd_ptd", JetFwd_ptd, &b_JetFwd_ptd);
   fChain->SetBranchAddress("JetFwd_axis2", JetFwd_axis2, &b_JetFwd_axis2);
   fChain->SetBranchAddress("JetFwd_mult", JetFwd_mult, &b_JetFwd_mult);
   fChain->SetBranchAddress("JetFwd_partonId", JetFwd_partonId, &b_JetFwd_partonId);
   fChain->SetBranchAddress("JetFwd_partonMotherId", JetFwd_partonMotherId, &b_JetFwd_partonMotherId);
   fChain->SetBranchAddress("JetFwd_nLeptons", JetFwd_nLeptons, &b_JetFwd_nLeptons);
   fChain->SetBranchAddress("JetFwd_id", JetFwd_id, &b_JetFwd_id);
   fChain->SetBranchAddress("JetFwd_puId", JetFwd_puId, &b_JetFwd_puId);
   fChain->SetBranchAddress("JetFwd_btagCSV", JetFwd_btagCSV, &b_JetFwd_btagCSV);
   fChain->SetBranchAddress("JetFwd_btagCMVA", JetFwd_btagCMVA, &b_JetFwd_btagCMVA);
   fChain->SetBranchAddress("JetFwd_rawPt", JetFwd_rawPt, &b_JetFwd_rawPt);
   fChain->SetBranchAddress("JetFwd_mcPt", JetFwd_mcPt, &b_JetFwd_mcPt);
   fChain->SetBranchAddress("JetFwd_mcFlavour", JetFwd_mcFlavour, &b_JetFwd_mcFlavour);
   fChain->SetBranchAddress("JetFwd_partonFlavour", JetFwd_partonFlavour, &b_JetFwd_partonFlavour);
   fChain->SetBranchAddress("JetFwd_hadronFlavour", JetFwd_hadronFlavour, &b_JetFwd_hadronFlavour);
   fChain->SetBranchAddress("JetFwd_mcMatchId", JetFwd_mcMatchId, &b_JetFwd_mcMatchId);
   fChain->SetBranchAddress("JetFwd_corr_JECUp", JetFwd_corr_JECUp, &b_JetFwd_corr_JECUp);
   fChain->SetBranchAddress("JetFwd_corr_JECDown", JetFwd_corr_JECDown, &b_JetFwd_corr_JECDown);
   fChain->SetBranchAddress("JetFwd_corr", JetFwd_corr, &b_JetFwd_corr);
   fChain->SetBranchAddress("JetFwd_corr_JERUp", JetFwd_corr_JERUp, &b_JetFwd_corr_JERUp);
   fChain->SetBranchAddress("JetFwd_corr_JERDown", JetFwd_corr_JERDown, &b_JetFwd_corr_JERDown);
   fChain->SetBranchAddress("JetFwd_corr_JER", JetFwd_corr_JER, &b_JetFwd_corr_JER);
   fChain->SetBranchAddress("JetFwd_pt", JetFwd_pt, &b_JetFwd_pt);
   fChain->SetBranchAddress("JetFwd_eta", JetFwd_eta, &b_JetFwd_eta);
   fChain->SetBranchAddress("JetFwd_phi", JetFwd_phi, &b_JetFwd_phi);
   fChain->SetBranchAddress("JetFwd_mass", JetFwd_mass, &b_JetFwd_mass);
   fChain->SetBranchAddress("JetFwd_prunedMass", JetFwd_prunedMass, &b_JetFwd_prunedMass);
   fChain->SetBranchAddress("JetFwd_mcNumPartons", JetFwd_mcNumPartons, &b_JetFwd_mcNumPartons);
   fChain->SetBranchAddress("JetFwd_mcNumLeptons", JetFwd_mcNumLeptons, &b_JetFwd_mcNumLeptons);
   fChain->SetBranchAddress("JetFwd_mcNumTaus", JetFwd_mcNumTaus, &b_JetFwd_mcNumTaus);
   fChain->SetBranchAddress("JetFwd_mcAnyPartonMass", JetFwd_mcAnyPartonMass, &b_JetFwd_mcAnyPartonMass);
   fChain->SetBranchAddress("JetFwd_nSubJets", JetFwd_nSubJets, &b_JetFwd_nSubJets);
   fChain->SetBranchAddress("JetFwd_nSubJets25", JetFwd_nSubJets25, &b_JetFwd_nSubJets25);
   fChain->SetBranchAddress("JetFwd_nSubJets30", JetFwd_nSubJets30, &b_JetFwd_nSubJets30);
   fChain->SetBranchAddress("JetFwd_nSubJets40", JetFwd_nSubJets40, &b_JetFwd_nSubJets40);
   fChain->SetBranchAddress("JetFwd_nSubJetsZ01", JetFwd_nSubJetsZ01, &b_JetFwd_nSubJetsZ01);
   fChain->SetBranchAddress("JetFwd_chHEF", JetFwd_chHEF, &b_JetFwd_chHEF);
   fChain->SetBranchAddress("JetFwd_neHEF", JetFwd_neHEF, &b_JetFwd_neHEF);
   fChain->SetBranchAddress("JetFwd_phEF", JetFwd_phEF, &b_JetFwd_phEF);
   fChain->SetBranchAddress("JetFwd_eEF", JetFwd_eEF, &b_JetFwd_eEF);
   fChain->SetBranchAddress("JetFwd_muEF", JetFwd_muEF, &b_JetFwd_muEF);
   fChain->SetBranchAddress("JetFwd_HFHEF", JetFwd_HFHEF, &b_JetFwd_HFHEF);
   fChain->SetBranchAddress("JetFwd_HFEMEF", JetFwd_HFEMEF, &b_JetFwd_HFEMEF);
   fChain->SetBranchAddress("JetFwd_chHMult", JetFwd_chHMult, &b_JetFwd_chHMult);
   fChain->SetBranchAddress("JetFwd_neHMult", JetFwd_neHMult, &b_JetFwd_neHMult);
   fChain->SetBranchAddress("JetFwd_phMult", JetFwd_phMult, &b_JetFwd_phMult);
   fChain->SetBranchAddress("JetFwd_eMult", JetFwd_eMult, &b_JetFwd_eMult);
   fChain->SetBranchAddress("JetFwd_muMult", JetFwd_muMult, &b_JetFwd_muMult);
   fChain->SetBranchAddress("JetFwd_HFHMult", JetFwd_HFHMult, &b_JetFwd_HFHMult);
   fChain->SetBranchAddress("JetFwd_HFEMMult", JetFwd_HFEMMult, &b_JetFwd_HFEMMult);
   fChain->SetBranchAddress("JetFwd_puMva", JetFwd_puMva, &b_JetFwd_puMva);
   fChain->SetBranchAddress("nGammaGood", &nGammaGood, &b_nGammaGood);
   fChain->SetBranchAddress("GammaGood_etaSc", GammaGood_etaSc, &b_GammaGood_etaSc);
   fChain->SetBranchAddress("GammaGood_idCutBased", GammaGood_idCutBased, &b_GammaGood_idCutBased);
   fChain->SetBranchAddress("GammaGood_hOverE", GammaGood_hOverE, &b_GammaGood_hOverE);
   fChain->SetBranchAddress("GammaGood_r9", GammaGood_r9, &b_GammaGood_r9);
   fChain->SetBranchAddress("GammaGood_sigmaIetaIeta", GammaGood_sigmaIetaIeta, &b_GammaGood_sigmaIetaIeta);
   fChain->SetBranchAddress("GammaGood_chHadIso04", GammaGood_chHadIso04, &b_GammaGood_chHadIso04);
   fChain->SetBranchAddress("GammaGood_chHadIso", GammaGood_chHadIso, &b_GammaGood_chHadIso);
   fChain->SetBranchAddress("GammaGood_phIso", GammaGood_phIso, &b_GammaGood_phIso);
   fChain->SetBranchAddress("GammaGood_neuHadIso", GammaGood_neuHadIso, &b_GammaGood_neuHadIso);
   fChain->SetBranchAddress("GammaGood_relIso", GammaGood_relIso, &b_GammaGood_relIso);
   fChain->SetBranchAddress("GammaGood_mcMatchId", GammaGood_mcMatchId, &b_GammaGood_mcMatchId);
   fChain->SetBranchAddress("GammaGood_mcPt", GammaGood_mcPt, &b_GammaGood_mcPt);
   fChain->SetBranchAddress("GammaGood_pdgId", GammaGood_pdgId, &b_GammaGood_pdgId);
   fChain->SetBranchAddress("GammaGood_pt", GammaGood_pt, &b_GammaGood_pt);
   fChain->SetBranchAddress("GammaGood_eta", GammaGood_eta, &b_GammaGood_eta);
   fChain->SetBranchAddress("GammaGood_phi", GammaGood_phi, &b_GammaGood_phi);
   fChain->SetBranchAddress("GammaGood_mass", GammaGood_mass, &b_GammaGood_mass);
   fChain->SetBranchAddress("GammaGood_genIso04", GammaGood_genIso04, &b_GammaGood_genIso04);
   fChain->SetBranchAddress("GammaGood_genIso03", GammaGood_genIso03, &b_GammaGood_genIso03);
   fChain->SetBranchAddress("GammaGood_chHadIsoRC04", GammaGood_chHadIsoRC04, &b_GammaGood_chHadIsoRC04);
   fChain->SetBranchAddress("GammaGood_chHadIsoRC", GammaGood_chHadIsoRC, &b_GammaGood_chHadIsoRC);
   fChain->SetBranchAddress("GammaGood_drMinParton", GammaGood_drMinParton, &b_GammaGood_drMinParton);
   Notify();
}

Bool_t rootNtupleClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void rootNtupleClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t rootNtupleClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef rootNtupleClass_cxx
