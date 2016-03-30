//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 21 09:11:44 2016 by ROOT version 6.02/05
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
   Int_t           nTrueInt;
   Float_t         genWeight;
   Float_t         json;
   Float_t         nPU0;
   Float_t         nPVs;
   Float_t         rho;
   Float_t         lheNj;
   Float_t         lheNb;
   Float_t         lheNc;
   Float_t         lheNg;
   Float_t         lheNl;
   Float_t         lheV_pt;
   Float_t         lheHT;
   Float_t         mhtJet30;
   Float_t         mhtPhiJet30;
   Float_t         htJet30;
   Float_t         met_rawpt;
   Float_t         metPuppi_pt;
   Float_t         metPuppi_phi;
   Float_t         metPuppi_rawpt;
   Float_t         metType1p2_pt;
   Float_t         metNoPU_pt;
   Float_t         metNoPU_phi;
   Float_t         tkMet_pt;
   Float_t         tkMet_phi;
   Float_t         tkMetPVchs_pt;
   Float_t         tkMetPVchs_phi;
   Float_t         Flag_hbheIsoFilter;
   Float_t         Flag_hbheFilterNew;
   Float_t         simPrimaryVertex_z;
   Float_t         genHiggsDecayMode;
   Float_t         bTagWeight_LFUp;
   Float_t         bTagWeight_Stats2Down;
   Float_t         bTagWeight_LFDown;
   Float_t         bTagWeight_HFUp;
   Float_t         bTagWeight_cErr1Down;
   Float_t         bTagWeight_JESDown;
   Float_t         bTagWeight_cErr1Up;
   Float_t         bTagWeight;
   Float_t         bTagWeight_HFDown;
   Float_t         bTagWeight_Stats2Up;
   Float_t         bTagWeight_cErr2Up;
   Float_t         bTagWeight_JESUp;
   Float_t         bTagWeight_Stats1Up;
   Float_t         bTagWeight_Stats1Down;
   Float_t         bTagWeight_cErr2Down;
   Float_t         Flag_EcalDeadCellTriggerPrimitiveFilter;
   Float_t         Flag_trkPOG_manystripclus53X;
   Float_t         Flag_ecalLaserCorrFilter;
   Float_t         Flag_trkPOG_toomanystripclus53X;
   Float_t         Flag_hcalLaserEventFilter;
   Float_t         Flag_trkPOG_logErrorTooManyClusters;
   Float_t         Flag_trkPOGFilters;
   Float_t         Flag_trackingFailureFilter;
   Float_t         Flag_CSCTightHaloFilter;
   Float_t         Flag_HBHENoiseFilter;
   Float_t         Flag_goodVertices;
   Float_t         Flag_METFilters;
   Float_t         Flag_eeBadScFilter;
   Float_t         HLT_BIT_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_BTagCSV0p7_v;
   Float_t         HLT_BIT_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_v;
   Float_t         HLT_BIT_HLT_PFMET90_PFMHT90_IDLoose_v;
   Float_t         HLT_BIT_HLT_PFMET100_PFMHT100_IDLoose_v;
   Float_t         HLT_BIT_HLT_PFMET110_PFMHT110_IDLoose_v;
   Float_t         HLT_BIT_HLT_PFMET120_PFMHT120_IDLoose_v;
   Float_t         HLT_BIT_HLT_PFMET170_NoiseCleaned_v;
   Float_t         HLT_BIT_HLT_DiCentralPFJet70_PFMET120_NoiseCleaned_v;
   Float_t         HLT_BIT_HLT_PFHT350_PFMET120_NoiseCleaned_v;
   Float_t         HLT_ZnnHbbAll;
   Float_t         HLT_BIT_HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v;
   Float_t         HLT_BIT_HLT_PFHT450_SixJet40_PFBTagCSV_v;
   Float_t         HLT_ttHhardonicLowLumi;
   Float_t         HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v;
   Float_t         HLT_WtaunHbbHighLumi;
   Float_t         HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v;
   Float_t         HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v;
   Float_t         HLT_VBFHbbLowLumi;
   Float_t         HLT_BIT_HLT_PFHT750_4Jet_v;
   Float_t         HLT_BIT_HLT_PFHT900_v;
   Float_t         HLT_BIT_HLT_PFJet40_v;
   Float_t         HLT_BIT_HLT_PFJet60_v;
   Float_t         HLT_BIT_HLT_PFJet80_v;
   Float_t         HLT_BIT_HLT_PFJet140_v;
   Float_t         HLT_BIT_HLT_PFJet200_v;
   Float_t         HLT_BIT_HLT_PFJet260_v;
   Float_t         HLT_BIT_HLT_PFJet320_v;
   Float_t         HLT_BIT_HLT_PFJet400_v;
   Float_t         HLT_BIT_HLT_PFJet450_v;
   Float_t         HLT_hadronic;
   Float_t         HLT_BIT_HLT_QuadJet45_TripleCSV0p5_v;
   Float_t         HLT_BIT_HLT_QuadJet45_DoubleCSV0p5_v;
   Float_t         HLT_BIT_HLT_DoubleJet90_Double30_TripleCSV0p5_v;
   Float_t         HLT_BIT_HLT_DoubleJet90_Double30_DoubleCSV0p5_v;
   Float_t         HLT_HH4bAll;
   Float_t         HLT_BIT_HLT_IsoMu24_eta2p1_v;
   Float_t         HLT_BIT_HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v;
   Float_t         HLT_BIT_HLT_Mu24_eta2p1_v;
   Float_t         HLT_BIT_HLT_TkMu24_eta2p1_v;
   Float_t         HLT_BIT_HLT_Mu24_v;
   Float_t         HLT_BIT_HLT_IsoMu27_v;
   Float_t         HLT_BIT_HLT_IsoTkMu27_v;
   Float_t         HLT_BIT_HLT_TkMu27_v;
   Float_t         HLT_BIT_HLT_Mu27_v;
   Float_t         HLT_BIT_HLT_IsoMu20_eta2p1_v;
   Float_t         HLT_BIT_HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v;
   Float_t         HLT_BIT_HLT_Mu20_v;
   Float_t         HLT_BIT_HLT_TkMu20_v;
   Float_t         HLT_BIT_HLT_IsoMu20_v;
   Float_t         HLT_BIT_HLT_IsoTkMu20_v;
   Float_t         HLT_BIT_HLT_Mu40_eta2p1_PFJet200_PFJet50_v;
   Float_t         HLT_BIT_HLT_IsoMu16_eta2p1_CaloMET30_v;
   Float_t         HLT_BIT_HLT_Mu16_eta2p1_CaloMET30_v;
   Float_t         HLT_BIT_HLT_PFMET120_NoiseCleaned_Mu5_v;
   Float_t         HLT_BIT_HLT_MonoCentralPFJet80_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v;
   Float_t         HLT_BIT_HLT_MonoCentralPFJet80_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v;
   Float_t         HLT_BIT_HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v;
   Float_t         HLT_BIT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v;
   Float_t         HLT_WmnHbbAll;
   Float_t         HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_v;
   Float_t         HLT_BIT_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v;
   Float_t         HLT_WenHbbHighLumi;
   Float_t         HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
   Float_t         HLT_BIT_HLT_DoubleEle24_22_eta2p1_WP75_Gsf_v;
   Float_t         HLT_ZeeHbbLowLumi;
   Float_t         HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v;
   Float_t         HLT_BIT_HLT_Ele27_WP85_Gsf_v;
   Float_t         HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_v;
   Float_t         HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v;
   Float_t         HLT_BIT_HLT_Ele105_CaloIdVT_GsfTrkIdT_v;
   Float_t         HLT_BIT_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v;
   Float_t         HLT_WenHbbAll;
   Float_t         HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v;
   Float_t         HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_v;
   Float_t         HLT_WtaunHbbAll;
   Float_t         HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
   Float_t         HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   Float_t         HLT_BIT_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v;
   Float_t         HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   Float_t         HLT_BIT_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   Float_t         HLT_ZeeHbbAll;
   Float_t         HLT_WtaunHbbLowLumi;
   Float_t         HLT_WmnHbbLowLumi;
   Float_t         HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   Float_t         HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v;
   Float_t         HLT_ttHleptonic;
   Float_t         HLT_ZeeHbbHighLumi;
   Float_t         HLT_BIT_HLT_PFHT450_SixJet40_v;
   Float_t         HLT_BIT_HLT_PFHT400_SixJet30_v;
   Float_t         HLT_BIT_HLT_PFHT350_v;
   Float_t         HLT_ttHhardonicAll;
   Float_t         HLT_BIT_HLT_Mu50_v;
   Float_t         HLT_HLT_Mu50;
   Float_t         HLT_BIT_HLT_Mu45_eta2p1_v;
   Float_t         HLT_HLT_Mu45_eta2p1;
   Float_t         HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v;
   Float_t         HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v;
   Float_t         HLT_ZmmHbbLowLumi;
   Float_t         HLT_WmnHbbHighLumi;
   Float_t         HLT_BIT_HLT_Mu17_TkMu8_DZ_v;
   Float_t         HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;
   Float_t         HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v;
   Float_t         HLT_BIT_HLT_DoubleIsoMu17_eta2p1_v;
   Float_t         HLT_ZmmHbbAll;
   Float_t         HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v;
   Float_t         HLT_HLT_Ele115_CaloIdVT_GsfTrkIdT;
   Float_t         HLT_WenHbbLowLumi;
   Float_t         HLT_ZnnHbbHighLumi;
   Float_t         HLT_HH4bLowLumi;
   Float_t         HLT_ZmmHbbHighLumi;
   Float_t         HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v;
   Float_t         HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v;
   Float_t         HLT_BIT_HLT_QuadPFJet_VBF_v;
   Float_t         HLT_BIT_HLT_L1_TripleJet_VBF_v;
   Float_t         HLT_BIT_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v;
   Float_t         HLT_BIT_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v;
   Float_t         HLT_VBFHbbAll;
   Float_t         HLT_VBFHbbHighLumi;
   Float_t         HLT_ttHhardonicHighLumi;
   Float_t         HLT_ZnnHbbLowLumi;
   Float_t         HLT_HH4bHighLumi;
   Float_t         met_pt;
   Float_t         met_eta;
   Float_t         met_phi;
   Float_t         met_mass;
   Float_t         met_sumEt;
   Float_t         met_genPt;
   Float_t         met_genPhi;
   Float_t         met_genEta;
   Float_t         met_shifted_UnclusteredEnUp_pt;
   Float_t         met_shifted_UnclusteredEnUp_phi;
   Float_t         met_shifted_UnclusteredEnUp_sumEt;
   Float_t         met_shifted_UnclusteredEnDown_pt;
   Float_t         met_shifted_UnclusteredEnDown_phi;
   Float_t         met_shifted_UnclusteredEnDown_sumEt;
   Float_t         met_shifted_JetResUp_pt;
   Float_t         met_shifted_JetResUp_phi;
   Float_t         met_shifted_JetResUp_sumEt;
   Float_t         met_shifted_JetResDown_pt;
   Float_t         met_shifted_JetResDown_phi;
   Float_t         met_shifted_JetResDown_sumEt;
   Float_t         met_shifted_JetEnUp_pt;
   Float_t         met_shifted_JetEnUp_phi;
   Float_t         met_shifted_JetEnUp_sumEt;
   Float_t         met_shifted_JetEnDown_pt;
   Float_t         met_shifted_JetEnDown_phi;
   Float_t         met_shifted_JetEnDown_sumEt;
   Float_t         met_shifted_MuonEnUp_pt;
   Float_t         met_shifted_MuonEnUp_phi;
   Float_t         met_shifted_MuonEnUp_sumEt;
   Float_t         met_shifted_MuonEnDown_pt;
   Float_t         met_shifted_MuonEnDown_phi;
   Float_t         met_shifted_MuonEnDown_sumEt;
   Float_t         met_shifted_ElectronEnUp_pt;
   Float_t         met_shifted_ElectronEnUp_phi;
   Float_t         met_shifted_ElectronEnUp_sumEt;
   Float_t         met_shifted_ElectronEnDown_pt;
   Float_t         met_shifted_ElectronEnDown_phi;
   Float_t         met_shifted_ElectronEnDown_sumEt;
   Float_t         met_shifted_TauEnUp_pt;
   Float_t         met_shifted_TauEnUp_phi;
   Float_t         met_shifted_TauEnUp_sumEt;
   Float_t         met_shifted_TauEnDown_pt;
   Float_t         met_shifted_TauEnDown_phi;
   Float_t         met_shifted_TauEnDown_sumEt;
   Int_t           nGenBQuarkFromHafterISR;
   Int_t           GenBQuarkFromHafterISR_pdgId[1];   //[nGenBQuarkFromHafterISR]
   Float_t         GenBQuarkFromHafterISR_pt[1];   //[nGenBQuarkFromHafterISR]
   Float_t         GenBQuarkFromHafterISR_eta[1];   //[nGenBQuarkFromHafterISR]
   Float_t         GenBQuarkFromHafterISR_phi[1];   //[nGenBQuarkFromHafterISR]
   Float_t         GenBQuarkFromHafterISR_mass[1];   //[nGenBQuarkFromHafterISR]
   Float_t         GenBQuarkFromHafterISR_charge[1];   //[nGenBQuarkFromHafterISR]
   Int_t           GenBQuarkFromHafterISR_status[1];   //[nGenBQuarkFromHafterISR]
   Int_t           npileUpVertex_ptHat;
   Float_t         pileUpVertex_ptHat[5];   //[npileUpVertex_ptHat]
   Int_t           nGenHiggsBoson;
   Int_t           GenHiggsBoson_pdgId[1];   //[nGenHiggsBoson]
   Float_t         GenHiggsBoson_pt[1];   //[nGenHiggsBoson]
   Float_t         GenHiggsBoson_eta[1];   //[nGenHiggsBoson]
   Float_t         GenHiggsBoson_phi[1];   //[nGenHiggsBoson]
   Float_t         GenHiggsBoson_mass[1];   //[nGenHiggsBoson]
   Float_t         GenHiggsBoson_charge[1];   //[nGenHiggsBoson]
   Int_t           GenHiggsBoson_status[1];   //[nGenHiggsBoson]
   Int_t           nGenLepFromTop;
   Int_t           GenLepFromTop_pdgId[1];   //[nGenLepFromTop]
   Float_t         GenLepFromTop_pt[1];   //[nGenLepFromTop]
   Float_t         GenLepFromTop_eta[1];   //[nGenLepFromTop]
   Float_t         GenLepFromTop_phi[1];   //[nGenLepFromTop]
   Float_t         GenLepFromTop_mass[1];   //[nGenLepFromTop]
   Float_t         GenLepFromTop_charge[1];   //[nGenLepFromTop]
   Int_t           GenLepFromTop_status[1];   //[nGenLepFromTop]
   Int_t           nGenVbosons;
   Int_t           GenVbosons_pdgId[2];   //[nGenVbosons]
   Float_t         GenVbosons_pt[2];   //[nGenVbosons]
   Float_t         GenVbosons_eta[2];   //[nGenVbosons]
   Float_t         GenVbosons_phi[2];   //[nGenVbosons]
   Float_t         GenVbosons_mass[2];   //[nGenVbosons]
   Float_t         GenVbosons_charge[2];   //[nGenVbosons]
   Int_t           GenVbosons_status[2];   //[nGenVbosons]
   Int_t           nFatjetCleanAK08prunedcalreg;
   Float_t         FatjetCleanAK08prunedcalreg_pt[5];   //[nFatjetCleanAK08prunedcalreg]
   Float_t         FatjetCleanAK08prunedcalreg_eta[5];   //[nFatjetCleanAK08prunedcalreg]
   Float_t         FatjetCleanAK08prunedcalreg_phi[5];   //[nFatjetCleanAK08prunedcalreg]
   Float_t         FatjetCleanAK08prunedcalreg_mass[5];   //[nFatjetCleanAK08prunedcalreg]
   Int_t           nFatjetAK08pruned;
   Float_t         FatjetAK08pruned_pt[6];   //[nFatjetAK08pruned]
   Float_t         FatjetAK08pruned_eta[6];   //[nFatjetAK08pruned]
   Float_t         FatjetAK08pruned_phi[6];   //[nFatjetAK08pruned]
   Float_t         FatjetAK08pruned_mass[6];   //[nFatjetAK08pruned]
   Int_t           nGenJet;
   Float_t         GenJet_charge[14];   //[nGenJet]
   Int_t           GenJet_status[14];   //[nGenJet]
   Int_t           GenJet_pdgId[14];   //[nGenJet]
   Float_t         GenJet_pt[14];   //[nGenJet]
   Float_t         GenJet_eta[14];   //[nGenJet]
   Float_t         GenJet_phi[14];   //[nGenJet]
   Float_t         GenJet_mass[14];   //[nGenJet]
   Int_t           GenJet_numBHadrons[14];   //[nGenJet]
   Int_t           GenJet_numCHadrons[14];   //[nGenJet]
   Int_t           GenJet_numBHadronsFromTop[14];   //[nGenJet]
   Int_t           GenJet_numCHadronsFromTop[14];   //[nGenJet]
   Int_t           GenJet_numBHadronsAfterTop[14];   //[nGenJet]
   Int_t           GenJet_numCHadronsAfterTop[14];   //[nGenJet]
   Float_t         GenJet_wNuPt[14];   //[nGenJet]
   Float_t         GenJet_wNuEta[14];   //[nGenJet]
   Float_t         GenJet_wNuPhi[14];   //[nGenJet]
   Float_t         GenJet_wNuM[14];   //[nGenJet]
   Int_t           nGenTau;
   Int_t           GenTau_pdgId[2];   //[nGenTau]
   Float_t         GenTau_pt[2];   //[nGenTau]
   Float_t         GenTau_eta[2];   //[nGenTau]
   Float_t         GenTau_phi[2];   //[nGenTau]
   Float_t         GenTau_mass[2];   //[nGenTau]
   Float_t         GenTau_charge[2];   //[nGenTau]
   Int_t           GenTau_status[2];   //[nGenTau]
   Int_t           nGenLep;
   Int_t           GenLep_pdgId[1];   //[nGenLep]
   Float_t         GenLep_pt[1];   //[nGenLep]
   Float_t         GenLep_eta[1];   //[nGenLep]
   Float_t         GenLep_phi[1];   //[nGenLep]
   Float_t         GenLep_mass[1];   //[nGenLep]
   Float_t         GenLep_charge[1];   //[nGenLep]
   Int_t           GenLep_status[1];   //[nGenLep]
   Int_t           nGenBQuarkFromH;
   Int_t           GenBQuarkFromH_pdgId[1];   //[nGenBQuarkFromH]
   Float_t         GenBQuarkFromH_pt[1];   //[nGenBQuarkFromH]
   Float_t         GenBQuarkFromH_eta[1];   //[nGenBQuarkFromH]
   Float_t         GenBQuarkFromH_phi[1];   //[nGenBQuarkFromH]
   Float_t         GenBQuarkFromH_mass[1];   //[nGenBQuarkFromH]
   Float_t         GenBQuarkFromH_charge[1];   //[nGenBQuarkFromH]
   Int_t           GenBQuarkFromH_status[1];   //[nGenBQuarkFromH]
   Int_t           nFatjetAK08prunedCal;
   Float_t         FatjetAK08prunedCal_pt[6];   //[nFatjetAK08prunedCal]
   Float_t         FatjetAK08prunedCal_eta[6];   //[nFatjetAK08prunedCal]
   Float_t         FatjetAK08prunedCal_phi[6];   //[nFatjetAK08prunedCal]
   Float_t         FatjetAK08prunedCal_mass[6];   //[nFatjetAK08prunedCal]
   Int_t           nincLeptons;
   Int_t           incLeptons_charge[8];   //[nincLeptons]
   Int_t           incLeptons_tightId[8];   //[nincLeptons]
   Int_t           incLeptons_eleCutIdCSA14_25ns_v1[8];   //[nincLeptons]
   Int_t           incLeptons_eleCutIdCSA14_50ns_v1[8];   //[nincLeptons]
   Int_t           incLeptons_eleCutIdSpring15_25ns_v1[8];   //[nincLeptons]
   Float_t         incLeptons_dxy[8];   //[nincLeptons]
   Float_t         incLeptons_dz[8];   //[nincLeptons]
   Float_t         incLeptons_edxy[8];   //[nincLeptons]
   Float_t         incLeptons_edz[8];   //[nincLeptons]
   Float_t         incLeptons_ip3d[8];   //[nincLeptons]
   Float_t         incLeptons_sip3d[8];   //[nincLeptons]
   Int_t           incLeptons_convVeto[8];   //[nincLeptons]
   Int_t           incLeptons_lostHits[8];   //[nincLeptons]
   Float_t         incLeptons_relIso03[8];   //[nincLeptons]
   Float_t         incLeptons_relIso04[8];   //[nincLeptons]
   Float_t         incLeptons_miniRelIso[8];   //[nincLeptons]
   Int_t           incLeptons_tightCharge[8];   //[nincLeptons]
   Int_t           incLeptons_mcMatchId[8];   //[nincLeptons]
   Int_t           incLeptons_mcMatchAny[8];   //[nincLeptons]
   Int_t           incLeptons_mcMatchTau[8];   //[nincLeptons]
   Float_t         incLeptons_mcPt[8];   //[nincLeptons]
   Int_t           incLeptons_mediumMuonId[8];   //[nincLeptons]
   Int_t           incLeptons_pdgId[8];   //[nincLeptons]
   Float_t         incLeptons_pt[8];   //[nincLeptons]
   Float_t         incLeptons_eta[8];   //[nincLeptons]
   Float_t         incLeptons_phi[8];   //[nincLeptons]
   Float_t         incLeptons_mass[8];   //[nincLeptons]
   Int_t           incLeptons_looseIdSusy[8];   //[nincLeptons]
   Int_t           incLeptons_looseIdPOG[8];   //[nincLeptons]
   Float_t         incLeptons_chargedHadRelIso03[8];   //[nincLeptons]
   Float_t         incLeptons_chargedHadRelIso04[8];   //[nincLeptons]
   Float_t         incLeptons_eleSieie[8];   //[nincLeptons]
   Float_t         incLeptons_e5x5[8];   //[nincLeptons]
   Float_t         incLeptons_e2x5Max[8];   //[nincLeptons]
   Float_t         incLeptons_e1x5[8];   //[nincLeptons]
   Float_t         incLeptons_isolTrkPt[8];   //[nincLeptons]
   Float_t         incLeptons_isolEmHadDepth1[8];   //[nincLeptons]
   Float_t         incLeptons_eleDEta[8];   //[nincLeptons]
   Float_t         incLeptons_eleDPhi[8];   //[nincLeptons]
   Float_t         incLeptons_eleHoE[8];   //[nincLeptons]
   Float_t         incLeptons_eleMissingHits[8];   //[nincLeptons]
   Float_t         incLeptons_eleChi2[8];   //[nincLeptons]
   Float_t         incLeptons_muonDX[8];   //[nincLeptons]
   Float_t         incLeptons_eleClusterEta[8];   //[nincLeptons]
   Float_t         incLeptons_eleClusterEnergy[8];   //[nincLeptons]
   Float_t         incLeptons_eleClusterDEta[8];   //[nincLeptons]
   Float_t         incLeptons_eleClusterDPhi[8];   //[nincLeptons]
   Float_t         incLeptons_nStations[8];   //[nincLeptons]
   Float_t         incLeptons_trkKink[8];   //[nincLeptons]
   Float_t         incLeptons_caloCompatibility[8];   //[nincLeptons]
   Float_t         incLeptons_globalTrackChi2[8];   //[nincLeptons]
   Float_t         incLeptons_nChamberHits[8];   //[nincLeptons]
   Int_t           incLeptons_isBarrelEle[8];   //[nincLeptons]
   Int_t           incLeptons_isEndCapEle[8];   //[nincLeptons]
   Int_t           incLeptons_isEcalDriven[8];   //[nincLeptons]
   Int_t           incLeptons_isMyGoodMuon[8];   //[nincLeptons]
   Int_t           incLeptons_isMyGoodMuon1[8];   //[nincLeptons]
   Int_t           incLeptons_isMyGoodElectron[8];   //[nincLeptons]
   Float_t         incLeptons_relPtError[8];   //[nincLeptons]
   Float_t         incLeptons_isPFMuon[8];   //[nincLeptons]
   Float_t         incLeptons_muon_dz[8];   //[nincLeptons]
   Int_t           incLeptons_isHighPtMuon[8];   //[nincLeptons]
   Float_t         incLeptons_muTrackIso[8];   //[nincLeptons]
   Float_t         incLeptons_isGlobalMuon[8];   //[nincLeptons]
   Float_t         incLeptons_isTrackerMuon[8];   //[nincLeptons]
   Float_t         incLeptons_pixelHits[8];   //[nincLeptons]
   Int_t           incLeptons_trackerLayers[8];   //[nincLeptons]
   Int_t           incLeptons_pixelLayers[8];   //[nincLeptons]
   Float_t         incLeptons_mvaTTH[8];   //[nincLeptons]
   Int_t           incLeptons_jetOverlapIdx[8];   //[nincLeptons]
   Float_t         incLeptons_jetPtRatio[8];   //[nincLeptons]
   Float_t         incLeptons_jetBTagCSV[8];   //[nincLeptons]
   Float_t         incLeptons_jetDR[8];   //[nincLeptons]
   Float_t         incLeptons_pfRelIso03[8];   //[nincLeptons]
   Float_t         incLeptons_pfRelIso04[8];   //[nincLeptons]
   Float_t         incLeptons_muonDB[8];   //[nincLeptons]
   Float_t         incLeptons_etaSc[8];   //[nincLeptons]
   Float_t         incLeptons_eleExpMissingInnerHits[8];   //[nincLeptons]
   Float_t         incLeptons_eleooEmooP[8];   //[nincLeptons]
   Float_t         incLeptons_muonTrackerLayers[8];   //[nincLeptons]
   Int_t           nGenTop;
   Int_t           GenTop_pdgId[1];   //[nGenTop]
   Float_t         GenTop_pt[1];   //[nGenTop]
   Float_t         GenTop_eta[1];   //[nGenTop]
   Float_t         GenTop_phi[1];   //[nGenTop]
   Float_t         GenTop_mass[1];   //[nGenTop]
   Float_t         GenTop_charge[1];   //[nGenTop]
   Int_t           GenTop_status[1];   //[nGenTop]
   Int_t           nGenLepFromTau;
   Int_t           GenLepFromTau_pdgId[1];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_pt[1];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_eta[1];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_phi[1];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_mass[1];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_charge[1];   //[nGenLepFromTau]
   Int_t           GenLepFromTau_status[1];   //[nGenLepFromTau]
   Int_t           nGenNuFromTop;
   Int_t           GenNuFromTop_pdgId[1];   //[nGenNuFromTop]
   Float_t         GenNuFromTop_pt[1];   //[nGenNuFromTop]
   Float_t         GenNuFromTop_eta[1];   //[nGenNuFromTop]
   Float_t         GenNuFromTop_phi[1];   //[nGenNuFromTop]
   Float_t         GenNuFromTop_mass[1];   //[nGenNuFromTop]
   Float_t         GenNuFromTop_charge[1];   //[nGenNuFromTop]
   Int_t           GenNuFromTop_status[1];   //[nGenNuFromTop]
   Int_t           nFatjetCleanAK08pruned;
   Float_t         FatjetCleanAK08pruned_pt[5];   //[nFatjetCleanAK08pruned]
   Float_t         FatjetCleanAK08pruned_eta[5];   //[nFatjetCleanAK08pruned]
   Float_t         FatjetCleanAK08pruned_phi[5];   //[nFatjetCleanAK08pruned]
   Float_t         FatjetCleanAK08pruned_mass[5];   //[nFatjetCleanAK08pruned]
   Int_t           nFatjetAK08ungroomed;
   Float_t         FatjetAK08ungroomed_pt[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_eta[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_phi[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_mass[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_chHEFrac[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_neHEFrac[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_chEmEFrac[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_neEmEFrac[6];   //[nFatjetAK08ungroomed]
   Int_t           FatjetAK08ungroomed_chMult[6];   //[nFatjetAK08ungroomed]
   Int_t           FatjetAK08ungroomed_neMult[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_muEFrac[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_tau1[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_tau2[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_tau3[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_msoftdrop[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_mpruned[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_mtrimmed[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_mfiltered[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_bbtag[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_PFLepton_ptrel[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_z_ratio[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_tau_dot[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_SV_mass_0[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_SV_EnergyRatio_0[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_SV_EnergyRatio_1[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_PFLepton_IP2D[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_tau_21[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_nSL[6];   //[nFatjetAK08ungroomed]
   Float_t         FatjetAK08ungroomed_vertexNTracks[6];   //[nFatjetAK08ungroomed]
   Int_t           nGenHiggsSisters;
   Int_t           GenHiggsSisters_pdgId[1];   //[nGenHiggsSisters]
   Float_t         GenHiggsSisters_pt[1];   //[nGenHiggsSisters]
   Float_t         GenHiggsSisters_eta[1];   //[nGenHiggsSisters]
   Float_t         GenHiggsSisters_phi[1];   //[nGenHiggsSisters]
   Float_t         GenHiggsSisters_mass[1];   //[nGenHiggsSisters]
   Float_t         GenHiggsSisters_charge[1];   //[nGenHiggsSisters]
   Int_t           GenHiggsSisters_status[1];   //[nGenHiggsSisters]
   Int_t           nSubjetCleanAK08pruned;
   Float_t         SubjetCleanAK08pruned_pt[10];   //[nSubjetCleanAK08pruned]
   Float_t         SubjetCleanAK08pruned_eta[10];   //[nSubjetCleanAK08pruned]
   Float_t         SubjetCleanAK08pruned_phi[10];   //[nSubjetCleanAK08pruned]
   Float_t         SubjetCleanAK08pruned_mass[10];   //[nSubjetCleanAK08pruned]
   Float_t         SubjetCleanAK08pruned_btag[10];   //[nSubjetCleanAK08pruned]
   Int_t           nGenNu;
   Int_t           GenNu_pdgId[4];   //[nGenNu]
   Float_t         GenNu_pt[4];   //[nGenNu]
   Float_t         GenNu_eta[4];   //[nGenNu]
   Float_t         GenNu_phi[4];   //[nGenNu]
   Float_t         GenNu_mass[4];   //[nGenNu]
   Float_t         GenNu_charge[4];   //[nGenNu]
   Int_t           GenNu_status[4];   //[nGenNu]
   Int_t           nselLeptons;
   Int_t           selLeptons_charge[8];   //[nselLeptons]
   Int_t           selLeptons_tightId[8];   //[nselLeptons]
   Int_t           selLeptons_eleCutIdCSA14_25ns_v1[8];   //[nselLeptons]
   Int_t           selLeptons_eleCutIdCSA14_50ns_v1[8];   //[nselLeptons]
   Int_t           selLeptons_eleCutIdSpring15_25ns_v1[8];   //[nselLeptons]
   Float_t         selLeptons_dxy[8];   //[nselLeptons]
   Float_t         selLeptons_dz[8];   //[nselLeptons]
   Float_t         selLeptons_edxy[8];   //[nselLeptons]
   Float_t         selLeptons_edz[8];   //[nselLeptons]
   Float_t         selLeptons_ip3d[8];   //[nselLeptons]
   Float_t         selLeptons_sip3d[8];   //[nselLeptons]
   Int_t           selLeptons_convVeto[8];   //[nselLeptons]
   Int_t           selLeptons_lostHits[8];   //[nselLeptons]
   Float_t         selLeptons_relIso03[8];   //[nselLeptons]
   Float_t         selLeptons_relIso04[8];   //[nselLeptons]
   Float_t         selLeptons_miniRelIso[8];   //[nselLeptons]
   Int_t           selLeptons_tightCharge[8];   //[nselLeptons]
   Int_t           selLeptons_mcMatchId[8];   //[nselLeptons]
   Int_t           selLeptons_mcMatchAny[8];   //[nselLeptons]
   Int_t           selLeptons_mcMatchTau[8];   //[nselLeptons]
   Float_t         selLeptons_mcPt[8];   //[nselLeptons]
   Int_t           selLeptons_mediumMuonId[8];   //[nselLeptons]
   Int_t           selLeptons_pdgId[8];   //[nselLeptons]
   Float_t         selLeptons_pt[8];   //[nselLeptons]
   Float_t         selLeptons_eta[8];   //[nselLeptons]
   Float_t         selLeptons_phi[8];   //[nselLeptons]
   Float_t         selLeptons_mass[8];   //[nselLeptons]
   Int_t           selLeptons_looseIdSusy[8];   //[nselLeptons]
   Int_t           selLeptons_looseIdPOG[8];   //[nselLeptons]
   Float_t         selLeptons_chargedHadRelIso03[8];   //[nselLeptons]
   Float_t         selLeptons_chargedHadRelIso04[8];   //[nselLeptons]
   Float_t         selLeptons_eleSieie[8];   //[nselLeptons]
   Float_t         selLeptons_e5x5[8];   //[nselLeptons]
   Float_t         selLeptons_e2x5Max[8];   //[nselLeptons]
   Float_t         selLeptons_e1x5[8];   //[nselLeptons]
   Float_t         selLeptons_isolTrkPt[8];   //[nselLeptons]
   Float_t         selLeptons_isolEmHadDepth1[8];   //[nselLeptons]
   Float_t         selLeptons_eleDEta[8];   //[nselLeptons]
   Float_t         selLeptons_eleDPhi[8];   //[nselLeptons]
   Float_t         selLeptons_eleHoE[8];   //[nselLeptons]
   Float_t         selLeptons_eleMissingHits[8];   //[nselLeptons]
   Float_t         selLeptons_eleChi2[8];   //[nselLeptons]
   Float_t         selLeptons_muonDX[8];   //[nselLeptons]
   Float_t         selLeptons_eleClusterEta[8];   //[nselLeptons]
   Float_t         selLeptons_eleClusterEnergy[8];   //[nselLeptons]
   Float_t         selLeptons_eleClusterDEta[8];   //[nselLeptons]
   Float_t         selLeptons_eleClusterDPhi[8];   //[nselLeptons]
   Float_t         selLeptons_nStations[8];   //[nselLeptons]
   Float_t         selLeptons_trkKink[8];   //[nselLeptons]
   Float_t         selLeptons_caloCompatibility[8];   //[nselLeptons]
   Float_t         selLeptons_globalTrackChi2[8];   //[nselLeptons]
   Float_t         selLeptons_nChamberHits[8];   //[nselLeptons]
   Int_t           selLeptons_isBarrelEle[8];   //[nselLeptons]
   Int_t           selLeptons_isEndCapEle[8];   //[nselLeptons]
   Int_t           selLeptons_isEcalDriven[8];   //[nselLeptons]
   Int_t           selLeptons_isMyGoodMuon[8];   //[nselLeptons]
   Int_t           selLeptons_isMyGoodMuon1[8];   //[nselLeptons]
   Int_t           selLeptons_isMyGoodElectron[8];   //[nselLeptons]
   Float_t         selLeptons_relPtError[8];   //[nselLeptons]
   Float_t         selLeptons_isPFMuon[8];   //[nselLeptons]
   Float_t         selLeptons_muon_dz[8];   //[nselLeptons]
   Int_t           selLeptons_isHighPtMuon[8];   //[nselLeptons]
   Float_t         selLeptons_muTrackIso[8];   //[nselLeptons]
   Float_t         selLeptons_isGlobalMuon[8];   //[nselLeptons]
   Float_t         selLeptons_isTrackerMuon[8];   //[nselLeptons]
   Float_t         selLeptons_pixelHits[8];   //[nselLeptons]
   Int_t           selLeptons_trackerLayers[8];   //[nselLeptons]
   Int_t           selLeptons_pixelLayers[8];   //[nselLeptons]
   Float_t         selLeptons_mvaTTH[8];   //[nselLeptons]
   Int_t           selLeptons_jetOverlapIdx[8];   //[nselLeptons]
   Float_t         selLeptons_jetPtRatio[8];   //[nselLeptons]
   Float_t         selLeptons_jetBTagCSV[8];   //[nselLeptons]
   Float_t         selLeptons_jetDR[8];   //[nselLeptons]
   Float_t         selLeptons_pfRelIso03[8];   //[nselLeptons]
   Float_t         selLeptons_pfRelIso04[8];   //[nselLeptons]
   Float_t         selLeptons_muonDB[8];   //[nselLeptons]
   Float_t         selLeptons_etaSc[8];   //[nselLeptons]
   Float_t         selLeptons_eleExpMissingInnerHits[8];   //[nselLeptons]
   Float_t         selLeptons_eleooEmooP[8];   //[nselLeptons]
   Float_t         selLeptons_muonTrackerLayers[8];   //[nselLeptons]
   Int_t           nFatjetCleanAK08ungroomed;
   Float_t         FatjetCleanAK08ungroomed_pt[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_eta[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_phi[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_mass[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_chHEFrac[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_neHEFrac[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_chEmEFrac[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_neEmEFrac[5];   //[nFatjetCleanAK08ungroomed]
   Int_t           FatjetCleanAK08ungroomed_chMult[5];   //[nFatjetCleanAK08ungroomed]
   Int_t           FatjetCleanAK08ungroomed_neMult[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_muEFrac[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_tau1[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_tau2[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_tau3[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_msoftdrop[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_mpruned[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_mtrimmed[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_mfiltered[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_bbtag[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_PFLepton_ptrel[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_z_ratio[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_tau_dot[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_SV_mass_0[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_SV_EnergyRatio_0[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_SV_EnergyRatio_1[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_PFLepton_IP2D[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_tau_21[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_nSL[5];   //[nFatjetCleanAK08ungroomed]
   Float_t         FatjetCleanAK08ungroomed_vertexNTracks[5];   //[nFatjetCleanAK08ungroomed]
   Int_t           nTauGood;
   Int_t           TauGood_charge[11];   //[nTauGood]
   Int_t           TauGood_decayMode[11];   //[nTauGood]
   Int_t           TauGood_idDecayMode[11];   //[nTauGood]
   Int_t           TauGood_idDecayModeNewDMs[11];   //[nTauGood]
   Float_t         TauGood_dxy[11];   //[nTauGood]
   Float_t         TauGood_dz[11];   //[nTauGood]
   Int_t           TauGood_idMVA[11];   //[nTauGood]
   Int_t           TauGood_idMVANewDM[11];   //[nTauGood]
   Int_t           TauGood_idCI3hit[11];   //[nTauGood]
   Int_t           TauGood_idAntiMu[11];   //[nTauGood]
   Int_t           TauGood_idAntiE[11];   //[nTauGood]
   Float_t         TauGood_isoCI3hit[11];   //[nTauGood]
   Int_t           TauGood_mcMatchId[11];   //[nTauGood]
   Int_t           TauGood_pdgId[11];   //[nTauGood]
   Float_t         TauGood_pt[11];   //[nTauGood]
   Float_t         TauGood_eta[11];   //[nTauGood]
   Float_t         TauGood_phi[11];   //[nTauGood]
   Float_t         TauGood_mass[11];   //[nTauGood]
   Int_t           nFatjetCleanAK08prunedcal;
   Float_t         FatjetCleanAK08prunedcal_pt[5];   //[nFatjetCleanAK08prunedcal]
   Float_t         FatjetCleanAK08prunedcal_eta[5];   //[nFatjetCleanAK08prunedcal]
   Float_t         FatjetCleanAK08prunedcal_phi[5];   //[nFatjetCleanAK08prunedcal]
   Float_t         FatjetCleanAK08prunedcal_mass[5];   //[nFatjetCleanAK08prunedcal]
   Int_t           nSubjetAK08pruned;
   Float_t         SubjetAK08pruned_pt[10];   //[nSubjetAK08pruned]
   Float_t         SubjetAK08pruned_eta[10];   //[nSubjetAK08pruned]
   Float_t         SubjetAK08pruned_phi[10];   //[nSubjetAK08pruned]
   Float_t         SubjetAK08pruned_mass[10];   //[nSubjetAK08pruned]
   Float_t         SubjetAK08pruned_btag[10];   //[nSubjetAK08pruned]
   Int_t           nGenStatus2bHad;
   Int_t           GenStatus2bHad_pdgId[11];   //[nGenStatus2bHad]
   Float_t         GenStatus2bHad_pt[11];   //[nGenStatus2bHad]
   Float_t         GenStatus2bHad_eta[11];   //[nGenStatus2bHad]
   Float_t         GenStatus2bHad_phi[11];   //[nGenStatus2bHad]
   Float_t         GenStatus2bHad_mass[11];   //[nGenStatus2bHad]
   Float_t         GenStatus2bHad_charge[11];   //[nGenStatus2bHad]
   Int_t           GenStatus2bHad_status[11];   //[nGenStatus2bHad]
   Int_t           nGenHadTaus;
   Float_t         GenHadTaus_charge[3];   //[nGenHadTaus]
   Int_t           GenHadTaus_status[3];   //[nGenHadTaus]
   Int_t           GenHadTaus_pdgId[3];   //[nGenHadTaus]
   Float_t         GenHadTaus_pt[3];   //[nGenHadTaus]
   Float_t         GenHadTaus_eta[3];   //[nGenHadTaus]
   Float_t         GenHadTaus_phi[3];   //[nGenHadTaus]
   Float_t         GenHadTaus_mass[3];   //[nGenHadTaus]
   Int_t           GenHadTaus_decayMode[3];   //[nGenHadTaus]
   Int_t           nJet;
   Int_t           Jet_id[23];   //[nJet]
   Int_t           Jet_puId[23];   //[nJet]
   Float_t         Jet_btagCSV[23];   //[nJet]
   Float_t         Jet_btagCMVA[23];   //[nJet]
   Float_t         Jet_rawPt[23];   //[nJet]
   Float_t         Jet_mcPt[23];   //[nJet]
   Int_t           Jet_mcFlavour[23];   //[nJet]
   Int_t           Jet_mcMatchId[23];   //[nJet]
   Float_t         Jet_corr_JECUp[23];   //[nJet]
   Float_t         Jet_corr_JECDown[23];   //[nJet]
   Float_t         Jet_corr[23];   //[nJet]
   Float_t         Jet_corr_JERUp[23];   //[nJet]
   Float_t         Jet_corr_JERDown[23];   //[nJet]
   Float_t         Jet_corr_JER[23];   //[nJet]
   Float_t         Jet_pt[23];   //[nJet]
   Float_t         Jet_eta[23];   //[nJet]
   Float_t         Jet_phi[23];   //[nJet]
   Float_t         Jet_mass[23];   //[nJet]
   Int_t           Jet_idxFirstTauMatch[23];   //[nJet]
   Int_t           Jet_heppyFlavour[23];   //[nJet]
   Int_t           Jet_hadronFlavour[23];   //[nJet]
   Float_t         Jet_btagBDT[23];   //[nJet]
   Float_t         Jet_btagProb[23];   //[nJet]
   Float_t         Jet_btagBProb[23];   //[nJet]
   Float_t         Jet_btagSoftEl[23];   //[nJet]
   Float_t         Jet_btagSoftMu[23];   //[nJet]
   Float_t         Jet_btagnew[23];   //[nJet]
   Float_t         Jet_btagCSVV0[23];   //[nJet]
   Int_t           Jet_isMyGoodLooseJet[23];   //[nJet]
   Int_t           Jet_isMyGoodTightJet[23];   //[nJet]
   Int_t           Jet_isMyGoodTightLepVetoJet[23];   //[nJet]
   Float_t         Jet_chHEF[23];   //[nJet]
   Float_t         Jet_neHEF[23];   //[nJet]
   Float_t         Jet_chEmEF[23];   //[nJet]
   Float_t         Jet_neEmEF[23];   //[nJet]
   Float_t         Jet_muEF[23];   //[nJet]
   Int_t           Jet_chMult[23];   //[nJet]
   Int_t           Jet_neMult[23];   //[nJet]
   Float_t         Jet_leadTrackPt[23];   //[nJet]
   Float_t         Jet_mcEta[23];   //[nJet]
   Float_t         Jet_mcPhi[23];   //[nJet]
   Float_t         Jet_mcM[23];   //[nJet]
   Float_t         Jet_leptonPdgId[23];   //[nJet]
   Float_t         Jet_leptonPt[23];   //[nJet]
   Float_t         Jet_leptonPtRel[23];   //[nJet]
   Float_t         Jet_leptonPtRelInv[23];   //[nJet]
   Float_t         Jet_leptonDeltaR[23];   //[nJet]
   Float_t         Jet_leptonDeltaPhi[23];   //[nJet]
   Float_t         Jet_leptonDeltaEta[23];   //[nJet]
   Float_t         Jet_vtxMass[23];   //[nJet]
   Float_t         Jet_vtxNtracks[23];   //[nJet]
   Float_t         Jet_vtxPt[23];   //[nJet]
   Float_t         Jet_vtx3DSig[23];   //[nJet]
   Float_t         Jet_vtx3DVal[23];   //[nJet]
   Float_t         Jet_vtxPosX[23];   //[nJet]
   Float_t         Jet_vtxPosY[23];   //[nJet]
   Float_t         Jet_vtxPosZ[23];   //[nJet]
   Float_t         Jet_pullVectorPhi[23];   //[nJet]
   Float_t         Jet_pullVectorMag[23];   //[nJet]
   Float_t         Jet_qgl[23];   //[nJet]
   Float_t         Jet_ptd[23];   //[nJet]
   Float_t         Jet_axis2[23];   //[nJet]
   Int_t           Jet_mult[23];   //[nJet]
   Int_t           Jet_numberOfDaughters[23];   //[nJet]
   Int_t           Jet_btagIdx[23];   //[nJet]
   Int_t           Jet_mcIdx[23];   //[nJet]
   Float_t         Jet_pt_reg[23];   //[nJet]
   Float_t         Jet_pt_regVBF[23];   //[nJet]
   Float_t         Jet_blike_VBF[23];   //[nJet]
   Float_t         Jet_bTagWeightJESUp[23];   //[nJet]
   Float_t         Jet_bTagWeightJESDown[23];   //[nJet]
   Float_t         Jet_bTagWeightLFUp[23];   //[nJet]
   Float_t         Jet_bTagWeightLFDown[23];   //[nJet]
   Float_t         Jet_bTagWeightHFUp[23];   //[nJet]
   Float_t         Jet_bTagWeightHFDown[23];   //[nJet]
   Float_t         Jet_bTagWeightStats1Up[23];   //[nJet]
   Float_t         Jet_bTagWeightStats1Down[23];   //[nJet]
   Float_t         Jet_bTagWeightStats2Up[23];   //[nJet]
   Float_t         Jet_bTagWeightStats2Down[23];   //[nJet]
   Float_t         Jet_bTagWeightcErr1Up[23];   //[nJet]
   Float_t         Jet_bTagWeightcErr1Down[23];   //[nJet]
   Float_t         Jet_bTagWeightcErr2Up[23];   //[nJet]
   Float_t         Jet_bTagWeightcErr2Down[23];   //[nJet]
   Float_t         Jet_bTagWeight[23];   //[nJet]
   Int_t           nFatjetGenAK08Pruned;
   Float_t         FatjetGenAK08Pruned_pt[8];   //[nFatjetGenAK08Pruned]
   Float_t         FatjetGenAK08Pruned_eta[8];   //[nFatjetGenAK08Pruned]
   Float_t         FatjetGenAK08Pruned_phi[8];   //[nFatjetGenAK08Pruned]
   Float_t         FatjetGenAK08Pruned_mass[8];   //[nFatjetGenAK08Pruned]
   Int_t           npileUpVertex_z;
   Float_t         pileUpVertex_z[5];   //[npileUpVertex_z]
   Int_t           nmyGenJetAK08;
   Float_t         myGenJetAK08_pt[10];   //[nmyGenJetAK08]
   Float_t         myGenJetAK08_eta[10];   //[nmyGenJetAK08]
   Float_t         myGenJetAK08_phi[10];   //[nmyGenJetAK08]
   Float_t         myGenJetAK08_mass[10];   //[nmyGenJetAK08]
   Float_t         myGenJetAK08_tau1[10];   //[nmyGenJetAK08]
   Float_t         myGenJetAK08_tau2[10];   //[nmyGenJetAK08]
   Float_t         myGenJetAK08_tau3[10];   //[nmyGenJetAK08]
   Int_t           nprimaryVertices;
   Float_t         primaryVertices_x[4];   //[nprimaryVertices]
   Float_t         primaryVertices_y[4];   //[nprimaryVertices]
   Float_t         primaryVertices_z[4];   //[nprimaryVertices]
   Float_t         primaryVertices_isFake[4];   //[nprimaryVertices]
   Float_t         primaryVertices_ndof[4];   //[nprimaryVertices]
   Float_t         primaryVertices_Rho[4];   //[nprimaryVertices]
   Int_t           nFatjetAK08prunedReg;
   Float_t         FatjetAK08prunedReg_pt[6];   //[nFatjetAK08prunedReg]
   Float_t         FatjetAK08prunedReg_eta[6];   //[nFatjetAK08prunedReg]
   Float_t         FatjetAK08prunedReg_phi[6];   //[nFatjetAK08prunedReg]
   Float_t         FatjetAK08prunedReg_mass[6];   //[nFatjetAK08prunedReg]
   Int_t           nGenWZQuark;
   Int_t           GenWZQuark_pdgId[2];   //[nGenWZQuark]
   Float_t         GenWZQuark_pt[2];   //[nGenWZQuark]
   Float_t         GenWZQuark_eta[2];   //[nGenWZQuark]
   Float_t         GenWZQuark_phi[2];   //[nGenWZQuark]
   Float_t         GenWZQuark_mass[2];   //[nGenWZQuark]
   Float_t         GenWZQuark_charge[2];   //[nGenWZQuark]
   Int_t           GenWZQuark_status[2];   //[nGenWZQuark]
   Int_t           nGenLepFromWPrime;
   Int_t           GenLepFromWPrime_pdgId[1];   //[nGenLepFromWPrime]
   Float_t         GenLepFromWPrime_pt[1];   //[nGenLepFromWPrime]
   Float_t         GenLepFromWPrime_eta[1];   //[nGenLepFromWPrime]
   Float_t         GenLepFromWPrime_phi[1];   //[nGenLepFromWPrime]
   Float_t         GenLepFromWPrime_mass[1];   //[nGenLepFromWPrime]
   Float_t         GenLepFromWPrime_charge[1];   //[nGenLepFromWPrime]
   Int_t           GenLepFromWPrime_status[1];   //[nGenLepFromWPrime]
   Int_t           nGenBQuarkFromTop;
   Int_t           GenBQuarkFromTop_pdgId[1];   //[nGenBQuarkFromTop]
   Float_t         GenBQuarkFromTop_pt[1];   //[nGenBQuarkFromTop]
   Float_t         GenBQuarkFromTop_eta[1];   //[nGenBQuarkFromTop]
   Float_t         GenBQuarkFromTop_phi[1];   //[nGenBQuarkFromTop]
   Float_t         GenBQuarkFromTop_mass[1];   //[nGenBQuarkFromTop]
   Float_t         GenBQuarkFromTop_charge[1];   //[nGenBQuarkFromTop]
   Int_t           GenBQuarkFromTop_status[1];   //[nGenBQuarkFromTop]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_xsec;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_json;   //!
   TBranch        *b_nPU0;   //!
   TBranch        *b_nPVs;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_lheNj;   //!
   TBranch        *b_lheNb;   //!
   TBranch        *b_lheNc;   //!
   TBranch        *b_lheNg;   //!
   TBranch        *b_lheNl;   //!
   TBranch        *b_lheV_pt;   //!
   TBranch        *b_lheHT;   //!
   TBranch        *b_mhtJet30;   //!
   TBranch        *b_mhtPhiJet30;   //!
   TBranch        *b_htJet30;   //!
   TBranch        *b_met_rawpt;   //!
   TBranch        *b_metPuppi_pt;   //!
   TBranch        *b_metPuppi_phi;   //!
   TBranch        *b_metPuppi_rawpt;   //!
   TBranch        *b_metType1p2_pt;   //!
   TBranch        *b_metNoPU_pt;   //!
   TBranch        *b_metNoPU_phi;   //!
   TBranch        *b_tkMet_pt;   //!
   TBranch        *b_tkMet_phi;   //!
   TBranch        *b_tkMetPVchs_pt;   //!
   TBranch        *b_tkMetPVchs_phi;   //!
   TBranch        *b_Flag_hbheIsoFilter;   //!
   TBranch        *b_Flag_hbheFilterNew;   //!
   TBranch        *b_simPrimaryVertex_z;   //!
   TBranch        *b_genHiggsDecayMode;   //!
   TBranch        *b_bTagWeight_LFUp;   //!
   TBranch        *b_bTagWeight_Stats2Down;   //!
   TBranch        *b_bTagWeight_LFDown;   //!
   TBranch        *b_bTagWeight_HFUp;   //!
   TBranch        *b_bTagWeight_cErr1Down;   //!
   TBranch        *b_bTagWeight_JESDown;   //!
   TBranch        *b_bTagWeight_cErr1Up;   //!
   TBranch        *b_bTagWeight;   //!
   TBranch        *b_bTagWeight_HFDown;   //!
   TBranch        *b_bTagWeight_Stats2Up;   //!
   TBranch        *b_bTagWeight_cErr2Up;   //!
   TBranch        *b_bTagWeight_JESUp;   //!
   TBranch        *b_bTagWeight_Stats1Up;   //!
   TBranch        *b_bTagWeight_Stats1Down;   //!
   TBranch        *b_bTagWeight_cErr2Down;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_trackingFailureFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_METFilters;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_HLT_BIT_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_BTagCSV0p7_v;   //!
   TBranch        *b_HLT_BIT_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET90_PFMHT90_IDLoose_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET100_PFMHT100_IDLoose_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET110_PFMHT110_IDLoose_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET120_PFMHT120_IDLoose_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET170_NoiseCleaned_v;   //!
   TBranch        *b_HLT_BIT_HLT_DiCentralPFJet70_PFMET120_NoiseCleaned_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT350_PFMET120_NoiseCleaned_v;   //!
   TBranch        *b_HLT_ZnnHbbAll;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT450_SixJet40_PFBTagCSV_v;   //!
   TBranch        *b_HLT_ttHhardonicLowLumi;   //!
   TBranch        *b_HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v;   //!
   TBranch        *b_HLT_WtaunHbbHighLumi;   //!
   TBranch        *b_HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v;   //!
   TBranch        *b_HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v;   //!
   TBranch        *b_HLT_VBFHbbLowLumi;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT750_4Jet_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT900_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet40_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet60_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet80_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet140_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet200_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet260_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet320_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet400_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFJet450_v;   //!
   TBranch        *b_HLT_hadronic;   //!
   TBranch        *b_HLT_BIT_HLT_QuadJet45_TripleCSV0p5_v;   //!
   TBranch        *b_HLT_BIT_HLT_QuadJet45_DoubleCSV0p5_v;   //!
   TBranch        *b_HLT_BIT_HLT_DoubleJet90_Double30_TripleCSV0p5_v;   //!
   TBranch        *b_HLT_BIT_HLT_DoubleJet90_Double30_DoubleCSV0p5_v;   //!
   TBranch        *b_HLT_HH4bAll;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu24_eta2p1_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v;   //!
   TBranch        *b_HLT_BIT_HLT_Mu24_eta2p1_v;   //!
   TBranch        *b_HLT_BIT_HLT_TkMu24_eta2p1_v;   //!
   TBranch        *b_HLT_BIT_HLT_Mu24_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu27_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoTkMu27_v;   //!
   TBranch        *b_HLT_BIT_HLT_TkMu27_v;   //!
   TBranch        *b_HLT_BIT_HLT_Mu27_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu20_eta2p1_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v;   //!
   TBranch        *b_HLT_BIT_HLT_Mu20_v;   //!
   TBranch        *b_HLT_BIT_HLT_TkMu20_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu20_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoTkMu20_v;   //!
   TBranch        *b_HLT_BIT_HLT_Mu40_eta2p1_PFJet200_PFJet50_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu16_eta2p1_CaloMET30_v;   //!
   TBranch        *b_HLT_BIT_HLT_Mu16_eta2p1_CaloMET30_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMET120_NoiseCleaned_Mu5_v;   //!
   TBranch        *b_HLT_BIT_HLT_MonoCentralPFJet80_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v;   //!
   TBranch        *b_HLT_BIT_HLT_MonoCentralPFJet80_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v;   //!
   TBranch        *b_HLT_WmnHbbAll;   //!
   TBranch        *b_HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v;   //!
   TBranch        *b_HLT_WenHbbHighLumi;   //!
   TBranch        *b_HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_HLT_BIT_HLT_DoubleEle24_22_eta2p1_WP75_Gsf_v;   //!
   TBranch        *b_HLT_ZeeHbbLowLumi;   //!
   TBranch        *b_HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele27_WP85_Gsf_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele105_CaloIdVT_GsfTrkIdT_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v;   //!
   TBranch        *b_HLT_WenHbbAll;   //!
   TBranch        *b_HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v;   //!
   TBranch        *b_HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_v;   //!
   TBranch        *b_HLT_WtaunHbbAll;   //!
   TBranch        *b_HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_ZeeHbbAll;   //!
   TBranch        *b_HLT_WtaunHbbLowLumi;   //!
   TBranch        *b_HLT_WmnHbbLowLumi;   //!
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_ttHleptonic;   //!
   TBranch        *b_HLT_ZeeHbbHighLumi;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT450_SixJet40_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT400_SixJet30_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT350_v;   //!
   TBranch        *b_HLT_ttHhardonicAll;   //!
   TBranch        *b_HLT_BIT_HLT_Mu50_v;   //!
   TBranch        *b_HLT_HLT_Mu50;   //!
   TBranch        *b_HLT_BIT_HLT_Mu45_eta2p1_v;   //!
   TBranch        *b_HLT_HLT_Mu45_eta2p1;   //!
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v;   //!
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v;   //!
   TBranch        *b_HLT_ZmmHbbLowLumi;   //!
   TBranch        *b_HLT_WmnHbbHighLumi;   //!
   TBranch        *b_HLT_BIT_HLT_Mu17_TkMu8_DZ_v;   //!
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;   //!
   TBranch        *b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v;   //!
   TBranch        *b_HLT_BIT_HLT_DoubleIsoMu17_eta2p1_v;   //!
   TBranch        *b_HLT_ZmmHbbAll;   //!
   TBranch        *b_HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v;   //!
   TBranch        *b_HLT_HLT_Ele115_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_WenHbbLowLumi;   //!
   TBranch        *b_HLT_ZnnHbbHighLumi;   //!
   TBranch        *b_HLT_HH4bLowLumi;   //!
   TBranch        *b_HLT_ZmmHbbHighLumi;   //!
   TBranch        *b_HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v;   //!
   TBranch        *b_HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v;   //!
   TBranch        *b_HLT_BIT_HLT_QuadPFJet_VBF_v;   //!
   TBranch        *b_HLT_BIT_HLT_L1_TripleJet_VBF_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v;   //!
   TBranch        *b_HLT_BIT_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v;   //!
   TBranch        *b_HLT_VBFHbbAll;   //!
   TBranch        *b_HLT_VBFHbbHighLumi;   //!
   TBranch        *b_HLT_ttHhardonicHighLumi;   //!
   TBranch        *b_HLT_ZnnHbbLowLumi;   //!
   TBranch        *b_HLT_HH4bHighLumi;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_mass;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_met_genPt;   //!
   TBranch        *b_met_genPhi;   //!
   TBranch        *b_met_genEta;   //!
   TBranch        *b_met_shifted_UnclusteredEnUp_pt;   //!
   TBranch        *b_met_shifted_UnclusteredEnUp_phi;   //!
   TBranch        *b_met_shifted_UnclusteredEnUp_sumEt;   //!
   TBranch        *b_met_shifted_UnclusteredEnDown_pt;   //!
   TBranch        *b_met_shifted_UnclusteredEnDown_phi;   //!
   TBranch        *b_met_shifted_UnclusteredEnDown_sumEt;   //!
   TBranch        *b_met_shifted_JetResUp_pt;   //!
   TBranch        *b_met_shifted_JetResUp_phi;   //!
   TBranch        *b_met_shifted_JetResUp_sumEt;   //!
   TBranch        *b_met_shifted_JetResDown_pt;   //!
   TBranch        *b_met_shifted_JetResDown_phi;   //!
   TBranch        *b_met_shifted_JetResDown_sumEt;   //!
   TBranch        *b_met_shifted_JetEnUp_pt;   //!
   TBranch        *b_met_shifted_JetEnUp_phi;   //!
   TBranch        *b_met_shifted_JetEnUp_sumEt;   //!
   TBranch        *b_met_shifted_JetEnDown_pt;   //!
   TBranch        *b_met_shifted_JetEnDown_phi;   //!
   TBranch        *b_met_shifted_JetEnDown_sumEt;   //!
   TBranch        *b_met_shifted_MuonEnUp_pt;   //!
   TBranch        *b_met_shifted_MuonEnUp_phi;   //!
   TBranch        *b_met_shifted_MuonEnUp_sumEt;   //!
   TBranch        *b_met_shifted_MuonEnDown_pt;   //!
   TBranch        *b_met_shifted_MuonEnDown_phi;   //!
   TBranch        *b_met_shifted_MuonEnDown_sumEt;   //!
   TBranch        *b_met_shifted_ElectronEnUp_pt;   //!
   TBranch        *b_met_shifted_ElectronEnUp_phi;   //!
   TBranch        *b_met_shifted_ElectronEnUp_sumEt;   //!
   TBranch        *b_met_shifted_ElectronEnDown_pt;   //!
   TBranch        *b_met_shifted_ElectronEnDown_phi;   //!
   TBranch        *b_met_shifted_ElectronEnDown_sumEt;   //!
   TBranch        *b_met_shifted_TauEnUp_pt;   //!
   TBranch        *b_met_shifted_TauEnUp_phi;   //!
   TBranch        *b_met_shifted_TauEnUp_sumEt;   //!
   TBranch        *b_met_shifted_TauEnDown_pt;   //!
   TBranch        *b_met_shifted_TauEnDown_phi;   //!
   TBranch        *b_met_shifted_TauEnDown_sumEt;   //!
   TBranch        *b_nGenBQuarkFromHafterISR;   //!
   TBranch        *b_GenBQuarkFromHafterISR_pdgId;   //!
   TBranch        *b_GenBQuarkFromHafterISR_pt;   //!
   TBranch        *b_GenBQuarkFromHafterISR_eta;   //!
   TBranch        *b_GenBQuarkFromHafterISR_phi;   //!
   TBranch        *b_GenBQuarkFromHafterISR_mass;   //!
   TBranch        *b_GenBQuarkFromHafterISR_charge;   //!
   TBranch        *b_GenBQuarkFromHafterISR_status;   //!
   TBranch        *b_npileUpVertex_ptHat;   //!
   TBranch        *b_pileUpVertex_ptHat;   //!
   TBranch        *b_nGenHiggsBoson;   //!
   TBranch        *b_GenHiggsBoson_pdgId;   //!
   TBranch        *b_GenHiggsBoson_pt;   //!
   TBranch        *b_GenHiggsBoson_eta;   //!
   TBranch        *b_GenHiggsBoson_phi;   //!
   TBranch        *b_GenHiggsBoson_mass;   //!
   TBranch        *b_GenHiggsBoson_charge;   //!
   TBranch        *b_GenHiggsBoson_status;   //!
   TBranch        *b_nGenLepFromTop;   //!
   TBranch        *b_GenLepFromTop_pdgId;   //!
   TBranch        *b_GenLepFromTop_pt;   //!
   TBranch        *b_GenLepFromTop_eta;   //!
   TBranch        *b_GenLepFromTop_phi;   //!
   TBranch        *b_GenLepFromTop_mass;   //!
   TBranch        *b_GenLepFromTop_charge;   //!
   TBranch        *b_GenLepFromTop_status;   //!
   TBranch        *b_nGenVbosons;   //!
   TBranch        *b_GenVbosons_pdgId;   //!
   TBranch        *b_GenVbosons_pt;   //!
   TBranch        *b_GenVbosons_eta;   //!
   TBranch        *b_GenVbosons_phi;   //!
   TBranch        *b_GenVbosons_mass;   //!
   TBranch        *b_GenVbosons_charge;   //!
   TBranch        *b_GenVbosons_status;   //!
   TBranch        *b_nFatjetCleanAK08prunedcalreg;   //!
   TBranch        *b_FatjetCleanAK08prunedcalreg_pt;   //!
   TBranch        *b_FatjetCleanAK08prunedcalreg_eta;   //!
   TBranch        *b_FatjetCleanAK08prunedcalreg_phi;   //!
   TBranch        *b_FatjetCleanAK08prunedcalreg_mass;   //!
   TBranch        *b_nFatjetAK08pruned;   //!
   TBranch        *b_FatjetAK08pruned_pt;   //!
   TBranch        *b_FatjetAK08pruned_eta;   //!
   TBranch        *b_FatjetAK08pruned_phi;   //!
   TBranch        *b_FatjetAK08pruned_mass;   //!
   TBranch        *b_nGenJet;   //!
   TBranch        *b_GenJet_charge;   //!
   TBranch        *b_GenJet_status;   //!
   TBranch        *b_GenJet_pdgId;   //!
   TBranch        *b_GenJet_pt;   //!
   TBranch        *b_GenJet_eta;   //!
   TBranch        *b_GenJet_phi;   //!
   TBranch        *b_GenJet_mass;   //!
   TBranch        *b_GenJet_numBHadrons;   //!
   TBranch        *b_GenJet_numCHadrons;   //!
   TBranch        *b_GenJet_numBHadronsFromTop;   //!
   TBranch        *b_GenJet_numCHadronsFromTop;   //!
   TBranch        *b_GenJet_numBHadronsAfterTop;   //!
   TBranch        *b_GenJet_numCHadronsAfterTop;   //!
   TBranch        *b_GenJet_wNuPt;   //!
   TBranch        *b_GenJet_wNuEta;   //!
   TBranch        *b_GenJet_wNuPhi;   //!
   TBranch        *b_GenJet_wNuM;   //!
   TBranch        *b_nGenTau;   //!
   TBranch        *b_GenTau_pdgId;   //!
   TBranch        *b_GenTau_pt;   //!
   TBranch        *b_GenTau_eta;   //!
   TBranch        *b_GenTau_phi;   //!
   TBranch        *b_GenTau_mass;   //!
   TBranch        *b_GenTau_charge;   //!
   TBranch        *b_GenTau_status;   //!
   TBranch        *b_nGenLep;   //!
   TBranch        *b_GenLep_pdgId;   //!
   TBranch        *b_GenLep_pt;   //!
   TBranch        *b_GenLep_eta;   //!
   TBranch        *b_GenLep_phi;   //!
   TBranch        *b_GenLep_mass;   //!
   TBranch        *b_GenLep_charge;   //!
   TBranch        *b_GenLep_status;   //!
   TBranch        *b_nGenBQuarkFromH;   //!
   TBranch        *b_GenBQuarkFromH_pdgId;   //!
   TBranch        *b_GenBQuarkFromH_pt;   //!
   TBranch        *b_GenBQuarkFromH_eta;   //!
   TBranch        *b_GenBQuarkFromH_phi;   //!
   TBranch        *b_GenBQuarkFromH_mass;   //!
   TBranch        *b_GenBQuarkFromH_charge;   //!
   TBranch        *b_GenBQuarkFromH_status;   //!
   TBranch        *b_nFatjetAK08prunedCal;   //!
   TBranch        *b_FatjetAK08prunedCal_pt;   //!
   TBranch        *b_FatjetAK08prunedCal_eta;   //!
   TBranch        *b_FatjetAK08prunedCal_phi;   //!
   TBranch        *b_FatjetAK08prunedCal_mass;   //!
   TBranch        *b_nincLeptons;   //!
   TBranch        *b_incLeptons_charge;   //!
   TBranch        *b_incLeptons_tightId;   //!
   TBranch        *b_incLeptons_eleCutIdCSA14_25ns_v1;   //!
   TBranch        *b_incLeptons_eleCutIdCSA14_50ns_v1;   //!
   TBranch        *b_incLeptons_eleCutIdSpring15_25ns_v1;   //!
   TBranch        *b_incLeptons_dxy;   //!
   TBranch        *b_incLeptons_dz;   //!
   TBranch        *b_incLeptons_edxy;   //!
   TBranch        *b_incLeptons_edz;   //!
   TBranch        *b_incLeptons_ip3d;   //!
   TBranch        *b_incLeptons_sip3d;   //!
   TBranch        *b_incLeptons_convVeto;   //!
   TBranch        *b_incLeptons_lostHits;   //!
   TBranch        *b_incLeptons_relIso03;   //!
   TBranch        *b_incLeptons_relIso04;   //!
   TBranch        *b_incLeptons_miniRelIso;   //!
   TBranch        *b_incLeptons_tightCharge;   //!
   TBranch        *b_incLeptons_mcMatchId;   //!
   TBranch        *b_incLeptons_mcMatchAny;   //!
   TBranch        *b_incLeptons_mcMatchTau;   //!
   TBranch        *b_incLeptons_mcPt;   //!
   TBranch        *b_incLeptons_mediumMuonId;   //!
   TBranch        *b_incLeptons_pdgId;   //!
   TBranch        *b_incLeptons_pt;   //!
   TBranch        *b_incLeptons_eta;   //!
   TBranch        *b_incLeptons_phi;   //!
   TBranch        *b_incLeptons_mass;   //!
   TBranch        *b_incLeptons_looseIdSusy;   //!
   TBranch        *b_incLeptons_looseIdPOG;   //!
   TBranch        *b_incLeptons_chargedHadRelIso03;   //!
   TBranch        *b_incLeptons_chargedHadRelIso04;   //!
   TBranch        *b_incLeptons_eleSieie;   //!
   TBranch        *b_incLeptons_e5x5;   //!
   TBranch        *b_incLeptons_e2x5Max;   //!
   TBranch        *b_incLeptons_e1x5;   //!
   TBranch        *b_incLeptons_isolTrkPt;   //!
   TBranch        *b_incLeptons_isolEmHadDepth1;   //!
   TBranch        *b_incLeptons_eleDEta;   //!
   TBranch        *b_incLeptons_eleDPhi;   //!
   TBranch        *b_incLeptons_eleHoE;   //!
   TBranch        *b_incLeptons_eleMissingHits;   //!
   TBranch        *b_incLeptons_eleChi2;   //!
   TBranch        *b_incLeptons_muonDX;   //!
   TBranch        *b_incLeptons_eleClusterEta;   //!
   TBranch        *b_incLeptons_eleClusterEnergy;   //!
   TBranch        *b_incLeptons_eleClusterDEta;   //!
   TBranch        *b_incLeptons_eleClusterDPhi;   //!
   TBranch        *b_incLeptons_nStations;   //!
   TBranch        *b_incLeptons_trkKink;   //!
   TBranch        *b_incLeptons_caloCompatibility;   //!
   TBranch        *b_incLeptons_globalTrackChi2;   //!
   TBranch        *b_incLeptons_nChamberHits;   //!
   TBranch        *b_incLeptons_isBarrelEle;   //!
   TBranch        *b_incLeptons_isEndCapEle;   //!
   TBranch        *b_incLeptons_isEcalDriven;   //!
   TBranch        *b_incLeptons_isMyGoodMuon;   //!
   TBranch        *b_incLeptons_isMyGoodMuon1;   //!
   TBranch        *b_incLeptons_isMyGoodElectron;   //!
   TBranch        *b_incLeptons_relPtError;   //!
   TBranch        *b_incLeptons_isPFMuon;   //!
   TBranch        *b_incLeptons_muon_dz;   //!
   TBranch        *b_incLeptons_isHighPtMuon;   //!
   TBranch        *b_incLeptons_muTrackIso;   //!
   TBranch        *b_incLeptons_isGlobalMuon;   //!
   TBranch        *b_incLeptons_isTrackerMuon;   //!
   TBranch        *b_incLeptons_pixelHits;   //!
   TBranch        *b_incLeptons_trackerLayers;   //!
   TBranch        *b_incLeptons_pixelLayers;   //!
   TBranch        *b_incLeptons_mvaTTH;   //!
   TBranch        *b_incLeptons_jetOverlapIdx;   //!
   TBranch        *b_incLeptons_jetPtRatio;   //!
   TBranch        *b_incLeptons_jetBTagCSV;   //!
   TBranch        *b_incLeptons_jetDR;   //!
   TBranch        *b_incLeptons_pfRelIso03;   //!
   TBranch        *b_incLeptons_pfRelIso04;   //!
   TBranch        *b_incLeptons_muonDB;   //!
   TBranch        *b_incLeptons_etaSc;   //!
   TBranch        *b_incLeptons_eleExpMissingInnerHits;   //!
   TBranch        *b_incLeptons_eleooEmooP;   //!
   TBranch        *b_incLeptons_muonTrackerLayers;   //!
   TBranch        *b_nGenTop;   //!
   TBranch        *b_GenTop_pdgId;   //!
   TBranch        *b_GenTop_pt;   //!
   TBranch        *b_GenTop_eta;   //!
   TBranch        *b_GenTop_phi;   //!
   TBranch        *b_GenTop_mass;   //!
   TBranch        *b_GenTop_charge;   //!
   TBranch        *b_GenTop_status;   //!
   TBranch        *b_nGenLepFromTau;   //!
   TBranch        *b_GenLepFromTau_pdgId;   //!
   TBranch        *b_GenLepFromTau_pt;   //!
   TBranch        *b_GenLepFromTau_eta;   //!
   TBranch        *b_GenLepFromTau_phi;   //!
   TBranch        *b_GenLepFromTau_mass;   //!
   TBranch        *b_GenLepFromTau_charge;   //!
   TBranch        *b_GenLepFromTau_status;   //!
   TBranch        *b_nGenNuFromTop;   //!
   TBranch        *b_GenNuFromTop_pdgId;   //!
   TBranch        *b_GenNuFromTop_pt;   //!
   TBranch        *b_GenNuFromTop_eta;   //!
   TBranch        *b_GenNuFromTop_phi;   //!
   TBranch        *b_GenNuFromTop_mass;   //!
   TBranch        *b_GenNuFromTop_charge;   //!
   TBranch        *b_GenNuFromTop_status;   //!
   TBranch        *b_nFatjetCleanAK08pruned;   //!
   TBranch        *b_FatjetCleanAK08pruned_pt;   //!
   TBranch        *b_FatjetCleanAK08pruned_eta;   //!
   TBranch        *b_FatjetCleanAK08pruned_phi;   //!
   TBranch        *b_FatjetCleanAK08pruned_mass;   //!
   TBranch        *b_nFatjetAK08ungroomed;   //!
   TBranch        *b_FatjetAK08ungroomed_pt;   //!
   TBranch        *b_FatjetAK08ungroomed_eta;   //!
   TBranch        *b_FatjetAK08ungroomed_phi;   //!
   TBranch        *b_FatjetAK08ungroomed_mass;   //!
   TBranch        *b_FatjetAK08ungroomed_chHEFrac;   //!
   TBranch        *b_FatjetAK08ungroomed_neHEFrac;   //!
   TBranch        *b_FatjetAK08ungroomed_chEmEFrac;   //!
   TBranch        *b_FatjetAK08ungroomed_neEmEFrac;   //!
   TBranch        *b_FatjetAK08ungroomed_chMult;   //!
   TBranch        *b_FatjetAK08ungroomed_neMult;   //!
   TBranch        *b_FatjetAK08ungroomed_muEFrac;   //!
   TBranch        *b_FatjetAK08ungroomed_tau1;   //!
   TBranch        *b_FatjetAK08ungroomed_tau2;   //!
   TBranch        *b_FatjetAK08ungroomed_tau3;   //!
   TBranch        *b_FatjetAK08ungroomed_msoftdrop;   //!
   TBranch        *b_FatjetAK08ungroomed_mpruned;   //!
   TBranch        *b_FatjetAK08ungroomed_mtrimmed;   //!
   TBranch        *b_FatjetAK08ungroomed_mfiltered;   //!
   TBranch        *b_FatjetAK08ungroomed_bbtag;   //!
   TBranch        *b_FatjetAK08ungroomed_PFLepton_ptrel;   //!
   TBranch        *b_FatjetAK08ungroomed_z_ratio;   //!
   TBranch        *b_FatjetAK08ungroomed_tau_dot;   //!
   TBranch        *b_FatjetAK08ungroomed_SV_mass_0;   //!
   TBranch        *b_FatjetAK08ungroomed_SV_EnergyRatio_0;   //!
   TBranch        *b_FatjetAK08ungroomed_SV_EnergyRatio_1;   //!
   TBranch        *b_FatjetAK08ungroomed_PFLepton_IP2D;   //!
   TBranch        *b_FatjetAK08ungroomed_tau_21;   //!
   TBranch        *b_FatjetAK08ungroomed_nSL;   //!
   TBranch        *b_FatjetAK08ungroomed_vertexNTracks;   //!
   TBranch        *b_nGenHiggsSisters;   //!
   TBranch        *b_GenHiggsSisters_pdgId;   //!
   TBranch        *b_GenHiggsSisters_pt;   //!
   TBranch        *b_GenHiggsSisters_eta;   //!
   TBranch        *b_GenHiggsSisters_phi;   //!
   TBranch        *b_GenHiggsSisters_mass;   //!
   TBranch        *b_GenHiggsSisters_charge;   //!
   TBranch        *b_GenHiggsSisters_status;   //!
   TBranch        *b_nSubjetCleanAK08pruned;   //!
   TBranch        *b_SubjetCleanAK08pruned_pt;   //!
   TBranch        *b_SubjetCleanAK08pruned_eta;   //!
   TBranch        *b_SubjetCleanAK08pruned_phi;   //!
   TBranch        *b_SubjetCleanAK08pruned_mass;   //!
   TBranch        *b_SubjetCleanAK08pruned_btag;   //!
   TBranch        *b_nGenNu;   //!
   TBranch        *b_GenNu_pdgId;   //!
   TBranch        *b_GenNu_pt;   //!
   TBranch        *b_GenNu_eta;   //!
   TBranch        *b_GenNu_phi;   //!
   TBranch        *b_GenNu_mass;   //!
   TBranch        *b_GenNu_charge;   //!
   TBranch        *b_GenNu_status;   //!
   TBranch        *b_nselLeptons;   //!
   TBranch        *b_selLeptons_charge;   //!
   TBranch        *b_selLeptons_tightId;   //!
   TBranch        *b_selLeptons_eleCutIdCSA14_25ns_v1;   //!
   TBranch        *b_selLeptons_eleCutIdCSA14_50ns_v1;   //!
   TBranch        *b_selLeptons_eleCutIdSpring15_25ns_v1;   //!
   TBranch        *b_selLeptons_dxy;   //!
   TBranch        *b_selLeptons_dz;   //!
   TBranch        *b_selLeptons_edxy;   //!
   TBranch        *b_selLeptons_edz;   //!
   TBranch        *b_selLeptons_ip3d;   //!
   TBranch        *b_selLeptons_sip3d;   //!
   TBranch        *b_selLeptons_convVeto;   //!
   TBranch        *b_selLeptons_lostHits;   //!
   TBranch        *b_selLeptons_relIso03;   //!
   TBranch        *b_selLeptons_relIso04;   //!
   TBranch        *b_selLeptons_miniRelIso;   //!
   TBranch        *b_selLeptons_tightCharge;   //!
   TBranch        *b_selLeptons_mcMatchId;   //!
   TBranch        *b_selLeptons_mcMatchAny;   //!
   TBranch        *b_selLeptons_mcMatchTau;   //!
   TBranch        *b_selLeptons_mcPt;   //!
   TBranch        *b_selLeptons_mediumMuonId;   //!
   TBranch        *b_selLeptons_pdgId;   //!
   TBranch        *b_selLeptons_pt;   //!
   TBranch        *b_selLeptons_eta;   //!
   TBranch        *b_selLeptons_phi;   //!
   TBranch        *b_selLeptons_mass;   //!
   TBranch        *b_selLeptons_looseIdSusy;   //!
   TBranch        *b_selLeptons_looseIdPOG;   //!
   TBranch        *b_selLeptons_chargedHadRelIso03;   //!
   TBranch        *b_selLeptons_chargedHadRelIso04;   //!
   TBranch        *b_selLeptons_eleSieie;   //!
   TBranch        *b_selLeptons_e5x5;   //!
   TBranch        *b_selLeptons_e2x5Max;   //!
   TBranch        *b_selLeptons_e1x5;   //!
   TBranch        *b_selLeptons_isolTrkPt;   //!
   TBranch        *b_selLeptons_isolEmHadDepth1;   //!
   TBranch        *b_selLeptons_eleDEta;   //!
   TBranch        *b_selLeptons_eleDPhi;   //!
   TBranch        *b_selLeptons_eleHoE;   //!
   TBranch        *b_selLeptons_eleMissingHits;   //!
   TBranch        *b_selLeptons_eleChi2;   //!
   TBranch        *b_selLeptons_muonDX;   //!
   TBranch        *b_selLeptons_eleClusterEta;   //!
   TBranch        *b_selLeptons_eleClusterEnergy;   //!
   TBranch        *b_selLeptons_eleClusterDEta;   //!
   TBranch        *b_selLeptons_eleClusterDPhi;   //!
   TBranch        *b_selLeptons_nStations;   //!
   TBranch        *b_selLeptons_trkKink;   //!
   TBranch        *b_selLeptons_caloCompatibility;   //!
   TBranch        *b_selLeptons_globalTrackChi2;   //!
   TBranch        *b_selLeptons_nChamberHits;   //!
   TBranch        *b_selLeptons_isBarrelEle;   //!
   TBranch        *b_selLeptons_isEndCapEle;   //!
   TBranch        *b_selLeptons_isEcalDriven;   //!
   TBranch        *b_selLeptons_isMyGoodMuon;   //!
   TBranch        *b_selLeptons_isMyGoodMuon1;   //!
   TBranch        *b_selLeptons_isMyGoodElectron;   //!
   TBranch        *b_selLeptons_relPtError;   //!
   TBranch        *b_selLeptons_isPFMuon;   //!
   TBranch        *b_selLeptons_muon_dz;   //!
   TBranch        *b_selLeptons_isHighPtMuon;   //!
   TBranch        *b_selLeptons_muTrackIso;   //!
   TBranch        *b_selLeptons_isGlobalMuon;   //!
   TBranch        *b_selLeptons_isTrackerMuon;   //!
   TBranch        *b_selLeptons_pixelHits;   //!
   TBranch        *b_selLeptons_trackerLayers;   //!
   TBranch        *b_selLeptons_pixelLayers;   //!
   TBranch        *b_selLeptons_mvaTTH;   //!
   TBranch        *b_selLeptons_jetOverlapIdx;   //!
   TBranch        *b_selLeptons_jetPtRatio;   //!
   TBranch        *b_selLeptons_jetBTagCSV;   //!
   TBranch        *b_selLeptons_jetDR;   //!
   TBranch        *b_selLeptons_pfRelIso03;   //!
   TBranch        *b_selLeptons_pfRelIso04;   //!
   TBranch        *b_selLeptons_muonDB;   //!
   TBranch        *b_selLeptons_etaSc;   //!
   TBranch        *b_selLeptons_eleExpMissingInnerHits;   //!
   TBranch        *b_selLeptons_eleooEmooP;   //!
   TBranch        *b_selLeptons_muonTrackerLayers;   //!
   TBranch        *b_nFatjetCleanAK08ungroomed;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_pt;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_eta;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_phi;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_mass;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_chHEFrac;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_neHEFrac;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_chEmEFrac;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_neEmEFrac;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_chMult;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_neMult;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_muEFrac;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_tau1;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_tau2;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_tau3;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_msoftdrop;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_mpruned;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_mtrimmed;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_mfiltered;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_bbtag;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_PFLepton_ptrel;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_z_ratio;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_tau_dot;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_SV_mass_0;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_SV_EnergyRatio_0;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_SV_EnergyRatio_1;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_PFLepton_IP2D;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_tau_21;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_nSL;   //!
   TBranch        *b_FatjetCleanAK08ungroomed_vertexNTracks;   //!
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
   TBranch        *b_nFatjetCleanAK08prunedcal;   //!
   TBranch        *b_FatjetCleanAK08prunedcal_pt;   //!
   TBranch        *b_FatjetCleanAK08prunedcal_eta;   //!
   TBranch        *b_FatjetCleanAK08prunedcal_phi;   //!
   TBranch        *b_FatjetCleanAK08prunedcal_mass;   //!
   TBranch        *b_nSubjetAK08pruned;   //!
   TBranch        *b_SubjetAK08pruned_pt;   //!
   TBranch        *b_SubjetAK08pruned_eta;   //!
   TBranch        *b_SubjetAK08pruned_phi;   //!
   TBranch        *b_SubjetAK08pruned_mass;   //!
   TBranch        *b_SubjetAK08pruned_btag;   //!
   TBranch        *b_nGenStatus2bHad;   //!
   TBranch        *b_GenStatus2bHad_pdgId;   //!
   TBranch        *b_GenStatus2bHad_pt;   //!
   TBranch        *b_GenStatus2bHad_eta;   //!
   TBranch        *b_GenStatus2bHad_phi;   //!
   TBranch        *b_GenStatus2bHad_mass;   //!
   TBranch        *b_GenStatus2bHad_charge;   //!
   TBranch        *b_GenStatus2bHad_status;   //!
   TBranch        *b_nGenHadTaus;   //!
   TBranch        *b_GenHadTaus_charge;   //!
   TBranch        *b_GenHadTaus_status;   //!
   TBranch        *b_GenHadTaus_pdgId;   //!
   TBranch        *b_GenHadTaus_pt;   //!
   TBranch        *b_GenHadTaus_eta;   //!
   TBranch        *b_GenHadTaus_phi;   //!
   TBranch        *b_GenHadTaus_mass;   //!
   TBranch        *b_GenHadTaus_decayMode;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_id;   //!
   TBranch        *b_Jet_puId;   //!
   TBranch        *b_Jet_btagCSV;   //!
   TBranch        *b_Jet_btagCMVA;   //!
   TBranch        *b_Jet_rawPt;   //!
   TBranch        *b_Jet_mcPt;   //!
   TBranch        *b_Jet_mcFlavour;   //!
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
   TBranch        *b_Jet_idxFirstTauMatch;   //!
   TBranch        *b_Jet_heppyFlavour;   //!
   TBranch        *b_Jet_hadronFlavour;   //!
   TBranch        *b_Jet_btagBDT;   //!
   TBranch        *b_Jet_btagProb;   //!
   TBranch        *b_Jet_btagBProb;   //!
   TBranch        *b_Jet_btagSoftEl;   //!
   TBranch        *b_Jet_btagSoftMu;   //!
   TBranch        *b_Jet_btagnew;   //!
   TBranch        *b_Jet_btagCSVV0;   //!
   TBranch        *b_Jet_isMyGoodLooseJet;   //!
   TBranch        *b_Jet_isMyGoodTightJet;   //!
   TBranch        *b_Jet_isMyGoodTightLepVetoJet;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_chMult;   //!
   TBranch        *b_Jet_neMult;   //!
   TBranch        *b_Jet_leadTrackPt;   //!
   TBranch        *b_Jet_mcEta;   //!
   TBranch        *b_Jet_mcPhi;   //!
   TBranch        *b_Jet_mcM;   //!
   TBranch        *b_Jet_leptonPdgId;   //!
   TBranch        *b_Jet_leptonPt;   //!
   TBranch        *b_Jet_leptonPtRel;   //!
   TBranch        *b_Jet_leptonPtRelInv;   //!
   TBranch        *b_Jet_leptonDeltaR;   //!
   TBranch        *b_Jet_leptonDeltaPhi;   //!
   TBranch        *b_Jet_leptonDeltaEta;   //!
   TBranch        *b_Jet_vtxMass;   //!
   TBranch        *b_Jet_vtxNtracks;   //!
   TBranch        *b_Jet_vtxPt;   //!
   TBranch        *b_Jet_vtx3DSig;   //!
   TBranch        *b_Jet_vtx3DVal;   //!
   TBranch        *b_Jet_vtxPosX;   //!
   TBranch        *b_Jet_vtxPosY;   //!
   TBranch        *b_Jet_vtxPosZ;   //!
   TBranch        *b_Jet_pullVectorPhi;   //!
   TBranch        *b_Jet_pullVectorMag;   //!
   TBranch        *b_Jet_qgl;   //!
   TBranch        *b_Jet_ptd;   //!
   TBranch        *b_Jet_axis2;   //!
   TBranch        *b_Jet_mult;   //!
   TBranch        *b_Jet_numberOfDaughters;   //!
   TBranch        *b_Jet_btagIdx;   //!
   TBranch        *b_Jet_mcIdx;   //!
   TBranch        *b_Jet_pt_reg;   //!
   TBranch        *b_Jet_pt_regVBF;   //!
   TBranch        *b_Jet_blike_VBF;   //!
   TBranch        *b_Jet_bTagWeightJESUp;   //!
   TBranch        *b_Jet_bTagWeightJESDown;   //!
   TBranch        *b_Jet_bTagWeightLFUp;   //!
   TBranch        *b_Jet_bTagWeightLFDown;   //!
   TBranch        *b_Jet_bTagWeightHFUp;   //!
   TBranch        *b_Jet_bTagWeightHFDown;   //!
   TBranch        *b_Jet_bTagWeightStats1Up;   //!
   TBranch        *b_Jet_bTagWeightStats1Down;   //!
   TBranch        *b_Jet_bTagWeightStats2Up;   //!
   TBranch        *b_Jet_bTagWeightStats2Down;   //!
   TBranch        *b_Jet_bTagWeightcErr1Up;   //!
   TBranch        *b_Jet_bTagWeightcErr1Down;   //!
   TBranch        *b_Jet_bTagWeightcErr2Up;   //!
   TBranch        *b_Jet_bTagWeightcErr2Down;   //!
   TBranch        *b_Jet_bTagWeight;   //!
   TBranch        *b_nFatjetGenAK08Pruned;   //!
   TBranch        *b_FatjetGenAK08Pruned_pt;   //!
   TBranch        *b_FatjetGenAK08Pruned_eta;   //!
   TBranch        *b_FatjetGenAK08Pruned_phi;   //!
   TBranch        *b_FatjetGenAK08Pruned_mass;   //!
   TBranch        *b_npileUpVertex_z;   //!
   TBranch        *b_pileUpVertex_z;   //!
   TBranch        *b_nmyGenJetAK08;   //!
   TBranch        *b_myGenJetAK08_pt;   //!
   TBranch        *b_myGenJetAK08_eta;   //!
   TBranch        *b_myGenJetAK08_phi;   //!
   TBranch        *b_myGenJetAK08_mass;   //!
   TBranch        *b_myGenJetAK08_tau1;   //!
   TBranch        *b_myGenJetAK08_tau2;   //!
   TBranch        *b_myGenJetAK08_tau3;   //!
   TBranch        *b_nprimaryVertices;   //!
   TBranch        *b_primaryVertices_x;   //!
   TBranch        *b_primaryVertices_y;   //!
   TBranch        *b_primaryVertices_z;   //!
   TBranch        *b_primaryVertices_isFake;   //!
   TBranch        *b_primaryVertices_ndof;   //!
   TBranch        *b_primaryVertices_Rho;   //!
   TBranch        *b_nFatjetAK08prunedReg;   //!
   TBranch        *b_FatjetAK08prunedReg_pt;   //!
   TBranch        *b_FatjetAK08prunedReg_eta;   //!
   TBranch        *b_FatjetAK08prunedReg_phi;   //!
   TBranch        *b_FatjetAK08prunedReg_mass;   //!
   TBranch        *b_nGenWZQuark;   //!
   TBranch        *b_GenWZQuark_pdgId;   //!
   TBranch        *b_GenWZQuark_pt;   //!
   TBranch        *b_GenWZQuark_eta;   //!
   TBranch        *b_GenWZQuark_phi;   //!
   TBranch        *b_GenWZQuark_mass;   //!
   TBranch        *b_GenWZQuark_charge;   //!
   TBranch        *b_GenWZQuark_status;   //!
   TBranch        *b_nGenLepFromWPrime;   //!
   TBranch        *b_GenLepFromWPrime_pdgId;   //!
   TBranch        *b_GenLepFromWPrime_pt;   //!
   TBranch        *b_GenLepFromWPrime_eta;   //!
   TBranch        *b_GenLepFromWPrime_phi;   //!
   TBranch        *b_GenLepFromWPrime_mass;   //!
   TBranch        *b_GenLepFromWPrime_charge;   //!
   TBranch        *b_GenLepFromWPrime_status;   //!
   TBranch        *b_nGenBQuarkFromTop;   //!
   TBranch        *b_GenBQuarkFromTop_pdgId;   //!
   TBranch        *b_GenBQuarkFromTop_pt;   //!
   TBranch        *b_GenBQuarkFromTop_eta;   //!
   TBranch        *b_GenBQuarkFromTop_phi;   //!
   TBranch        *b_GenBQuarkFromTop_mass;   //!
   TBranch        *b_GenBQuarkFromTop_charge;   //!
   TBranch        *b_GenBQuarkFromTop_status;   //!

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
      chain->Add("/pnfs/roma1.infn.it/data/cms/store/user/sgelli/BulkGravToWWToWlepWhad_narrow_M-1000_13TeV-madgraph/tt_14_BulkGravToWWToWlepWhad_narrow_M-1000_13TeV-madgraph_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/160315_142153/0000/tree_1.root/tree");
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
   fChain->SetBranchAddress("json", &json, &b_json);
   fChain->SetBranchAddress("nPU0", &nPU0, &b_nPU0);
   fChain->SetBranchAddress("nPVs", &nPVs, &b_nPVs);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("lheNj", &lheNj, &b_lheNj);
   fChain->SetBranchAddress("lheNb", &lheNb, &b_lheNb);
   fChain->SetBranchAddress("lheNc", &lheNc, &b_lheNc);
   fChain->SetBranchAddress("lheNg", &lheNg, &b_lheNg);
   fChain->SetBranchAddress("lheNl", &lheNl, &b_lheNl);
   fChain->SetBranchAddress("lheV_pt", &lheV_pt, &b_lheV_pt);
   fChain->SetBranchAddress("lheHT", &lheHT, &b_lheHT);
   fChain->SetBranchAddress("mhtJet30", &mhtJet30, &b_mhtJet30);
   fChain->SetBranchAddress("mhtPhiJet30", &mhtPhiJet30, &b_mhtPhiJet30);
   fChain->SetBranchAddress("htJet30", &htJet30, &b_htJet30);
   fChain->SetBranchAddress("met_rawpt", &met_rawpt, &b_met_rawpt);
   fChain->SetBranchAddress("metPuppi_pt", &metPuppi_pt, &b_metPuppi_pt);
   fChain->SetBranchAddress("metPuppi_phi", &metPuppi_phi, &b_metPuppi_phi);
   fChain->SetBranchAddress("metPuppi_rawpt", &metPuppi_rawpt, &b_metPuppi_rawpt);
   fChain->SetBranchAddress("metType1p2_pt", &metType1p2_pt, &b_metType1p2_pt);
   fChain->SetBranchAddress("metNoPU_pt", &metNoPU_pt, &b_metNoPU_pt);
   fChain->SetBranchAddress("metNoPU_phi", &metNoPU_phi, &b_metNoPU_phi);
   fChain->SetBranchAddress("tkMet_pt", &tkMet_pt, &b_tkMet_pt);
   fChain->SetBranchAddress("tkMet_phi", &tkMet_phi, &b_tkMet_phi);
   fChain->SetBranchAddress("tkMetPVchs_pt", &tkMetPVchs_pt, &b_tkMetPVchs_pt);
   fChain->SetBranchAddress("tkMetPVchs_phi", &tkMetPVchs_phi, &b_tkMetPVchs_phi);
   fChain->SetBranchAddress("Flag_hbheIsoFilter", &Flag_hbheIsoFilter, &b_Flag_hbheIsoFilter);
   fChain->SetBranchAddress("Flag_hbheFilterNew", &Flag_hbheFilterNew, &b_Flag_hbheFilterNew);
   fChain->SetBranchAddress("simPrimaryVertex_z", &simPrimaryVertex_z, &b_simPrimaryVertex_z);
   fChain->SetBranchAddress("genHiggsDecayMode", &genHiggsDecayMode, &b_genHiggsDecayMode);
   fChain->SetBranchAddress("bTagWeight_LFUp", &bTagWeight_LFUp, &b_bTagWeight_LFUp);
   fChain->SetBranchAddress("bTagWeight_Stats2Down", &bTagWeight_Stats2Down, &b_bTagWeight_Stats2Down);
   fChain->SetBranchAddress("bTagWeight_LFDown", &bTagWeight_LFDown, &b_bTagWeight_LFDown);
   fChain->SetBranchAddress("bTagWeight_HFUp", &bTagWeight_HFUp, &b_bTagWeight_HFUp);
   fChain->SetBranchAddress("bTagWeight_cErr1Down", &bTagWeight_cErr1Down, &b_bTagWeight_cErr1Down);
   fChain->SetBranchAddress("bTagWeight_JESDown", &bTagWeight_JESDown, &b_bTagWeight_JESDown);
   fChain->SetBranchAddress("bTagWeight_cErr1Up", &bTagWeight_cErr1Up, &b_bTagWeight_cErr1Up);
   fChain->SetBranchAddress("bTagWeight", &bTagWeight, &b_bTagWeight);
   fChain->SetBranchAddress("bTagWeight_HFDown", &bTagWeight_HFDown, &b_bTagWeight_HFDown);
   fChain->SetBranchAddress("bTagWeight_Stats2Up", &bTagWeight_Stats2Up, &b_bTagWeight_Stats2Up);
   fChain->SetBranchAddress("bTagWeight_cErr2Up", &bTagWeight_cErr2Up, &b_bTagWeight_cErr2Up);
   fChain->SetBranchAddress("bTagWeight_JESUp", &bTagWeight_JESUp, &b_bTagWeight_JESUp);
   fChain->SetBranchAddress("bTagWeight_Stats1Up", &bTagWeight_Stats1Up, &b_bTagWeight_Stats1Up);
   fChain->SetBranchAddress("bTagWeight_Stats1Down", &bTagWeight_Stats1Down, &b_bTagWeight_Stats1Down);
   fChain->SetBranchAddress("bTagWeight_cErr2Down", &bTagWeight_cErr2Down, &b_bTagWeight_cErr2Down);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, &b_Flag_trackingFailureFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("HLT_BIT_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_BTagCSV0p7_v", &HLT_BIT_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_BTagCSV0p7_v, &b_HLT_BIT_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_BTagCSV0p7_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_v", &HLT_BIT_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_v, &b_HLT_BIT_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET90_PFMHT90_IDLoose_v", &HLT_BIT_HLT_PFMET90_PFMHT90_IDLoose_v, &b_HLT_BIT_HLT_PFMET90_PFMHT90_IDLoose_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET100_PFMHT100_IDLoose_v", &HLT_BIT_HLT_PFMET100_PFMHT100_IDLoose_v, &b_HLT_BIT_HLT_PFMET100_PFMHT100_IDLoose_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET110_PFMHT110_IDLoose_v", &HLT_BIT_HLT_PFMET110_PFMHT110_IDLoose_v, &b_HLT_BIT_HLT_PFMET110_PFMHT110_IDLoose_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET120_PFMHT120_IDLoose_v", &HLT_BIT_HLT_PFMET120_PFMHT120_IDLoose_v, &b_HLT_BIT_HLT_PFMET120_PFMHT120_IDLoose_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET170_NoiseCleaned_v", &HLT_BIT_HLT_PFMET170_NoiseCleaned_v, &b_HLT_BIT_HLT_PFMET170_NoiseCleaned_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_DiCentralPFJet70_PFMET120_NoiseCleaned_v", &HLT_BIT_HLT_DiCentralPFJet70_PFMET120_NoiseCleaned_v, &b_HLT_BIT_HLT_DiCentralPFJet70_PFMET120_NoiseCleaned_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT350_PFMET120_NoiseCleaned_v", &HLT_BIT_HLT_PFHT350_PFMET120_NoiseCleaned_v, &b_HLT_BIT_HLT_PFHT350_PFMET120_NoiseCleaned_v);
   fChain->SetBranchAddress("HLT_ZnnHbbAll", &HLT_ZnnHbbAll, &b_HLT_ZnnHbbAll);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v", &HLT_BIT_HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v, &b_HLT_BIT_HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT450_SixJet40_PFBTagCSV_v", &HLT_BIT_HLT_PFHT450_SixJet40_PFBTagCSV_v, &b_HLT_BIT_HLT_PFHT450_SixJet40_PFBTagCSV_v);
   fChain->SetBranchAddress("HLT_ttHhardonicLowLumi", &HLT_ttHhardonicLowLumi, &b_HLT_ttHhardonicLowLumi);
   fChain->SetBranchAddress("HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v", &HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v, &b_HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v);
   fChain->SetBranchAddress("HLT_WtaunHbbHighLumi", &HLT_WtaunHbbHighLumi, &b_HLT_WtaunHbbHighLumi);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v", &HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v, &b_HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v", &HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v, &b_HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v);
   fChain->SetBranchAddress("HLT_VBFHbbLowLumi", &HLT_VBFHbbLowLumi, &b_HLT_VBFHbbLowLumi);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT750_4Jet_v", &HLT_BIT_HLT_PFHT750_4Jet_v, &b_HLT_BIT_HLT_PFHT750_4Jet_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT900_v", &HLT_BIT_HLT_PFHT900_v, &b_HLT_BIT_HLT_PFHT900_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet40_v", &HLT_BIT_HLT_PFJet40_v, &b_HLT_BIT_HLT_PFJet40_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet60_v", &HLT_BIT_HLT_PFJet60_v, &b_HLT_BIT_HLT_PFJet60_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet80_v", &HLT_BIT_HLT_PFJet80_v, &b_HLT_BIT_HLT_PFJet80_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet140_v", &HLT_BIT_HLT_PFJet140_v, &b_HLT_BIT_HLT_PFJet140_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet200_v", &HLT_BIT_HLT_PFJet200_v, &b_HLT_BIT_HLT_PFJet200_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet260_v", &HLT_BIT_HLT_PFJet260_v, &b_HLT_BIT_HLT_PFJet260_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet320_v", &HLT_BIT_HLT_PFJet320_v, &b_HLT_BIT_HLT_PFJet320_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet400_v", &HLT_BIT_HLT_PFJet400_v, &b_HLT_BIT_HLT_PFJet400_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFJet450_v", &HLT_BIT_HLT_PFJet450_v, &b_HLT_BIT_HLT_PFJet450_v);
   fChain->SetBranchAddress("HLT_hadronic", &HLT_hadronic, &b_HLT_hadronic);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadJet45_TripleCSV0p5_v", &HLT_BIT_HLT_QuadJet45_TripleCSV0p5_v, &b_HLT_BIT_HLT_QuadJet45_TripleCSV0p5_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadJet45_DoubleCSV0p5_v", &HLT_BIT_HLT_QuadJet45_DoubleCSV0p5_v, &b_HLT_BIT_HLT_QuadJet45_DoubleCSV0p5_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_TripleCSV0p5_v", &HLT_BIT_HLT_DoubleJet90_Double30_TripleCSV0p5_v, &b_HLT_BIT_HLT_DoubleJet90_Double30_TripleCSV0p5_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_DoubleCSV0p5_v", &HLT_BIT_HLT_DoubleJet90_Double30_DoubleCSV0p5_v, &b_HLT_BIT_HLT_DoubleJet90_Double30_DoubleCSV0p5_v);
   fChain->SetBranchAddress("HLT_HH4bAll", &HLT_HH4bAll, &b_HLT_HH4bAll);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu24_eta2p1_v", &HLT_BIT_HLT_IsoMu24_eta2p1_v, &b_HLT_BIT_HLT_IsoMu24_eta2p1_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v", &HLT_BIT_HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v, &b_HLT_BIT_HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu24_eta2p1_v", &HLT_BIT_HLT_Mu24_eta2p1_v, &b_HLT_BIT_HLT_Mu24_eta2p1_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_TkMu24_eta2p1_v", &HLT_BIT_HLT_TkMu24_eta2p1_v, &b_HLT_BIT_HLT_TkMu24_eta2p1_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu24_v", &HLT_BIT_HLT_Mu24_v, &b_HLT_BIT_HLT_Mu24_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu27_v", &HLT_BIT_HLT_IsoMu27_v, &b_HLT_BIT_HLT_IsoMu27_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoTkMu27_v", &HLT_BIT_HLT_IsoTkMu27_v, &b_HLT_BIT_HLT_IsoTkMu27_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_TkMu27_v", &HLT_BIT_HLT_TkMu27_v, &b_HLT_BIT_HLT_TkMu27_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu27_v", &HLT_BIT_HLT_Mu27_v, &b_HLT_BIT_HLT_Mu27_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu20_eta2p1_v", &HLT_BIT_HLT_IsoMu20_eta2p1_v, &b_HLT_BIT_HLT_IsoMu20_eta2p1_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v", &HLT_BIT_HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v, &b_HLT_BIT_HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu20_v", &HLT_BIT_HLT_Mu20_v, &b_HLT_BIT_HLT_Mu20_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_TkMu20_v", &HLT_BIT_HLT_TkMu20_v, &b_HLT_BIT_HLT_TkMu20_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu20_v", &HLT_BIT_HLT_IsoMu20_v, &b_HLT_BIT_HLT_IsoMu20_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoTkMu20_v", &HLT_BIT_HLT_IsoTkMu20_v, &b_HLT_BIT_HLT_IsoTkMu20_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu40_eta2p1_PFJet200_PFJet50_v", &HLT_BIT_HLT_Mu40_eta2p1_PFJet200_PFJet50_v, &b_HLT_BIT_HLT_Mu40_eta2p1_PFJet200_PFJet50_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu16_eta2p1_CaloMET30_v", &HLT_BIT_HLT_IsoMu16_eta2p1_CaloMET30_v, &b_HLT_BIT_HLT_IsoMu16_eta2p1_CaloMET30_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu16_eta2p1_CaloMET30_v", &HLT_BIT_HLT_Mu16_eta2p1_CaloMET30_v, &b_HLT_BIT_HLT_Mu16_eta2p1_CaloMET30_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMET120_NoiseCleaned_Mu5_v", &HLT_BIT_HLT_PFMET120_NoiseCleaned_Mu5_v, &b_HLT_BIT_HLT_PFMET120_NoiseCleaned_Mu5_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_MonoCentralPFJet80_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v", &HLT_BIT_HLT_MonoCentralPFJet80_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v, &b_HLT_BIT_HLT_MonoCentralPFJet80_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_MonoCentralPFJet80_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v", &HLT_BIT_HLT_MonoCentralPFJet80_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v, &b_HLT_BIT_HLT_MonoCentralPFJet80_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v", &HLT_BIT_HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v, &b_HLT_BIT_HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v", &HLT_BIT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v, &b_HLT_BIT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v);
   fChain->SetBranchAddress("HLT_WmnHbbAll", &HLT_WmnHbbAll, &b_HLT_WmnHbbAll);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_v", &HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_v, &b_HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v", &HLT_BIT_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v, &b_HLT_BIT_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v);
   fChain->SetBranchAddress("HLT_WenHbbHighLumi", &HLT_WenHbbHighLumi, &b_HLT_WenHbbHighLumi);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_DoubleEle24_22_eta2p1_WP75_Gsf_v", &HLT_BIT_HLT_DoubleEle24_22_eta2p1_WP75_Gsf_v, &b_HLT_BIT_HLT_DoubleEle24_22_eta2p1_WP75_Gsf_v);
   fChain->SetBranchAddress("HLT_ZeeHbbLowLumi", &HLT_ZeeHbbLowLumi, &b_HLT_ZeeHbbLowLumi);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v", &HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v, &b_HLT_BIT_HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele27_WP85_Gsf_v", &HLT_BIT_HLT_Ele27_WP85_Gsf_v, &b_HLT_BIT_HLT_Ele27_WP85_Gsf_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_v", &HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_v, &b_HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v", &HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v, &b_HLT_BIT_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele105_CaloIdVT_GsfTrkIdT_v", &HLT_BIT_HLT_Ele105_CaloIdVT_GsfTrkIdT_v, &b_HLT_BIT_HLT_Ele105_CaloIdVT_GsfTrkIdT_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v", &HLT_BIT_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v, &b_HLT_BIT_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v);
   fChain->SetBranchAddress("HLT_WenHbbAll", &HLT_WenHbbAll, &b_HLT_WenHbbAll);
   fChain->SetBranchAddress("HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v", &HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v, &b_HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_v", &HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_v, &b_HLT_BIT_HLT_LooseIsoPFTau50_Trk30_eta2p1_v);
   fChain->SetBranchAddress("HLT_WtaunHbbAll", &HLT_WtaunHbbAll, &b_HLT_WtaunHbbAll);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", &HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v", &HLT_BIT_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_BIT_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v", &HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v", &HLT_BIT_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_BIT_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_ZeeHbbAll", &HLT_ZeeHbbAll, &b_HLT_ZeeHbbAll);
   fChain->SetBranchAddress("HLT_WtaunHbbLowLumi", &HLT_WtaunHbbLowLumi, &b_HLT_WtaunHbbLowLumi);
   fChain->SetBranchAddress("HLT_WmnHbbLowLumi", &HLT_WmnHbbLowLumi, &b_HLT_WmnHbbLowLumi);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v", &HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_ttHleptonic", &HLT_ttHleptonic, &b_HLT_ttHleptonic);
   fChain->SetBranchAddress("HLT_ZeeHbbHighLumi", &HLT_ZeeHbbHighLumi, &b_HLT_ZeeHbbHighLumi);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT450_SixJet40_v", &HLT_BIT_HLT_PFHT450_SixJet40_v, &b_HLT_BIT_HLT_PFHT450_SixJet40_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT400_SixJet30_v", &HLT_BIT_HLT_PFHT400_SixJet30_v, &b_HLT_BIT_HLT_PFHT400_SixJet30_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT350_v", &HLT_BIT_HLT_PFHT350_v, &b_HLT_BIT_HLT_PFHT350_v);
   fChain->SetBranchAddress("HLT_ttHhardonicAll", &HLT_ttHhardonicAll, &b_HLT_ttHhardonicAll);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu50_v", &HLT_BIT_HLT_Mu50_v, &b_HLT_BIT_HLT_Mu50_v);
   fChain->SetBranchAddress("HLT_HLT_Mu50", &HLT_HLT_Mu50, &b_HLT_HLT_Mu50);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu45_eta2p1_v", &HLT_BIT_HLT_Mu45_eta2p1_v, &b_HLT_BIT_HLT_Mu45_eta2p1_v);
   fChain->SetBranchAddress("HLT_HLT_Mu45_eta2p1", &HLT_HLT_Mu45_eta2p1, &b_HLT_HLT_Mu45_eta2p1);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", &HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v);
   fChain->SetBranchAddress("HLT_ZmmHbbLowLumi", &HLT_ZmmHbbLowLumi, &b_HLT_ZmmHbbLowLumi);
   fChain->SetBranchAddress("HLT_WmnHbbHighLumi", &HLT_WmnHbbHighLumi, &b_HLT_WmnHbbHighLumi);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TkMu8_DZ_v", &HLT_BIT_HLT_Mu17_TkMu8_DZ_v, &b_HLT_BIT_HLT_Mu17_TkMu8_DZ_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v", &HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v, &b_HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_DoubleIsoMu17_eta2p1_v", &HLT_BIT_HLT_DoubleIsoMu17_eta2p1_v, &b_HLT_BIT_HLT_DoubleIsoMu17_eta2p1_v);
   fChain->SetBranchAddress("HLT_ZmmHbbAll", &HLT_ZmmHbbAll, &b_HLT_ZmmHbbAll);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v", &HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v, &b_HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v);
   fChain->SetBranchAddress("HLT_HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_HLT_Ele115_CaloIdVT_GsfTrkIdT, &b_HLT_HLT_Ele115_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_WenHbbLowLumi", &HLT_WenHbbLowLumi, &b_HLT_WenHbbLowLumi);
   fChain->SetBranchAddress("HLT_ZnnHbbHighLumi", &HLT_ZnnHbbHighLumi, &b_HLT_ZnnHbbHighLumi);
   fChain->SetBranchAddress("HLT_HH4bLowLumi", &HLT_HH4bLowLumi, &b_HLT_HH4bLowLumi);
   fChain->SetBranchAddress("HLT_ZmmHbbHighLumi", &HLT_ZmmHbbHighLumi, &b_HLT_ZmmHbbHighLumi);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v", &HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v, &b_HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v", &HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v, &b_HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_VBF_v", &HLT_BIT_HLT_QuadPFJet_VBF_v, &b_HLT_BIT_HLT_QuadPFJet_VBF_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_L1_TripleJet_VBF_v", &HLT_BIT_HLT_L1_TripleJet_VBF_v, &b_HLT_BIT_HLT_L1_TripleJet_VBF_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v", &HLT_BIT_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v, &b_HLT_BIT_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v", &HLT_BIT_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v, &b_HLT_BIT_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v);
   fChain->SetBranchAddress("HLT_VBFHbbAll", &HLT_VBFHbbAll, &b_HLT_VBFHbbAll);
   fChain->SetBranchAddress("HLT_VBFHbbHighLumi", &HLT_VBFHbbHighLumi, &b_HLT_VBFHbbHighLumi);
   fChain->SetBranchAddress("HLT_ttHhardonicHighLumi", &HLT_ttHhardonicHighLumi, &b_HLT_ttHhardonicHighLumi);
   fChain->SetBranchAddress("HLT_ZnnHbbLowLumi", &HLT_ZnnHbbLowLumi, &b_HLT_ZnnHbbLowLumi);
   fChain->SetBranchAddress("HLT_HH4bHighLumi", &HLT_HH4bHighLumi, &b_HLT_HH4bHighLumi);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_mass", &met_mass, &b_met_mass);
   fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   fChain->SetBranchAddress("met_genPt", &met_genPt, &b_met_genPt);
   fChain->SetBranchAddress("met_genPhi", &met_genPhi, &b_met_genPhi);
   fChain->SetBranchAddress("met_genEta", &met_genEta, &b_met_genEta);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnUp_pt", &met_shifted_UnclusteredEnUp_pt, &b_met_shifted_UnclusteredEnUp_pt);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnUp_phi", &met_shifted_UnclusteredEnUp_phi, &b_met_shifted_UnclusteredEnUp_phi);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnUp_sumEt", &met_shifted_UnclusteredEnUp_sumEt, &b_met_shifted_UnclusteredEnUp_sumEt);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnDown_pt", &met_shifted_UnclusteredEnDown_pt, &b_met_shifted_UnclusteredEnDown_pt);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnDown_phi", &met_shifted_UnclusteredEnDown_phi, &b_met_shifted_UnclusteredEnDown_phi);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnDown_sumEt", &met_shifted_UnclusteredEnDown_sumEt, &b_met_shifted_UnclusteredEnDown_sumEt);
   fChain->SetBranchAddress("met_shifted_JetResUp_pt", &met_shifted_JetResUp_pt, &b_met_shifted_JetResUp_pt);
   fChain->SetBranchAddress("met_shifted_JetResUp_phi", &met_shifted_JetResUp_phi, &b_met_shifted_JetResUp_phi);
   fChain->SetBranchAddress("met_shifted_JetResUp_sumEt", &met_shifted_JetResUp_sumEt, &b_met_shifted_JetResUp_sumEt);
   fChain->SetBranchAddress("met_shifted_JetResDown_pt", &met_shifted_JetResDown_pt, &b_met_shifted_JetResDown_pt);
   fChain->SetBranchAddress("met_shifted_JetResDown_phi", &met_shifted_JetResDown_phi, &b_met_shifted_JetResDown_phi);
   fChain->SetBranchAddress("met_shifted_JetResDown_sumEt", &met_shifted_JetResDown_sumEt, &b_met_shifted_JetResDown_sumEt);
   fChain->SetBranchAddress("met_shifted_JetEnUp_pt", &met_shifted_JetEnUp_pt, &b_met_shifted_JetEnUp_pt);
   fChain->SetBranchAddress("met_shifted_JetEnUp_phi", &met_shifted_JetEnUp_phi, &b_met_shifted_JetEnUp_phi);
   fChain->SetBranchAddress("met_shifted_JetEnUp_sumEt", &met_shifted_JetEnUp_sumEt, &b_met_shifted_JetEnUp_sumEt);
   fChain->SetBranchAddress("met_shifted_JetEnDown_pt", &met_shifted_JetEnDown_pt, &b_met_shifted_JetEnDown_pt);
   fChain->SetBranchAddress("met_shifted_JetEnDown_phi", &met_shifted_JetEnDown_phi, &b_met_shifted_JetEnDown_phi);
   fChain->SetBranchAddress("met_shifted_JetEnDown_sumEt", &met_shifted_JetEnDown_sumEt, &b_met_shifted_JetEnDown_sumEt);
   fChain->SetBranchAddress("met_shifted_MuonEnUp_pt", &met_shifted_MuonEnUp_pt, &b_met_shifted_MuonEnUp_pt);
   fChain->SetBranchAddress("met_shifted_MuonEnUp_phi", &met_shifted_MuonEnUp_phi, &b_met_shifted_MuonEnUp_phi);
   fChain->SetBranchAddress("met_shifted_MuonEnUp_sumEt", &met_shifted_MuonEnUp_sumEt, &b_met_shifted_MuonEnUp_sumEt);
   fChain->SetBranchAddress("met_shifted_MuonEnDown_pt", &met_shifted_MuonEnDown_pt, &b_met_shifted_MuonEnDown_pt);
   fChain->SetBranchAddress("met_shifted_MuonEnDown_phi", &met_shifted_MuonEnDown_phi, &b_met_shifted_MuonEnDown_phi);
   fChain->SetBranchAddress("met_shifted_MuonEnDown_sumEt", &met_shifted_MuonEnDown_sumEt, &b_met_shifted_MuonEnDown_sumEt);
   fChain->SetBranchAddress("met_shifted_ElectronEnUp_pt", &met_shifted_ElectronEnUp_pt, &b_met_shifted_ElectronEnUp_pt);
   fChain->SetBranchAddress("met_shifted_ElectronEnUp_phi", &met_shifted_ElectronEnUp_phi, &b_met_shifted_ElectronEnUp_phi);
   fChain->SetBranchAddress("met_shifted_ElectronEnUp_sumEt", &met_shifted_ElectronEnUp_sumEt, &b_met_shifted_ElectronEnUp_sumEt);
   fChain->SetBranchAddress("met_shifted_ElectronEnDown_pt", &met_shifted_ElectronEnDown_pt, &b_met_shifted_ElectronEnDown_pt);
   fChain->SetBranchAddress("met_shifted_ElectronEnDown_phi", &met_shifted_ElectronEnDown_phi, &b_met_shifted_ElectronEnDown_phi);
   fChain->SetBranchAddress("met_shifted_ElectronEnDown_sumEt", &met_shifted_ElectronEnDown_sumEt, &b_met_shifted_ElectronEnDown_sumEt);
   fChain->SetBranchAddress("met_shifted_TauEnUp_pt", &met_shifted_TauEnUp_pt, &b_met_shifted_TauEnUp_pt);
   fChain->SetBranchAddress("met_shifted_TauEnUp_phi", &met_shifted_TauEnUp_phi, &b_met_shifted_TauEnUp_phi);
   fChain->SetBranchAddress("met_shifted_TauEnUp_sumEt", &met_shifted_TauEnUp_sumEt, &b_met_shifted_TauEnUp_sumEt);
   fChain->SetBranchAddress("met_shifted_TauEnDown_pt", &met_shifted_TauEnDown_pt, &b_met_shifted_TauEnDown_pt);
   fChain->SetBranchAddress("met_shifted_TauEnDown_phi", &met_shifted_TauEnDown_phi, &b_met_shifted_TauEnDown_phi);
   fChain->SetBranchAddress("met_shifted_TauEnDown_sumEt", &met_shifted_TauEnDown_sumEt, &b_met_shifted_TauEnDown_sumEt);
   fChain->SetBranchAddress("nGenBQuarkFromHafterISR", &nGenBQuarkFromHafterISR, &b_nGenBQuarkFromHafterISR);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_pdgId", &GenBQuarkFromHafterISR_pdgId, &b_GenBQuarkFromHafterISR_pdgId);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_pt", &GenBQuarkFromHafterISR_pt, &b_GenBQuarkFromHafterISR_pt);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_eta", &GenBQuarkFromHafterISR_eta, &b_GenBQuarkFromHafterISR_eta);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_phi", &GenBQuarkFromHafterISR_phi, &b_GenBQuarkFromHafterISR_phi);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_mass", &GenBQuarkFromHafterISR_mass, &b_GenBQuarkFromHafterISR_mass);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_charge", &GenBQuarkFromHafterISR_charge, &b_GenBQuarkFromHafterISR_charge);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_status", &GenBQuarkFromHafterISR_status, &b_GenBQuarkFromHafterISR_status);
   fChain->SetBranchAddress("npileUpVertex_ptHat", &npileUpVertex_ptHat, &b_npileUpVertex_ptHat);
   fChain->SetBranchAddress("pileUpVertex_ptHat", pileUpVertex_ptHat, &b_pileUpVertex_ptHat);
   fChain->SetBranchAddress("nGenHiggsBoson", &nGenHiggsBoson, &b_nGenHiggsBoson);
   fChain->SetBranchAddress("GenHiggsBoson_pdgId", &GenHiggsBoson_pdgId, &b_GenHiggsBoson_pdgId);
   fChain->SetBranchAddress("GenHiggsBoson_pt", &GenHiggsBoson_pt, &b_GenHiggsBoson_pt);
   fChain->SetBranchAddress("GenHiggsBoson_eta", &GenHiggsBoson_eta, &b_GenHiggsBoson_eta);
   fChain->SetBranchAddress("GenHiggsBoson_phi", &GenHiggsBoson_phi, &b_GenHiggsBoson_phi);
   fChain->SetBranchAddress("GenHiggsBoson_mass", &GenHiggsBoson_mass, &b_GenHiggsBoson_mass);
   fChain->SetBranchAddress("GenHiggsBoson_charge", &GenHiggsBoson_charge, &b_GenHiggsBoson_charge);
   fChain->SetBranchAddress("GenHiggsBoson_status", &GenHiggsBoson_status, &b_GenHiggsBoson_status);
   fChain->SetBranchAddress("nGenLepFromTop", &nGenLepFromTop, &b_nGenLepFromTop);
   fChain->SetBranchAddress("GenLepFromTop_pdgId", &GenLepFromTop_pdgId, &b_GenLepFromTop_pdgId);
   fChain->SetBranchAddress("GenLepFromTop_pt", &GenLepFromTop_pt, &b_GenLepFromTop_pt);
   fChain->SetBranchAddress("GenLepFromTop_eta", &GenLepFromTop_eta, &b_GenLepFromTop_eta);
   fChain->SetBranchAddress("GenLepFromTop_phi", &GenLepFromTop_phi, &b_GenLepFromTop_phi);
   fChain->SetBranchAddress("GenLepFromTop_mass", &GenLepFromTop_mass, &b_GenLepFromTop_mass);
   fChain->SetBranchAddress("GenLepFromTop_charge", &GenLepFromTop_charge, &b_GenLepFromTop_charge);
   fChain->SetBranchAddress("GenLepFromTop_status", &GenLepFromTop_status, &b_GenLepFromTop_status);
   fChain->SetBranchAddress("nGenVbosons", &nGenVbosons, &b_nGenVbosons);
   fChain->SetBranchAddress("GenVbosons_pdgId", GenVbosons_pdgId, &b_GenVbosons_pdgId);
   fChain->SetBranchAddress("GenVbosons_pt", GenVbosons_pt, &b_GenVbosons_pt);
   fChain->SetBranchAddress("GenVbosons_eta", GenVbosons_eta, &b_GenVbosons_eta);
   fChain->SetBranchAddress("GenVbosons_phi", GenVbosons_phi, &b_GenVbosons_phi);
   fChain->SetBranchAddress("GenVbosons_mass", GenVbosons_mass, &b_GenVbosons_mass);
   fChain->SetBranchAddress("GenVbosons_charge", GenVbosons_charge, &b_GenVbosons_charge);
   fChain->SetBranchAddress("GenVbosons_status", GenVbosons_status, &b_GenVbosons_status);
   fChain->SetBranchAddress("nFatjetCleanAK08prunedcalreg", &nFatjetCleanAK08prunedcalreg, &b_nFatjetCleanAK08prunedcalreg);
   fChain->SetBranchAddress("FatjetCleanAK08prunedcalreg_pt", FatjetCleanAK08prunedcalreg_pt, &b_FatjetCleanAK08prunedcalreg_pt);
   fChain->SetBranchAddress("FatjetCleanAK08prunedcalreg_eta", FatjetCleanAK08prunedcalreg_eta, &b_FatjetCleanAK08prunedcalreg_eta);
   fChain->SetBranchAddress("FatjetCleanAK08prunedcalreg_phi", FatjetCleanAK08prunedcalreg_phi, &b_FatjetCleanAK08prunedcalreg_phi);
   fChain->SetBranchAddress("FatjetCleanAK08prunedcalreg_mass", FatjetCleanAK08prunedcalreg_mass, &b_FatjetCleanAK08prunedcalreg_mass);
   fChain->SetBranchAddress("nFatjetAK08pruned", &nFatjetAK08pruned, &b_nFatjetAK08pruned);
   fChain->SetBranchAddress("FatjetAK08pruned_pt", FatjetAK08pruned_pt, &b_FatjetAK08pruned_pt);
   fChain->SetBranchAddress("FatjetAK08pruned_eta", FatjetAK08pruned_eta, &b_FatjetAK08pruned_eta);
   fChain->SetBranchAddress("FatjetAK08pruned_phi", FatjetAK08pruned_phi, &b_FatjetAK08pruned_phi);
   fChain->SetBranchAddress("FatjetAK08pruned_mass", FatjetAK08pruned_mass, &b_FatjetAK08pruned_mass);
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   fChain->SetBranchAddress("GenJet_charge", GenJet_charge, &b_GenJet_charge);
   fChain->SetBranchAddress("GenJet_status", GenJet_status, &b_GenJet_status);
   fChain->SetBranchAddress("GenJet_pdgId", GenJet_pdgId, &b_GenJet_pdgId);
   fChain->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
   fChain->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
   fChain->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
   fChain->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
   fChain->SetBranchAddress("GenJet_numBHadrons", GenJet_numBHadrons, &b_GenJet_numBHadrons);
   fChain->SetBranchAddress("GenJet_numCHadrons", GenJet_numCHadrons, &b_GenJet_numCHadrons);
   fChain->SetBranchAddress("GenJet_numBHadronsFromTop", GenJet_numBHadronsFromTop, &b_GenJet_numBHadronsFromTop);
   fChain->SetBranchAddress("GenJet_numCHadronsFromTop", GenJet_numCHadronsFromTop, &b_GenJet_numCHadronsFromTop);
   fChain->SetBranchAddress("GenJet_numBHadronsAfterTop", GenJet_numBHadronsAfterTop, &b_GenJet_numBHadronsAfterTop);
   fChain->SetBranchAddress("GenJet_numCHadronsAfterTop", GenJet_numCHadronsAfterTop, &b_GenJet_numCHadronsAfterTop);
   fChain->SetBranchAddress("GenJet_wNuPt", GenJet_wNuPt, &b_GenJet_wNuPt);
   fChain->SetBranchAddress("GenJet_wNuEta", GenJet_wNuEta, &b_GenJet_wNuEta);
   fChain->SetBranchAddress("GenJet_wNuPhi", GenJet_wNuPhi, &b_GenJet_wNuPhi);
   fChain->SetBranchAddress("GenJet_wNuM", GenJet_wNuM, &b_GenJet_wNuM);
   fChain->SetBranchAddress("nGenTau", &nGenTau, &b_nGenTau);
   fChain->SetBranchAddress("GenTau_pdgId", GenTau_pdgId, &b_GenTau_pdgId);
   fChain->SetBranchAddress("GenTau_pt", GenTau_pt, &b_GenTau_pt);
   fChain->SetBranchAddress("GenTau_eta", GenTau_eta, &b_GenTau_eta);
   fChain->SetBranchAddress("GenTau_phi", GenTau_phi, &b_GenTau_phi);
   fChain->SetBranchAddress("GenTau_mass", GenTau_mass, &b_GenTau_mass);
   fChain->SetBranchAddress("GenTau_charge", GenTau_charge, &b_GenTau_charge);
   fChain->SetBranchAddress("GenTau_status", GenTau_status, &b_GenTau_status);
   fChain->SetBranchAddress("nGenLep", &nGenLep, &b_nGenLep);
   fChain->SetBranchAddress("GenLep_pdgId", GenLep_pdgId, &b_GenLep_pdgId);
   fChain->SetBranchAddress("GenLep_pt", GenLep_pt, &b_GenLep_pt);
   fChain->SetBranchAddress("GenLep_eta", GenLep_eta, &b_GenLep_eta);
   fChain->SetBranchAddress("GenLep_phi", GenLep_phi, &b_GenLep_phi);
   fChain->SetBranchAddress("GenLep_mass", GenLep_mass, &b_GenLep_mass);
   fChain->SetBranchAddress("GenLep_charge", GenLep_charge, &b_GenLep_charge);
   fChain->SetBranchAddress("GenLep_status", GenLep_status, &b_GenLep_status);
   fChain->SetBranchAddress("nGenBQuarkFromH", &nGenBQuarkFromH, &b_nGenBQuarkFromH);
   fChain->SetBranchAddress("GenBQuarkFromH_pdgId", &GenBQuarkFromH_pdgId, &b_GenBQuarkFromH_pdgId);
   fChain->SetBranchAddress("GenBQuarkFromH_pt", &GenBQuarkFromH_pt, &b_GenBQuarkFromH_pt);
   fChain->SetBranchAddress("GenBQuarkFromH_eta", &GenBQuarkFromH_eta, &b_GenBQuarkFromH_eta);
   fChain->SetBranchAddress("GenBQuarkFromH_phi", &GenBQuarkFromH_phi, &b_GenBQuarkFromH_phi);
   fChain->SetBranchAddress("GenBQuarkFromH_mass", &GenBQuarkFromH_mass, &b_GenBQuarkFromH_mass);
   fChain->SetBranchAddress("GenBQuarkFromH_charge", &GenBQuarkFromH_charge, &b_GenBQuarkFromH_charge);
   fChain->SetBranchAddress("GenBQuarkFromH_status", &GenBQuarkFromH_status, &b_GenBQuarkFromH_status);
   fChain->SetBranchAddress("nFatjetAK08prunedCal", &nFatjetAK08prunedCal, &b_nFatjetAK08prunedCal);
   fChain->SetBranchAddress("FatjetAK08prunedCal_pt", FatjetAK08prunedCal_pt, &b_FatjetAK08prunedCal_pt);
   fChain->SetBranchAddress("FatjetAK08prunedCal_eta", FatjetAK08prunedCal_eta, &b_FatjetAK08prunedCal_eta);
   fChain->SetBranchAddress("FatjetAK08prunedCal_phi", FatjetAK08prunedCal_phi, &b_FatjetAK08prunedCal_phi);
   fChain->SetBranchAddress("FatjetAK08prunedCal_mass", FatjetAK08prunedCal_mass, &b_FatjetAK08prunedCal_mass);
   fChain->SetBranchAddress("nincLeptons", &nincLeptons, &b_nincLeptons);
   fChain->SetBranchAddress("incLeptons_charge", incLeptons_charge, &b_incLeptons_charge);
   fChain->SetBranchAddress("incLeptons_tightId", incLeptons_tightId, &b_incLeptons_tightId);
   fChain->SetBranchAddress("incLeptons_eleCutIdCSA14_25ns_v1", incLeptons_eleCutIdCSA14_25ns_v1, &b_incLeptons_eleCutIdCSA14_25ns_v1);
   fChain->SetBranchAddress("incLeptons_eleCutIdCSA14_50ns_v1", incLeptons_eleCutIdCSA14_50ns_v1, &b_incLeptons_eleCutIdCSA14_50ns_v1);
   fChain->SetBranchAddress("incLeptons_eleCutIdSpring15_25ns_v1", incLeptons_eleCutIdSpring15_25ns_v1, &b_incLeptons_eleCutIdSpring15_25ns_v1);
   fChain->SetBranchAddress("incLeptons_dxy", incLeptons_dxy, &b_incLeptons_dxy);
   fChain->SetBranchAddress("incLeptons_dz", incLeptons_dz, &b_incLeptons_dz);
   fChain->SetBranchAddress("incLeptons_edxy", incLeptons_edxy, &b_incLeptons_edxy);
   fChain->SetBranchAddress("incLeptons_edz", incLeptons_edz, &b_incLeptons_edz);
   fChain->SetBranchAddress("incLeptons_ip3d", incLeptons_ip3d, &b_incLeptons_ip3d);
   fChain->SetBranchAddress("incLeptons_sip3d", incLeptons_sip3d, &b_incLeptons_sip3d);
   fChain->SetBranchAddress("incLeptons_convVeto", incLeptons_convVeto, &b_incLeptons_convVeto);
   fChain->SetBranchAddress("incLeptons_lostHits", incLeptons_lostHits, &b_incLeptons_lostHits);
   fChain->SetBranchAddress("incLeptons_relIso03", incLeptons_relIso03, &b_incLeptons_relIso03);
   fChain->SetBranchAddress("incLeptons_relIso04", incLeptons_relIso04, &b_incLeptons_relIso04);
   fChain->SetBranchAddress("incLeptons_miniRelIso", incLeptons_miniRelIso, &b_incLeptons_miniRelIso);
   fChain->SetBranchAddress("incLeptons_tightCharge", incLeptons_tightCharge, &b_incLeptons_tightCharge);
   fChain->SetBranchAddress("incLeptons_mcMatchId", incLeptons_mcMatchId, &b_incLeptons_mcMatchId);
   fChain->SetBranchAddress("incLeptons_mcMatchAny", incLeptons_mcMatchAny, &b_incLeptons_mcMatchAny);
   fChain->SetBranchAddress("incLeptons_mcMatchTau", incLeptons_mcMatchTau, &b_incLeptons_mcMatchTau);
   fChain->SetBranchAddress("incLeptons_mcPt", incLeptons_mcPt, &b_incLeptons_mcPt);
   fChain->SetBranchAddress("incLeptons_mediumMuonId", incLeptons_mediumMuonId, &b_incLeptons_mediumMuonId);
   fChain->SetBranchAddress("incLeptons_pdgId", incLeptons_pdgId, &b_incLeptons_pdgId);
   fChain->SetBranchAddress("incLeptons_pt", incLeptons_pt, &b_incLeptons_pt);
   fChain->SetBranchAddress("incLeptons_eta", incLeptons_eta, &b_incLeptons_eta);
   fChain->SetBranchAddress("incLeptons_phi", incLeptons_phi, &b_incLeptons_phi);
   fChain->SetBranchAddress("incLeptons_mass", incLeptons_mass, &b_incLeptons_mass);
   fChain->SetBranchAddress("incLeptons_looseIdSusy", incLeptons_looseIdSusy, &b_incLeptons_looseIdSusy);
   fChain->SetBranchAddress("incLeptons_looseIdPOG", incLeptons_looseIdPOG, &b_incLeptons_looseIdPOG);
   fChain->SetBranchAddress("incLeptons_chargedHadRelIso03", incLeptons_chargedHadRelIso03, &b_incLeptons_chargedHadRelIso03);
   fChain->SetBranchAddress("incLeptons_chargedHadRelIso04", incLeptons_chargedHadRelIso04, &b_incLeptons_chargedHadRelIso04);
   fChain->SetBranchAddress("incLeptons_eleSieie", incLeptons_eleSieie, &b_incLeptons_eleSieie);
   fChain->SetBranchAddress("incLeptons_e5x5", incLeptons_e5x5, &b_incLeptons_e5x5);
   fChain->SetBranchAddress("incLeptons_e2x5Max", incLeptons_e2x5Max, &b_incLeptons_e2x5Max);
   fChain->SetBranchAddress("incLeptons_e1x5", incLeptons_e1x5, &b_incLeptons_e1x5);
   fChain->SetBranchAddress("incLeptons_isolTrkPt", incLeptons_isolTrkPt, &b_incLeptons_isolTrkPt);
   fChain->SetBranchAddress("incLeptons_isolEmHadDepth1", incLeptons_isolEmHadDepth1, &b_incLeptons_isolEmHadDepth1);
   fChain->SetBranchAddress("incLeptons_eleDEta", incLeptons_eleDEta, &b_incLeptons_eleDEta);
   fChain->SetBranchAddress("incLeptons_eleDPhi", incLeptons_eleDPhi, &b_incLeptons_eleDPhi);
   fChain->SetBranchAddress("incLeptons_eleHoE", incLeptons_eleHoE, &b_incLeptons_eleHoE);
   fChain->SetBranchAddress("incLeptons_eleMissingHits", incLeptons_eleMissingHits, &b_incLeptons_eleMissingHits);
   fChain->SetBranchAddress("incLeptons_eleChi2", incLeptons_eleChi2, &b_incLeptons_eleChi2);
   fChain->SetBranchAddress("incLeptons_muonDX", incLeptons_muonDX, &b_incLeptons_muonDX);
   fChain->SetBranchAddress("incLeptons_eleClusterEta", incLeptons_eleClusterEta, &b_incLeptons_eleClusterEta);
   fChain->SetBranchAddress("incLeptons_eleClusterEnergy", incLeptons_eleClusterEnergy, &b_incLeptons_eleClusterEnergy);
   fChain->SetBranchAddress("incLeptons_eleClusterDEta", incLeptons_eleClusterDEta, &b_incLeptons_eleClusterDEta);
   fChain->SetBranchAddress("incLeptons_eleClusterDPhi", incLeptons_eleClusterDPhi, &b_incLeptons_eleClusterDPhi);
   fChain->SetBranchAddress("incLeptons_nStations", incLeptons_nStations, &b_incLeptons_nStations);
   fChain->SetBranchAddress("incLeptons_trkKink", incLeptons_trkKink, &b_incLeptons_trkKink);
   fChain->SetBranchAddress("incLeptons_caloCompatibility", incLeptons_caloCompatibility, &b_incLeptons_caloCompatibility);
   fChain->SetBranchAddress("incLeptons_globalTrackChi2", incLeptons_globalTrackChi2, &b_incLeptons_globalTrackChi2);
   fChain->SetBranchAddress("incLeptons_nChamberHits", incLeptons_nChamberHits, &b_incLeptons_nChamberHits);
   fChain->SetBranchAddress("incLeptons_isBarrelEle", incLeptons_isBarrelEle, &b_incLeptons_isBarrelEle);
   fChain->SetBranchAddress("incLeptons_isEndCapEle", incLeptons_isEndCapEle, &b_incLeptons_isEndCapEle);
   fChain->SetBranchAddress("incLeptons_isEcalDriven", incLeptons_isEcalDriven, &b_incLeptons_isEcalDriven);
   fChain->SetBranchAddress("incLeptons_isMyGoodMuon", incLeptons_isMyGoodMuon, &b_incLeptons_isMyGoodMuon);
   fChain->SetBranchAddress("incLeptons_isMyGoodMuon1", incLeptons_isMyGoodMuon1, &b_incLeptons_isMyGoodMuon1);
   fChain->SetBranchAddress("incLeptons_isMyGoodElectron", incLeptons_isMyGoodElectron, &b_incLeptons_isMyGoodElectron);
   fChain->SetBranchAddress("incLeptons_relPtError", incLeptons_relPtError, &b_incLeptons_relPtError);
   fChain->SetBranchAddress("incLeptons_isPFMuon", incLeptons_isPFMuon, &b_incLeptons_isPFMuon);
   fChain->SetBranchAddress("incLeptons_muon_dz", incLeptons_muon_dz, &b_incLeptons_muon_dz);
   fChain->SetBranchAddress("incLeptons_isHighPtMuon", incLeptons_isHighPtMuon, &b_incLeptons_isHighPtMuon);
   fChain->SetBranchAddress("incLeptons_muTrackIso", incLeptons_muTrackIso, &b_incLeptons_muTrackIso);
   fChain->SetBranchAddress("incLeptons_isGlobalMuon", incLeptons_isGlobalMuon, &b_incLeptons_isGlobalMuon);
   fChain->SetBranchAddress("incLeptons_isTrackerMuon", incLeptons_isTrackerMuon, &b_incLeptons_isTrackerMuon);
   fChain->SetBranchAddress("incLeptons_pixelHits", incLeptons_pixelHits, &b_incLeptons_pixelHits);
   fChain->SetBranchAddress("incLeptons_trackerLayers", incLeptons_trackerLayers, &b_incLeptons_trackerLayers);
   fChain->SetBranchAddress("incLeptons_pixelLayers", incLeptons_pixelLayers, &b_incLeptons_pixelLayers);
   fChain->SetBranchAddress("incLeptons_mvaTTH", incLeptons_mvaTTH, &b_incLeptons_mvaTTH);
   fChain->SetBranchAddress("incLeptons_jetOverlapIdx", incLeptons_jetOverlapIdx, &b_incLeptons_jetOverlapIdx);
   fChain->SetBranchAddress("incLeptons_jetPtRatio", incLeptons_jetPtRatio, &b_incLeptons_jetPtRatio);
   fChain->SetBranchAddress("incLeptons_jetBTagCSV", incLeptons_jetBTagCSV, &b_incLeptons_jetBTagCSV);
   fChain->SetBranchAddress("incLeptons_jetDR", incLeptons_jetDR, &b_incLeptons_jetDR);
   fChain->SetBranchAddress("incLeptons_pfRelIso03", incLeptons_pfRelIso03, &b_incLeptons_pfRelIso03);
   fChain->SetBranchAddress("incLeptons_pfRelIso04", incLeptons_pfRelIso04, &b_incLeptons_pfRelIso04);
   fChain->SetBranchAddress("incLeptons_muonDB", incLeptons_muonDB, &b_incLeptons_muonDB);
   fChain->SetBranchAddress("incLeptons_etaSc", incLeptons_etaSc, &b_incLeptons_etaSc);
   fChain->SetBranchAddress("incLeptons_eleExpMissingInnerHits", incLeptons_eleExpMissingInnerHits, &b_incLeptons_eleExpMissingInnerHits);
   fChain->SetBranchAddress("incLeptons_eleooEmooP", incLeptons_eleooEmooP, &b_incLeptons_eleooEmooP);
   fChain->SetBranchAddress("incLeptons_muonTrackerLayers", incLeptons_muonTrackerLayers, &b_incLeptons_muonTrackerLayers);
   fChain->SetBranchAddress("nGenTop", &nGenTop, &b_nGenTop);
   fChain->SetBranchAddress("GenTop_pdgId", &GenTop_pdgId, &b_GenTop_pdgId);
   fChain->SetBranchAddress("GenTop_pt", &GenTop_pt, &b_GenTop_pt);
   fChain->SetBranchAddress("GenTop_eta", &GenTop_eta, &b_GenTop_eta);
   fChain->SetBranchAddress("GenTop_phi", &GenTop_phi, &b_GenTop_phi);
   fChain->SetBranchAddress("GenTop_mass", &GenTop_mass, &b_GenTop_mass);
   fChain->SetBranchAddress("GenTop_charge", &GenTop_charge, &b_GenTop_charge);
   fChain->SetBranchAddress("GenTop_status", &GenTop_status, &b_GenTop_status);
   fChain->SetBranchAddress("nGenLepFromTau", &nGenLepFromTau, &b_nGenLepFromTau);
   fChain->SetBranchAddress("GenLepFromTau_pdgId", &GenLepFromTau_pdgId, &b_GenLepFromTau_pdgId);
   fChain->SetBranchAddress("GenLepFromTau_pt", &GenLepFromTau_pt, &b_GenLepFromTau_pt);
   fChain->SetBranchAddress("GenLepFromTau_eta", &GenLepFromTau_eta, &b_GenLepFromTau_eta);
   fChain->SetBranchAddress("GenLepFromTau_phi", &GenLepFromTau_phi, &b_GenLepFromTau_phi);
   fChain->SetBranchAddress("GenLepFromTau_mass", &GenLepFromTau_mass, &b_GenLepFromTau_mass);
   fChain->SetBranchAddress("GenLepFromTau_charge", &GenLepFromTau_charge, &b_GenLepFromTau_charge);
   fChain->SetBranchAddress("GenLepFromTau_status", &GenLepFromTau_status, &b_GenLepFromTau_status);
   fChain->SetBranchAddress("nGenNuFromTop", &nGenNuFromTop, &b_nGenNuFromTop);
   fChain->SetBranchAddress("GenNuFromTop_pdgId", &GenNuFromTop_pdgId, &b_GenNuFromTop_pdgId);
   fChain->SetBranchAddress("GenNuFromTop_pt", &GenNuFromTop_pt, &b_GenNuFromTop_pt);
   fChain->SetBranchAddress("GenNuFromTop_eta", &GenNuFromTop_eta, &b_GenNuFromTop_eta);
   fChain->SetBranchAddress("GenNuFromTop_phi", &GenNuFromTop_phi, &b_GenNuFromTop_phi);
   fChain->SetBranchAddress("GenNuFromTop_mass", &GenNuFromTop_mass, &b_GenNuFromTop_mass);
   fChain->SetBranchAddress("GenNuFromTop_charge", &GenNuFromTop_charge, &b_GenNuFromTop_charge);
   fChain->SetBranchAddress("GenNuFromTop_status", &GenNuFromTop_status, &b_GenNuFromTop_status);
   fChain->SetBranchAddress("nFatjetCleanAK08pruned", &nFatjetCleanAK08pruned, &b_nFatjetCleanAK08pruned);
   fChain->SetBranchAddress("FatjetCleanAK08pruned_pt", FatjetCleanAK08pruned_pt, &b_FatjetCleanAK08pruned_pt);
   fChain->SetBranchAddress("FatjetCleanAK08pruned_eta", FatjetCleanAK08pruned_eta, &b_FatjetCleanAK08pruned_eta);
   fChain->SetBranchAddress("FatjetCleanAK08pruned_phi", FatjetCleanAK08pruned_phi, &b_FatjetCleanAK08pruned_phi);
   fChain->SetBranchAddress("FatjetCleanAK08pruned_mass", FatjetCleanAK08pruned_mass, &b_FatjetCleanAK08pruned_mass);
   fChain->SetBranchAddress("nFatjetAK08ungroomed", &nFatjetAK08ungroomed, &b_nFatjetAK08ungroomed);
   fChain->SetBranchAddress("FatjetAK08ungroomed_pt", FatjetAK08ungroomed_pt, &b_FatjetAK08ungroomed_pt);
   fChain->SetBranchAddress("FatjetAK08ungroomed_eta", FatjetAK08ungroomed_eta, &b_FatjetAK08ungroomed_eta);
   fChain->SetBranchAddress("FatjetAK08ungroomed_phi", FatjetAK08ungroomed_phi, &b_FatjetAK08ungroomed_phi);
   fChain->SetBranchAddress("FatjetAK08ungroomed_mass", FatjetAK08ungroomed_mass, &b_FatjetAK08ungroomed_mass);
   fChain->SetBranchAddress("FatjetAK08ungroomed_chHEFrac", FatjetAK08ungroomed_chHEFrac, &b_FatjetAK08ungroomed_chHEFrac);
   fChain->SetBranchAddress("FatjetAK08ungroomed_neHEFrac", FatjetAK08ungroomed_neHEFrac, &b_FatjetAK08ungroomed_neHEFrac);
   fChain->SetBranchAddress("FatjetAK08ungroomed_chEmEFrac", FatjetAK08ungroomed_chEmEFrac, &b_FatjetAK08ungroomed_chEmEFrac);
   fChain->SetBranchAddress("FatjetAK08ungroomed_neEmEFrac", FatjetAK08ungroomed_neEmEFrac, &b_FatjetAK08ungroomed_neEmEFrac);
   fChain->SetBranchAddress("FatjetAK08ungroomed_chMult", FatjetAK08ungroomed_chMult, &b_FatjetAK08ungroomed_chMult);
   fChain->SetBranchAddress("FatjetAK08ungroomed_neMult", FatjetAK08ungroomed_neMult, &b_FatjetAK08ungroomed_neMult);
   fChain->SetBranchAddress("FatjetAK08ungroomed_muEFrac", FatjetAK08ungroomed_muEFrac, &b_FatjetAK08ungroomed_muEFrac);
   fChain->SetBranchAddress("FatjetAK08ungroomed_tau1", FatjetAK08ungroomed_tau1, &b_FatjetAK08ungroomed_tau1);
   fChain->SetBranchAddress("FatjetAK08ungroomed_tau2", FatjetAK08ungroomed_tau2, &b_FatjetAK08ungroomed_tau2);
   fChain->SetBranchAddress("FatjetAK08ungroomed_tau3", FatjetAK08ungroomed_tau3, &b_FatjetAK08ungroomed_tau3);
   fChain->SetBranchAddress("FatjetAK08ungroomed_msoftdrop", FatjetAK08ungroomed_msoftdrop, &b_FatjetAK08ungroomed_msoftdrop);
   fChain->SetBranchAddress("FatjetAK08ungroomed_mpruned", FatjetAK08ungroomed_mpruned, &b_FatjetAK08ungroomed_mpruned);
   fChain->SetBranchAddress("FatjetAK08ungroomed_mtrimmed", FatjetAK08ungroomed_mtrimmed, &b_FatjetAK08ungroomed_mtrimmed);
   fChain->SetBranchAddress("FatjetAK08ungroomed_mfiltered", FatjetAK08ungroomed_mfiltered, &b_FatjetAK08ungroomed_mfiltered);
   fChain->SetBranchAddress("FatjetAK08ungroomed_bbtag", FatjetAK08ungroomed_bbtag, &b_FatjetAK08ungroomed_bbtag);
   fChain->SetBranchAddress("FatjetAK08ungroomed_PFLepton_ptrel", FatjetAK08ungroomed_PFLepton_ptrel, &b_FatjetAK08ungroomed_PFLepton_ptrel);
   fChain->SetBranchAddress("FatjetAK08ungroomed_z_ratio", FatjetAK08ungroomed_z_ratio, &b_FatjetAK08ungroomed_z_ratio);
   fChain->SetBranchAddress("FatjetAK08ungroomed_tau_dot", FatjetAK08ungroomed_tau_dot, &b_FatjetAK08ungroomed_tau_dot);
   fChain->SetBranchAddress("FatjetAK08ungroomed_SV_mass_0", FatjetAK08ungroomed_SV_mass_0, &b_FatjetAK08ungroomed_SV_mass_0);
   fChain->SetBranchAddress("FatjetAK08ungroomed_SV_EnergyRatio_0", FatjetAK08ungroomed_SV_EnergyRatio_0, &b_FatjetAK08ungroomed_SV_EnergyRatio_0);
   fChain->SetBranchAddress("FatjetAK08ungroomed_SV_EnergyRatio_1", FatjetAK08ungroomed_SV_EnergyRatio_1, &b_FatjetAK08ungroomed_SV_EnergyRatio_1);
   fChain->SetBranchAddress("FatjetAK08ungroomed_PFLepton_IP2D", FatjetAK08ungroomed_PFLepton_IP2D, &b_FatjetAK08ungroomed_PFLepton_IP2D);
   fChain->SetBranchAddress("FatjetAK08ungroomed_tau_21", FatjetAK08ungroomed_tau_21, &b_FatjetAK08ungroomed_tau_21);
   fChain->SetBranchAddress("FatjetAK08ungroomed_nSL", FatjetAK08ungroomed_nSL, &b_FatjetAK08ungroomed_nSL);
   fChain->SetBranchAddress("FatjetAK08ungroomed_vertexNTracks", FatjetAK08ungroomed_vertexNTracks, &b_FatjetAK08ungroomed_vertexNTracks);
   fChain->SetBranchAddress("nGenHiggsSisters", &nGenHiggsSisters, &b_nGenHiggsSisters);
   fChain->SetBranchAddress("GenHiggsSisters_pdgId", &GenHiggsSisters_pdgId, &b_GenHiggsSisters_pdgId);
   fChain->SetBranchAddress("GenHiggsSisters_pt", &GenHiggsSisters_pt, &b_GenHiggsSisters_pt);
   fChain->SetBranchAddress("GenHiggsSisters_eta", &GenHiggsSisters_eta, &b_GenHiggsSisters_eta);
   fChain->SetBranchAddress("GenHiggsSisters_phi", &GenHiggsSisters_phi, &b_GenHiggsSisters_phi);
   fChain->SetBranchAddress("GenHiggsSisters_mass", &GenHiggsSisters_mass, &b_GenHiggsSisters_mass);
   fChain->SetBranchAddress("GenHiggsSisters_charge", &GenHiggsSisters_charge, &b_GenHiggsSisters_charge);
   fChain->SetBranchAddress("GenHiggsSisters_status", &GenHiggsSisters_status, &b_GenHiggsSisters_status);
   fChain->SetBranchAddress("nSubjetCleanAK08pruned", &nSubjetCleanAK08pruned, &b_nSubjetCleanAK08pruned);
   fChain->SetBranchAddress("SubjetCleanAK08pruned_pt", SubjetCleanAK08pruned_pt, &b_SubjetCleanAK08pruned_pt);
   fChain->SetBranchAddress("SubjetCleanAK08pruned_eta", SubjetCleanAK08pruned_eta, &b_SubjetCleanAK08pruned_eta);
   fChain->SetBranchAddress("SubjetCleanAK08pruned_phi", SubjetCleanAK08pruned_phi, &b_SubjetCleanAK08pruned_phi);
   fChain->SetBranchAddress("SubjetCleanAK08pruned_mass", SubjetCleanAK08pruned_mass, &b_SubjetCleanAK08pruned_mass);
   fChain->SetBranchAddress("SubjetCleanAK08pruned_btag", SubjetCleanAK08pruned_btag, &b_SubjetCleanAK08pruned_btag);
   fChain->SetBranchAddress("nGenNu", &nGenNu, &b_nGenNu);
   fChain->SetBranchAddress("GenNu_pdgId", GenNu_pdgId, &b_GenNu_pdgId);
   fChain->SetBranchAddress("GenNu_pt", GenNu_pt, &b_GenNu_pt);
   fChain->SetBranchAddress("GenNu_eta", GenNu_eta, &b_GenNu_eta);
   fChain->SetBranchAddress("GenNu_phi", GenNu_phi, &b_GenNu_phi);
   fChain->SetBranchAddress("GenNu_mass", GenNu_mass, &b_GenNu_mass);
   fChain->SetBranchAddress("GenNu_charge", GenNu_charge, &b_GenNu_charge);
   fChain->SetBranchAddress("GenNu_status", GenNu_status, &b_GenNu_status);
   fChain->SetBranchAddress("nselLeptons", &nselLeptons, &b_nselLeptons);
   fChain->SetBranchAddress("selLeptons_charge", selLeptons_charge, &b_selLeptons_charge);
   fChain->SetBranchAddress("selLeptons_tightId", selLeptons_tightId, &b_selLeptons_tightId);
   fChain->SetBranchAddress("selLeptons_eleCutIdCSA14_25ns_v1", selLeptons_eleCutIdCSA14_25ns_v1, &b_selLeptons_eleCutIdCSA14_25ns_v1);
   fChain->SetBranchAddress("selLeptons_eleCutIdCSA14_50ns_v1", selLeptons_eleCutIdCSA14_50ns_v1, &b_selLeptons_eleCutIdCSA14_50ns_v1);
   fChain->SetBranchAddress("selLeptons_eleCutIdSpring15_25ns_v1", selLeptons_eleCutIdSpring15_25ns_v1, &b_selLeptons_eleCutIdSpring15_25ns_v1);
   fChain->SetBranchAddress("selLeptons_dxy", selLeptons_dxy, &b_selLeptons_dxy);
   fChain->SetBranchAddress("selLeptons_dz", selLeptons_dz, &b_selLeptons_dz);
   fChain->SetBranchAddress("selLeptons_edxy", selLeptons_edxy, &b_selLeptons_edxy);
   fChain->SetBranchAddress("selLeptons_edz", selLeptons_edz, &b_selLeptons_edz);
   fChain->SetBranchAddress("selLeptons_ip3d", selLeptons_ip3d, &b_selLeptons_ip3d);
   fChain->SetBranchAddress("selLeptons_sip3d", selLeptons_sip3d, &b_selLeptons_sip3d);
   fChain->SetBranchAddress("selLeptons_convVeto", selLeptons_convVeto, &b_selLeptons_convVeto);
   fChain->SetBranchAddress("selLeptons_lostHits", selLeptons_lostHits, &b_selLeptons_lostHits);
   fChain->SetBranchAddress("selLeptons_relIso03", selLeptons_relIso03, &b_selLeptons_relIso03);
   fChain->SetBranchAddress("selLeptons_relIso04", selLeptons_relIso04, &b_selLeptons_relIso04);
   fChain->SetBranchAddress("selLeptons_miniRelIso", selLeptons_miniRelIso, &b_selLeptons_miniRelIso);
   fChain->SetBranchAddress("selLeptons_tightCharge", selLeptons_tightCharge, &b_selLeptons_tightCharge);
   fChain->SetBranchAddress("selLeptons_mcMatchId", selLeptons_mcMatchId, &b_selLeptons_mcMatchId);
   fChain->SetBranchAddress("selLeptons_mcMatchAny", selLeptons_mcMatchAny, &b_selLeptons_mcMatchAny);
   fChain->SetBranchAddress("selLeptons_mcMatchTau", selLeptons_mcMatchTau, &b_selLeptons_mcMatchTau);
   fChain->SetBranchAddress("selLeptons_mcPt", selLeptons_mcPt, &b_selLeptons_mcPt);
   fChain->SetBranchAddress("selLeptons_mediumMuonId", selLeptons_mediumMuonId, &b_selLeptons_mediumMuonId);
   fChain->SetBranchAddress("selLeptons_pdgId", selLeptons_pdgId, &b_selLeptons_pdgId);
   fChain->SetBranchAddress("selLeptons_pt", selLeptons_pt, &b_selLeptons_pt);
   fChain->SetBranchAddress("selLeptons_eta", selLeptons_eta, &b_selLeptons_eta);
   fChain->SetBranchAddress("selLeptons_phi", selLeptons_phi, &b_selLeptons_phi);
   fChain->SetBranchAddress("selLeptons_mass", selLeptons_mass, &b_selLeptons_mass);
   fChain->SetBranchAddress("selLeptons_looseIdSusy", selLeptons_looseIdSusy, &b_selLeptons_looseIdSusy);
   fChain->SetBranchAddress("selLeptons_looseIdPOG", selLeptons_looseIdPOG, &b_selLeptons_looseIdPOG);
   fChain->SetBranchAddress("selLeptons_chargedHadRelIso03", selLeptons_chargedHadRelIso03, &b_selLeptons_chargedHadRelIso03);
   fChain->SetBranchAddress("selLeptons_chargedHadRelIso04", selLeptons_chargedHadRelIso04, &b_selLeptons_chargedHadRelIso04);
   fChain->SetBranchAddress("selLeptons_eleSieie", selLeptons_eleSieie, &b_selLeptons_eleSieie);
   fChain->SetBranchAddress("selLeptons_e5x5", selLeptons_e5x5, &b_selLeptons_e5x5);
   fChain->SetBranchAddress("selLeptons_e2x5Max", selLeptons_e2x5Max, &b_selLeptons_e2x5Max);
   fChain->SetBranchAddress("selLeptons_e1x5", selLeptons_e1x5, &b_selLeptons_e1x5);
   fChain->SetBranchAddress("selLeptons_isolTrkPt", selLeptons_isolTrkPt, &b_selLeptons_isolTrkPt);
   fChain->SetBranchAddress("selLeptons_isolEmHadDepth1", selLeptons_isolEmHadDepth1, &b_selLeptons_isolEmHadDepth1);
   fChain->SetBranchAddress("selLeptons_eleDEta", selLeptons_eleDEta, &b_selLeptons_eleDEta);
   fChain->SetBranchAddress("selLeptons_eleDPhi", selLeptons_eleDPhi, &b_selLeptons_eleDPhi);
   fChain->SetBranchAddress("selLeptons_eleHoE", selLeptons_eleHoE, &b_selLeptons_eleHoE);
   fChain->SetBranchAddress("selLeptons_eleMissingHits", selLeptons_eleMissingHits, &b_selLeptons_eleMissingHits);
   fChain->SetBranchAddress("selLeptons_eleChi2", selLeptons_eleChi2, &b_selLeptons_eleChi2);
   fChain->SetBranchAddress("selLeptons_muonDX", selLeptons_muonDX, &b_selLeptons_muonDX);
   fChain->SetBranchAddress("selLeptons_eleClusterEta", selLeptons_eleClusterEta, &b_selLeptons_eleClusterEta);
   fChain->SetBranchAddress("selLeptons_eleClusterEnergy", selLeptons_eleClusterEnergy, &b_selLeptons_eleClusterEnergy);
   fChain->SetBranchAddress("selLeptons_eleClusterDEta", selLeptons_eleClusterDEta, &b_selLeptons_eleClusterDEta);
   fChain->SetBranchAddress("selLeptons_eleClusterDPhi", selLeptons_eleClusterDPhi, &b_selLeptons_eleClusterDPhi);
   fChain->SetBranchAddress("selLeptons_nStations", selLeptons_nStations, &b_selLeptons_nStations);
   fChain->SetBranchAddress("selLeptons_trkKink", selLeptons_trkKink, &b_selLeptons_trkKink);
   fChain->SetBranchAddress("selLeptons_caloCompatibility", selLeptons_caloCompatibility, &b_selLeptons_caloCompatibility);
   fChain->SetBranchAddress("selLeptons_globalTrackChi2", selLeptons_globalTrackChi2, &b_selLeptons_globalTrackChi2);
   fChain->SetBranchAddress("selLeptons_nChamberHits", selLeptons_nChamberHits, &b_selLeptons_nChamberHits);
   fChain->SetBranchAddress("selLeptons_isBarrelEle", selLeptons_isBarrelEle, &b_selLeptons_isBarrelEle);
   fChain->SetBranchAddress("selLeptons_isEndCapEle", selLeptons_isEndCapEle, &b_selLeptons_isEndCapEle);
   fChain->SetBranchAddress("selLeptons_isEcalDriven", selLeptons_isEcalDriven, &b_selLeptons_isEcalDriven);
   fChain->SetBranchAddress("selLeptons_isMyGoodMuon", selLeptons_isMyGoodMuon, &b_selLeptons_isMyGoodMuon);
   fChain->SetBranchAddress("selLeptons_isMyGoodMuon1", selLeptons_isMyGoodMuon1, &b_selLeptons_isMyGoodMuon1);
   fChain->SetBranchAddress("selLeptons_isMyGoodElectron", selLeptons_isMyGoodElectron, &b_selLeptons_isMyGoodElectron);
   fChain->SetBranchAddress("selLeptons_relPtError", selLeptons_relPtError, &b_selLeptons_relPtError);
   fChain->SetBranchAddress("selLeptons_isPFMuon", selLeptons_isPFMuon, &b_selLeptons_isPFMuon);
   fChain->SetBranchAddress("selLeptons_muon_dz", selLeptons_muon_dz, &b_selLeptons_muon_dz);
   fChain->SetBranchAddress("selLeptons_isHighPtMuon", selLeptons_isHighPtMuon, &b_selLeptons_isHighPtMuon);
   fChain->SetBranchAddress("selLeptons_muTrackIso", selLeptons_muTrackIso, &b_selLeptons_muTrackIso);
   fChain->SetBranchAddress("selLeptons_isGlobalMuon", selLeptons_isGlobalMuon, &b_selLeptons_isGlobalMuon);
   fChain->SetBranchAddress("selLeptons_isTrackerMuon", selLeptons_isTrackerMuon, &b_selLeptons_isTrackerMuon);
   fChain->SetBranchAddress("selLeptons_pixelHits", selLeptons_pixelHits, &b_selLeptons_pixelHits);
   fChain->SetBranchAddress("selLeptons_trackerLayers", selLeptons_trackerLayers, &b_selLeptons_trackerLayers);
   fChain->SetBranchAddress("selLeptons_pixelLayers", selLeptons_pixelLayers, &b_selLeptons_pixelLayers);
   fChain->SetBranchAddress("selLeptons_mvaTTH", selLeptons_mvaTTH, &b_selLeptons_mvaTTH);
   fChain->SetBranchAddress("selLeptons_jetOverlapIdx", selLeptons_jetOverlapIdx, &b_selLeptons_jetOverlapIdx);
   fChain->SetBranchAddress("selLeptons_jetPtRatio", selLeptons_jetPtRatio, &b_selLeptons_jetPtRatio);
   fChain->SetBranchAddress("selLeptons_jetBTagCSV", selLeptons_jetBTagCSV, &b_selLeptons_jetBTagCSV);
   fChain->SetBranchAddress("selLeptons_jetDR", selLeptons_jetDR, &b_selLeptons_jetDR);
   fChain->SetBranchAddress("selLeptons_pfRelIso03", selLeptons_pfRelIso03, &b_selLeptons_pfRelIso03);
   fChain->SetBranchAddress("selLeptons_pfRelIso04", selLeptons_pfRelIso04, &b_selLeptons_pfRelIso04);
   fChain->SetBranchAddress("selLeptons_muonDB", selLeptons_muonDB, &b_selLeptons_muonDB);
   fChain->SetBranchAddress("selLeptons_etaSc", selLeptons_etaSc, &b_selLeptons_etaSc);
   fChain->SetBranchAddress("selLeptons_eleExpMissingInnerHits", selLeptons_eleExpMissingInnerHits, &b_selLeptons_eleExpMissingInnerHits);
   fChain->SetBranchAddress("selLeptons_eleooEmooP", selLeptons_eleooEmooP, &b_selLeptons_eleooEmooP);
   fChain->SetBranchAddress("selLeptons_muonTrackerLayers", selLeptons_muonTrackerLayers, &b_selLeptons_muonTrackerLayers);
   fChain->SetBranchAddress("nFatjetCleanAK08ungroomed", &nFatjetCleanAK08ungroomed, &b_nFatjetCleanAK08ungroomed);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_pt", FatjetCleanAK08ungroomed_pt, &b_FatjetCleanAK08ungroomed_pt);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_eta", FatjetCleanAK08ungroomed_eta, &b_FatjetCleanAK08ungroomed_eta);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_phi", FatjetCleanAK08ungroomed_phi, &b_FatjetCleanAK08ungroomed_phi);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_mass", FatjetCleanAK08ungroomed_mass, &b_FatjetCleanAK08ungroomed_mass);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_chHEFrac", FatjetCleanAK08ungroomed_chHEFrac, &b_FatjetCleanAK08ungroomed_chHEFrac);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_neHEFrac", FatjetCleanAK08ungroomed_neHEFrac, &b_FatjetCleanAK08ungroomed_neHEFrac);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_chEmEFrac", FatjetCleanAK08ungroomed_chEmEFrac, &b_FatjetCleanAK08ungroomed_chEmEFrac);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_neEmEFrac", FatjetCleanAK08ungroomed_neEmEFrac, &b_FatjetCleanAK08ungroomed_neEmEFrac);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_chMult", FatjetCleanAK08ungroomed_chMult, &b_FatjetCleanAK08ungroomed_chMult);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_neMult", FatjetCleanAK08ungroomed_neMult, &b_FatjetCleanAK08ungroomed_neMult);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_muEFrac", FatjetCleanAK08ungroomed_muEFrac, &b_FatjetCleanAK08ungroomed_muEFrac);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_tau1", FatjetCleanAK08ungroomed_tau1, &b_FatjetCleanAK08ungroomed_tau1);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_tau2", FatjetCleanAK08ungroomed_tau2, &b_FatjetCleanAK08ungroomed_tau2);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_tau3", FatjetCleanAK08ungroomed_tau3, &b_FatjetCleanAK08ungroomed_tau3);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_msoftdrop", FatjetCleanAK08ungroomed_msoftdrop, &b_FatjetCleanAK08ungroomed_msoftdrop);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_mpruned", FatjetCleanAK08ungroomed_mpruned, &b_FatjetCleanAK08ungroomed_mpruned);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_mtrimmed", FatjetCleanAK08ungroomed_mtrimmed, &b_FatjetCleanAK08ungroomed_mtrimmed);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_mfiltered", FatjetCleanAK08ungroomed_mfiltered, &b_FatjetCleanAK08ungroomed_mfiltered);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_bbtag", FatjetCleanAK08ungroomed_bbtag, &b_FatjetCleanAK08ungroomed_bbtag);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_PFLepton_ptrel", FatjetCleanAK08ungroomed_PFLepton_ptrel, &b_FatjetCleanAK08ungroomed_PFLepton_ptrel);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_z_ratio", FatjetCleanAK08ungroomed_z_ratio, &b_FatjetCleanAK08ungroomed_z_ratio);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_tau_dot", FatjetCleanAK08ungroomed_tau_dot, &b_FatjetCleanAK08ungroomed_tau_dot);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_SV_mass_0", FatjetCleanAK08ungroomed_SV_mass_0, &b_FatjetCleanAK08ungroomed_SV_mass_0);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_SV_EnergyRatio_0", FatjetCleanAK08ungroomed_SV_EnergyRatio_0, &b_FatjetCleanAK08ungroomed_SV_EnergyRatio_0);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_SV_EnergyRatio_1", FatjetCleanAK08ungroomed_SV_EnergyRatio_1, &b_FatjetCleanAK08ungroomed_SV_EnergyRatio_1);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_PFLepton_IP2D", FatjetCleanAK08ungroomed_PFLepton_IP2D, &b_FatjetCleanAK08ungroomed_PFLepton_IP2D);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_tau_21", FatjetCleanAK08ungroomed_tau_21, &b_FatjetCleanAK08ungroomed_tau_21);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_nSL", FatjetCleanAK08ungroomed_nSL, &b_FatjetCleanAK08ungroomed_nSL);
   fChain->SetBranchAddress("FatjetCleanAK08ungroomed_vertexNTracks", FatjetCleanAK08ungroomed_vertexNTracks, &b_FatjetCleanAK08ungroomed_vertexNTracks);
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
   fChain->SetBranchAddress("nFatjetCleanAK08prunedcal", &nFatjetCleanAK08prunedcal, &b_nFatjetCleanAK08prunedcal);
   fChain->SetBranchAddress("FatjetCleanAK08prunedcal_pt", FatjetCleanAK08prunedcal_pt, &b_FatjetCleanAK08prunedcal_pt);
   fChain->SetBranchAddress("FatjetCleanAK08prunedcal_eta", FatjetCleanAK08prunedcal_eta, &b_FatjetCleanAK08prunedcal_eta);
   fChain->SetBranchAddress("FatjetCleanAK08prunedcal_phi", FatjetCleanAK08prunedcal_phi, &b_FatjetCleanAK08prunedcal_phi);
   fChain->SetBranchAddress("FatjetCleanAK08prunedcal_mass", FatjetCleanAK08prunedcal_mass, &b_FatjetCleanAK08prunedcal_mass);
   fChain->SetBranchAddress("nSubjetAK08pruned", &nSubjetAK08pruned, &b_nSubjetAK08pruned);
   fChain->SetBranchAddress("SubjetAK08pruned_pt", SubjetAK08pruned_pt, &b_SubjetAK08pruned_pt);
   fChain->SetBranchAddress("SubjetAK08pruned_eta", SubjetAK08pruned_eta, &b_SubjetAK08pruned_eta);
   fChain->SetBranchAddress("SubjetAK08pruned_phi", SubjetAK08pruned_phi, &b_SubjetAK08pruned_phi);
   fChain->SetBranchAddress("SubjetAK08pruned_mass", SubjetAK08pruned_mass, &b_SubjetAK08pruned_mass);
   fChain->SetBranchAddress("SubjetAK08pruned_btag", SubjetAK08pruned_btag, &b_SubjetAK08pruned_btag);
   fChain->SetBranchAddress("nGenStatus2bHad", &nGenStatus2bHad, &b_nGenStatus2bHad);
   fChain->SetBranchAddress("GenStatus2bHad_pdgId", GenStatus2bHad_pdgId, &b_GenStatus2bHad_pdgId);
   fChain->SetBranchAddress("GenStatus2bHad_pt", GenStatus2bHad_pt, &b_GenStatus2bHad_pt);
   fChain->SetBranchAddress("GenStatus2bHad_eta", GenStatus2bHad_eta, &b_GenStatus2bHad_eta);
   fChain->SetBranchAddress("GenStatus2bHad_phi", GenStatus2bHad_phi, &b_GenStatus2bHad_phi);
   fChain->SetBranchAddress("GenStatus2bHad_mass", GenStatus2bHad_mass, &b_GenStatus2bHad_mass);
   fChain->SetBranchAddress("GenStatus2bHad_charge", GenStatus2bHad_charge, &b_GenStatus2bHad_charge);
   fChain->SetBranchAddress("GenStatus2bHad_status", GenStatus2bHad_status, &b_GenStatus2bHad_status);
   fChain->SetBranchAddress("nGenHadTaus", &nGenHadTaus, &b_nGenHadTaus);
   fChain->SetBranchAddress("GenHadTaus_charge", GenHadTaus_charge, &b_GenHadTaus_charge);
   fChain->SetBranchAddress("GenHadTaus_status", GenHadTaus_status, &b_GenHadTaus_status);
   fChain->SetBranchAddress("GenHadTaus_pdgId", GenHadTaus_pdgId, &b_GenHadTaus_pdgId);
   fChain->SetBranchAddress("GenHadTaus_pt", GenHadTaus_pt, &b_GenHadTaus_pt);
   fChain->SetBranchAddress("GenHadTaus_eta", GenHadTaus_eta, &b_GenHadTaus_eta);
   fChain->SetBranchAddress("GenHadTaus_phi", GenHadTaus_phi, &b_GenHadTaus_phi);
   fChain->SetBranchAddress("GenHadTaus_mass", GenHadTaus_mass, &b_GenHadTaus_mass);
   fChain->SetBranchAddress("GenHadTaus_decayMode", GenHadTaus_decayMode, &b_GenHadTaus_decayMode);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_id", Jet_id, &b_Jet_id);
   fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   fChain->SetBranchAddress("Jet_btagCSV", Jet_btagCSV, &b_Jet_btagCSV);
   fChain->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA, &b_Jet_btagCMVA);
   fChain->SetBranchAddress("Jet_rawPt", Jet_rawPt, &b_Jet_rawPt);
   fChain->SetBranchAddress("Jet_mcPt", Jet_mcPt, &b_Jet_mcPt);
   fChain->SetBranchAddress("Jet_mcFlavour", Jet_mcFlavour, &b_Jet_mcFlavour);
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
   fChain->SetBranchAddress("Jet_idxFirstTauMatch", Jet_idxFirstTauMatch, &b_Jet_idxFirstTauMatch);
   fChain->SetBranchAddress("Jet_heppyFlavour", Jet_heppyFlavour, &b_Jet_heppyFlavour);
   fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   fChain->SetBranchAddress("Jet_btagBDT", Jet_btagBDT, &b_Jet_btagBDT);
   fChain->SetBranchAddress("Jet_btagProb", Jet_btagProb, &b_Jet_btagProb);
   fChain->SetBranchAddress("Jet_btagBProb", Jet_btagBProb, &b_Jet_btagBProb);
   fChain->SetBranchAddress("Jet_btagSoftEl", Jet_btagSoftEl, &b_Jet_btagSoftEl);
   fChain->SetBranchAddress("Jet_btagSoftMu", Jet_btagSoftMu, &b_Jet_btagSoftMu);
   fChain->SetBranchAddress("Jet_btagnew", Jet_btagnew, &b_Jet_btagnew);
   fChain->SetBranchAddress("Jet_btagCSVV0", Jet_btagCSVV0, &b_Jet_btagCSVV0);
   fChain->SetBranchAddress("Jet_isMyGoodLooseJet", Jet_isMyGoodLooseJet, &b_Jet_isMyGoodLooseJet);
   fChain->SetBranchAddress("Jet_isMyGoodTightJet", Jet_isMyGoodTightJet, &b_Jet_isMyGoodTightJet);
   fChain->SetBranchAddress("Jet_isMyGoodTightLepVetoJet", Jet_isMyGoodTightLepVetoJet, &b_Jet_isMyGoodTightLepVetoJet);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_chMult", Jet_chMult, &b_Jet_chMult);
   fChain->SetBranchAddress("Jet_neMult", Jet_neMult, &b_Jet_neMult);
   fChain->SetBranchAddress("Jet_leadTrackPt", Jet_leadTrackPt, &b_Jet_leadTrackPt);
   fChain->SetBranchAddress("Jet_mcEta", Jet_mcEta, &b_Jet_mcEta);
   fChain->SetBranchAddress("Jet_mcPhi", Jet_mcPhi, &b_Jet_mcPhi);
   fChain->SetBranchAddress("Jet_mcM", Jet_mcM, &b_Jet_mcM);
   fChain->SetBranchAddress("Jet_leptonPdgId", Jet_leptonPdgId, &b_Jet_leptonPdgId);
   fChain->SetBranchAddress("Jet_leptonPt", Jet_leptonPt, &b_Jet_leptonPt);
   fChain->SetBranchAddress("Jet_leptonPtRel", Jet_leptonPtRel, &b_Jet_leptonPtRel);
   fChain->SetBranchAddress("Jet_leptonPtRelInv", Jet_leptonPtRelInv, &b_Jet_leptonPtRelInv);
   fChain->SetBranchAddress("Jet_leptonDeltaR", Jet_leptonDeltaR, &b_Jet_leptonDeltaR);
   fChain->SetBranchAddress("Jet_leptonDeltaPhi", Jet_leptonDeltaPhi, &b_Jet_leptonDeltaPhi);
   fChain->SetBranchAddress("Jet_leptonDeltaEta", Jet_leptonDeltaEta, &b_Jet_leptonDeltaEta);
   fChain->SetBranchAddress("Jet_vtxMass", Jet_vtxMass, &b_Jet_vtxMass);
   fChain->SetBranchAddress("Jet_vtxNtracks", Jet_vtxNtracks, &b_Jet_vtxNtracks);
   fChain->SetBranchAddress("Jet_vtxPt", Jet_vtxPt, &b_Jet_vtxPt);
   fChain->SetBranchAddress("Jet_vtx3DSig", Jet_vtx3DSig, &b_Jet_vtx3DSig);
   fChain->SetBranchAddress("Jet_vtx3DVal", Jet_vtx3DVal, &b_Jet_vtx3DVal);
   fChain->SetBranchAddress("Jet_vtxPosX", Jet_vtxPosX, &b_Jet_vtxPosX);
   fChain->SetBranchAddress("Jet_vtxPosY", Jet_vtxPosY, &b_Jet_vtxPosY);
   fChain->SetBranchAddress("Jet_vtxPosZ", Jet_vtxPosZ, &b_Jet_vtxPosZ);
   fChain->SetBranchAddress("Jet_pullVectorPhi", Jet_pullVectorPhi, &b_Jet_pullVectorPhi);
   fChain->SetBranchAddress("Jet_pullVectorMag", Jet_pullVectorMag, &b_Jet_pullVectorMag);
   fChain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   fChain->SetBranchAddress("Jet_ptd", Jet_ptd, &b_Jet_ptd);
   fChain->SetBranchAddress("Jet_axis2", Jet_axis2, &b_Jet_axis2);
   fChain->SetBranchAddress("Jet_mult", Jet_mult, &b_Jet_mult);
   fChain->SetBranchAddress("Jet_numberOfDaughters", Jet_numberOfDaughters, &b_Jet_numberOfDaughters);
   fChain->SetBranchAddress("Jet_btagIdx", Jet_btagIdx, &b_Jet_btagIdx);
   fChain->SetBranchAddress("Jet_mcIdx", Jet_mcIdx, &b_Jet_mcIdx);
   fChain->SetBranchAddress("Jet_pt_reg", Jet_pt_reg, &b_Jet_pt_reg);
   fChain->SetBranchAddress("Jet_pt_regVBF", Jet_pt_regVBF, &b_Jet_pt_regVBF);
   fChain->SetBranchAddress("Jet_blike_VBF", Jet_blike_VBF, &b_Jet_blike_VBF);
   fChain->SetBranchAddress("Jet_bTagWeightJESUp", Jet_bTagWeightJESUp, &b_Jet_bTagWeightJESUp);
   fChain->SetBranchAddress("Jet_bTagWeightJESDown", Jet_bTagWeightJESDown, &b_Jet_bTagWeightJESDown);
   fChain->SetBranchAddress("Jet_bTagWeightLFUp", Jet_bTagWeightLFUp, &b_Jet_bTagWeightLFUp);
   fChain->SetBranchAddress("Jet_bTagWeightLFDown", Jet_bTagWeightLFDown, &b_Jet_bTagWeightLFDown);
   fChain->SetBranchAddress("Jet_bTagWeightHFUp", Jet_bTagWeightHFUp, &b_Jet_bTagWeightHFUp);
   fChain->SetBranchAddress("Jet_bTagWeightHFDown", Jet_bTagWeightHFDown, &b_Jet_bTagWeightHFDown);
   fChain->SetBranchAddress("Jet_bTagWeightStats1Up", Jet_bTagWeightStats1Up, &b_Jet_bTagWeightStats1Up);
   fChain->SetBranchAddress("Jet_bTagWeightStats1Down", Jet_bTagWeightStats1Down, &b_Jet_bTagWeightStats1Down);
   fChain->SetBranchAddress("Jet_bTagWeightStats2Up", Jet_bTagWeightStats2Up, &b_Jet_bTagWeightStats2Up);
   fChain->SetBranchAddress("Jet_bTagWeightStats2Down", Jet_bTagWeightStats2Down, &b_Jet_bTagWeightStats2Down);
   fChain->SetBranchAddress("Jet_bTagWeightcErr1Up", Jet_bTagWeightcErr1Up, &b_Jet_bTagWeightcErr1Up);
   fChain->SetBranchAddress("Jet_bTagWeightcErr1Down", Jet_bTagWeightcErr1Down, &b_Jet_bTagWeightcErr1Down);
   fChain->SetBranchAddress("Jet_bTagWeightcErr2Up", Jet_bTagWeightcErr2Up, &b_Jet_bTagWeightcErr2Up);
   fChain->SetBranchAddress("Jet_bTagWeightcErr2Down", Jet_bTagWeightcErr2Down, &b_Jet_bTagWeightcErr2Down);
   fChain->SetBranchAddress("Jet_bTagWeight", Jet_bTagWeight, &b_Jet_bTagWeight);
   fChain->SetBranchAddress("nFatjetGenAK08Pruned", &nFatjetGenAK08Pruned, &b_nFatjetGenAK08Pruned);
   fChain->SetBranchAddress("FatjetGenAK08Pruned_pt", FatjetGenAK08Pruned_pt, &b_FatjetGenAK08Pruned_pt);
   fChain->SetBranchAddress("FatjetGenAK08Pruned_eta", FatjetGenAK08Pruned_eta, &b_FatjetGenAK08Pruned_eta);
   fChain->SetBranchAddress("FatjetGenAK08Pruned_phi", FatjetGenAK08Pruned_phi, &b_FatjetGenAK08Pruned_phi);
   fChain->SetBranchAddress("FatjetGenAK08Pruned_mass", FatjetGenAK08Pruned_mass, &b_FatjetGenAK08Pruned_mass);
   fChain->SetBranchAddress("npileUpVertex_z", &npileUpVertex_z, &b_npileUpVertex_z);
   fChain->SetBranchAddress("pileUpVertex_z", pileUpVertex_z, &b_pileUpVertex_z);
   fChain->SetBranchAddress("nmyGenJetAK08", &nmyGenJetAK08, &b_nmyGenJetAK08);
   fChain->SetBranchAddress("myGenJetAK08_pt", myGenJetAK08_pt, &b_myGenJetAK08_pt);
   fChain->SetBranchAddress("myGenJetAK08_eta", myGenJetAK08_eta, &b_myGenJetAK08_eta);
   fChain->SetBranchAddress("myGenJetAK08_phi", myGenJetAK08_phi, &b_myGenJetAK08_phi);
   fChain->SetBranchAddress("myGenJetAK08_mass", myGenJetAK08_mass, &b_myGenJetAK08_mass);
   fChain->SetBranchAddress("myGenJetAK08_tau1", myGenJetAK08_tau1, &b_myGenJetAK08_tau1);
   fChain->SetBranchAddress("myGenJetAK08_tau2", myGenJetAK08_tau2, &b_myGenJetAK08_tau2);
   fChain->SetBranchAddress("myGenJetAK08_tau3", myGenJetAK08_tau3, &b_myGenJetAK08_tau3);
   fChain->SetBranchAddress("nprimaryVertices", &nprimaryVertices, &b_nprimaryVertices);
   fChain->SetBranchAddress("primaryVertices_x", primaryVertices_x, &b_primaryVertices_x);
   fChain->SetBranchAddress("primaryVertices_y", primaryVertices_y, &b_primaryVertices_y);
   fChain->SetBranchAddress("primaryVertices_z", primaryVertices_z, &b_primaryVertices_z);
   fChain->SetBranchAddress("primaryVertices_isFake", primaryVertices_isFake, &b_primaryVertices_isFake);
   fChain->SetBranchAddress("primaryVertices_ndof", primaryVertices_ndof, &b_primaryVertices_ndof);
   fChain->SetBranchAddress("primaryVertices_Rho", primaryVertices_Rho, &b_primaryVertices_Rho);
   fChain->SetBranchAddress("nFatjetAK08prunedReg", &nFatjetAK08prunedReg, &b_nFatjetAK08prunedReg);
   fChain->SetBranchAddress("FatjetAK08prunedReg_pt", FatjetAK08prunedReg_pt, &b_FatjetAK08prunedReg_pt);
   fChain->SetBranchAddress("FatjetAK08prunedReg_eta", FatjetAK08prunedReg_eta, &b_FatjetAK08prunedReg_eta);
   fChain->SetBranchAddress("FatjetAK08prunedReg_phi", FatjetAK08prunedReg_phi, &b_FatjetAK08prunedReg_phi);
   fChain->SetBranchAddress("FatjetAK08prunedReg_mass", FatjetAK08prunedReg_mass, &b_FatjetAK08prunedReg_mass);
   fChain->SetBranchAddress("nGenWZQuark", &nGenWZQuark, &b_nGenWZQuark);
   fChain->SetBranchAddress("GenWZQuark_pdgId", GenWZQuark_pdgId, &b_GenWZQuark_pdgId);
   fChain->SetBranchAddress("GenWZQuark_pt", GenWZQuark_pt, &b_GenWZQuark_pt);
   fChain->SetBranchAddress("GenWZQuark_eta", GenWZQuark_eta, &b_GenWZQuark_eta);
   fChain->SetBranchAddress("GenWZQuark_phi", GenWZQuark_phi, &b_GenWZQuark_phi);
   fChain->SetBranchAddress("GenWZQuark_mass", GenWZQuark_mass, &b_GenWZQuark_mass);
   fChain->SetBranchAddress("GenWZQuark_charge", GenWZQuark_charge, &b_GenWZQuark_charge);
   fChain->SetBranchAddress("GenWZQuark_status", GenWZQuark_status, &b_GenWZQuark_status);
   fChain->SetBranchAddress("nGenLepFromWPrime", &nGenLepFromWPrime, &b_nGenLepFromWPrime);
   fChain->SetBranchAddress("GenLepFromWPrime_pdgId", &GenLepFromWPrime_pdgId, &b_GenLepFromWPrime_pdgId);
   fChain->SetBranchAddress("GenLepFromWPrime_pt", &GenLepFromWPrime_pt, &b_GenLepFromWPrime_pt);
   fChain->SetBranchAddress("GenLepFromWPrime_eta", &GenLepFromWPrime_eta, &b_GenLepFromWPrime_eta);
   fChain->SetBranchAddress("GenLepFromWPrime_phi", &GenLepFromWPrime_phi, &b_GenLepFromWPrime_phi);
   fChain->SetBranchAddress("GenLepFromWPrime_mass", &GenLepFromWPrime_mass, &b_GenLepFromWPrime_mass);
   fChain->SetBranchAddress("GenLepFromWPrime_charge", &GenLepFromWPrime_charge, &b_GenLepFromWPrime_charge);
   fChain->SetBranchAddress("GenLepFromWPrime_status", &GenLepFromWPrime_status, &b_GenLepFromWPrime_status);
   fChain->SetBranchAddress("nGenBQuarkFromTop", &nGenBQuarkFromTop, &b_nGenBQuarkFromTop);
   fChain->SetBranchAddress("GenBQuarkFromTop_pdgId", &GenBQuarkFromTop_pdgId, &b_GenBQuarkFromTop_pdgId);
   fChain->SetBranchAddress("GenBQuarkFromTop_pt", &GenBQuarkFromTop_pt, &b_GenBQuarkFromTop_pt);
   fChain->SetBranchAddress("GenBQuarkFromTop_eta", &GenBQuarkFromTop_eta, &b_GenBQuarkFromTop_eta);
   fChain->SetBranchAddress("GenBQuarkFromTop_phi", &GenBQuarkFromTop_phi, &b_GenBQuarkFromTop_phi);
   fChain->SetBranchAddress("GenBQuarkFromTop_mass", &GenBQuarkFromTop_mass, &b_GenBQuarkFromTop_mass);
   fChain->SetBranchAddress("GenBQuarkFromTop_charge", &GenBQuarkFromTop_charge, &b_GenBQuarkFromTop_charge);
   fChain->SetBranchAddress("GenBQuarkFromTop_status", &GenBQuarkFromTop_status, &b_GenBQuarkFromTop_status);
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
