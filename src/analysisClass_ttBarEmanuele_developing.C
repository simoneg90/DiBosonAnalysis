#define analysisClass_cxx
#include "analysisClass.h"
#include "setTDRStyle.h"
#include "utility.h"

#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TMath.h>
#include "TStopwatch.h"

#define e_mass 0.0005 //GeV
#define mu_mass 0.1 //GeV
#define W_mass 80 //GeV

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

/*  string jetAlgo = getPreCutString1("jetAlgo");
  double rParam = getPreCutValue1("DeltaR");

  if( jetAlgo == "AntiKt" )
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::antikt_algorithm, rParam) );
  else if( jetAlgo == "Kt" )
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::kt_algorithm, rParam) );
  else 
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::cambridge_algorithm, rParam) );
*/
  std::cout << "analysisClass::analysisClass(): ends " << std::endl;
}

analysisClass::~analysisClass()
{
  std::cout << "analysisClass::~analysisClass(): begins " << std::endl;

  std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   
   TStopwatch time;
   frame("Starting the analysis");
   time.Start(true);
   setTDRStyle();
   if (fChain == 0) return;
   
   //////////book histos here

   // TH1F *h_nJetFinal = new TH1F ("h_nJetFinal","",10,0,10);
   // h_nJetFinal->Sumw2();      
   // TH1F *h_nVtx = new TH1F ("h_nVtx","",30,0,30);
   // h_nVtx->Sumw2(); 
   // TH1F *h_trueVtx = new TH1F ("h_trueVtx","",40,0,40);
   // h_trueVtx->Sumw2();  
   // TH1F *h_pT1stJet = new TH1F ("h_pT1stJet","",100,0,3000);
   // h_pT1stJet->Sumw2();
   // TH1F *h_pT2ndJet = new TH1F ("h_pT2ndJet","",100,0,3000);
   // h_pT2ndJet->Sumw2();
   // TH1F *h_eta1stJet = new TH1F ("h_eta1stJet","",5,-2.5,2.5);
   // h_eta1stJet->Sumw2();
   // TH1F *h_eta2ndJet = new TH1F ("h_eta2ndJet","",5,-2.5,2.5);
   // h_eta2ndJet->Sumw2();
   // TH1F *h_DijetMass = new TH1F ("h_DijetMass","",600,0,6000);
   // h_DijetMass->Sumw2();
   // TH1F *h_DeltaETAjj = new TH1F ("h_DeltaETAjj","",120,0,3.);
   // h_DeltaETAjj->Sumw2();

   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   double lepton_mass;
   std::vector <int> goodLepton, looseLepton, goodAK08, goodPuppiAK04, goodPuppiAK04_lep, goodAK04, goodAK04_lep, goodEle, goodMuon, looseEle, looseMuon;
   TLorentzVector electron, muon, genW, ak04, ak08, ak08Pruned, lepton, leptonLoose, W, MET, wGenQ1, wGenQ2, wGenSumi, subjet1, subjet2, subjetSum, wGenSum, bGen1, bGen2, puppiAK4;
   int btag_ak04_loose, btag_ak04_medium, btag_ak04_tight,btag_ak04Puppi_loose, btag_ak04Puppi_medium, btag_ak04Puppi_tight, b_counter;
   int subjet_index1, subjet_index2;
   double drMin;
   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<1000;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     //if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     if(jentry < 10 || jentry%10000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;
     if(jentry%1000 == 0 && jentry%10000 !=0) std::cout<<"."<<std::flush;
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done for every event
     //std::cout<<"Initializing variables"<<std::endl;
     b_counter=0;
     btag_ak04_loose=0;
     btag_ak04_medium=0;
     btag_ak04_tight=0;
     goodLepton.clear();
     looseLepton.clear();
     goodAK04.clear();
     goodPuppiAK04.clear();
     goodAK08.clear();
     goodAK04_lep.clear();
     goodPuppiAK04_lep.clear();
     goodEle.clear();
     goodMuon.clear();
     looseEle.clear();
     looseMuon.clear();
     //std::cout<<"resetCuts()"<<std::endl;
     resetCuts();
   
     //std::cout<<"reset cuts"<<std::endl;
//=== Forse da inserire ===

/*     fillVariableWithValue("run",runNo);     
     fillVariableWithValue("event",evtNo);     
     fillVariableWithValue("lumi",lumi);     
     fillVariableWithValue("nVtx",nvtx);     
     fillVariableWithValue("nJet",widejets.size());
*/
    
    //std::cout<<"Looping over nLepGood"<<std::endl;
    for(int i=0; i<nLepGood; ++i){
      //loop for muons
      //std::cout<<"Loop muon"<<std::endl;
      //std::cout<<"PDGID: "<<abs(LepGood_pdgId[i])<<std::endl;
      if((abs(LepGood_pdgId[i])==13 && LepGood_isGlobalMuon[i] && LepGood_nChamberHits[i] > 0 && LepGood_nStations[i] > 1 && LepGood_muonPtRatio[i]<0.3 && LepGood_muonDB[i]< 0.2 && LepGood_pixelHits[i]>0 && LepGood_trackerLayers[i]>5&&LepGood_pt[i]>30 && abs(LepGood_eta[i])<2.1)){
      
        goodMuon.push_back(i);
        goodLepton.push_back(i);
      
      }
    }//end loop for muons

    for(int i=0; i<nLepGood; ++i){
      //std::cout<<"Loop ele"<<std::endl;
      //loop for electrons
      //if((abs(LepGood_pdgId[i])==11)){
      if((abs(LepGood_pdgId[i])==11&&(abs(LepGood_etaSc[i])<1.4442 && LepGood_isEcalDriven[i] && abs(LepGood_eleClusterDEta[i])<0.004 && abs(LepGood_dPhiScTrkIn[i])<0.06 && (LepGood_hadronicOverEm[i]<1.0/LepGood_superCluster_energy[i]+0.05) && (LepGood_e2x5Max[i]/LepGood_e5x5[i]>0.94 || LepGood_e1x5[i]/LepGood_e5x5[i]>0.83) && abs(LepGood_dxy[i])<0.02) || (abs(LepGood_etaSc[i])>1.566 && abs(LepGood_etaSc[i])<2.5 && LepGood_isEcalDriven[i] && abs(LepGood_eleClusterDEta[i])<0.006 && abs(LepGood_dPhiScTrkIn[i])< 0.06 && (LepGood_hadronicOverEm[i]<5.0/LepGood_superCluster_energy[i]+0.05) && abs(LepGood_dxy[i])<0.05 && LepGood_full5x5_sigmaIetaIeta[i]<0.03) )&&  LepGood_pt[i]>30 && (abs(LepGood_eta[i])<1.442 || (abs(LepGood_eta[i])>1.56 && abs(LepGood_eta[i])<2.5))){

        //std::cout<<"elettrone"<<std::endl;
        electron.SetPtEtaPhiM(LepGood_pt[i],LepGood_eta[i],LepGood_phi[i],e_mass);

        if(goodMuon.size()==0){
          //std::cout<<"Eccolo"<<std::endl;
          goodEle.push_back(i);
          goodLepton.push_back(i);
        }
        for(int j=0; j<goodMuon.size(); ++j){
          muon.SetPtEtaPhiM(LepGood_pt[goodMuon[j]],LepGood_eta[goodMuon[j]],LepGood_phi[goodMuon[j]],LepGood_mass[goodMuon[j]]);
          //std::cout<<"DR muon ele "<<muon.DeltaR(electron)<<std::endl;
          fillVariableWithValue("muon_ele_DR",muon.DeltaR(electron));
          if(muon.DeltaR(electron)>.1){
            //std::cout<<"deltaR cut passed"<<std::endl;
            goodEle.push_back(i);
            goodLepton.push_back(i);
          }
        }//electron cleaning
      }
    }//end loop for electrons
    //std::cout<<"Electrons: "<<goodEle.size()<<" muons: "<<goodMuon.size()<<std::endl;
    drMin=100.;
    for(int i=0; i<nFatJet;++i){
      if((((FatJet_neHEF[i]<0.99 && FatJet_neEmEF[i]<0.99 && (FatJet_chMult[i]+FatJet_neMult[i])>1) && ((abs(FatJet_eta[i])<=2.4 && FatJet_chHEF[i]>0 && FatJet_chMult[i]>0 && FatJet_chEmEF[i]<0.99) || abs(FatJet_eta[i])>2.4) && abs(FatJet_eta[i])<=3.0)||((FatJet_neEmEF[i]<0.90 && FatJet_neMult[i]>10 && abs(FatJet_eta[i])>3.0 ))) && FatJet_pt[i]>100 && abs(FatJet_eta[i])<2.4){
        ak08.SetPtEtaPhiM(FatJet_pt[i], FatJet_eta[i], FatJet_phi[i], FatJet_mass[i]);
        for(int j=0; j<goodLepton.size(); ++j){
          lepton.SetPtEtaPhiM(LepGood_pt[goodLepton[j]],LepGood_eta[goodLepton[j]],LepGood_phi[goodLepton[j]],LepGood_mass[goodLepton[j]]);
          if(ak08.DeltaR(lepton)<drMin) drMin=ak08.DeltaR(lepton);
          fillVariableWithValue("ak08_lepton_DR", ak08.DeltaR(lepton));
        }
        if(drMin>.1) goodAK08.push_back(i);
      }//end if good AK08
    }//end loop over nFatJet

    drMin=100.;
    for(int i=0; i<nJet; ++i){
      if((((Jet_neHEF[i]<0.99 && Jet_phEF[i]<0.99 && (Jet_chHMult[i]+Jet_neHMult[i]+Jet_phMult[i]+Jet_eMult[i])>1) && ((abs(Jet_eta[i])<=2.4 && Jet_chHEF[i]>0 && (Jet_chHMult[i]+Jet_eMult[i])>0 && Jet_eEF[i]<0.99) || abs(Jet_eta[i])>2.4) && abs(Jet_eta[i])<=3.0)||((Jet_phEF[i]<0.90 && (Jet_neHMult[i]+Jet_phMult[i])>10 && abs(Jet_eta[i])>3.0 ))) && Jet_pt[i]>30 && abs(Jet_eta[i])<2.4){
        //CreateAndFillUserTH1D("goodAk04LooseSelection", 2,-.5,1.5, 1);
        ak04.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
        for(int j=0; j<goodLepton.size(); ++j){
          lepton.SetPtEtaPhiM(LepGood_pt[goodLepton[j]],LepGood_eta[goodLepton[j]],LepGood_phi[goodLepton[j]],LepGood_mass[goodLepton[j]]);
          if(ak04.DeltaR(lepton)<drMin)drMin=ak04.DeltaR(lepton);
        }
        if(drMin>.3) goodAK04.push_back(i);
      }
    }//end loop over nJet

    drMin=100.;
    //std::cout<<"PUPPI"<<std::endl;
    for(int i=0; i<nPuppiJet; ++i){
      //std::cout<<"In Puppi"<<std::endl;
      if((((PuppiJet_neHEF[i]<0.99 && PuppiJet_phEF[i]<0.99 && (PuppiJet_chHMult[i]+PuppiJet_neHMult[i]+PuppiJet_phMult[i]+PuppiJet_eMult[i])>1) && ((abs(PuppiJet_eta[i])<=2.4 && PuppiJet_chHEF[i]>0 && (PuppiJet_chHMult[i]+PuppiJet_eMult[i])>0 && PuppiJet_eEF[i]<0.99) || abs(PuppiJet_eta[i])>2.4) && abs(PuppiJet_eta[i])<=3.0)||((PuppiJet_phEF[i]<0.90 && (PuppiJet_neHMult[i]+PuppiJet_phMult[i])>10 && abs(PuppiJet_eta[i])>3.0 ))) && PuppiJet_pt[i]>30 && abs(PuppiJet_eta[i])<2.4){
        //CreateAndFillUserTH1D("goodAk04LooseSelection", 2,-.5,1.5, 1);
        puppiAK4.SetPtEtaPhiM(PuppiJet_pt[i], PuppiJet_eta[i], PuppiJet_phi[i], PuppiJet_mass[i]);
        for(int j=0; j<goodLepton.size(); ++j){
          lepton.SetPtEtaPhiM(LepGood_pt[goodLepton[j]],LepGood_eta[goodLepton[j]],LepGood_phi[goodLepton[j]],LepGood_mass[goodLepton[j]]);
          if(puppiAK4.DeltaR(lepton)<drMin)drMin=puppiAK4.DeltaR(lepton);
        }
        if(drMin>.3) goodPuppiAK04.push_back(i);
      }
    }//end loop over nPuppiJet

    fillVariableWithValue("lepton_goodNumber", goodLepton.size());
   // fillVariableWithValue("lepton_looseNumber", looseLepton.size());
    fillVariableWithValue("ak08_goodNumber", goodAK08.size());
    fillVariableWithValue("nInitialLepton",nLepGood);
    fillVariableWithValue("nGoodEle", goodEle.size());
  //  fillVariableWithValue("nLooseEle", looseEle.size());
    fillVariableWithValue("nGoodMuon", goodMuon.size());
    fillVariableWithValue("nGoodLepton",goodEle.size()+goodMuon.size());
  //  fillVariableWithValue("nLooseMuon", looseMuon.size());
    if(goodEle.size()>=1){
      fillVariableWithValue("electron_1_pt", LepGood_pt[goodEle[0]]);
      fillVariableWithValue("electron_1_eta", LepGood_eta[goodEle[0]]);
      fillVariableWithValue("electron_1_phi", LepGood_phi[goodEle[0]]);

    }
    if(goodMuon.size()>=1){
      fillVariableWithValue("muon_1_pt", LepGood_pt[goodMuon[0]]);
      fillVariableWithValue("muon_1_eta", LepGood_eta[goodMuon[0]]);
      fillVariableWithValue("muon_1_phi", LepGood_phi[goodMuon[0]]);
      
    }

    if(goodAK08.size()>=1){
      fillVariableWithValue("ak08Ungroomed_1_pt", FatJet_pt[goodAK08[0]]);
      fillVariableWithValue("ak08Ungroomed_1_eta", FatJet_eta[goodAK08[0]]);
      fillVariableWithValue("ak08Ungroomed_1_phi", FatJet_phi[goodAK08[0]]);
      fillVariableWithValue("ak08Ungroomed_1_mass", FatJet_mass[goodAK08[0]]);
      fillVariableWithValue("ak08Ungroomed_1_tau21", FatJet_tau2[goodAK08[0]]/FatJet_tau1[goodAK08[0]]);
      fillVariableWithValue("ak08Pruned_1_mass", FatJet_prunedMass[goodAK08[0]]);
      fillVariableWithValue("ak08Ungroomed_1_CSV", FatJet_btagCSV[goodAK08[0]]);
      fillVariableWithValue("ak08Puppi_1_pt", FatJet_puppiPt[goodAK08[0]]);
      fillVariableWithValue("ak08Puppi_1_eta", FatJet_puppiEta[goodAK08[0]]);
      fillVariableWithValue("ak08Puppi_1_phi", FatJet_puppiPhi[goodAK08[0]]);
      fillVariableWithValue("ak08Puppi_1_mass", FatJet_puppiMass[goodAK08[0]]);
      fillVariableWithValue("ak08Puppi_1_tau21", FatJet_puppiTau2[goodAK08[0]]/FatJet_puppiTau1[goodAK08[0]]);
      
      ak08.SetPtEtaPhiM(FatJet_pt[goodAK08[0]], FatJet_eta[goodAK08[0]], FatJet_phi[goodAK08[0]], FatJet_mass[goodAK08[0]]);
    }
    
    //std::cout<<"Prima del custom "<<ncustomPuppiAK8<<std::endl;
    if(ncustomPuppiAK8>0){
      //std::cout<<" in customPuppi"<<std::endl;
      fillVariableWithValue("ak08CustomPuppi_Pt", customPuppiAK8_pt[0]);
      fillVariableWithValue("ak08CustomPuppi_Eta", customPuppiAK8_eta[0]);
      fillVariableWithValue("ak08CustomPuppi_Phi", customPuppiAK8_phi[0]);
      fillVariableWithValue("ak08CustomPuppi_Mass", customPuppiAK8_mass[0]);
  /*    if(isData==0){
        fillVariableWithValue("ak08CustomPuppi_MCorr", customPuppiAK8_mass[0]);
      }else{
    */    fillVariableWithValue("ak08CustomPuppi_MCorr", customPuppiAK8_massCorrected[0]);
      //}
    }

    double minDR_W=999;
    int w_counter=0;
    //std::cout<<"gen part"<<std::endl;
    if(isData==0){
      for (int w=0; w<nGenPart; ++w){
        if(abs(GenPart_pdgId[w])==24){
          genW.SetPtEtaPhiM(GenPart_pt[w], GenPart_eta[w], GenPart_phi[w], W_mass);
          if(goodLepton.size()>=1){
            if(abs(LepGood_pdgId[goodLepton[0]])==11) lepton_mass=e_mass;
            if(abs(LepGood_pdgId[goodLepton[0]])==13) lepton_mass=mu_mass;
            lepton.SetPtEtaPhiM(LepGood_pt[goodLepton[0]],LepGood_eta[goodLepton[0]],LepGood_phi[goodLepton[0]], lepton_mass);
          } 
          if(ak08.DeltaR(genW)<minDR_W){
            minDR_W=ak08.DeltaR(genW);
            w_counter=w;
          }
        }//if gen Boson==W
      }
    
      //std::cout<<"Second call for gen part"<<std::endl;
      genW.SetPtEtaPhiM(GenPart_pt[w_counter], GenPart_eta[w_counter], GenPart_phi[w_counter], W_mass);
      fillVariableWithValue("ak08Ungroomed_WGen_DR",ak08.DeltaR(genW));
      fillVariableWithValue("WGen_quark_DR", minDR_W);
      fillVariableWithValue("W_Gen_pt", GenPart_pt[w_counter]);
      fillVariableWithValue("W_Gen_eta", GenPart_eta[w_counter]);
      fillVariableWithValue("W_Gen_phi", GenPart_phi[w_counter]);
      fillVariableWithValue("W_Gen_mass", GenPart_mass[w_counter]);
      fillVariableWithValue("lepton_WGen_DR", lepton.DeltaR(genW));
      b_counter=0;
      //std::cout<<"nGenPart"<<std::endl;
      for(int p=0; p<nGenPart; ++p){
        if((abs(GenPart_pdgId[p])==5) && (abs(GenPart_motherId[p])==6)){//b found
          if(b_counter==0) bGen1.SetPtEtaPhiM(GenPart_pt[p],GenPart_eta[p],GenPart_phi[p],GenPart_mass[p]);
          if(b_counter==1) bGen2.SetPtEtaPhiM(GenPart_pt[p],GenPart_eta[p],GenPart_phi[p],GenPart_mass[p]);
          ++b_counter;
        }//end if b quark
      }//end loop over genParticles to find the good bGenJets
      fillVariableWithValue("genW_genBquark1_DR", genW.DeltaR(bGen1));
      fillVariableWithValue("genW_genBquarkMin_DR", TMath::Min(genW.DeltaR(bGen1), genW.DeltaR(bGen2)));
      fillVariableWithValue("genW_genBquarkMax_DR", TMath::Max(genW.DeltaR(bGen1), genW.DeltaR(bGen2)));
      //std::cout<<"1. "<<genW.DeltaR(bGen1)<<std::endl;
      fillVariableWithValue("genW_genBquark2_DR", genW.DeltaR(bGen2));
    }else{
      fillVariableWithValue("ak08Ungroomed_WGen_DR",ak08.DeltaR(genW));
      fillVariableWithValue("WGen_quark_DR", -999.);
      fillVariableWithValue("W_Gen_pt", -999.);
      fillVariableWithValue("W_Gen_eta", -999.);
      fillVariableWithValue("W_Gen_phi", -999.);
      fillVariableWithValue("lepton_WGen_DR", -999.);
      fillVariableWithValue("genW_genBquark1_DR", -999.);
      fillVariableWithValue("genW_genBquarkMin_DR", -999.);
      fillVariableWithValue("genW_genBquarkMax_DR", -999.);
      fillVariableWithValue("genW_genBquark2_DR", -999.);


    }//end if isData

    //std::cout<<"end if IsData"<<std::endl;
    if(goodLepton.size()>=1){
      if(abs(LepGood_pdgId[goodLepton[0]])==11) lepton_mass=e_mass;
      if(abs(LepGood_pdgId[goodLepton[0]])==13) lepton_mass=mu_mass;
      lepton.SetPtEtaPhiM(LepGood_pt[goodLepton[0]],LepGood_eta[goodLepton[0]],LepGood_phi[goodLepton[0]], lepton_mass);
      fillVariableWithValue("lepton_pt", LepGood_pt[goodLepton[0]]);
      fillVariableWithValue("lepton_eta", LepGood_eta[goodLepton[0]]);
      fillVariableWithValue("lepton_phi", LepGood_phi[goodLepton[0]]);
      fillVariableWithValue("lepton_pdgID", LepGood_pdgId[goodLepton[0]]);
      fillVariableWithValue("muonRelIso03", LepGood_relIso03[0]);
      fillVariableWithValue("muontrkIso", LepGood_muTrackIso[0]);
      MET.SetPtEtaPhiM(met_pt, met_eta, met_phi, 0);
      W=lepton+MET;
      fillVariableWithValue("WType1_pt", W.Pt()); //cambiati
      fillVariableWithValue("WType1_eta",W.Eta());
      fillVariableWithValue("WType1_phi",W.Phi());
      fillVariableWithValue("WType1_mT", 2*abs(MET.Pt())*abs(lepton.Pt())*(1-cos(lepton.DeltaPhi(MET))));
    }
    if(goodAK08.size()>=1) fillVariableWithValue("ak08Ungroomed_lepton_DR", lepton.DeltaR(ak08)); //DO NOT USE IT for the analysis cut!!!!!
    fillVariableWithValue("nAK04", goodAK04.size());
    fillVariableWithValue("nPuppiAK04", goodPuppiAK04.size());
    for(int j=0; j<goodAK04.size(); ++j){
      ak04.SetPtEtaPhiM(Jet_pt[goodAK04[j]], Jet_eta[goodAK04[j]], Jet_phi[goodAK04[j]], Jet_mass[goodAK04[j]]);
      if(ak04.DeltaR(ak08)>.8 ){//&& ak04.DeltaR(lepton)>.3){
        //CreateAndFillUserTH1D("Ak04_lepton&AK08_DRCut", 2,-.5,1.5, 1);
        goodAK04_lep.push_back(goodAK04[j]);
        if(Jet_btagCSV[goodAK04[j]]>0.605){
          ++btag_ak04_loose;
        }
        if(Jet_btagCSV[goodAK04[j]]>0.890){
          ++btag_ak04_medium;
        }
        if(Jet_btagCSV[goodAK04[j]]>0.97){
          ++btag_ak04_tight;
        }//end if for CSV medium working point
      }//end if ak04 without leptons and ak08 nearby
    }//end loop over goodAK04.size()
    //std::cout<<"before loop over goodPuppiAK04"<<std::endl;
    for(int j=0; j<goodPuppiAK04.size(); ++j){
      puppiAK4.SetPtEtaPhiM(PuppiJet_pt[goodPuppiAK04[j]], PuppiJet_eta[goodPuppiAK04[j]], PuppiJet_phi[goodPuppiAK04[j]], PuppiJet_mass[goodPuppiAK04[j]]);
      if(puppiAK4.DeltaR(ak08)>.8 ){
        goodPuppiAK04_lep.push_back(goodPuppiAK04[j]);
        if(PuppiJet_btagCSV[goodPuppiAK04[j]]>0.605){
          ++btag_ak04Puppi_loose;
        }
        if(PuppiJet_btagCSV[goodPuppiAK04[j]]>0.890){
          ++btag_ak04Puppi_medium;
        }
        if(PuppiJet_btagCSV[goodPuppiAK04[j]]>0.97){
          ++btag_ak04Puppi_tight;
        }
      }//end if NOT nearby ak08
    }//end loop over goodPuppiAK04
   
      

    if(goodLepton.size()>=2){
      fillVariableWithValue("lepton2_pt", LepGood_pt[goodLepton[1]]);
      fillVariableWithValue("lepton2_eta", LepGood_eta[goodLepton[1]]);
      fillVariableWithValue("lepton2_phi", LepGood_phi[goodLepton[1]]);
      fillVariableWithValue("lepton2_pdgID", LepGood_pdgId[goodLepton[1]]);
    }
    //std::cout<<"after goodLepton.size()>=2"<<std::endl;
    fillVariableWithValue("btag_loose",btag_ak04_loose);
    fillVariableWithValue("btag_medium",btag_ak04_medium);
    fillVariableWithValue("btag_tight",btag_ak04_tight);
    //std::cout<<"after filling btag"<<std::endl;
    if(goodAK04_lep.size()>=1){
      fillVariableWithValue("ak04_1_pt", Jet_pt[goodAK04_lep[0]]);
      fillVariableWithValue("ak04_1_eta", Jet_eta[goodAK04_lep[0]]);
      fillVariableWithValue("ak04_1_phi", Jet_phi[goodAK04_lep[0]]);
      fillVariableWithValue("ak04_1_mass", Jet_mass[goodAK04_lep[0]]);
      //=== TH1D to check the fillReduceSkim procedure ===
      CreateAndFillUserTH1D("ak04_first_pt", 1000,0,500, Jet_pt[goodAK04_lep[0]]);
    }
    if(goodPuppiAK04_lep.size()>=1){
      fillVariableWithValue("ak04Puppi_1_pt", PuppiJet_pt[goodPuppiAK04_lep[0]]);
      fillVariableWithValue("ak04Puppi_1_eta", PuppiJet_eta[goodPuppiAK04_lep[0]]);
      fillVariableWithValue("ak04Puppi_1_phi", PuppiJet_phi[goodPuppiAK04_lep[0]]);
      fillVariableWithValue("ak04Puppi_1_mass", PuppiJet_mass[goodPuppiAK04_lep[0]]);
      fillVariableWithValue("ak04Puppi_1_MCorr", PuppiJet_massCorrected[goodPuppiAK04_lep[0]]);
    }
    if(goodAK04_lep.size()>=2){
      fillVariableWithValue("ak04_2_pt", Jet_pt[goodAK04_lep[1]]);
      fillVariableWithValue("ak04_2_eta", Jet_eta[goodAK04_lep[1]]);
      fillVariableWithValue("ak04_2_phi", Jet_phi[goodAK04_lep[1]]);
      fillVariableWithValue("ak04_2_mass", Jet_mass[goodAK04_lep[1]]);
    }
    fillVariableWithValue("metPuppi", metPuppi_pt);
    fillVariableWithValue("met",met_pt);
    fillVariableWithValue("nPrimaryVertexes", nVert);
    fillVariableWithValue("CSC_filter",Flag_CSCTightHaloFilter);
    fillVariableWithValue("eeBADFilter", Flag_eeBadScFilter);
    fillVariableWithValue("nTrueInteractions", nTrueInt);
    //fillVariableWithValue("run", run);


    //std::cout<<"End of loop"<<std::endl;

     // Evaluate cuts (but do not apply them)
     evaluateCuts();
     //fillReducedSkimTree(); 

     // optional call to fill a skim with the full content of the input roottuple
     //if( passedCut("nJetFinal") ) fillSkimTree();
     
     // optional call to fill a skim with a subset of the variables defined in the cutFile (use flag SAVE)
     //if( passedAllPreviousCuts("mjj") && passedCut("mjj") ) 
     if( passedCut("met"))//passedAllPreviousCuts("eeBADFilter") && passedCut("eeBADFilter"))
       {
         //frame("Beware! This part can be set outside the if -Passed cuts-");
	 fillReducedSkimTree();
	 
	 // ===== Take a look at this =====
	 // //Example on how to investigate quickly the data
 	 // if(getVariableValue("mjj")>4000)
	 //   {
	 //     //fast creation and filling of histograms
	 //     CreateAndFillUserTH1D("h_dphijj_mjjgt4000", 100, 0, 3.15, getVariableValue("deltaPHIjj"));
	 //     CreateAndFillUserTH1D("h_htak4_mjjgt4000", 1000, 0, 10000, getVariableValue("HTAK4"));
	 //     CreateAndFillUserTH1D("h_nvtx_mjjgt4000", 31, -0.5, 30.5, getVariableValue("nVtx"));
	 //   }

       }

     // ===== Example of mjj spectrum after HLT selection =====
     // if( passedAllPreviousCuts("mjj") )
     //   {
     // 	 if(getVariableValue("passHLT")>0)
     // 	   {
     // 	     //fast creation and filling of histograms
     // 	     CreateAndFillUserTH1D("h_mjj_passHLT", getHistoNBins("mjj"), getHistoMin("mjj"), getHistoMax("mjj"), getVariableValue("mjj"));
     // 	   }
     //   }

     // reject events that did not pass level 0 cuts
     //if( !passedCut("0") ) continue;
     // ......
     
     // reject events that did not pass level 1 cuts
     //if( !passedCut("1") ) continue;
     // ......

     // reject events that did not pass the full cut list
     //if( !passedCut("all") ) continue;
     // ......

     // if( widejets.size() >= 2) {
     //  h_nJetFinal->Fill(widejets.size());
     //  h_DijetMass->Fill(wdijet.M());
     //  h_pT1stJet->Fill(widejets[0].Pt());
     //  h_pT2ndJet->Fill(widejets[1].Pt());
     //  h_eta1stJet->Fill(widejets[0].Eta());
     //  h_eta2ndJet->Fill(widejets[1].Eta());
     // }
     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   //////////write histos 

   // h_nVtx->Write();
   // h_trueVtx->Write();
   // h_nJetFinal->Write();
   // h_pT1stJet->Write();
   // h_pT2ndJet->Write();
   // h_DijetMass->Write();
   // h_eta1stJet->Write();
   // h_eta2ndJet->Write();

   // //pT of both jets, to be built using the histograms produced automatically by baseClass
   // TH1F * h_pTJets = new TH1F ("h_pTJets","", getHistoNBins("pT1stJet"), getHistoMin("pT1stJet"), getHistoMax("pT1stJet"));
   // h_pTJets->Add( & getHisto_noCuts_or_skim("pT1stJet") ); // all histos can be retrieved, see other getHisto_xxxx methods in baseClass.h
   // h_pTJets->Add( & getHisto_noCuts_or_skim("pT2ndJet") );
   // //one could also do:  *h_pTJets = getHisto_noCuts_or_skim("pT1stJet") + getHisto_noCuts_or_skim("pT2ndJet");
   // h_pTJets->Write();
   // //one could also do:   const TH1F& h = getHisto_noCuts_or_skim// and use h

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
   std::cout<<""<<std::endl;
   std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << time.RealTime() << " s" <<std::endl;
   std::cout<<""<<std::endl;
}
