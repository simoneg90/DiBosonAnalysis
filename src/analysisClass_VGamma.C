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

#define el_Mass 0.0005 //GeV
#define mu_Mass 0.1 //GeV
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

   /////////initialize and define variables
   TLorentzVector photon, electron, muon, ak04, ak08, ak04_photon, ak08_photon;
   std::vector <int> goodPhoton, goodElectron, goodMuon, goodAk04, goodAk08;
   bool mu_jet_DR, el_jet_DR;
   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<1000;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     if(jentry%100 == 0 && jentry%1000 !=0) std::cout<<"."<<std::flush; 


     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done for every event
     goodPhoton.clear();
     goodElectron.clear();
     goodMuon.clear();
     goodAk04.clear();
     goodAk08.clear();

     resetCuts();

     ////Start coding

     for(int i=0; i<ph_N; ++i){

        if(ph_passLooseId->at(i)==1 && ph_pt->at(i)>180 && abs(ph_eta->at(i))<2.4){
          std::cout<<"Good photon"<<std::endl;
          goodPhoton.push_back(i);
        }

     }//end loop over # of photons


/////     for(int i=0; i<mu_N; ++i){
/////
/////       if(mu_isLooseMuon->at(i)==1 && mu_pt->at(i)>30 && abs(mu_eta->at(i))<2.1){
/////         std::cout<<"Good muon"<<std::endl; 
/////         goodMuon.push_back(i);
/////       }
/////
/////     }//end loop over # of muons
/////
/////     for(int i=0; i<el_N; ++i){
/////       
/////       if(el_isLooseElectron->at(i)==1 && el_pt->at(i)>30 && abs(el_eta->at(i))<2.5){
/////         electron.SetPtEtaPhiM(el_pt->at(i), el_eta->at(i), el_phi->at(i), el_Mass);
/////         for(int j=0; j<goodMuon.size(); ++j){
/////            muon.SetPtEtaPhiM(mu_pt->at(j), mu_eta->at(j), mu_phi->at(j), mu_Mass);
/////            if(muon.DeltaR(electron)>.1) std::cout<<"Good electron"<<std::endl; goodElectron.push_back(i);
/////         }//muon cleaning
/////        
/////       }//end 'if' selection for electrons
/////
/////     }//end loop over # of muons


     for(int i=0; i<jetAK8_N; ++i){
    
        if(jetAK8_pt->at(i)>180 && abs(jetAK8_eta->at(i))<2. && jetAK8_IDLoose->at(i)==1){
            ak08.SetPtEtaPhiM(jetAK8_pt->at(i), jetAK8_eta->at(i), jetAK8_phi->at(i), jetAK8_mass->at(i));
/////            for(int j=0; j<goodMuon.size(); ++j){
/////              muon.SetPtEtaPhiM(mu_pt->at(j), mu_eta->at(j), mu_phi->at(j), mu_Mass);
/////              if(muon.DeltaR(ak08)>.8){
/////                std::cout<<"Good ak08"<<std::endl; 
/////                //goodAk08.push_back(i);
/////                mu_jet_DR=1;
/////                }else{
/////                  mu_jet_DR=0;
/////              }
/////            }//end for 'muon cleaning'
/////            for(int j=0; j<goodElectron.size(); ++j){
/////              electron.SetPtEtaPhiM(el_pt->at(j), el_eta->at(j), el_phi->at(j), el_Mass);
/////              if(electron.DeltaR(ak08)>.8){
/////                std::cout<<"Good ak08"<<std::endl; 
/////                //goodAk08.push_back(i);
/////                el_jet_DR=1;
/////              }else{
/////                el_jet_DR=0;
/////              }
/////            }//end for 'electron cleaning'
/////            if((goodElectron.size()==0 && goodMuon.size()==0) || (mu_jet_DR==1 && el_jet_DR==1)){
              std::cout<<"Good ak08"<<std::endl; 
              goodAk08.push_back(i);
/////            }
        }//end if for good jet ak08

     }//end loop over # of ak8

     for(int i=0; i<jetAK4_N; ++i){
    
        if(jetAK4_pt->at(i)>180 && abs(jetAK4_eta->at(i))<2. && jetAK4_IDLoose->at(i)==1){
            ak04.SetPtEtaPhiM(jetAK4_pt->at(i), jetAK4_eta->at(i), jetAK4_phi->at(i), jetAK4_mass->at(i));
/////            for(int j=0; j<goodMuon.size(); ++j){
/////              muon.SetPtEtaPhiM(mu_pt->at(j), mu_eta->at(j), mu_phi->at(j), mu_Mass);
/////              if(muon.DeltaR(ak04)>.4){
/////                std::cout<<"Good ak04"<<std::endl; 
/////                //goodAk04.push_back(i);
/////                mu_jet_DR=1;
/////              }else{
/////                mu_jet_DR=0;
/////              }
/////            }//end for 'muon cleaning'
/////            for(int j=0; j<goodElectron.size(); ++j){
/////              electron.SetPtEtaPhiM(el_pt->at(j), el_eta->at(j), el_phi->at(j), el_Mass);
/////              if(electron.DeltaR(ak04)>.4){
/////                std::cout<<"Good ak04"<<std::endl; 
/////                //goodAk04.push_back(i);
/////                el_jet_DR=1;
/////              }else{
/////                el_jet_DR=0;
/////              }
/////            }//end for 'electron cleaning'
/////            if((goodElectron.size()==0 && goodMuon.size()==0) || (mu_jet_DR==1 && el_jet_DR==1)){
              std::cout<<"Good ak04"<<std::endl;
              goodAk04.push_back(i);
/////            }
        }//end if for good jet ak04

     }//end loop over # of ak4

     ////Filling outputTree
     if(goodPhoton.size()>0){
       fillVariableWithValue("photon_pt", ph_pt->at(goodPhoton[0]));
       fillVariableWithValue("photon_eta", ph_eta->at(goodPhoton[0]));
       fillVariableWithValue("photon_phi", ph_phi->at(goodPhoton[0]));
       std::cout<<"++++++++++++ "<<ph_pt->at(goodPhoton[0])<<std::endl;
     } 

/////     if(goodElectron.size()>0){
/////       fillVariableWithValue("electron_pt", el_pt->at(goodElectron[0]));
/////       fillVariableWithValue("electron_eta", el_eta->at(goodElectron[0]));
/////       fillVariableWithValue("electron_phi", el_phi->at(goodElectron[0]));
/////     }
/////
/////     if(goodMuon.size()>0){
/////       fillVariableWithValue("muon_pt", mu_pt->at(goodMuon[0]));
/////       fillVariableWithValue("muon_eta", mu_eta->at(goodMuon[0]));
/////       fillVariableWithValue("muon_phi", mu_phi->at(goodMuon[0]));
/////     }

     if(goodAk08.size()>0){
       fillVariableWithValue("ak08_pt", jetAK8_pt->at(goodAk08[0]));
       fillVariableWithValue("ak08_eta", jetAK8_eta->at(goodAk08[0]));
       fillVariableWithValue("ak08_phi", jetAK8_phi->at(goodAk08[0]));
       fillVariableWithValue("ak08_mass", jetAK8_mass->at(goodAk08[0]));
       fillVariableWithValue("ak08_prunedMass", jetAK8_pruned_mass->at(goodAk08[0]));
       fillVariableWithValue("ak08_CSV", jetAK8_csv->at(goodAk08[0]));
       fillVariableWithValue("ak08_tau21", (jetAK8_tau2->at(goodAk08[0]))/(jetAK8_tau1->at(goodAk08[0])));
     }

     if(goodAk04.size()>0){
       fillVariableWithValue("ak04_pt", jetAK4_pt->at(goodAk04[0]));
       fillVariableWithValue("ak04_eta", jetAK4_eta->at(goodAk04[0]));
       fillVariableWithValue("ak04_phi", jetAK4_phi->at(goodAk04[0]));
       fillVariableWithValue("ak04_mass", jetAK4_mass->at(goodAk04[0]));
     }

     if(goodPhoton.size()>0 && goodAk08.size()>0){
       ak08.SetPtEtaPhiM(jetAK8_pt->at(goodAk08[0]), jetAK8_eta->at(goodAk08[0]), jetAK8_phi->at(goodAk08[0]), jetAK8_mass->at(goodAk08[0]));
       photon.SetPtEtaPhiM(ph_pt->at(goodPhoton[0]), ph_eta->at(goodPhoton[0]), ph_phi->at(goodPhoton[0]), 0);
       ak08_photon = ak08 + photon;
       fillVariableWithValue("ak08_photon_pt", ak08_photon.Pt());
       fillVariableWithValue("ak08_photon_eta", ak08_photon.Eta());
       fillVariableWithValue("ak08_photon_phi", ak08_photon.Phi());
       fillVariableWithValue("ak08_photon_mass", ak08_photon.M());
     }
     if(goodPhoton.size()>0 && goodAk04.size()>0){
       ak04.SetPtEtaPhiM(jetAK4_pt->at(goodAk04[0]), jetAK4_eta->at(goodAk04[0]), jetAK4_phi->at(goodAk04[0]), jetAK4_mass->at(goodAk04[0]));
       photon.SetPtEtaPhiM(ph_pt->at(goodPhoton[0]), ph_eta->at(goodPhoton[0]), ph_phi->at(goodPhoton[0]), 0);
       ak04_photon = ak04 + photon;
       fillVariableWithValue("ak04_photon_pt", ak04_photon.Pt());
       fillVariableWithValue("ak04_photon_eta", ak04_photon.Eta());
       fillVariableWithValue("ak04_photon_phi", ak04_photon.Phi());
       fillVariableWithValue("ak04_photon_mass", ak04_photon.M());
     }
     fillVariableWithValue("metRaw", METraw_sumEt->at(0));
     fillVariableWithValue("met", MET_sumEt->at(0));

     evaluateCuts();
     fillReducedSkimTree(); 
     if( passedCut("ak04_photon_mass")){
        //fillReducedSkimTree();  
     }
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
