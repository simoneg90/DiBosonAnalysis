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
#include "TStopwatch.h"

#define e_mass 0.0005 //GeV
#define mu_mass 0.1 //GeV

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
   std::vector <int> goodLepton, looseLepton, goodAK08, goodAK04, goodAK08_lep, goodAK04_lep, goodAK08Pruned;
   TLorentzVector ak04, ak08, ak08Pruned, lepton, leptonLoose, W, MET;
   int btag_ak04;
   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<1000;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done for every event
     btag_ak04=0;
     goodLepton.clear();
     looseLepton.clear();
     goodAK04.clear();
     goodAK08.clear();
     goodAK08_lep.clear();
     goodAK04_lep.clear();
     goodAK08Pruned.clear();

     resetCuts();
   

//=== Forse da inserire ===

/*     fillVariableWithValue("run",runNo);     
     fillVariableWithValue("event",evtNo);     
     fillVariableWithValue("lumi",lumi);     
     fillVariableWithValue("nVtx",nvtx);     
     fillVariableWithValue("nJet",widejets.size());
*/
    

    if(nselLeptons>0){
      if((abs(selLeptons_pdgId[0])==11 && selLeptons_isMyGoodElectron[0]==1 && selLeptons_pt[0]>30 && (abs(selLeptons_eta[0])<1.442 || (abs(selLeptons_eta[0])>1.56 && abs(selLeptons_eta[0])<2.5))) || (abs(selLeptons_pdgId[0])==13 && selLeptons_isMyGoodMuon[0]==1 && selLeptons_pt[0]>20 && abs(selLeptons_eta[0])<2.1 && selLeptons_relIso03[0]<0.1)){
        goodLepton.push_back(0);
        std::cout<<"Passing tight selection"<<std::endl;
        CreateAndFillUserTH1D("goodEleTightSelection", 2,-.5,1.5, 1);
      }//end if 'good' lepton cuts
    }//end if nselLeptons

    if(nFatjetAK08ungroomed>0){
      if((((FatjetAK08ungroomed_neHEFrac[0]<0.99 && FatjetAK08ungroomed_neEmEFrac[0]<0.99 && (FatjetAK08ungroomed_chMult[0]+FatjetAK08ungroomed_neMult[0])>1) && ((abs(FatjetAK08ungroomed_eta[0])<=2.4 && FatjetAK08ungroomed_chHEFrac[0]>0 && FatjetAK08ungroomed_chMult[0]>0 && FatjetAK08ungroomed_chEmEFrac[0]<0.99) || abs(FatjetAK08ungroomed_eta[0])>2.4) && abs(FatjetAK08ungroomed_eta[0])<=3.0)||((FatjetAK08ungroomed_neEmEFrac[0]<0.90 && FatjetAK08ungroomed_neMult[0]>10 && abs(FatjetAK08ungroomed_eta[0])>3.0 ))) && FatjetAK08ungroomed_pt[0]>100 && abs(FatjetAK08ungroomed_eta[0])<2.4){
        goodAK08.push_back(0);
      }//end if good AK08
    }//end loop over nFatjetAK08ungroomed

    for(int i=0; i<nJet; ++i){
      if((((Jet_neHEF[i]<0.99 && Jet_neEmEF[i]<0.99 && (Jet_chMult[i]+Jet_neMult[i])>1) && ((abs(Jet_eta[i])<=2.4 && Jet_chHEF[i]>0 && Jet_chMult[i]>0 && Jet_chEmEF[i]<0.99) || abs(Jet_eta[i])>2.4) && abs(Jet_eta[i])<=3.0)||((Jet_neEmEF[i]<0.90 && Jet_neMult[i]>10 && abs(Jet_eta[i])>3.0 ))) && Jet_pt[i]>30 && abs(Jet_eta[i])<2.4){
        CreateAndFillUserTH1D("goodAk04LooseSelection", 2,-.5,1.5, 1);
        goodAK04.push_back(i);
      }
    }//end loop over nJet

    fillVariableWithValue("lepton_goodNumber", goodLepton.size());
    fillVariableWithValue("ak08_goodNumber", goodAK08.size());
      if(goodAK08.size()==1){
        fillVariableWithValue("ak08Ungroomed_1_pt", FatjetAK08ungroomed_pt[goodAK08[0]]);
        fillVariableWithValue("ak08Ungroomed_1_eta", FatjetAK08ungroomed_eta[goodAK08[0]]);
        fillVariableWithValue("ak08Ungroomed_1_phi", FatjetAK08ungroomed_phi[goodAK08[0]]);
        fillVariableWithValue("ak08Ungroomed_1_mass", FatjetAK08ungroomed_mass[goodAK08[0]]);
        fillVariableWithValue("ak08Ungroomed_1_tau21", FatjetAK08ungroomed_tau2[goodAK08[0]]/FatjetAK08ungroomed_tau1[goodAK08[0]]);
        ak08.SetPtEtaPhiM(FatjetAK08ungroomed_pt[goodAK08[0]], FatjetAK08ungroomed_eta[goodAK08[0]], FatjetAK08ungroomed_phi[goodAK08[0]], FatjetAK08ungroomed_mass[goodAK08[0]]);
        double dr_tmp=99;
        int index=0;
        for(int ii=0; ii<nFatjetAK08pruned; ++ii){
          ak08Pruned.SetPtEtaPhiM(FatjetAK08pruned_pt[ii],FatjetAK08pruned_eta[ii], FatjetAK08pruned_phi[ii], FatjetAK08pruned_mass[ii]);
          if(ak08.DeltaR(ak08Pruned)<dr_tmp){
            dr_tmp=ak08.DeltaR(ak08Pruned);
            index=ii;
          }//matched good pruned AK08
          goodAK08Pruned.push_back(index);
        }//end loop over nFatjetAK08pruned

      }

      if(goodLepton.size()==1){
        if(abs(selLeptons_pdgId[goodLepton[0]])==11) lepton_mass=e_mass;
        if(abs(selLeptons_pdgId[goodLepton[0]])==13) lepton_mass=mu_mass;
        lepton.SetPtEtaPhiM(selLeptons_pt[goodLepton[0]],selLeptons_eta[goodLepton[0]],selLeptons_phi[goodLepton[0]], lepton_mass);
        fillVariableWithValue("lepton_pt", selLeptons_pt[goodLepton[0]]);
        fillVariableWithValue("lepton_eta", selLeptons_eta[goodLepton[0]]);
        fillVariableWithValue("lepton_phi", selLeptons_phi[goodLepton[0]]);
        fillVariableWithValue("lepton_pdgID", selLeptons_pdgId[goodLepton[0]]);
        MET.SetPtEtaPhiM(met_pt, met_eta, met_phi, 0);
        W=lepton+MET;
        fillVariableWithValue("W_pt", W.Pt());
        fillVariableWithValue("W_eta",W.Eta());
        fillVariableWithValue("W_phi",W.Phi());
        fillVariableWithValue("W_mT",W.Mt());
        MET.SetPtEtaPhiM(metType1p2_pt, met_eta, met_phi, 0);
        W=lepton+MET;
        fillVariableWithValue("WType1_pt", W.Pt());
        fillVariableWithValue("WType1_eta",W.Eta());
        fillVariableWithValue("WType1_phi",W.Phi());
        fillVariableWithValue("WType1_mT",W.Mt());
        MET.SetPtEtaPhiM(metPuppi_pt, met_eta, met_phi, 0); 
        W=lepton+MET;
        fillVariableWithValue("WPuppi_pt", W.Pt());
        fillVariableWithValue("WPuppi_eta",W.Eta());
        fillVariableWithValue("WPuppi_phi",W.Phi());
        fillVariableWithValue("WPuppi_mT",W.Mt());

        if(goodAK08.size()==1) fillVariableWithValue("ak08Ungroomed_lepton_DR", lepton.DeltaR(ak08));
          for(int j=0; j<goodAK04.size(); ++j){
            ak04.SetPtEtaPhiM(Jet_pt[goodAK04[j]], Jet_eta[goodAK04[j]], Jet_phi[goodAK04[j]], Jet_mass[goodAK04[j]]);
            if(ak04.DeltaR(ak08)>.8 && ak04.DeltaR(lepton)>.3){
              CreateAndFillUserTH1D("Ak04_lepton&AK08_DRCut", 2,-.5,1.5, 1);
              goodAK04_lep.push_back(goodAK04[j]);
              if(Jet_btagCSV[goodAK04[j]]>0.9){
                ++btag_ak04;
              }//end if for CSV medium working point
            }//end if ak04 without leptons and ak08 nearby
          }//end loop over goodAK04.size()

        }//end loop over goodAK08.size()
        fillVariableWithValue("btag",btag_ak04);
        if(goodAK08Pruned.size()>=1){
          fillVariableWithValue("ak08Pruned_1_pt", FatjetAK08pruned_pt[goodAK08Pruned[0]]);
          fillVariableWithValue("ak08Pruned_1_eta", FatjetAK08pruned_eta[goodAK08Pruned[0]]);
          fillVariableWithValue("ak08Pruned_1_phi", FatjetAK08pruned_phi[goodAK08Pruned[0]]);
          fillVariableWithValue("ak08Pruned_1_mass", FatjetAK08pruned_mass[goodAK08Pruned[0]]);
        }
        if(goodAK04_lep.size()>=1){
          fillVariableWithValue("ak04_1_pt", Jet_pt[goodAK04_lep[0]]);
          fillVariableWithValue("ak04_1_eta", Jet_eta[goodAK04_lep[0]]);
          fillVariableWithValue("ak04_1_phi", Jet_phi[goodAK04_lep[0]]);
          fillVariableWithValue("ak04_1_mass", Jet_mass[goodAK04_lep[0]]);
          //=== TH1D to check the fillReduceSkim procedure ===
          CreateAndFillUserTH1D("ak04_first_pt", 1000,0,500, Jet_pt[goodAK04_lep[0]]);
        }
        if(goodAK04_lep.size()>=2){
          fillVariableWithValue("ak04_2_pt", Jet_pt[goodAK04_lep[1]]);
          fillVariableWithValue("ak04_2_eta", Jet_eta[goodAK04_lep[1]]);
          fillVariableWithValue("ak04_2_phi", Jet_phi[goodAK04_lep[1]]);
          fillVariableWithValue("ak04_2_mass", Jet_mass[goodAK04_lep[1]]);
        }
        fillVariableWithValue("metPuppi",metPuppi_pt);
        fillVariableWithValue("metType1", metType1p2_pt);
        fillVariableWithValue("met",met_pt);
        fillVariableWithValue("nPrimaryVertexes", nprimaryVertices);

        



     // Evaluate cuts (but do not apply them)
     evaluateCuts();
     //fillReducedSkimTree(); 

     // optional call to fill a skim with the full content of the input roottuple
     //if( passedCut("nJetFinal") ) fillSkimTree();
     
     // optional call to fill a skim with a subset of the variables defined in the cutFile (use flag SAVE)
     //if( passedAllPreviousCuts("mjj") && passedCut("mjj") ) 
     if( passedAllPreviousCuts("btag") && passedCut("btag"))
       {
         frame("Beware! This part can be set outside the if -Passed cuts-");
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
