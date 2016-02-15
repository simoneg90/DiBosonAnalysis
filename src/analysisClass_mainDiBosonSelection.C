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

#define b_mass 4.70 //GeV
#define e_mass 0.0005 //GeV
#define mu_mass 0.1 //GeV

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

 
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
  int problema=0; 
  if (fChain == 0) return;
 
  TStopwatch time;
  time.Start(true);
  
  setTDRStyle();
  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   
  TLorentzVector b, anti_b, ak4, jet1, jet2, mu, nu_mu, e, nu_e, tmp, nu_tmp, H, W_plus, W_minus, W, tau, nu_tau;
  bool isB, isAntiB, isE, isNuE, isMu, isNuMu, isTau, isNuTau;
  isB=isAntiB=isE=isNuE=isMu=isNuMu=isTau=isNuTau=0;
  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////    
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry == 0 || jentry%100 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
    // if (Cut(ientry) < 0) continue;
    
    ////////////////////// User's code starts here ///////////////////////
    
    ///Stuff to be done for every event
    
    resetCuts();

    frame("Hadronic part");
    std::cout<<"We have "<<nGenBQuarkFromH<<" b"<<std::endl;
    for(int i=0; i<nGenBQuarkFromH;++i){
      if(GenBQuarkFromH_pdgId[i]==5){
        CreateAndFillUserTH1F("bGen_pt",440,0,2200,GenBQuarkFromH_pt[i]);
        CreateAndFillUserTH1F("bGen_eta",100,-5,5,GenBQuarkFromH_eta[i]);
        CreateAndFillUserTH1F("bGen_phi",100, -3.14, 3.14,GenBQuarkFromH_phi[i]);
        b.SetPtEtaPhiM(GenBQuarkFromH_pt[i],GenBQuarkFromH_eta[i],GenBQuarkFromH_phi[i],GenBQuarkFromH_mass[i]);//b_mass);
        isB=1;
      }else if(GenBQuarkFromH_pdgId[i]==-5){
        CreateAndFillUserTH1F("AntibGen_pt",440,0,2200,GenBQuarkFromH_pt[i]);
        CreateAndFillUserTH1F("AntibGen_eta",100,-5,5,GenBQuarkFromH_eta[i]);
        CreateAndFillUserTH1F("AntibGen_phi",100, -3.14, 3.14,GenBQuarkFromH_phi[i]);
        anti_b.SetPtEtaPhiM(GenBQuarkFromH_pt[i],GenBQuarkFromH_eta[i],GenBQuarkFromH_phi[i],GenBQuarkFromH_mass[i]);//b_mass);
        isAntiB=1;
      }//end if b or anti-b
    }//loop over GenB
    if(isB>0 && isAntiB>0){
      CreateAndFillUserTH1F("Higgs_pt_GenLevel",500,0,2500,(b+anti_b).Pt());
      CreateAndFillUserTH1F("Higgs_M_GenLevel",8000,100,200,(b+anti_b).M());
      CreateAndFillUserTH1F("Higgs_eta_GenLevel",100,-5,5,(b+anti_b).Eta());
      CreateAndFillUserTH1F("Higgs_phi_GenLevel",100, -3.14, 3.14,(b+anti_b).Phi());
    }
    std::cout<<"We have "<<nGenHiggsBoson<<" Higgs Bosons O_o"<<std::endl;
    for(int i=0; i<nGenHiggsBoson;++i){
      H.SetPtEtaPhiM(GenHiggsBoson_pt[i], GenHiggsBoson_eta[i], GenHiggsBoson_phi[i], GenHiggsBoson_mass[i]);
      CreateAndFillUserTH1F("Higgs_pt",500,0,2500, GenHiggsBoson_pt[i]);
      CreateAndFillUserTH1F("Higgs_M",8000,100,200,GenHiggsBoson_mass[i]);
      CreateAndFillUserTH1F("Higgs_eta",100,-5,5,GenHiggsBoson_eta[i]);
      CreateAndFillUserTH1F("Higgs_phi",100, -3.14, 3.14,GenHiggsBoson_phi[i]);
      //=== Comparison part ===
      CreateAndFillUserTH2F("H_Comparison_pt",500,0,2500,500,0,2500,(b+anti_b).Pt(),GenHiggsBoson_pt[i]);
      CreateAndFillUserTH2F("H_Comparison_eta",100,-5,5,100,-5,5,(b+anti_b).Eta(),GenHiggsBoson_eta[i]);
      CreateAndFillUserTH2F("H_Comparison_phi",100, -3.14, 3.14,100, -3.14, 3.14,(b+anti_b).Phi(),GenHiggsBoson_phi[i]);

    }//end loop over # of Higgs...
    std::cout<<"We have "<<nGenJet<<" jets"<<std::endl;
    for(int i=0; i<nGenJet; ++i){
      CreateAndFillUserTH1F("GenJet_pt",320,0,1600, GenJet_pt[i]);
      CreateAndFillUserTH1F("GenJet_eta",100,-5,5,GenJet_eta[i]);
      CreateAndFillUserTH1F("GenJet_phi",100, -3.14, 3.14,GenJet_phi[i]);
    }//end loop over GenJet

    frame("Lepton part");
    std::cout<<"We have "<<nGenHiggsSisters<<" Higgs' Sisters"<<std::endl;   
    //W and Lepton part
    
    std::cout<<"We have "<<nGenVbosons<<" VBosons"<<std::endl;
    for(int i=0; i<nGenVbosons; ++i){
      CreateAndFillUserTH1F("Boson_PdgId",100,-50.5,50.5,GenVbosons_pdgId[i]);
      CreateAndFillUserTH1F("Boson_pt",400,0,2000,GenHiggsSisters_pt[i]);
      CreateAndFillUserTH1F("Boson_eta",100,-5,5,GenHiggsSisters_eta[i]);
      CreateAndFillUserTH1F("Boson_phi",100, -3.14, 3.14,GenHiggsSisters_phi[i]);
      
      if(abs(GenVbosons_pdgId[i])==24){
        if(GenVbosons_pdgId[i]==24) {
          W_plus.SetPtEtaPhiM(GenVbosons_pt[i], GenVbosons_eta[i], GenVbosons_phi[i], GenVbosons_mass[i]);
          W.SetPtEtaPhiM(GenVbosons_pt[i], GenVbosons_eta[i], GenVbosons_phi[i], GenVbosons_mass[i]);
          frame("X+ Study");

          CreateAndFillUserTH1F("XPlus_pt",400,0,2000,(H+W_plus).Pt());
          CreateAndFillUserTH1F("XPlus_eta",100,-5,5,(H+W_plus).Eta());
          CreateAndFillUserTH1F("XPlus_phi",100, -3.14, 3.14,(H+W_plus).Phi());
          CreateAndFillUserTH1F("XPlus_M",4000, 975, 2000,(H+W_plus).M());

          //std::cout<<"------------- "<<(H+W_plus).M()<<"------------------"<<std::endl;
        }

        if(GenVbosons_pdgId[i]==-24){
          W_minus.SetPtEtaPhiM(GenVbosons_pt[i], GenVbosons_eta[i], GenVbosons_phi[i], GenVbosons_mass[i]);
          W.SetPtEtaPhiM(GenVbosons_pt[i], GenVbosons_eta[i], GenVbosons_phi[i], GenVbosons_mass[i]);
          frame("X- Study");

          CreateAndFillUserTH1F("XMinus_pt",800,0,4000,(H+W_minus).Pt());
          CreateAndFillUserTH1F("XMinus_eta",100,-5,5,(H+W_minus).Eta());
          CreateAndFillUserTH1F("XMinus_phi",100, -3.14, 3.14,(H+W_minus).Phi());
          CreateAndFillUserTH1F("XMinus_M",4000, 975, 2000,(H+W_minus).M());
        
        }
        CreateAndFillUserTH1F("X_pt",800,0,4000,(H+W).Pt());
        CreateAndFillUserTH1F("X_eta",100,-5,5,(H+W).Eta());
        CreateAndFillUserTH1F("X_phi",100, -3.14, 3.14,(H+W).Phi());
        CreateAndFillUserTH1F("X_M",4000, 975, 2000,(H+W).M());
 
      }//end if for W+-
      
    }//end loop over nGenVBoson 

    //=== Starting loop over taus ===

    for(int i=0; i<nGenTau; ++i){
        tau.SetPtEtaPhiM(GenTau_pt[i], GenTau_eta[i], GenTau_phi[i], GenTau_mass[i]);
        nu_tau.SetPtEtaPhiM(GenNu_pt[i], GenNu_eta[i], GenNu_phi[i], 0);
        CreateAndFillUserTH1F("Tau_pt",440,0,2200,GenTau_pt[i]);
        CreateAndFillUserTH1F("Tau_eta",100,-5,5,GenTau_eta[i]);
        CreateAndFillUserTH1F("Tau_phi",100, -3.14, 3.14,GenTau_phi[i]);
        CreateAndFillUserTH1F("NuTau_pt",440,0,2200,GenNu_pt[i]);
        CreateAndFillUserTH1F("NuTau_eta",100,-5,5,GenNu_eta[i]);
        CreateAndFillUserTH1F("NuTau_phi",100, -3.14, 3.14,GenNu_phi[i]);
        tmp=tau;
        nu_tmp=nu_tau;
        isTau=1;
        isNuTau=1;
    }
    if(isTau==1 && isNuTau==1){
        CreateAndFillUserTH1F("Boson_pt_Tau",500,0,2500,(tau+nu_tau).Pt());
        CreateAndFillUserTH1F("Boson_eta_Tau",100,-5,5,(tau+nu_tau).Eta());
        CreateAndFillUserTH1F("Boson_phi_Tau",100, -3.14, 3.14,(tau+nu_tau).Phi());
        CreateAndFillUserTH1F("Boson_M_Tau",200, 60, 100,(tau+nu_tau).M());
        //isTau=isNuTau=0;
    }

    for(int i=0; i<nGenLep; ++i){
       if(abs(GenLep_pdgId[i])==11 && (GenLep_pdgId[i]*GenVbosons_pdgId[i]<0)){
         e.SetPtEtaPhiM(GenLep_pt[i],GenLep_eta[i],GenLep_phi[i],e_mass);
         CreateAndFillUserTH1F("Electron_pt",440,0,2200,GenLep_pt[i]);
         CreateAndFillUserTH1F("Electron_eta",100,-5,5,GenLep_eta[i]);
         CreateAndFillUserTH1F("Electron_phi",100, -3.14, 3.14,GenLep_phi[i]);
         tmp=e;
         isE=1;
         CreateAndFillUserTH1F("NuE_pt_Met",440,0,2200,met_genPt);
         CreateAndFillUserTH1F("NuE_eta_Met",100,-5,5,met_genEta);
         CreateAndFillUserTH1F("NuE_phi_Met",100, -3.14, 3.14,met_genPhi);
         //attenzione
         if(i<nGenNu && abs(GenNu_pdgId[i])==12 && (GenNu_pdgId[i]*GenLep_pdgId[i])<0){
           std::cout<<"I have a e neutrino!"<<std::endl;
           CreateAndFillUserTH1F("NuE_pt",440,0,2200,GenNu_pt[i]);
           CreateAndFillUserTH1F("NuE_eta",100,-5,5,GenNu_eta[i]);
           CreateAndFillUserTH1F("NuE_phi",100, -3.14, 3.14,GenNu_phi[i]);
           nu_e.SetPtEtaPhiM(GenNu_pt[i],GenNu_eta[i],GenNu_phi[i],0);
           nu_tmp=nu_e;
           isNuE=1;
         }
       }else if (abs(GenLep_pdgId[i])==13 && (GenLep_pdgId[i]*GenVbosons_pdgId[i]<0)){
         mu.SetPtEtaPhiM(GenLep_pt[i],GenLep_eta[i],GenLep_phi[i],mu_mass);
         CreateAndFillUserTH1F("Muon_pt",440,0,2200,GenLep_pt[i]);
         CreateAndFillUserTH1F("Muon_eta",100,-5,5,GenLep_eta[i]);
         CreateAndFillUserTH1F("Muon_phi",100, -3.14, 3.14,GenLep_phi[i]);
         tmp=mu;
         isMu=1;
         CreateAndFillUserTH1F("NuMu_pt_Met",440,0,2200,met_genPt);
         CreateAndFillUserTH1F("NuMu_eta_Met",100,-5,5,met_genEta);
         CreateAndFillUserTH1F("NuMu_phi_Met",100, -3.14, 3.14,met_genPhi);
         if(i<nGenNu && abs(GenNu_pdgId[i])==14 && (GenNu_pdgId[i]*GenLep_pdgId[i])<0){
           CreateAndFillUserTH1F("NuMu_pt",440,0,2200,GenNu_pt[i]);
           CreateAndFillUserTH1F("NuMu_eta",100,-5,5,GenNu_eta[i]);
           CreateAndFillUserTH1F("NuMu_phi",100, -3.14, 3.14,GenNu_phi[i]);
           nu_mu.SetPtEtaPhiM(GenNu_pt[i],GenNu_eta[i],GenNu_phi[i],0);
           nu_tmp=nu_mu;
           isNuMu=1;
         }

       }//end else if muons...

     
       if(isE==1 && isNuE==1){
         CreateAndFillUserTH1F("Boson_pt_El",500,0,2500,(e+nu_e).Pt());
         CreateAndFillUserTH1F("Boson_eta_El",100,-5,5,(e+nu_e).Eta());
         CreateAndFillUserTH1F("Boson_phi_El",100, -3.14, 3.14,(e+nu_e).Phi());
         CreateAndFillUserTH1F("Boson_M_El",200, 60, 100,(e+nu_e).M());
         
       }
       if(isMu==1 && isNuMu==1){
         CreateAndFillUserTH1F("Boson_pt_Mu",500,0,2500,(mu+nu_mu).Pt());
         CreateAndFillUserTH1F("Boson_eta_Mu",100,-5,5,(mu+nu_mu).Eta());
         CreateAndFillUserTH1F("Boson_phi_Mu",100, -3.14, 3.14,(mu+nu_mu).Phi());
         CreateAndFillUserTH1F("Boson_M_Mu",200, 60, 100,(mu+nu_mu).M());
       }
       
     }//end loop over nGenLep
     if((isE==1 && isNuE==1) || (isMu==1 && isNuMu==1) || (isTau==1 && isNuTau==1)){
       CreateAndFillUserTH1F("Boson_pt_TauMuEl",500,0,2500,(tmp+nu_tmp).Pt());
       CreateAndFillUserTH1F("Boson_eta_TauMuEl",100,-5,5,(tmp+nu_tmp).Eta());
       CreateAndFillUserTH1F("Boson_phi_TauMuEl",100, -3.14, 3.14,(tmp+nu_tmp).Phi());
       CreateAndFillUserTH1F("Boson_M_TauMuEl",200, 60, 100,(tmp+nu_tmp).M());
       isE=isNuE=0;
       isMu=isNuMu=0;
       isTau=isNuTau=0;
     }
     for(int i=0; i<nGenHiggsSisters;++i){
        CreateAndFillUserTH1F("HiggsSister_pt",500,0,2500,GenHiggsSisters_pt[i]);
        CreateAndFillUserTH1F("HiggsSister_eta",100,-5,5,GenHiggsSisters_eta[i]);
        CreateAndFillUserTH1F("HiggsSister_phi",100, -3.14, 3.14,GenHiggsSisters_phi[i]);
        //=== Comparison Part ===
        if((GenVbosons_pt[i]/(tmp+nu_tmp).Pt())> 1.001 || (GenVbosons_pt[i]/(tmp+nu_tmp).Pt())<0.999){
          std::cout<<"Ratio: "<<(GenVbosons_pt[i]/(tmp+nu_tmp).Pt())<<std::endl;
          ++problema;
        }
        CreateAndFillUserTH2F("HiggsSister_Comparison_pt",800,0,4000,2000,0,4000,(tmp+nu_tmp).Pt(),GenVbosons_pt[i]);//GenHiggsSisters_pt[i]);
        CreateAndFillUserTH2F("HiggsSister_Comparison_phi",100, -3.14, 3.14,100, -3.14, 3.14,(tmp+nu_tmp).Phi(),GenHiggsSisters_phi[i]);
        CreateAndFillUserTH2F("HiggsSister_Comparison_eta",100,-5,5,100,-5,5,(tmp+nu_tmp).Eta(),GenHiggsSisters_eta[i]);
     }
  }//end loop over entries
  std::cout<<"problema: "<<problema<<std::endl;
  std::cout << "analysisClass::Loop() ends" <<std::endl;   
  std::cout<<""<<std::endl;
  std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << time.RealTime() << " s" <<std::endl;
  std::cout<<""<<std::endl;

} 
