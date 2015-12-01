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

#define b_mass 4.20 //GeV

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
  
  if (fChain == 0) return;
  
  
  setTDRStyle();
  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   
  
  TLorentzVector b, anti_b, bReco, anti_bReco, jet1, jet2, mu, nu_mu;
  bool isB, isAntiB;
  isB=isAntiB=0;
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
    for(int i=0; i<sizeof(GenBQuarkFromH_pdgId);++i){
      std::cout<<"We have "<<sizeof(GenBQuarkFromH_pdgId)<<" b"<<std::endl;
      if(GenBQuarkFromH_pdgId[i]==5){
        CreateAndFillUserTH1F("bGen_pt",2000,0,4000,GenBQuarkFromH_pt[i]);
        CreateAndFillUserTH1F("bGen_eta",100,-5,5,GenBQuarkFromH_eta[i]);
        CreateAndFillUserTH1F("bGen_phi",100, -3.14, 3.14,GenBQuarkFromH_phi[i]);
        b.SetPtEtaPhiM(GenBQuarkFromH_pt[i],GenBQuarkFromH_eta[i],GenBQuarkFromH_phi[i],b_mass);
        isB=1;
      }else if(GenBQuarkFromH_pdgId[i]==-5){
        CreateAndFillUserTH1F("AntibGen_pt",2000,0,4000,GenBQuarkFromH_pt[i]);
        CreateAndFillUserTH1F("AntibGen_eta",100,-5,5,GenBQuarkFromH_eta[i]);
        CreateAndFillUserTH1F("AntibGen_phi",100, -3.14, 3.14,GenBQuarkFromH_phi[i]);
        anti_b.SetPtEtaPhiM(GenBQuarkFromH_pt[i],GenBQuarkFromH_eta[i],GenBQuarkFromH_phi[i],b_mass);
        isAntiB=1;
      }//end if b or anti-b
    }//loop over GenB
    if(isB>0 && isAntiB>0){
      CreateAndFillUserTH1F("Higgs_pt_GenLevel",1000,0,2000,(b+anti_b).Pt());
      CreateAndFillUserTH1F("Higgs_M_GenLevel",8000,100,200,(b+anti_b).M());
      CreateAndFillUserTH1F("Higgs_eta_GenLevel",100,-5,5,(b+anti_b).Eta());
    }
    for(int i=0; i<sizeof(GenHiggsBoson_pdgId);++i){
      std::cout<<"We have "<<sizeof(GenHiggsBoson_pdgId)<<" Higgs Bosons O_o"<<std::endl;
      CreateAndFillUserTH1F("Higgs_pt",1000,0,2000, GenHiggsBoson_pt[i]);
      CreateAndFillUserTH1F("Higgs_M",8000,100,200,GenHiggsBoson_mass[i]);
      CreateAndFillUserTH1F("Higgs_eta",100,-5,5,GenHiggsBoson_eta[i]);
      CreateAndFillUserTH1F("Higgs_phi",100, -3.14, 3.14,GenHiggsBoson_phi[i]);
    }//end loop over # of Higgs...


    //W and Lepton part
    for(int i=0; i<nGenHiggsSisters;++i){
      std::cout<<"We have "<<nGenHiggsSisters<<" Higgs' Sisters"<<std::endl;
      CreateAndFillUserTH1F("HiggsSister_pt",1000,0,2000,GenHiggsSisters_pt[i]);
      CreateAndFillUserTH1F("HiggsSister_eta",100,-5,5,GenHiggsSisters_eta[i]);
      CreateAndFillUserTH1F("HiggsSister_phi",100, -3.14, 3.14,GenHiggsSisters_phi[i]);
    }

    for(int i=0; i<nGenVbosons; ++i){
      std::cout<<"We have "<<nGenVbosons<<" VBosons"<<std::endl;
      CreateAndFillUserTH1F("Boson_PdgId",100,-50.5,50.5,GenVbosons_pdgId[i]);
      CreateAndFillUserTH1F("Boson_pt",1000,0,2000,GenHiggsSisters_pt[i]);
      CreateAndFillUserTH1F("Boson_eta",100,-5,5,GenHiggsSisters_eta[i]);
      CreateAndFillUserTH1F("Boson_phi",100, -3.14, 3.14,GenHiggsSisters_phi[i]);
    }







  }
  
  
  
  std::cout << "analysisClass::Loop() ends" <<std::endl;   
} 
