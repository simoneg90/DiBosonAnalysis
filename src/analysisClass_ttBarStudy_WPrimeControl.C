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
#define mu_mass 0.105 //GeV

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;
  std::cout<<"List: "<<inputList_<<std::endl; 
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
 
  TStopwatch time;
  frame("Starting the analysis");
  time.Start(true);
  
  setTDRStyle();
  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////    
  
  Long64_t nbytes = 0, nb = 0;

  std::cout<<"Prova counter: "<<getCrabProcessedEvents()<<std::endl;
  int ele, ele_new, mu;
  ele=ele_new=mu=0;
  std::vector<int> goodLepton, goodJet; //vectors containing the indexes of 'good' particles
  //double varBin[]= {-.5, 49.5, 399.5, 599.5, 799.5, 999.5, 1299.5, 1449.5, 1749.5, 1949.5, 2149.5, 2399.5, 2599.5, 2799.5, 3049.5};
  double varBin[]={-.5,249.5, 489.5, 649.5, 849.5, 1149.5, 1349.5, 1489.5, 1599.5, 1849.5, 2049.5, 2249.5, 2449.5,2699.5, 2949.5,3149.5, 3249.5 };
  //double varBin[15]={-.5,29.5,49.5, 149.5, 349.5, 549.5, 749.5, 949.5,1149.5, 1349.5, 1459.5, 1649.5,1849.5, 1949.5, 2149.5};
  //double varBin1[25]={-.5,29.5,49.5, 149.5, 249.5, 349.5, 449.5, 549.5, 649.5, 749.5, 849.5, 949.5, 1049.5, 1149.5, 1249.5, 1349.5, 1449.5, 1459.5, 1549.5, 1649.5, 1749.5, 1849.5, 1859.5, 1949.5, 2049.5};


  TLorentzVector genEle, recoEle, goodRecoEle, genMuon, recoMuon, goodRecoMuon; //vectors for efficiences test
  int n_eleGen, n_eleReco, n_eleRecoGood;
  n_eleGen=n_eleReco=n_eleRecoGood=0;
  double tmp;
  std::vector<int> goodElectrons, fakeElectrons, goodMuons, fakeMuons, goodGenElectrons, goodGenMuons;
  double dr_tmp; 
  int index_tmp;
  index_tmp=0;
  //std::cout<<"List: "<<inputList_->c_str()<<std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry == 0 || jentry%100 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
    // if (Cut(ientry) < 0) continue;
    
    ////////////////////// User's code starts here ///////////////////////
    
    ///Stuff to be done for every event
    goodElectrons.clear();
    fakeElectrons.clear();
    goodMuons.clear();
    fakeMuons.clear();
    dr_tmp=100;
    resetCuts();
    //pre-cut study
    for(int i=0; i<nselLeptons; ++i){
      
      if(selLeptons_pt[i]>30){
        //=== Implementate solo le ID+Iso degli elettroni ===  => da fare quelle per i muoni
        if((abs(selLeptons_eta[i])<2.4 && abs(selLeptons_pdgId[i])==11) || (abs(selLeptons_eta[i])<2.1 && abs(selLeptons_pdgId[i])==13)){
          CreateAndFillUserTH1D("isEcalDriven", 2,0,1,selLeptons_isEcalDriven[i]);
          CreateAndFillUserTH1D("EleDEta", 100, -1,1, selLeptons_eleDEta[i]);
          CreateAndFillUserTH1D("EleDPhi", 100, -1,1, selLeptons_eleDPhi[i]);
          CreateAndFillUserTH1D("EleHoE", 100, -2,2, selLeptons_eleHoE[i]);
          CreateAndFillUserTH1D("EleFull5x5", 100, -2,2, selLeptons_e2x5Max[i]/selLeptons_e5x5[i]);
          CreateAndFillUserTH1D("EleHadDepth", 100, -2,2, selLeptons_isolEmHadDepth1[i]);
          CreateAndFillUserTH1D("EleTrkPt", 1000, 0,600, selLeptons_isolTrkPt[i]);
          CreateAndFillUserTH1D("EleLostHits", 100, -.5,99.5, selLeptons_eleMissingHits[i]);
          CreateAndFillUserTH1D("EleDxy", 100, -2,2, selLeptons_dxy[i]);
        }//end if eta Cut
      }//end if pt>30GeV

    }


    //efficiency test 

    for(int k=0; k<=((sizeof(varBin)/sizeof(double))-1); ++k){
      tmp=varBin[k]+1;
      for(int i=0; i<nGenLepFromWPrime; ++i){
        //electrons
        if(abs(GenLepFromWPrime_pdgId[i])==11 && abs(GenLepFromWPrime_eta[i])<2.5){
          if((GenLepFromWPrime_pt[i]>varBin[k] && GenLepFromWPrime_pt[i]<varBin[k+1]) || ((k==(sizeof(varBin)/sizeof(double))-1) && GenLepFromWPrime_pt[i]>varBin[k])){
            //good gen electrons
            CreateAndFillUserVariableTH1D("ele_N_", sizeof(varBin)/sizeof(double)-1, varBin, GenLepFromWPrime_pt[i]);
            CreateAndFillUserTH1D("eleGen_pt", 1000,0,2000,GenLepFromWPrime_pt[i]);
            CreateAndFillUserTH1D("eleGen_eta", 100,-5,5,GenLepFromWPrime_eta[i]);
            CreateAndFillUserTH1D("eleGen_phi", 100,-3.16,3.16,GenLepFromWPrime_phi[i]);
          }
          genEle.SetPtEtaPhiM(GenLepFromWPrime_pt[i],GenLepFromWPrime_eta[i],GenLepFromWPrime_phi[i], e_mass);
          for(int j=0; j<nselLeptons; ++j){
            if(abs(selLeptons_pdgId[j])==11){
              recoEle.SetPtEtaPhiM(selLeptons_pt[j],selLeptons_eta[j],selLeptons_phi[j], e_mass);
              CreateAndFillUserTH1D("GenReco_ele_DR", 1100, -.5, 10, genEle.DeltaR(recoEle));
              if(genEle.DeltaR(recoEle)<dr_tmp){
                dr_tmp=genEle.DeltaR(recoEle); 
                index_tmp=j;
              }
              if(genEle.DeltaR(recoEle)>2.) fakeElectrons.push_back(j); //fake electrons->jet
            }
          }//end loop over SelLeptons
          CreateAndFillUserTH1D("eleDRMin", 1100,-.5,10,dr_tmp);
          if(dr_tmp<0.15){
            goodElectrons.push_back(index_tmp);
            goodGenElectrons.push_back(i);
          }
          dr_tmp=100;
        }//end if good gen electrons
        //muons
        if(abs(GenLepFromWPrime_pdgId[i])==13 && abs(GenLepFromWPrime_eta[i])<2.5){
          if((GenLepFromWPrime_pt[i]>varBin[k] && GenLepFromWPrime_pt[i]<varBin[k+1]) || ((k==(sizeof(varBin)/sizeof(double))-1) && GenLepFromWPrime_pt[i]>varBin[k])){
            //good gen muon
            CreateAndFillUserVariableTH1D("muon_N_", sizeof(varBin)/sizeof(double)-1, varBin, GenLepFromWPrime_pt[i]);
            CreateAndFillUserTH1D("muonGen_pt", 1000,0,2000,GenLepFromWPrime_pt[i]);
            CreateAndFillUserTH1D("muonGen_eta", 100,-5,5,GenLepFromWPrime_eta[i]);
            CreateAndFillUserTH1D("muonGen_phi", 100,-3.16,3.16,GenLepFromWPrime_phi[i]);
          }
          genMuon.SetPtEtaPhiM(GenLepFromWPrime_pt[i],GenLepFromWPrime_eta[i],GenLepFromWPrime_phi[i], mu_mass);
          for(int j=0; j<nselLeptons; ++j){
            if(abs(selLeptons_pdgId[j])==13){
              recoMuon.SetPtEtaPhiM(selLeptons_pt[j],selLeptons_eta[j],selLeptons_phi[j], mu_mass);
              CreateAndFillUserTH1D("GenReco_muon_DR", 1100, -.5, 10, genMuon.DeltaR(recoMuon));
              if(genMuon.DeltaR(recoMuon)<dr_tmp){
                dr_tmp=genMuon.DeltaR(recoMuon);
                index_tmp=j;
              }
              if(genMuon.DeltaR(recoMuon)>2.) fakeMuons.push_back(j); //fake muons->jet
            }
          }//end loop over SelLeptons
          if(dr_tmp<0.15){
            goodMuons.push_back(index_tmp);
            goodGenMuons.push_back(i);
          }
          CreateAndFillUserTH1D("muonDRMin", 1100,-.5,10,dr_tmp);
          dr_tmp=100;

        }//end if good gen muon
      }//end GenLep loop
      for(int j=0; j<goodElectrons.size(); ++j){
        CreateAndFillUserTH1D("eleGoodNumber", 10,-.5,9.5, goodElectrons.size());
        //if((GenLep_pt[goodGenElectrons[j]]>varBin[k] && GenLep_pt[goodGenElectrons[j]]<varBin[k+1]) || ((k==(sizeof(varBin)/sizeof(double))-1) && GenLep_pt[goodGenElectrons[j]]>varBin[k])){ 
        if((selLeptons_pt[goodElectrons[j]]>varBin[k] && selLeptons_pt[goodElectrons[j]]<varBin[k+1]) || ((k==(sizeof(varBin)/sizeof(double))-1)&&selLeptons_pt[goodElectrons[j]]>varBin[k])){
          CreateAndFillUserVariableTH1D("ele_RecoPrimo_", sizeof(varBin)/sizeof(double)-1, varBin, selLeptons_pt[goodElectrons[j]]);//GenLep_pt[goodGenElectrons[j]]);
          CreateAndFillUserTH1D("ele_RECO",14,0,3150, selLeptons_pt[goodElectrons[j]]);
          if(selLeptons_isMyGoodElectron[goodElectrons[j]]==1){// && selLeptons_isolEmHadDepth1[goodElectrons[j]]<(2+0.03*selLeptons_eleClusterEnergy[goodElectrons[j]]+0.28*rho)){
            CreateAndFillUserVariableTH1D("ele_recoGoodPrimo_", sizeof(varBin)/sizeof(double)-1, varBin, selLeptons_pt[goodElectrons[j]]);//GenLep_pt[goodGenElectrons[j]]);
            CreateAndFillUserTH1D("ele_GOOD",14,0,3150, selLeptons_pt[goodElectrons[j]]);
            CreateAndFillUserTH1D("ele_ResoPt", 1000,-10,10,GenLep_pt[goodGenElectrons[j]]-selLeptons_pt[goodElectrons[j]]);
            CreateAndFillUserTH1D("ele_Resolution", 1000,-10,10,(selLeptons_pt[goodElectrons[j]]-GenLep_pt[goodGenElectrons[j]])/GenLep_pt[goodGenElectrons[j]]);
          }
        }
        if((selLeptons_pt[goodElectrons[j]]>varBin[k] && selLeptons_pt[goodElectrons[j]]<varBin[k+1]) || ((k==(sizeof(varBin)/sizeof(double))-1) && selLeptons_pt[goodElectrons[j]]>varBin[k])){
          //CreateAndFillUserVariableTH1D("ele_Reco_", sizeof(varBin)/sizeof(double)-1, varBin, tmp);
          
          CreateAndFillUserVariableTH1D("ele_Reco_", sizeof(varBin)/sizeof(double)-1, varBin, selLeptons_pt[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_pt", 1000,0,2000,selLeptons_pt[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_eta", 100,-5,5,selLeptons_eta[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_phi", 100,-3.16,3.16,selLeptons_phi[goodElectrons[j]]);

          //Qui mettere i TH1D con le variabili su cui taglio!
          CreateAndFillUserTH1D("eleReco_ClusterDEta",1000,-5,5,abs(selLeptons_eleClusterDEta[goodElectrons[j]]));
          CreateAndFillUserTH1D("eleReco_ClusterDPhi",1000,-3.16,3.16,abs(selLeptons_eleClusterDPhi[goodElectrons[j]]));
          CreateAndFillUserTH1D("eleReco_ClusterEnergy",1000,0,2000,selLeptons_eleClusterEnergy[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_HoE",100,-2,2,selLeptons_eleHoE[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_e2x5", 1000,0,2000,selLeptons_e2x5Max[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_e1x5", 1000,0,2000,selLeptons_e1x5[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_e5x5", 1000,0,2000,selLeptons_e5x5[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_ClusterEta",100,-5,5,selLeptons_eleClusterEta[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_dxy", 1000,-1,1,abs(selLeptons_dxy[goodElectrons[j]]));
          CreateAndFillUserTH1D("eleReco_Sieie",1000,-1,1,selLeptons_eleSieie[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_isEcalDriven",2,-.5,1.5,selLeptons_isEcalDriven[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_e2x5_e5x5", 1000,-10,10,selLeptons_e2x5Max[goodElectrons[j]]/selLeptons_e5x5[goodElectrons[j]]);
          CreateAndFillUserTH1D("eleReco_relIso03", 1000,0,10,selLeptons_relIso03[goodElectrons[j]]);



          if(selLeptons_isMyGoodElectron[goodElectrons[j]]==1){
            CreateAndFillUserVariableTH1D("ele_recoGood_", sizeof(varBin)/sizeof(double)-1, varBin, selLeptons_pt[goodElectrons[j]]);
            CreateAndFillUserTH1D("eleRecoGood_pt", 1000,0,2000,selLeptons_pt[goodElectrons[j]]);
            CreateAndFillUserTH1D("eleRecoGood_eta", 100,-5,5,selLeptons_eta[goodElectrons[j]]);
            CreateAndFillUserTH1D("eleRecoGood_phi", 100,-3.16,3.16,selLeptons_phi[goodElectrons[j]]);
          }//end good reco electron

          if((abs(selLeptons_eleClusterEta[j])<1.4442 && selLeptons_isEcalDriven[j]==1 && abs(selLeptons_eleClusterDEta[j])<0.004 && abs(selLeptons_eleClusterDPhi[j])< 0.06 && (selLeptons_eleHoE[j]<1.0/selLeptons_eleClusterEnergy[j]+0.05) && (selLeptons_e2x5Max[j]/selLeptons_e5x5[j]>0.94 || selLeptons_e1x5[j]/selLeptons_e5x5[j]>0.83) && abs(selLeptons_dxy[j])<0.02) || (abs(selLeptons_eleClusterEta[j])>1.566 && abs(selLeptons_eleClusterEta[j])<2.5 && selLeptons_isEcalDriven[j] && abs(selLeptons_eleClusterDEta[j])<0.006 && abs(selLeptons_eleClusterDPhi[j])< 0.06 && (selLeptons_eleHoE[j]<5.0/selLeptons_eleClusterEnergy[j]+0.05) && abs(selLeptons_dxy[j])<0.05 && selLeptons_eleSieie[j]<0.03)){
            CreateAndFillUserVariableTH1D("ele_recoGood_primo_", sizeof(varBin)/sizeof(double)-1, varBin, tmp);//just to control the offline and online cuts!
          }
        }//end if passed pt bins
      }//end for reco ele
      for(int j=0; j<fakeElectrons.size(); ++j){
        //qui i TH1D con le variabili per i fake  
        CreateAndFillUserTH1D("eleFake_pt", 1000,0,2000,selLeptons_pt[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_eta", 100,-5,5,selLeptons_eta[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_phi", 100,-3.16,3.16,selLeptons_phi[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_ClusterDEta",1000,-5,5,abs(selLeptons_eleClusterDEta[fakeElectrons[j]]));
        CreateAndFillUserTH1D("eleFake_ClusterDPhi",1000,-3.16,3.16,abs(selLeptons_eleClusterDPhi[fakeElectrons[j]]));
        CreateAndFillUserTH1D("eleFake_ClusterEnergy",1000,0,2000,selLeptons_eleClusterEnergy[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_HoE",100,2,2,selLeptons_eleHoE[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_e2x5", 1000,0,2000,selLeptons_e2x5Max[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_e1x5", 1000,0,2000,selLeptons_e1x5[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_e5x5", 1000,0,2000,selLeptons_e5x5[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_e2x5_e5x5", 1000,-10,10,selLeptons_e2x5Max[fakeElectrons[j]]/selLeptons_e5x5[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_ClusterEta",100,-5,5,selLeptons_eleClusterEta[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_dxy", 1000,-1,1,abs(selLeptons_dxy[fakeElectrons[j]]));
        CreateAndFillUserTH1D("eleFake_Sieie",1000,-1,1,selLeptons_eleSieie[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_isEcalDriven",2,-.5,1.5,selLeptons_isEcalDriven[fakeElectrons[j]]);
        CreateAndFillUserTH1D("eleFake_relIso03", 1000,0,10,selLeptons_relIso03[fakeElectrons[j]]);
      }//end for Fake Electrons

      for(int w=0; w<nselLeptons;++w){
        std::cout<<"in sel"<<std::endl;
        if(abs(selLeptons_pdgId[w])==11 ){
          if((selLeptons_pt[w]>varBin[k] && selLeptons_pt[w]<varBin[k+1])){//|| ((k==(sizeof(varBin)/sizeof(double))-1) && selLeptons_pt[w]>varBin[k])){ 
            //CreateAndFillUserTH1D("muon_RECO_1",18,0,4500, selLeptons_pt[w]);
            CreateAndFillUserVariableTH1D("ele_RECO_1",sizeof(varBin)/sizeof(double)-1, varBin, selLeptons_pt[w]);
            if(selLeptons_isMyGoodElectron[w]==1){
              std::cout<<"in good"<<std::endl;
              //CreateAndFillUserTH1D("muon_GOOD_1",18,0,4500, selLeptons_pt[w]);
              CreateAndFillUserVariableTH1D("ele_GOOD_1",sizeof(varBin)/sizeof(double)-1, varBin,selLeptons_pt[w]); 
          
            }
          }
        }
      }
      for(int w=0; w<nselLeptons;++w){
        //std::cout<<"in sel"<<std::endl;
        if(abs(selLeptons_pdgId[w])==13 ){
          if((selLeptons_pt[w]>varBin[k] && selLeptons_pt[w]<varBin[k+1])){//|| ((k==(sizeof(varBin)/sizeof(double))-1) && selLeptons_pt[w]>varBin[k])){ 
            //std::cout<<"in reco"<<std::endl;
            //CreateAndFillUserTH1D("muon_RECO_1",18,0,4500, selLeptons_pt[w]);
            CreateAndFillUserVariableTH1D("muon_RECO_1",sizeof(varBin)/sizeof(double)-1, varBin, selLeptons_pt[w]);
            if(selLeptons_isMyGoodMuon[w]==1){
              //std::cout<<"in good"<<std::endl;
              //CreateAndFillUserTH1D("muon_GOOD_1",18,0,4500, selLeptons_pt[w]);
              CreateAndFillUserVariableTH1D("muon_GOOD_1",sizeof(varBin)/sizeof(double)-1, varBin,selLeptons_pt[w]); 
          
            }
          }
        }
      }
       // std::cout<<"after"<<std::endl;
      for(int j=0; j<goodMuons.size(); ++j){
        CreateAndFillUserTH1D("muonGoodNumber", 10,-.5,9.5, goodMuons.size());
        if((selLeptons_pt[goodMuons[j]]>varBin[k] && selLeptons_pt[goodMuons[j]]<varBin[k+1])|| ((k==(sizeof(varBin)/sizeof(double))-1) && selLeptons_pt[goodMuons[j]]>varBin[k])){    
          CreateAndFillUserVariableTH1D("muon_RecoPrimo_", sizeof(varBin)/sizeof(double)-1, varBin, selLeptons_pt[goodMuons[j]]);
          std::cout<<"IN"<<std::endl;
          CreateAndFillUserTH1D("muon_RECO",14,0,3150, selLeptons_pt[goodMuons[j]]);
          
          if(selLeptons_isMyGoodMuon[goodMuons[j]]==1){
            CreateAndFillUserVariableTH1D("muon_recoGoodPrimo_", sizeof(varBin)/sizeof(double)-1, varBin, selLeptons_pt[goodMuons[j]]);
            CreateAndFillUserTH1D("muon_GOOD",14,0,3150, selLeptons_pt[goodMuons[j]]);
            CreateAndFillUserTH1D("muon_ResoPt", 1000,-10,10,GenLep_pt[goodGenMuons[j]]-selLeptons_pt[goodMuons[j]]);
            CreateAndFillUserTH1D("muon_Resolution", 1000,-10,10,(selLeptons_pt[goodMuons[j]]-GenLep_pt[goodGenMuons[j]])/GenLep_pt[goodGenMuons[j]]);
          }
        }
        if((selLeptons_pt[goodMuons[j]]>varBin[k] && selLeptons_pt[goodMuons[j]]<varBin[k+1]) || ((k==(sizeof(varBin)/sizeof(double))-1) && selLeptons_pt[goodMuons[j]]>varBin[k])){
          CreateAndFillUserVariableTH1D("muon_Reco_", sizeof(varBin)/sizeof(double)-1, varBin, tmp);
          CreateAndFillUserTH1D("muonReco_pt", 1000,0,2000,selLeptons_pt[goodMuons[j]]);
          CreateAndFillUserTH1D("muonReco_eta", 100,-5,5,selLeptons_eta[goodMuons[j]]);
          CreateAndFillUserTH1D("muonReco_phi", 100,-3.16,3.16,selLeptons_phi[goodMuons[j]]);

          //TH1D for fake/non fake muons
          CreateAndFillUserTH1D("muonReco_nChamberHits", 100,-.5,99.5, selLeptons_nChamberHits[goodMuons[j]]);
          CreateAndFillUserTH1D("muonReco_isGlobalMuon", 2,-.5,1.5, selLeptons_isGlobalMuon[goodMuons[j]]);
          CreateAndFillUserTH1D("muonReco_nStations", 100,-.5,99.5, selLeptons_nStations[goodMuons[j]]);
          CreateAndFillUserTH1D("muonReco_relPtError", 100,-2,2, selLeptons_relPtError[goodMuons[j]]);
          CreateAndFillUserTH1D("muonReco_muonDB", 100,-2,2, selLeptons_muonDB[goodMuons[j]]);
          CreateAndFillUserTH1D("muonReco_pixelHits", 100,-.5,99.5, selLeptons_pixelHits[goodMuons[j]]);
          CreateAndFillUserTH1D("muonReco_muonTrackerLayers", 100,-.5,99.5, selLeptons_muonTrackerLayers[goodMuons[j]]);
          CreateAndFillUserTH1D("muonReco_relIso04", 1000,0,10,selLeptons_relIso04[goodMuons[j]]);



          if(selLeptons_isMyGoodMuon[goodMuons[j]]==1){
            CreateAndFillUserVariableTH1D("muon_recoGood_", sizeof(varBin)/sizeof(double)-1, varBin, tmp);
            CreateAndFillUserTH1D("muonRecoGood_pt", 1000,0,2000,selLeptons_pt[goodMuons[j]]);
            CreateAndFillUserTH1D("muonRecoGood_eta", 100,-5,5,selLeptons_eta[goodMuons[j]]);
            CreateAndFillUserTH1D("muonRecoGood_phi", 100,-3.16,3.16,selLeptons_phi[goodMuons[j]]);
          }//end good reco muon
        }//end if passed pt bins
      }//end for reco muon
      for(int j=0; j<fakeMuons.size(); ++j){
      
        CreateAndFillUserTH1D("muonFake_pt", 1000,0,2000,selLeptons_pt[fakeMuons[j]]);
        CreateAndFillUserTH1D("muonFake_eta", 100,-5,5,selLeptons_eta[fakeMuons[j]]);
        CreateAndFillUserTH1D("muonFake_phi", 100,-3.16,3.16,selLeptons_phi[fakeMuons[j]]);
        CreateAndFillUserTH1D("muonFake_nChamberHits", 100,-.5,99.5, selLeptons_nChamberHits[fakeMuons[j]]);
        CreateAndFillUserTH1D("muonFake_isGlobalMuon", 2,-.5,1.5, selLeptons_isGlobalMuon[fakeMuons[j]]);
        CreateAndFillUserTH1D("muonFake_nStations", 100,-.5,99.5, selLeptons_nStations[fakeMuons[j]]);
        CreateAndFillUserTH1D("muonFake_relPtError", 100,-2,2, selLeptons_relPtError[fakeMuons[j]]);
        CreateAndFillUserTH1D("muonFake_muonDB", 100,-2,2, selLeptons_muonDB[fakeMuons[j]]);
        CreateAndFillUserTH1D("muonFake_pixelHits", 100,-.5,99.5, selLeptons_pixelHits[fakeMuons[j]]);
        CreateAndFillUserTH1D("muonFake_muonTrackerLayers", 100,-.5,99.5, selLeptons_muonTrackerLayers[fakeMuons[j]]);
        CreateAndFillUserTH1D("muonFake_relIso04", 1000,0,10,selLeptons_relIso04[fakeMuons[j]]);

      }//end for Fake Muons
      goodElectrons.clear();
      goodMuons.clear();
    }//end for pt bins
















    goodLepton.clear();
    goodJet.clear();
    //=== Filling variables ===
    for(int i=0; i<nincLeptons; ++i){
      if((incLeptons_isMyGoodMuon[i]==1 && abs(incLeptons_pdgId[i])==13) || (incLeptons_isMyGoodElectron[i]==1 && abs(incLeptons_pdgId[i])==11)){
        goodLepton.push_back(i);
      }
      //put here isEcalDriven, and dynamic cuts
      //if(abs(incLeptons_pdgId[i])==11) ++ele;
      //if(incLeptons_isMyGoodElectron[i]==1 && abs(incLeptons_eta[i])<1.2) ++ele;
      if((abs(incLeptons_pdgId[i])==11)&&(incLeptons_pt[i]>35)&&(abs(incLeptons_eta[i])<1.4442)&&(incLeptons_isEcalDriven[i]==1)&&(abs(incLeptons_eleDEta[i])<0.004)&&(abs(incLeptons_eleDPhi[i])<0.06)&&(incLeptons_eleHoE[i]<0.05+1/incLeptons_pt[i])&&(incLeptons_e1x5[i]/incLeptons_e5x5[i]>0.83)&&(incLeptons_isolEmHadDepth1[i]<(2+0.03*incLeptons_pt[i]+0.28*rho))&&(incLeptons_isolTrkPt[i]<5)&&(incLeptons_lostHits[i]<=1)&&(abs(incLeptons_dxy[i])<0.02)) ++ele_new;
      //study of cut variables
      if(abs(incLeptons_pdgId[i])==13){
        fillVariableWithValue("GlobalMuon", incLeptons_isGlobalMuon[i]);
        fillVariableWithValue("nChamberHits", incLeptons_nChamberHits[i]);
        fillVariableWithValue("nStations", incLeptons_nStations[i]);
        fillVariableWithValue("relPtErr", incLeptons_relPtError[i]);
        fillVariableWithValue("pixelHits", incLeptons_pixelHits[i]);
        fillVariableWithValue("trackerLayers", incLeptons_trackerLayers[i]);
      }
      if((incLeptons_isEcalDriven[i]==1) && (incLeptons_eleHoE[i]<0.05+1/incLeptons_pt[i]) && (incLeptons_isolEmHadDepth1[i]<(2+0.03*incLeptons_pt[i]+0.28*rho))){
        fillVariableWithValue("lepton_pdgId", abs(incLeptons_pdgId[i]));
        fillVariableWithValue("lepton_pt", incLeptons_pt[i]);
        fillVariableWithValue("lepton_eta", incLeptons_eta[i]); //remember->its cut are applied on the %!
        fillVariableWithValue("lepton_DeltaEtaSeed", incLeptons_eleDEta[i]); //remember->its cut are applied on the %!
        fillVariableWithValue("lepton_DeltaPhi", incLeptons_eleDPhi[i]); //remember->its cut are applied on the %!
        fillVariableWithValue("lepton_5x5", incLeptons_e2x5Max[i]/incLeptons_e5x5[i]);
        fillVariableWithValue("lepton_TrkPt", incLeptons_isolTrkPt[i]);
        fillVariableWithValue("lepton_lostHits", incLeptons_lostHits[i]);
        fillVariableWithValue("lepton_dxy", incLeptons_dxy[i]); //remember->its cut are applied on the %!
      }//end if for isEcalDriven, HoE, Em and Hadr isolation
      //end of study of cut variables


    
    }//end loop over 'nincLeptons'

    //=== JET CATEGORY ===
    for(int i=0; i<nJet; ++i){
      if( ((Jet_neHEF[i]<0.99) && (Jet_neEmEF[i]<0.99) && ((Jet_chMult[i]+Jet_neMult[i])>1) && ((abs(Jet_eta[i])<=2.4 && (Jet_chHEF[i])>0 && (Jet_chMult[i]>0) && Jet_chEmEF[i]<0.99) || abs(Jet_eta[i])>2.4) && abs(Jet_eta[i])<=3.0) || (Jet_neEmEF[i]<0.90 && Jet_neMult[i]>10 && abs(Jet_eta[i])>3.0 )){

	goodJet.push_back(i);     
            
      }else if(Jet_neEmEF[i]<0.90 && Jet_neMult[i]>10 && abs(Jet_eta[i])>3.0 ) {
	goodJet.push_back(i);
      }
    }//end for nJet
    
    
    
    //fillVariableWithValue("nLepton", nincLeptons);
    //fillVariableWithValue("leadingAk08_pt",FatjetAK08pruned_pt[0]);
    //evaluate cuts without applying them
    evaluateCuts();
    
    for(int size=0; size<goodLepton.size(); ++size){
      if(abs(incLeptons_pdgId[size])==11){
        CreateAndFillUserTH1D("goodEle_pt", 2000,0,1000, incLeptons_pt[goodLepton[size]]);
        CreateAndFillUserTH1D("goodEle_eta", 100,-5,5, incLeptons_eta[goodLepton[size]]);
        CreateAndFillUserTH1D("goodEle_phi", 100,-3.16,3.16, incLeptons_phi[goodLepton[size]]);
      }else if(abs(incLeptons_pdgId[size])==13){
        CreateAndFillUserTH1D("goodMuon_pt", 2000,0,1000, incLeptons_pt[goodLepton[size]]);
        CreateAndFillUserTH1D("goodMuon_eta", 100,-5,5, incLeptons_eta[goodLepton[size]]);
        CreateAndFillUserTH1D("goodMuon_phi", 100,-3.16,3.16, incLeptons_phi[goodLepton[size]]);
      }
    }//end for over goodLeptons
    for(int size=0; size<goodJet.size(); ++size){
      CreateAndFillUserTH1D("goodJet_pt", 2000,0,1000, Jet_pt[goodJet[size]]);
      CreateAndFillUserTH1D("goodJet_eta", 100,-5,5, Jet_eta[goodJet[size]]);
      CreateAndFillUserTH1D("goodJet_phi", 100,-3.16,3.16, Jet_phi[goodJet[size]]);
    }//end for over goodJet


    if(passedAllPreviousCuts("lepton_pt") && passedCut("lepton_pt")){
      fillReducedSkimTree();
    }//end pass 'if'
  }//end loop over entries
  std::cout << "analysisClass::Loop() ends" <<std::endl;   
  std::cout<<"Electrons: "<<ele<<std::endl;
  std::cout<<"Muons: "<<mu<<std::endl;
  std::cout<<"Ele new: "<<ele_new<<std::endl; 
  writeTree();
  std::cout<<""<<std::endl;
  std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << time.RealTime() << " s" <<std::endl;
  std::cout<<""<<std::endl;




  std::cout<<n_eleGen<<" "<< n_eleReco<<" "<< n_eleRecoGood<<std::endl;
} 
