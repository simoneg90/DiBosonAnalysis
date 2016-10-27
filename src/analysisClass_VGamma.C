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
#include "TGraphAsymmErrors.h"
#include "TF1.h"

#define el_Mass 0.0005 //GeV
#define mu_Mass 0.1 //GeV
#define W_mass 80 //GeV

class Chebyshev {
   public:
      Chebyshev(int n, double xmin, double xmax) :
         fA(xmin), fB(xmax),
         fT(std::vector<double>(n) )  {}
    
      double operator() (const double * xx, const double *p) {
          double x = (xx[0] - fA -fB)/(fB-fA);
          int order = fT.size();
          if (order == 1) return p[0];
          if (order == 2) return p[0] + x*p[1];
          // build the polynomials
          fT[0] = 1;
          fT[1] = x;
          for (int i = 1; i< order; ++i) {
              fT[i+1] =  2 *x * fT[i] - fT[i-1];
          }
          double sum = p[0]*fT[0];
          for (int i = 1; i<= order; ++i) {
              sum += p[i] * fT[i];
          }
          return sum;
     }
     
   private:
      double fA;
      double fB;
      std::vector<double> fT; // polynomial
      std::vector<double> fC; // coefficients
};


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
   std::cout<<"InputList: "<<inputList_->c_str()<<std::endl;
   int isGammaJets=0;
   if(system(Form("grep 'GJets' %s", inputList_->c_str()))==0){
     frame("Going to analyse a GJets sample");
     isGammaJets=1;
   }
   std::cout << "analysisClass::Loop() begins" <<std::endl;   
   TStopwatch time;
   frame("Starting the analysis");
   time.Start(true);
   setTDRStyle();
   if (fChain == 0) return;
   

   /////////initialize and define variables
   TLorentzVector photon, genPhoton, electron, muon, ak04, ak08, ak04_photon, ak08_photon, ak08_puppi, ak08_puppiSoftdrop;
   std::vector <int> goodPhoton, goodElectron, goodMuon, goodAk04, goodAk08;
   bool mu_jet_DR, el_jet_DR;
   double drMin;
   int matchSoftDrop;
   std::vector <int> jetID, photonID;
   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   
   TFile* gJet_correction = TFile::Open("/cmshome/gellisim/CMSSW_VGamma/src/DiBosonAnalysis/uncertainties_EWK_24bins.root");
   if (!gJet_correction ){
     frame("GammaJets corrections file not found! ERROR!");
     exit(-1);
   }
   TFile* puwFile = TFile::Open("/cmshome/gellisim/CMSSW_VGamma/src/DiBosonAnalysis/pileUp_17_2fb.root");//puw_2016_13fb_200.root");//pileup_profile_runs_271036_279931.root");//puw_2016_13fb_200.root");
   if (!puwFile){
     frame("Pile Up file not found! ERROR!");
     exit(-1);
   }
   TH1F* gCorrNominal = (TH1F*)gJet_correction->Get("GJets_1j_NLO/nominal_G");
   TH1F* gCorrInv = (TH1F*)gJet_correction->Get("GJets_LO/inv_pt_G");
   TH1D* puw = (TH1D*)puwFile->Get("pileup");
   Chebyshev * cheb = new Chebyshev(4,200,1230);
   TF1 * f1_ = new TF1("f1_",cheb,200,1230,5,"Chebyshev");
   TGraphAsymmErrors *grCorr = new TGraphAsymmErrors(0);
   grCorr->Divide(gCorrNominal,gCorrInv, "pois");
   grCorr->Fit(f1_, "");//, "R");
   double pu_weight=1.;
   double gJets_correction=1.;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<1000;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%10000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     if(jentry%1000 == 0 && jentry%10000 !=0) std::cout<<"."<<std::flush; 
    
     if(!isData && nTrueInt>0) pu_weight=puw->GetBinContent(puw->FindBin(nTrueInt));

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done for every event
     goodPhoton.clear();
     goodElectron.clear();
     goodMuon.clear();
     goodAk04.clear();
     goodAk08.clear();
     jetID.clear();
     photonID.clear();

     resetCuts();

     ////Start coding

     //std::cout<<"Initial good photons "<<nGammaGood<<std::endl;
     for(int i=0; i<nGammaGood; ++i){
       //std::cout<<"pt: "<<GammaGood_pt[i]<<" eta: "<<abs(GammaGood_eta[i])<<std::endl;
      //  if(GammaGood_pt[i]>80 && abs(GammaGood_eta[i])<2.4){
          //std::cout<<"Good photon"<<std::endl;
          goodPhoton.push_back(i);
          if((abs(GammaGood_eta[i])<1.4442&&GammaGood_hOverE[i]<.05&&GammaGood_sigmaIetaIeta[i]<0.0102&&GammaGood_chHadIso[i]<3.32&&GammaGood_neuHadIso[i]<(1.92+0.014*GammaGood_pt[i]+0.000019*GammaGood_pt[i]*GammaGood_pt[i])&&GammaGood_phIso[i]<(.81+0.0053*GammaGood_pt[i]))|| ((abs(GammaGood_eta[i])<2.5&&abs(GammaGood_eta[i])>1.566)&&GammaGood_hOverE[i]<.05&&GammaGood_sigmaIetaIeta[i]<0.0274&&GammaGood_chHadIso[i]<1.97&&GammaGood_neuHadIso[i]<(11.86+0.0139*GammaGood_pt[i]+0.000025*GammaGood_pt[i]*GammaGood_pt[i])&&GammaGood_phIso[i]<(.83+0.0034*GammaGood_pt[i]))){
          //if((abs(GammaGood_eta[i])<1.4442&&GammaGood_hOverE[i]<.05&&GammaGood_sigmaIetaIeta[i]<0.0102&&GammaGood_chHadIso[i]<3.32&&GammaGood_neuHadIso[i]<(11.86+0.0139*GammaGood_pt[i]+0.000025*GammaGood_pt[i]*GammaGood_pt[i])&&GammaGood_phIso[i]<(.83+0.0034*GammaGood_pt[i]))|| ((abs(GammaGood_eta[i])<2.5||abs(GammaGood_eta[i])>1.566)&&GammaGood_hOverE[i]<.05&&GammaGood_sigmaIetaIeta[i]<0.0274&&GammaGood_chHadIso[i]<1.97&&GammaGood_neuHadIso[i]<(11.86+0.0139*GammaGood_pt[i]+0.000025*GammaGood_pt[i]*GammaGood_pt[i])&&GammaGood_phIso[i]<(.83+0.0034*GammaGood_pt[i]))){
            photonID.push_back(1);
          }else{
            photonID.push_back(0);
          }

        //}//end if for pt and eta requirements

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


     for(int i=0; i<nFatJet; ++i){
    
        if(FatJet_pt[i]>80 && abs(FatJet_eta[i])<2.5){//&& FatJet_id[i]==1){ //FatJet_id[i]==1 corresponds to the Loose ID (https://github.com/CERN-PH-CMG/cmg-cmssw/blob/fc6c8f8fa537a417ae622797bbf5d61a7a90a175/PhysicsTools/Heppy/python/physicsobjects/Jet.py)
            ak08.SetPtEtaPhiM(FatJet_pt[i], FatJet_eta[i], FatJet_phi[i], FatJet_mass[i]);
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
              //std::cout<<"Good ak08"<<std::endl; 
            drMin=100.;
            for(int j=0; j<goodPhoton.size();++j){
              //std::cout<<"Photon loop in jet"<<std::endl;
              photon.SetPtEtaPhiM(GammaGood_pt[goodPhoton[j]],GammaGood_eta[goodPhoton[j]],GammaGood_phi[goodPhoton[j]],0);
              //std::cout<<"-------------- "<<drMin<<std::endl;
              if(photon.DeltaR(ak08)<drMin){
                drMin=photon.DeltaR(ak08);
                //std::cout<<"-------------- "<<drMin<<std::endl;
              }
              //std::cout<<"-------------- "<<drMin<<std::endl;
            }
            fillVariableWithValue("ak08_Photon_DR", drMin);
            if(drMin>=.9){
              goodAk08.push_back(i);
              //std::cout<<"-- 1++"<<i<<std::endl;
              if(((FatJet_neHEF[i]<.9 && FatJet_neEmEF[i]<.9 && (FatJet_neMult[i]+FatJet_chMult[i])>1) && ((abs(FatJet_eta[i])<=2.4 && FatJet_chHEF[i]>0 && FatJet_chMult[i]>0 && FatJet_chEmEF[i]<.9) || abs(FatJet_eta[i])>2.4) && abs(FatJet_eta[i])<=3.) || ((FatJet_neEmEF[i]<.9 && FatJet_neMult[i]>10 && abs(FatJet_eta[i])>3.) )){
                jetID.push_back(1);
                //std::cout<<"++ 1++"<<i<<std::endl;
              }else{
                //std::cout<<"++ 0++"<<i<<std::endl;
                jetID.push_back(0);
              }
            }
/////            }
        }//end if for good jet ak08

     }//end loop over # of ak8

     for(int i=0; i<nJet; ++i){
    
        if(Jet_pt[i]>180 && abs(Jet_eta[i])<2. && Jet_id[i]==1){
            ak04.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
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
              //std::cout<<"Good ak04"<<std::endl;
              goodAk04.push_back(i);
/////            }
        }//end if for good jet ak04

     }//end loop over # of ak4

     ////Filling outputTree
     //std::cout<<"n of good photons: "<<goodPhoton.size()<<std::endl;
     //std::cout<<"Diff nphotons: "<<nGammaGood-goodPhoton.size()<<std::endl;
     
     if(goodPhoton.size()>0){
       int closestGen=-1;
       if(isGammaJets==1){
         drMin=100.;
         photon.SetPtEtaPhiM(GammaGood_pt[goodPhoton[0]],GammaGood_eta[goodPhoton[0]],GammaGood_phi[goodPhoton[0]],0);
         for(int j=0; j<nGenPart; ++j){
           genPhoton.SetPtEtaPhiM(GenPart_pt[j], GenPart_eta[j], GenPart_phi[j], 0);
           if(photon.DeltaR(genPhoton)<drMin){
             drMin=photon.DeltaR(genPhoton);
             closestGen=j;
           }
         }
         if(closestGen>=0) {
           //std::cout<<"deriving correction"<<std::endl;
           //std::cout<<"Pt: "<<GenPart_pt[closestGen]<<std::endl;

           gJets_correction=f1_->Eval(GenPart_pt[closestGen]);//(gCorrNominal->GetBinContent(GenPart_pt[closestGen])/gCorrInv->GetBinContent(GenPart_pt[closestGen]));
         //  std::cout<<"---> "<<f1_->Eval(GenPart_pt[closestGen])<<" "<<std::endl;
         }
       } 
       //std::cout<<"---> "<< gJets_correction<<" "<<std::endl;      

       fillVariableWithValue("pw_correction",pu_weight);
       fillVariableWithValue("photonCorrection", gJets_correction);
       fillVariableWithValue("photon_pt", GammaGood_pt[goodPhoton[0]]);
       fillVariableWithValue("photon_eta", GammaGood_eta[goodPhoton[0]]);
       fillVariableWithValue("photon_phi", GammaGood_phi[goodPhoton[0]]);
       fillVariableWithValue("photonLooseID", photonID[0]);
       //std::cout<<"++++++++++++ "<<GammaGood_pt[goodPhoton[0]]<<std::endl;
     } 

/////     if(goodElectron.size()>0){
/////       fillVariableWithValue("electron_pt", el_pt[goodElectron[0]));
/////       fillVariableWithValue("electron_eta", el_eta[goodElectron[0]));
/////       fillVariableWithValue("electron_phi", el_phi[goodElectron[0]));
/////     }
/////
/////     if(goodMuon.size()>0){
/////       fillVariableWithValue("muon_pt", mu_pt->at(goodMuon[0]));
/////       fillVariableWithValue("muon_eta", mu_eta->at(goodMuon[0]));
/////       fillVariableWithValue("muon_phi", mu_phi->at(goodMuon[0]));
/////     }

     if(goodAk08.size()>0){
       fillVariableWithValue("ak08_pt", FatJet_pt[goodAk08[0]] );
       fillVariableWithValue("ak08_eta", FatJet_eta[goodAk08[0]] );
       fillVariableWithValue("ak08_phi", FatJet_phi[goodAk08[0]] );
       fillVariableWithValue("ak08_mass", FatJet_mass[goodAk08[0]] );
       fillVariableWithValue("ak08_prunedMass", FatJet_prunedMass[goodAk08[0]] );
       fillVariableWithValue("ak08_CSV", FatJet_btagCSV[goodAk08[0]] );
       fillVariableWithValue("ak08_tau21", (FatJet_tau2[goodAk08[0]])/(FatJet_tau1[goodAk08[0]]) );
       fillVariableWithValue("ak08_puppiTau21", (FatJet_puppiTau2[goodAk08[0]])/(FatJet_puppiTau1[goodAk08[0]]) );
       fillVariableWithValue("ak08_puppiPt",FatJet_puppiPt[goodAk08[0]]);
       fillVariableWithValue("ak08_puppiEta",FatJet_puppiEta[goodAk08[0]]);
       fillVariableWithValue("ak08_puppiPhi",FatJet_puppiPhi[goodAk08[0]]);
       fillVariableWithValue("ak08_neEmHF", FatJet_neEmEF[goodAk08[0]]);
       fillVariableWithValue("jetTightLepVetoID", jetID[0]);
     }

     if(goodAk04.size()>0){
       fillVariableWithValue("ak04_pt", Jet_pt[goodAk04[0]] );
       fillVariableWithValue("ak04_eta", Jet_eta[goodAk04[0]] );
       fillVariableWithValue("ak04_phi", Jet_phi[goodAk04[0]] );
       fillVariableWithValue("ak04_mass", Jet_mass[goodAk04[0]] );
     }

     if(goodPhoton.size()>0 && goodAk08.size()>0){
       ak08.SetPtEtaPhiM(FatJet_pt[goodAk08[0]], FatJet_eta[goodAk08[0]], FatJet_phi[goodAk08[0]], FatJet_mass[goodAk08[0]]);
       photon.SetPtEtaPhiM(GammaGood_pt[goodPhoton[0]], GammaGood_eta[goodPhoton[0]], GammaGood_phi[goodPhoton[0]], 0);
       ak08_photon = ak08 + photon;
       fillVariableWithValue("ak08_photon_pt", ak08_photon.Pt() );
       fillVariableWithValue("ak08_photon_eta", ak08_photon.Eta() );
       fillVariableWithValue("ak08_photon_phi", ak08_photon.Phi() );
       fillVariableWithValue("ak08_photon_mass", ak08_photon.M() );;
       fillVariableWithValue("photonPt_over_ak08PhotMass", GammaGood_pt[goodPhoton[0]]/ak08_photon.M());
       fillVariableWithValue("cosThetaStar_ak08Photon",TMath::TanH((GammaGood_eta[goodPhoton[0]]-FatJet_eta[goodAk08[0]])/2));
       fillVariableWithValue("DeltaEta_photonAk08", GammaGood_eta[goodPhoton[0]]-FatJet_eta[goodAk08[0]]);
       drMin=100.;
       for(int i=0; i<ncustomPuppiSoftDropAK8; ++i){
         ak08_puppiSoftdrop.SetPtEtaPhiM(customPuppiSoftDropAK8_pt[i],customPuppiSoftDropAK8_eta[i],customPuppiSoftDropAK8_phi[i],customPuppiSoftDropAK8_mass[i]);
         if(ak08.DeltaR(ak08_puppiSoftdrop)<drMin) {
           drMin=ak08.DeltaR(ak08_puppiSoftdrop);
           matchSoftDrop=i;
         }
       }//end loop over puppi/softDrop Quantities
       if(drMin<.3){

         ak08_puppiSoftdrop.SetPtEtaPhiM(customPuppiSoftDropAK8_pt[matchSoftDrop],customPuppiSoftDropAK8_eta[matchSoftDrop],customPuppiSoftDropAK8_phi[matchSoftDrop],customPuppiSoftDropAK8_massCorrected[matchSoftDrop]);
         ak08_photon = ak08_puppiSoftdrop + photon;
         fillVariableWithValue("ak08_photon_massSoftDropCorr", ak08_photon.M());
         ak08_puppiSoftdrop.SetPtEtaPhiM(customPuppiSoftDropAK8_pt[matchSoftDrop],customPuppiSoftDropAK8_pt[matchSoftDrop],customPuppiSoftDropAK8_pt[matchSoftDrop],customPuppiSoftDropAK8_mass[matchSoftDrop]);
         ak08_photon = ak08_puppiSoftdrop + photon;
         fillVariableWithValue("ak08_photon_massSoftDrop", ak08_photon.M());

       }

       
     }
     if(goodPhoton.size()>0 && goodAk04.size()>0){
       ak04.SetPtEtaPhiM(Jet_pt[goodAk04[0]], Jet_eta[goodAk04[0]], Jet_phi[goodAk04[0]], Jet_mass[goodAk04[0]]);
       photon.SetPtEtaPhiM(GammaGood_pt[goodPhoton[0]], GammaGood_eta[goodPhoton[0]], GammaGood_phi[goodPhoton[0]], 0);
       ak04_photon = ak04 + photon;
       fillVariableWithValue("ak04_photon_pt", ak04_photon.Pt() );
       fillVariableWithValue("ak04_photon_eta", ak04_photon.Eta() );
       fillVariableWithValue("ak04_photon_phi", ak04_photon.Phi() );
       fillVariableWithValue("ak04_photon_mass", ak04_photon.M() );
     }
     fillVariableWithValue("nPhotons", goodPhoton.size());
     fillVariableWithValue("nAk04", goodAk04.size());
     fillVariableWithValue("nAk08", goodAk08.size());
     fillVariableWithValue("nTrueInteractions", nTrueInt);
     fillVariableWithValue("nVert", nVert);
     fillVariableWithValue("HLT_BIT_HLT_PFJet40_v",HLT_BIT_HLT_PFJet40_v);
     fillVariableWithValue("HLT_BIT_HLT_PFHT800_v", HLT_BIT_HLT_PFHT800_v);
     fillVariableWithValue("HLT_BIT_HLT_Photon165_HE10_v",HLT_BIT_HLT_Photon165_HE10_v);
     fillVariableWithValue("HLT_BIT_HLT_Photon175_v",HLT_BIT_HLT_Photon175_v);
     fillVariableWithValue("HLT_BIT_HLT_PFJet400_v",HLT_BIT_HLT_PFJet400_v);
     fillVariableWithValue("HLT_BIT_HLT_PFJet140_v",HLT_BIT_HLT_PFJet140_v);
     fillVariableWithValue("HLT_BIT_HLT_PFJet500_v",HLT_BIT_HLT_PFJet500_v);
     fillVariableWithValue("HLT_BIT_HLT_PFJet260_v",HLT_BIT_HLT_PFJet260_v);
     fillVariableWithValue("HLT_BIT_HLT_PFJet200_v",HLT_BIT_HLT_PFJet200_v);
     fillVariableWithValue("metRaw", met_rawPt );
     fillVariableWithValue("met", met_pt );

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
