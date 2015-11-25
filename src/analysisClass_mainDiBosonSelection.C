#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

  /*std::string jetAlgo = getPreCutString1("jetAlgo");
  double rParam = getPreCutValue1("DeltaR");

  if( jetAlgo == "AntiKt" )
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::antikt_algorithm, rParam) );
  else if( jetAlgo == "Kt" )
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::kt_algorithm, rParam) );
  else 
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::cambridge_algorithm, rParam) );
*/
  // For JECs  ---  removed for the DiBoson analysis
/*  if( int(getPreCutValue1("useJECs"))==1 )
  {
    std::cout << "Reapplying JECs on the fly" << std::endl;
    //std::string L1Path = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_6_patch6/src/CMSDIJET/DijetRootTreeMaker/data/Summer15_V5/Summer15_V5_MC_L1FastJet_AK4PFchs.txt";
    //std::string L2Path = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_6_patch6/src/CMSDIJET/DijetRootTreeMaker/data/Summer15_V5/Summer15_V5_MC_L2Relative_AK4PFchs.txt"; 
    //std::string L3Path = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_6_patch6/src/CMSDIJET/DijetRootTreeMaker/data/Summer15_V5/Summer15_V5_MC_L3Absolute_AK4PFchs.txt";
    //std::string L1Path = "data/Summer15_50nsV2/Summer15_50nsV2_MC_L1FastJet_AK4PFchs.txt";
    //std::string L2Path = "data/Summer15_50nsV2/Summer15_50nsV2_MC_L2Relative_AK4PFchs.txt"; 
    //std::string L3Path = "data/Summer15_50nsV2/Summer15_50nsV2_MC_L3Absolute_AK4PFchs.txt";
    //std::string L1DATAPath = "data/Summer15_50nsV2/Summer15_50nsV2_MC_L1FastJet_AK4PFchs.txt";
    //std::string L2DATAPath = "data/Summer15_50nsV2/Summer15_50nsV2_MC_L2Relative_AK4PFchs.txt"; 
    //std::string L3DATAPath = "data/Summer15_50nsV2/Summer15_50nsV2_MC_L3Absolute_AK4PFchs.txt";
    //
    //std::string L1Path = "data/Summer15_50nsV4/Summer15_50nsV4_MC_L1FastJet_AK4PFchs.txt";
    //std::string L2Path = "data/Summer15_50nsV4/Summer15_50nsV4_MC_L2Relative_AK4PFchs.txt"; 
    //std::string L3Path = "data/Summer15_50nsV4/Summer15_50nsV4_MC_L3Absolute_AK4PFchs.txt";
    //std::string L1DATAPath = "data/Summer15_50nsV4/Summer15_50nsV4_DATA_L1FastJet_AK4PFchs.txt";
    //std::string L2DATAPath = "data/Summer15_50nsV4/Summer15_50nsV4_DATA_L2Relative_AK4PFchs.txt"; 
    //std::string L3DATAPath = "data/Summer15_50nsV4/Summer15_50nsV4_DATA_L3Absolute_AK4PFchs.txt";
    //std::string L2L3ResidualPath = "data/Summer15_50nsV4/Summer15_50nsV4_DATA_L2L3Residual_AK4PFchs.txt" ;
    //
    //std::string L1Path = "data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L1FastJet_AK4PFchs.txt";
    //std::string L2Path = "data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L2Relative_AK4PFchs.txt"; 
    //std::string L3Path = "data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L3Absolute_AK4PFchs.txt";
    //std::string L1DATAPath = "data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L1FastJet_AK4PFchs.txt";
    //std::string L2DATAPath = "data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2Relative_AK4PFchs.txt"; 
    //std::string L3DATAPath = "data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L3Absolute_AK4PFchs.txt";
    //std::string L2L3ResidualPath = "data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2L3Residual_AK4PFchs.txt" ;
    //
    //std::string L1Path = "data/Summer15_25nsV5_MC/Summer15_25nsV5_MC_L1FastJet_AK4PFchs.txt";
    //std::string L2Path = "data/Summer15_25nsV5_MC/Summer15_25nsV5_MC_L2Relative_AK4PFchs.txt"; 
    //std::string L3Path = "data/Summer15_25nsV5_MC/Summer15_25nsV5_MC_L3Absolute_AK4PFchs.txt";
    //std::string L1DATAPath = "data/Summer15_25nsV5_DATA/Summer15_25nsV5_DATA_L1FastJet_AK4PFchs.txt";
    //std::string L2DATAPath = "data/Summer15_25nsV5_DATA/Summer15_25nsV5_DATA_L2Relative_AK4PFchs.txt"; 
    //std::string L3DATAPath = "data/Summer15_25nsV5_DATA/Summer15_25nsV5_DATA_L3Absolute_AK4PFchs.txt";
    //std::string L2L3ResidualPath = "data/Summer15_25nsV5_DATA/Summer15_25nsV5_DATA_L2L3Residual_AK4PFchs.txt" ;
    //
    std::string L1Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt";
    std::string L2Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt"; 
    std::string L3Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt";
    std::string L1DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt";
    std::string L2DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt"; 
    std::string L3DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt";
    std::string L2L3ResidualPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt" ;

    
    L1Par = new JetCorrectorParameters(L1Path);
    L2Par = new JetCorrectorParameters(L2Path);
    L3Par = new JetCorrectorParameters(L3Path);
    L1DATAPar = new JetCorrectorParameters(L1DATAPath);
    L2DATAPar = new JetCorrectorParameters(L2DATAPath);
    L3DATAPar = new JetCorrectorParameters(L3DATAPath);
    L2L3Residual = new JetCorrectorParameters(L2L3ResidualPath);

    std::vector<JetCorrectorParameters> vPar;
    std::vector<JetCorrectorParameters> vPar_data;
    vPar.push_back(*L1Par);
    vPar.push_back(*L2Par);
    vPar.push_back(*L3Par);
   
    //residuals are applied only to data
    vPar_data.push_back(*L1DATAPar);
    vPar_data.push_back(*L2DATAPar);
    vPar_data.push_back(*L3DATAPar);
    vPar_data.push_back(*L2L3Residual);

    JetCorrector = new FactorizedJetCorrector(vPar);
    JetCorrector_data = new FactorizedJetCorrector(vPar_data);

    //uncertainty
    //unc = new JetCorrectionUncertainty("data/Summer15_50nsV5_DATA/Summer15_50nsV5_DATA_Uncertainty_AK4PFchs.txt");
    //unc = new JetCorrectionUncertainty("data/Summer15_25nsV5_DATA/Summer15_25nsV5_DATA_Uncertainty_AK4PFchs.txt");
    unc = new JetCorrectionUncertainty("data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_Uncertainty_AK4PFchs.txt");
  }*/
  
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

   // variable binning for mjj trigger efficiency plots
   /*const int nMassBins = 103;

   double massBoundaries[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325,
     354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,
     1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509,
     4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 
     10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000};


   char* HLTname[50] = {"noTrig","PFHT475","PFHT800","PFHT650MJJ900","PFHT800_OR_PFHT650MJJ900","PFHT800_noPFHT475", 
                        "Mu45Eta2p1", "PFHT800AndMu45Eta2p1"};
   TH1F* h_mjj_HLTpass[8];
   char name_histoHLT[50];
   for (int i=0; i<8; i++){  
     sprintf(name_histoHLT,"h_mjj_HLTpass_%s",HLTname[i]);
     h_mjj_HLTpass[i]= new TH1F(name_histoHLT,"",103,massBoundaries);
   }
  */

   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<2000;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done for every event
     
     resetCuts();

     double Jet=Jet_pt[0];
     fillVariableWithValue("Jet",triggerResult->at(1));// HLT_QuadPFJet_VBF_v*

     CreateAndFillUserTH1D("First_Jet",100,-10.,10.,Jet_pt[0]);

     // Evaluate cuts (but do not apply them)
     evaluateCuts();
     
     // optional call to fill a skim with the full content of the input roottuple
     //if( passedCut("nJetFinal") ) fillSkimTree();
     /*if( passedCut("PassJSON")
	 && passedCut("nVtx") 
	 && passedCut("IdTight_j1")
	 && passedCut("IdTight_j2")
	 && passedCut("nJet")
	 && passedCut("pTWJ_j1")
	 && passedCut("etaWJ_j1")
	 && passedCut("pTWJ_j2")
	 && passedCut("etaWJ_j2")
	 && getVariableValue("deltaETAjj") <  getPreCutValue1("DetaJJforTrig") )*/

     // optional call to fill a skim with a subset of the variables defined in the cutFile (use flag SAVE)
     //if( passedAllPreviousCuts("mjj") && passedCut("mjj") )  //ADD IT AGAIN!!!!!!!!!!!!! 
     //  {
	// fillReducedSkimTree();

	 // ===== Take a look at this =====
	 // //Example on how to investigate quickly the data
	 // if(getVariableValue("mjj")>4000)
	 //   {
	 //     //fast creation and filling of histograms
	 //     CreateAndFillUserTH1D("h_dphijj_mjjgt4000", 100, 0, 3.15, getVariableValue("deltaPHIjj"));
	 //     CreateAndFillUserTH1D("h_htak4_mjjgt4000", 1000, 0, 10000, getVariableValue("HTAK4"));
	 //     CreateAndFillUserTH1D("h_nvtx_mjjgt4000", 31, -0.5, 30.5, getVariableValue("nVtx"));
	 //   }

      // }

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

   //Writing histos (different from TH1DUserHisto...
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
   }
