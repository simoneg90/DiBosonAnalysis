#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "RooArgusBG.h"
#include "RooDCBShape.cxx"
//#include "RooDCBShape_cxx.so"
#include <TF1.h>
#include <TMath.h>
using namespace RooFit ;

#define nbins 18 //24
#define xmin 40  //30
#define xMax 130 //150
#define Lumi 3.99*1000 //[pb-1]

void fit_Signal_new()
{
  bool doWeight = 0;  //if you want to have reweighted dataset->1 otherwise use 0
  gSystem->AddIncludePath("-I$ROOFITSYS/include");
  gROOT->GetInterpreter()->AddIncludePath("$ROOFITSYS/include");
  gSystem->SetIncludePath("-I$ROOFITSYS/include");
  gSystem->Load("$ROOFITSYS/lib/libRooFit.so") ;
  gSystem->Load("$ROOFITSYS/lib/libRooFitCore.so") ;
  gSystem->Load("$ROOFITSYS/lib/libRooStats.so") ;
  using namespace RooFit ;
  gSystem->Load("libRooFit");
  //gSystem->Load("RooDCBShape_cxx.so");
  //Take variable from tree
  //model_ttbar.fitTo(*dh_ttbarMatch) ;
  //ak08_subjetDR<.1&&subjet1_btagLoose==0&&subjet2_btagLoose==0&&

  std::string match= "ak08Ungroomed_1_tau21<.6&&genW_genBquark2_DR>.8&&genW_genBquark1_DR>.8&&ak08Ungroomed_WGen_DR<.1";
  std::string unmatch= "ak08Ungroomed_1_tau21<.6&&((genW_genBquark1_DR<.8&&ak08Ungroomed_WGen_DR<.1)||(genW_genBquark2_DR<.8&&ak08Ungroomed_WGen_DR<.1)||(genW_genBquark2_DR<.8&&genW_genBquark1_DR<.8&&ak08Ungroomed_WGen_DR<.1)||ak08Ungroomed_WGen_DR>.1)";
  std::string pass= "ak08Ungroomed_1_tau21<.6";




  std::string cut_matched_pass = "((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&(btag_medium>0)&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&lepton_goodNumber==1&&ak08Ungroomed_1_tau21<6&&genW_genBquark2_DR>.8&&genW_genBquark1_DR>.8&&ak08Ungroomed_WGen_DR<.1";//((genW_genBquark2_DR<1.&&genW_genBquark1_DR<1.&&ak08Ungroomed_WGen_DR<.1)||ak08Ungroomed_WGen_DR>1.)";
  //std::string cut_matched_pass = "((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&btag_medium>0&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&lepton_goodNumber==1&&ak08Ungroomed_WGen_DR<.1&&genW_genBquark1_DR>1.&&genW_genBquark2_DR>1.";//&&lepton_WGen_DR>1.";//genW_genBquark1_DR>1.&&genW_genBquark2_DR>1.&&ak08Ungroomed_1_tau21<2";//&&lepton_WGen_DR>1.";//&&genW_genBquark1_DR>1.&&genW_genBquark2_DR>1.";//&&genW_genBquarkMin_DR>3.&&genW_genBquarkMax_DR>2";//subjet1_qGen1_DR<0.1&&subjet2_qGen2_DR<0.1&&ak08_subjetDR<0.1&&subjet1_btagLoose==0&&subjet2_btagLoose==0";
  std::string cut_unmatched_pass = "((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&(btag_medium>0)&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&lepton_goodNumber==1&&ak08Ungroomed_1_tau21>.6&&((genW_genBquark1_DR<.8&&ak08Ungroomed_WGen_DR<.1)||(genW_genBquark2_DR<.8&&ak08Ungroomed_WGen_DR<.1)||(genW_genBquark2_DR<.8&&genW_genBquark1_DR<.8&&ak08Ungroomed_WGen_DR<.1)||ak08Ungroomed_WGen_DR>.1)";//&&subjet1_btagTight==0&&subjet2_btagTight==0&&ak08_subjetDR<.3";
//"((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&(btag_tight>0)&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Ungroomed_1_tau21<0.6&&ak08Ungroomed_WGen_DR>.1&&subjet1_qGen1_DR<0.1&&subjet2_qGen2_DR<0.1&&ak08_subjetDR<0.1&&lepton_goodNumber==1";//&&lepton_WGen_DR>1";
  std::string cut_pass = "((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&(btag_loose>0&&btag_medium>0)&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Ungroomed_1_tau21>.6&lepton_goodNumber==1";//&&subjet1_btagLoose==0&&subjet2_btagLoose==0&&ak08_subjetDR<0.1";//"((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1&&lepton_goodNumber==1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&btag_medium>0&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Ungroomed_1_tau21<0.6";
  std::string cut_fail = "((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1&&lepton_goodNumber==1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&btag_medium>0&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Ungroomed_1_tau21>0.6";
  std::string cut = "abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&btag_medium>0&&lepton_pt>40&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Pruned_1_mass>40&&ak08Pruned_1_mass<150&&ak08Ungroomed_1_tau21<0.6&&ak08Ungroomed_WGen_DR<0.1";
  std::string variable="ak08Pruned_1_mass";

  std::cout<<"String 1: "<<cut_matched_pass.c_str()<<std::endl;
  std::cout<<"String 2: "<<cut.c_str()<<std::endl;
 
  std::cout<<"+++++++++++++++++++"<<std::endl;
  std::cout<<"      TopTop       "<<std::endl;
  std::cout<<"+++++++++++++++++++"<<std::endl;


  TFile *file_ttbar;
  file_ttbar=TFile::Open("../ttbar_output_new/TT_Inclusive_20160703_190901/rootfile_tt_inclusive__20160703_190901_0_reduced_skim.root");//"../ttbar_output_new/TT_Inclusive_20160623_191418/rootfile_inputListTest__20160623_191418_0_reduced_skim.root");//"/afs/cern.ch/work/s/sgelli/private/CMSSW_Analysis/src/DiBosonAnalysis/rootFolder/TTBar_Inclusive.root");//"/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/TTBar_Inclusive.root");///cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/output/bulkGraviton1000_reduced_skim.root");//");
  std::cout<<"opened file"<<std::endl;
  double nEvt_ttbar= ((TH1D *)(file_ttbar->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  std::cout<<"+++++++++++++++++"<<nEvt_ttbar<<std::endl;
  double xSec_ttbar = 831.76;

  RooRealVar ak08Pruned_1_mass("ak08Pruned_1_mass","ak08Pruned_1_mass", 40,130);
  RooRealVar ak08Ungroomed_1_tau21("ak08Ungroomed_1_tau21","ak08Ungroomed_1_tau21", 0,1);
  RooRealVar genW_genBquark1_DR("genW_genBquark1_DR","genW_genBquark1_DR",-10,10);
  RooRealVar genW_genBquark2_DR("genW_genBquark2_DR","genW_genBquark2_DR",-10,10);
  RooRealVar ak08Ungroomed_WGen_DR("ak08Ungroomed_WGen_DR","ak08Ungroomed_WGen_DR",-10,10);

  TFile *file_ttbar2;
  file_ttbar2=TFile::Open("../ttbar_output_new/TT_Inclusive_20160703_190901/rootfile_tt_inclusive__20160703_190901_0_reduced_skim_redu.root");//"/afs/cern.ch/work/s/sgelli/private/CMSSW_Analysis/src/DiBosonAnalysis/ttbar_output_new/TT_Inclusive_redu.root");//ttbar_redu.root");
  std::cout<<"After tree assigned"<<std::endl;
  TTree *tree_ttbar=(TTree *)file_ttbar2->Get("mio");//tree_ttbar->CopyTree(cut_matched_pass.c_str());
  
  RooDataSet *ds_ttbar = new RooDataSet("ds_ttbar","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_ttbar));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  std::cout<<"RooDataset defined"<<std::endl;
  ds_ttbar->Print();

  RooDataSet* dsC_match = (RooDataSet*) ds_ttbar->reduce(match.c_str()) ;
  RooDataSet* dsC_unmatch = (RooDataSet*) ds_ttbar->reduce(unmatch.c_str()) ;
  dsC_match->SetName("dsC_match");
  dsC_match->Print();
  dsC_unmatch->SetName("dsC_unmatch");
  dsC_unmatch->Print();
  double weight_ttbar=xSec_ttbar*Lumi/nEvt_ttbar;
  if(!doWeight) weight_ttbar=1;
  //RooConstVar *mioPeso[10];
  //mioPeso[0]= new RooConstVar("mioPeso0","mioPeso0",weight_ttbar);
  //RooFormulaVar *wFunc_mio[10];
  //std::cout<<"++++++++++++*"<< mioPeso[0]->getValV()<<std::endl;
  std::cout<<"++*"<<std::endl; 
  RooRealVar w_ttbar ("w_ttbar", "w_ttbar", weight_ttbar);
  //w_ttbar.setVal(2);
  RooFormulaVar wFunc_ttbar("wFunc_ttbar","event weight","w_ttbar",RooArgList(w_ttbar,ak08Pruned_1_mass)) ;
  RooFormulaVar *wFunc_ttbar1[10];
  std::cout<<"++*"<<std::endl;
  wFunc_ttbar1[0]= new RooFormulaVar("wFunc_ttbar1","event weight","w_ttbar",RooArgList(w_ttbar,ak08Pruned_1_mass)) ;
  std::cout<<"++*"<<std::endl;
  //wFunc_ttbar1[0]= new RooFormulaVar("wFunc_ttbar1_0","event weight","w_ttbar",RooArgList(w_ttbar,ak08Pruned_1_mass)) ;
  //wFunc_mio[0]=new RooFormulaVar("wFunc_mio0","wFunc_mio0","mioPeso0",RooArgList((const RooAbsArg)mioPeso[0],ak08Pruned_1_mass)) ;
  //RooRealVar* W_ttbar = (RooRealVar*) dsC_match->addColumn(((RooAbsArg &)*wFunc_ttbar1[0]));//wFunc_ttbar) ;
  RooRealVar* W_ttbar[2];
  W_ttbar[0]= (RooRealVar*) dsC_match->addColumn(((RooAbsArg &)*wFunc_ttbar1[0]));//wFunc_ttbar) ;
  std::cout<<"+++++++++++++++++++"<<((RooAbsArg &)*wFunc_ttbar1[0]).GetName()<<std::endl;
  //return;
  RooDataSet* wdata_ttbar_match= new RooDataSet(dsC_match->GetName(),dsC_match->GetTitle(),dsC_match,*dsC_match->get(),0,W_ttbar[0]->GetName()) ;
  RooRealVar* w1 = (RooRealVar*) dsC_unmatch->addColumn(wFunc_ttbar) ;
  RooDataSet wdata_ttbar_unmatch(dsC_unmatch->GetName(),dsC_unmatch->GetTitle(),dsC_unmatch,*dsC_unmatch->get(),0,w1->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
  std::cout<<"Fit def"<<std::endl;
  wdata_ttbar_match->Print();
  //wdata.append(wdata1);
  wdata_ttbar_unmatch.Print();

  //w_ttbar.setVal(2);
  //std::cout<<w_ttbar.getValV();
  //return;
  RooRealVar Nttbar_true("Nttbar_true","Nttbar_true", 1000,0,10000);
  RooRealVar mean_ttbar("mean_ttbar","mean",78,74,85) ;
  RooRealVar sigma_ttbar("sigma_ttbar","sigma",4,0,20) ;
  RooGaussian gx_ttbar("gx_ttbar","gx",ak08Pruned_1_mass,mean_ttbar,sigma_ttbar) ;
  RooExtendPdf gx_ttbar_norm("gx_ttbar_norm", "gx_ttbar_norm", gx_ttbar, Nttbar_true);
  RooRealVar NttbarBKG_true("NttbarBKG_true","NttbarBKG_true", 100,0,10000);
  RooRealVar mean_ttbar1("mean_ttbar1","mean",80,60,100) ;
  RooRealVar sigma_ttbar1("sigma_ttbar1","sigma",50,0,150) ;
  RooGaussian gx_ttbar1("gx_ttbar1","gx",ak08Pruned_1_mass,mean_ttbar1,sigma_ttbar1) ;
  RooExtendPdf gx_ttbar1_norm("gx_ttbar1_norm", "gx_ttbar1_norm", gx_ttbar1, NttbarBKG_true);
  RooRealVar k("k", "k", 100, 0,10000); 



  RooRealVar dCBPowL_ttbar("dCBPowL_ttbar","a_ttbar", 1., 0., 50.);
  RooRealVar dCBCutL_ttbar("dCBCutL_ttbar","a1_ttbar", 1., 0., 50.);//0.1,-1,1);
  //RooRealVar a2_ttbar("a2_ttbar","a2_ttbar", -.1,-10,10);
  RooRealVar a_ttbar("a_ttbar","a_pass_sim_MC", .1,-10,10);//1,-100,100);
  RooRealVar a1_ttbar("a1_ttbar","a1_pass_sim_MC", .1,-1,1);//0.1,-1,1);
  RooRealVar a2_ttbar("a2_ttbar","a2_pass_sim_MC", .1,-1,1);//0.1,-1,1);
  /*RooCBShape test_("test_","", ak08Pruned_1_mass_ttbar, mean_ttbar, sigma_ttbar, a1_ttbar, a_ttbar); 
  */RooChebychev cheby_ttbar("cheby_ttbar","cheby_ttbar",ak08Pruned_1_mass,RooArgSet(a_ttbar,a1_ttbar, a2_ttbar)) ;
  /*RooPolynomial p2_ttbar("p2_ttbar","p2",ak08Pruned_1_mass,RooArgList(a_ttbar,a1_ttbar, a2_ttbar),0) ;
  RooAddPdf model_ttbar("model_ttbar","model",RooArgList(gx_ttbar, cheby_ttbar), RooArgList(Nttbar_true,k));//, px_ctl)) ;  
  */RooRealVar dCBPowR_ttbar("dCBPowR_ttbar","a_ttbar",.3, 0.2, 15.);
  RooRealVar dCBCutR_ttbar("dCBCutR_ttbar","a1_ttbar", 0.5, 0.2, 2.);//0.1,-1,1);
  RooRealVar const_mean_pass_sim_MC("const_mean_pass_sim_MC", "Double CB mean", 75.5, 75,76); //85
  RooRealVar const_sigma_pass_sim_MC("const_sigma_pass_sim_MC", "Double CB Width", 10.,9,11); //working points 10, 0, 40 for tau21<.4
  RooRealVar dCBCutL_pass_sim_MC("dCBCutL", "Double CB Cut left", 0.48,0.3,0.52); //for tau21<0.4 1., 0.1, 50. works
  RooRealVar dCBCutR_pass_sim_MC("dCBCutR", "Double CB Cut right", 1.,0.9,1.1);
  RooRealVar dCBPowerL_pass_sim_MC("dCBPowerL", "Double CB Power left",  1.,0.9,1.1);
  RooRealVar dCBPowerR_pass_sim_MC("dCBPowerR", "Double CB Power right", 9.5,9.,10.); 
  RooDCBShape dcb("dcb", "double crystal ball", ak08Pruned_1_mass, mean_ttbar, sigma_ttbar, dCBCutL_ttbar, dCBPowL_ttbar, dCBCutR_ttbar,dCBPowR_ttbar); //const_mean_pass_sim_MC, const_sigma_pass_sim_MC,// const_dCBCutL_pass_sim_MC, const_dCBCutR_pass_sim_MC, const_dCBPowerL_pass_sim_MC, const_dCBPowerR_pass_sim_MC);
  

  RooRealVar dCBPowL_ttbar_unmatch("dCBPowL_ttbar_unmatch","a_ttbar", 1., 0., 50.);
  RooRealVar dCBCutL_ttbar_unmatch("dCBCutL_ttbar_unmatch","a1_ttbar", 1., 0., 50.);
  RooRealVar dCBPowR_ttbar_unmatch("dCBPowR_ttbar_unmatch","a_ttbar", 2., 0., 50.);
  RooRealVar dCBCutR_ttbar_unmatch("dCBCutR_ttbar_unmatch","a1_ttbar", 2., 0., 50.);
  RooRealVar mean_ttbar_unmatch("mean_ttbar_unmatch","mean",80,75,90) ;
  RooRealVar sigma_ttbar_unmatch("sigma_ttbar_unmatch","sigma",10,0,50) ;
  RooDCBShape dcb_unmatch("dcb_unmatch", "double crystal ball", ak08Pruned_1_mass, mean_ttbar_unmatch, sigma_ttbar_unmatch, dCBCutL_ttbar_unmatch, dCBPowL_ttbar_unmatch, dCBCutR_ttbar_unmatch,dCBPowR_ttbar_unmatch);
  //model_ttbar.fitTo(*dh_ttbarMatch) ;
  //test_.fitTo(*dh_ttbarMatch) ;
  RooRealVar n2_ttbar("n2_ttbar","n2_ttbar", 20,-100,100);
  //My_double_CB chib= My_double_CB("chib", "chib", ak08Pruned_1_mass, mean_ttbar, sigma_ttbar, a1_ttbar, a_ttbar, a2_ttbar, n2_ttbar);
  //dcb.fitTo(*dh_ttbarMatch) ;
  //cheby_ttbar.fitTo(*dh_ttbarMatch) ;
  //dcb_unmatch.fitTo(*dh_ttbarUnMatch);
  //dcb.fitTo(*dsC_match, RooFit::SumW2Error(kTRUE));
  dcb.fitTo(*wdata_ttbar_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2)) ;
  dcb.fitTo(*wdata_ttbar_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2)) ;
  dcb.fitTo(*wdata_ttbar_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2)) ;

  RooPlot* frame1_ttbar = ak08Pruned_1_mass.frame(Bins(nbins),Title("TopTop Pass/Match sample")) ;
  //dh_ttbarMatch->plotOn(frame1_ttbar, LineColor(2), DataError(RooAbsData::SumW2)) ;
  ////////
  //dsC->plotOn(frame1_ttbar);//,DataError(RooAbsData::SumW2)); 
  
  wdata_ttbar_match/*dsC_match*/->plotOn(frame1_ttbar,DataError(RooAbsData::SumW2));
  dcb.plotOn(frame1_ttbar);
  wdata_ttbar_unmatch.plotOn(frame1_ttbar,DataError(RooAbsData::SumW2), LineColor(2));
  

  //test_.plotOn(frame1_ttbar, LineStyle(8));
  //model_ttbar.plotOn(frame1_ttbar);//, LineStyle(8));
  //model_ttbar.plotOn(frame1_ttbar, Components(gx_ttbar), LineColor(kRed), Normalization(1.0,RooAbsReal::RelativeExpected), LineStyle(8));
  //model_ttbar.plotOn(frame1_ttbar, Components(cheby_ttbar), LineColor(kGreen), Normalization(1.0,RooAbsReal::RelativeExpected), LineStyle(8));
  //model_ttbar.paramOn(frame1_ttbar);
  RooPlot* frame2_ttbar = ak08Pruned_1_mass.frame(Bins(nbins),Title("TopTop Pass/UnMatch sample")) ;

  TCanvas* c_ttbar = new TCanvas("c_ttbar","c_ttbar",1);//,800,400) ;
  //c_ttbar->Divide(2) ;
  /*c_ttbar->cd(1) ;*/ gPad->SetLeftMargin(0.15) ; frame1_ttbar->GetYaxis()->SetTitleOffset(1.4) ; frame1_ttbar->Draw() ;
  //c_ttbar->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_ttbar->GetYaxis()->SetTitleOffset(1.4) ; frame2_ttbar->Draw() ;
  TCanvas* c_ttbar_fail = new TCanvas("c_ttbar_fail","c_ttbar_fail",1);//,800,400) ;
  gPad->SetLeftMargin(0.15) ; frame2_ttbar->GetYaxis()->SetTitleOffset(1.4) ; frame2_ttbar->Draw() ;

  //return;



  std::cout<<"+++++++++++++++++++"<<std::endl;
  std::cout<<"    Single Top     "<<std::endl;
  std::cout<<"+++++++++++++++++++"<<std::endl;
  TFile *file_Stop1;
  file_Stop1=TFile::Open("../ttbar_output_new/SingleTop_sChannel_20160703_141827/rootfile_SingleTop_sChannel__20160703_141827_0_reduced_skim.root");//"/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_s-channel_4f.root");//testMuonEff/firstTest_reduced_skim.root");
  double nEvt_Stop1= ((TH1D *)(file_Stop1->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_Stop1 = 47.13;
  TFile *file_Stop1_primo;
  file_Stop1_primo=TFile::Open("../ttbar_output_new/SingleTop_sChannel_20160703_141827/rootfile_SingleTop_sChannel__20160703_141827_0_reduced_skim_redu.root");//"/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_s-channel_4f_redu.root");
  TTree *tree_Stop1=(TTree *)file_Stop1_primo->Get("mio");
  
  RooDataSet ds_Stop1("ds_Stop1","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_Stop1));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  ds_ttbar->Print();
  
  RooDataSet* dsC_Stop1_match = (RooDataSet*) ds_Stop1.reduce(match.c_str()) ;
  RooDataSet* dsC_Stop1_unmatch = (RooDataSet*) ds_Stop1.reduce(unmatch.c_str()) ;
  dsC_Stop1_match->SetName("dsC_Stop1_match");
  dsC_Stop1_match->Print();
  dsC_Stop1_unmatch->SetName("dsC_Stop1_unmatch");
  dsC_Stop1_unmatch->Print();
  double weight_Stop1=xSec_Stop1*Lumi/nEvt_Stop1;
  if(!doWeight) weight_Stop1=1;
  RooConstVar w_Stop1 ("w_Stop1", "w_Stop1", weight_Stop1);
  RooFormulaVar wFunc_Stop1("wFunc_Stop1","event weight","w_Stop1",RooArgList(w_Stop1,ak08Pruned_1_mass)) ;
  RooRealVar* W_Stop1_match = (RooRealVar*) dsC_Stop1_match->addColumn(wFunc_Stop1) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_Stop1.getValV()<<std::endl;
  RooDataSet wdata_Stop1_match(dsC_Stop1_match->GetName(),dsC_Stop1_match->GetTitle(),dsC_Stop1_match,*dsC_Stop1_match->get(),0,W_Stop1_match->GetName()) ;
  RooRealVar* w1_Stop1 = (RooRealVar*) dsC_Stop1_unmatch->addColumn(wFunc_Stop1) ;
  RooDataSet wdata_Stop1_unmatch(dsC_Stop1_unmatch->GetName(),dsC_Stop1_unmatch->GetTitle(),dsC_Stop1_unmatch,*dsC_Stop1_unmatch->get(),0,w1_Stop1->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;



  TFile *file_Stop2;
  file_Stop2=TFile::Open("../ttbar_output_new/SingleTop_tW_20160703_141828/rootfile_SingleTop_tW__20160703_141828_0_reduced_skim.root");//"/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_tW-channel_top.root");//testMuonEff/firstTest_reduced_skim.root");
  double nEvt_Stop2= ((TH1D *)(file_Stop2->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_Stop2 = 35.6;
  TFile *file_Stop2_primo;
  file_Stop2_primo=TFile::Open("../ttbar_output_new/SingleTop_tW_20160703_141828/rootfile_SingleTop_tW__20160703_141828_0_reduced_skim_redu.root");
  TTree *tree_Stop2=(TTree *)file_Stop2_primo->Get("mio");
  
  RooDataSet ds_Stop2("ds_Stop2","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_Stop2));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_Stop2_match = (RooDataSet*) ds_Stop2.reduce(match.c_str()) ;
  RooDataSet* dsC_Stop2_unmatch = (RooDataSet*) ds_Stop2.reduce(unmatch.c_str()) ;
  dsC_Stop2_match->SetName("dsC_Stop2_match");
  dsC_Stop2_match->Print();
  dsC_Stop2_unmatch->SetName("dsC_Stop2_unmatch");
  dsC_Stop2_unmatch->Print();
  double weight_Stop2=xSec_Stop2*Lumi/nEvt_Stop2;
  if(!doWeight) weight_Stop2=1;
  RooConstVar w_Stop2 ("w_Stop2", "w_Stop2", weight_Stop2);
  RooFormulaVar wFunc_Stop2("wFunc_Stop2","event weight","w_Stop2",RooArgList(w_Stop2,ak08Pruned_1_mass)) ;
  RooRealVar* W_Stop2_match = (RooRealVar*) dsC_Stop2_match->addColumn(wFunc_Stop2) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_Stop2.getValV()<<std::endl;
  RooDataSet wdata_Stop2_match(dsC_Stop2_match->GetName(),dsC_Stop2_match->GetTitle(),dsC_Stop2_match,*dsC_Stop2_match->get(),0,W_Stop2_match->GetName()) ;
  RooRealVar* w1_Stop2 = (RooRealVar*) dsC_Stop2_unmatch->addColumn(wFunc_Stop2) ;
  RooDataSet wdata_Stop2_unmatch(dsC_Stop2_unmatch->GetName(),dsC_Stop2_unmatch->GetTitle(),dsC_Stop2_unmatch,*dsC_Stop2_unmatch->get(),0,w1_Stop2->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;

  TFile *file_Stop3; 
  file_Stop3=TFile::Open("../ttbar_output_new/SingleTop_tbarW_20160703_141828/rootfile_SingleTop_tbarW__20160703_141828_0_reduced_skim.root");//testMuonEff/firstTest_reduced_skim.root");
  double nEvt_Stop3= ((TH1D *)(file_Stop3->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_Stop3 = 35.6;
  TFile *file_Stop3_primo;
  file_Stop3_primo=TFile::Open("../ttbar_output_new/SingleTop_tbarW_20160703_141828/rootfile_SingleTop_tbarW__20160703_141828_0_reduced_skim_redu.root");
  TTree *tree_Stop3=(TTree *)file_Stop3_primo->Get("mio");
  
  RooDataSet ds_Stop3("ds_Stop3","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_Stop3));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_Stop3_match = (RooDataSet*) ds_Stop3.reduce(match.c_str()) ;
  RooDataSet* dsC_Stop3_unmatch = (RooDataSet*) ds_Stop3.reduce(unmatch.c_str()) ;
  dsC_Stop3_match->SetName("dsC_Stop3_match");
  dsC_Stop3_match->Print();
  dsC_Stop3_unmatch->SetName("dsC_Stop3_unmatch");
  dsC_Stop3_unmatch->Print();
  double weight_Stop3=xSec_Stop3*Lumi/nEvt_Stop3;
  if(!doWeight) weight_Stop3=1;
  RooConstVar w_Stop3 ("w_Stop3", "w_Stop3", weight_Stop3);
  RooFormulaVar wFunc_Stop3("wFunc_Stop3","event weight","w_Stop3",RooArgList(w_Stop3,ak08Pruned_1_mass)) ;
  RooRealVar* W_Stop3_match = (RooRealVar*) dsC_Stop3_match->addColumn(wFunc_Stop3) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_Stop3.getValV()<<std::endl;
  RooDataSet wdata_Stop3_match(dsC_Stop3_match->GetName(),dsC_Stop3_match->GetTitle(),dsC_Stop3_match,*dsC_Stop3_match->get(),0,W_Stop3_match->GetName()) ;
  RooRealVar* w1_Stop3 = (RooRealVar*) dsC_Stop3_unmatch->addColumn(wFunc_Stop3) ;
  RooDataSet wdata_Stop3_unmatch(dsC_Stop3_unmatch->GetName(),dsC_Stop3_unmatch->GetTitle(),dsC_Stop3_unmatch,*dsC_Stop3_unmatch->get(),0,w1_Stop3->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;



///////  TFile *file_Stop4; 
///////  file_Stop4=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_t-channel_antitop.root");//testMuonEff/firstTest_reduced_skim.root");
///////  double nEvt_Stop4= ((TH1D *)(file_Stop4->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
///////  double xSec_Stop4 = 26;
///////  TFile *file_Stop4_primo;
///////  file_Stop4_primo=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_t-channel_antitop_redu.root");
///////  TTree *tree_Stop4=(TTree *)file_Stop4_primo->Get("mio");
///////  
///////  RooDataSet ds_Stop4("ds_Stop4","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_Stop4));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
///////      
///////  RooDataSet* dsC_Stop4_match = (RooDataSet*) ds_Stop4.reduce(match.c_str()) ; 
///////  RooDataSet* dsC_Stop4_unmatch = (RooDataSet*) ds_Stop4.reduce(unmatch.c_str()) ;
///////  dsC_Stop4_match->SetName("dsC_Stop4_match");
///////  dsC_Stop4_match->Print();
///////  dsC_Stop4_unmatch->SetName("dsC_Stop4_unmatch");
///////  dsC_Stop4_unmatch->Print();
///////  double weight_Stop4=xSec_Stop4*Lumi/nEvt_Stop4;
///////  RooConstVar w_Stop4 ("w_Stop4", "w_Stop4", weight_Stop4);
///////  RooFormulaVar wFunc_Stop4("wFunc_Stop4","event weight","w_Stop4",RooArgList(w_Stop4,ak08Pruned_1_mass)) ;
///////  RooRealVar* W_Stop4_match = (RooRealVar*) dsC_Stop4_match->addColumn(wFunc_Stop4) ;
///////  std::cout<<"+++++++++++++++++++"<<wFunc_Stop4.getValV()<<std::endl;
///////  RooDataSet wdata_Stop4_match(dsC_Stop4_match->GetName(),dsC_Stop4_match->GetTitle(),dsC_Stop4_match,*dsC_Stop4_match->get(),0,W_Stop4_match->GetName()) ;
///////  RooRealVar* w1_Stop4 = (RooRealVar*) dsC_Stop4_unmatch->addColumn(wFunc_Stop4) ;
///////  RooDataSet wdata_Stop4_unmatch(dsC_Stop4_unmatch->GetName(),dsC_Stop4_unmatch->GetTitle(),dsC_Stop4_unmatch,*dsC_Stop4_unmatch->get(),0,w1_Stop4->GetName()) ;
///////  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
///////  
///////  TFile *file_Stop5; 
///////  file_Stop5=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_t-channel_top.root");//testMuonEff/firstTest_reduced_skim.root");
///////  double nEvt_Stop5= ((TH1D *)(file_Stop5->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
///////  double xSec_Stop5 = 43.8;
///////  TFile *file_Stop5_primo;
///////  file_Stop5_primo=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_t-channel_top_redu.root");
///////  TTree *tree_Stop5=(TTree *)file_Stop5_primo->Get("mio");
///////  
///////  RooDataSet ds_Stop5("ds_Stop5","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_Stop5));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
///////      
///////  RooDataSet* dsC_Stop5_match = (RooDataSet*) ds_Stop5.reduce(match.c_str()) ; 
///////  RooDataSet* dsC_Stop5_unmatch = (RooDataSet*) ds_Stop5.reduce(unmatch.c_str()) ;
///////  dsC_Stop5_match->SetName("dsC_Stop5_match");
///////  dsC_Stop5_match->Print();
///////  dsC_Stop5_unmatch->SetName("dsC_Stop5_unmatch");
///////  dsC_Stop5_unmatch->Print();
///////  double weight_Stop5=xSec_Stop5*Lumi/nEvt_Stop5;
///////  RooConstVar w_Stop5 ("w_Stop5", "w_Stop5", weight_Stop5);
///////  RooFormulaVar wFunc_Stop5("wFunc_Stop5","event weight","w_Stop5",RooArgList(w_Stop5,ak08Pruned_1_mass)) ;
///////  RooRealVar* W_Stop5_match = (RooRealVar*) dsC_Stop5_match->addColumn(wFunc_Stop5) ;
///////  std::cout<<"+++++++++++++++++++"<<wFunc_Stop5.getValV()<<std::endl;
///////  RooDataSet wdata_Stop5_match(dsC_Stop5_match->GetName(),dsC_Stop5_match->GetTitle(),dsC_Stop5_match,*dsC_Stop5_match->get(),0,W_Stop5_match->GetName()) ;
///////  RooRealVar* w1_Stop5 = (RooRealVar*) dsC_Stop5_unmatch->addColumn(wFunc_Stop5) ;
///////  RooDataSet wdata_Stop5_unmatch(dsC_Stop5_unmatch->GetName(),dsC_Stop5_unmatch->GetTitle(),dsC_Stop5_unmatch,*dsC_Stop5_unmatch->get(),0,w1_Stop5->GetName()) ;
///////  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;

  
  
  
////  RooDataSet *nuovo;
////  nuovo= new RooDataSet(dsC_Stop5_unmatch->GetName(),dsC_Stop5_unmatch->GetTitle(),dsC_Stop5_unmatch,*dsC_Stop5_unmatch->get(),0,w1_Stop5->GetName()) ;
////  nuovo->append(*wdata_ttbar_match); 
  //return;
  dcb.fitTo(*wdata_ttbar_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2)) ;
  dcb.fitTo(*wdata_ttbar_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2)) ;
  dcb.fitTo(*wdata_ttbar_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2)) ;
  RooPlot* testMio= ak08Pruned_1_mass.frame(Bins(nbins),Title("total Match sample")) ;
  //ds_ttbarMatch->plotOn(testMio, DataError(RooAbsData::SumW2));
  wdata_ttbar_match->plotOn(testMio,DataError(RooAbsData::SumW2));
  //nuovo->plotOn(testMio,LineColor(2),DataError(RooAbsData::SumW2));
  dcb.plotOn(testMio);
  dcb.paramOn(testMio);
  TCanvas* c_new = new TCanvas("c_new","all matched",1);//,800,400) ;
  gPad->SetLeftMargin(0.15) ; testMio->GetYaxis()->SetTitleOffset(1.4) ; testMio->Draw() ;

  //end new addition
  //return;
  RooRealVar mean_Stop("mean_Stop","mean",80,60,100) ;
  RooRealVar sigma_Stop("sigma_Stop","sigma",2,0,8) ;
  RooGaussian gx_Stop("gx_Stop","gx",ak08Pruned_1_mass,mean_Stop,sigma_Stop) ;
  RooRealVar NStop_true("NStop_true","NStop_true", 100,0,10000);
  RooExtendPdf gx_Stop_norm("gx_Stop_norm", "gx_Stop_norm", gx_Stop, NStop_true);
  RooRealVar mean_Stop1("mean_Stop1","mean",80,60,100) ;
  RooRealVar sigma_Stop1("sigma_Stop1","sigma",50,0,100) ;
  RooGaussian gx_Stop1("gx_Stop1","gx",ak08Pruned_1_mass,mean_Stop1,sigma_Stop1) ;
  RooRealVar k_Stop("k_Stop", "k", 100, 0,10000);
  RooRealVar a_Stop("a_Stop","a_Stop", 1,-100,100);
  RooRealVar a1_Stop("a1_Stop","a1_Stop", 0.1,-1,1);
  RooRealVar a2_Stop("a2_Stop","a2_Stop", -0.1,-1,1);
  //RooPolynomial p2_ttbar("p2_ttbar","p2",ak08Pruned_1_mass,RooArgList(a_ttbar,a1_ttbar, a2_ttbar),0) ;
  RooChebychev cheby_Stop("cheby_Stop","cheby_Stop",ak08Pruned_1_mass,RooArgSet(a_Stop,a1_Stop, a2_Stop)) ;

  RooAddPdf model_Stop("model_Stop","model",RooArgList(gx_Stop, cheby_Stop/*gx_Stop1*/), RooArgList(NStop_true,k_Stop));//, px_ctl)) ;
  wdata_Stop1_match.append(wdata_Stop2_match);
  wdata_Stop1_match.append(wdata_Stop3_match);
///  wdata_Stop1_match.append(wdata_Stop4_match);
  model_Stop.fitTo(wdata_Stop1_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2));
  model_Stop.fitTo(wdata_Stop1_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2));
  model_Stop.fitTo(wdata_Stop1_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2));

  //model_Stop.fitTo(*dh_StopMatch) ;
  //gx_Stop_norm.fitTo(*dh_StopMatch) ;

  //Plotting
  RooPlot* frame1_Stop = ak08Pruned_1_mass.frame(Bins(nbins),Title("Single Top Pass/Match sample")) ;
  //gx_Stop.plotOn(frame1_Stop, LineColor(kRed)) ;
  //cheby_Stop.plotOn(frame1_Stop, LineColor(kGreen)) ;
  //gx_Stop1.plotOn(frame1_Stop, LineColor(kGreen)) ;
  wdata_Stop1_match.plotOn(frame1_Stop,DataError(RooAbsData::SumW2));
  wdata_Stop2_match.plotOn(frame1_Stop,DataError(RooAbsData::SumW2),  LineColor(kRed), MarkerColor(kRed));
  model_Stop.plotOn(frame1_Stop) ;
  model_Stop.plotOn(frame1_Stop, Components(gx_Stop), LineColor(kRed), Normalization(1.0,RooAbsReal::RelativeExpected));
  model_Stop.plotOn(frame1_Stop, Components(cheby_Stop), LineColor(kGreen), Normalization(1.0,RooAbsReal::RelativeExpected));
  model_Stop.paramOn(frame1_Stop);

//  gx_Stop1.plotOn(frame1_Stop, LineColor(kGreen));
//  model_Stop.plotOn(frame1_Stop) ;
////////  combData.plotOn(frame1,Cut("sample==sample::physics")) ;
////  simPdf_Stop.plotOn(frame1_Stop,Slice(sampleStop,"pass"),ProjWData(sampleStop,combDataStop)) ; 
//////  //gx.plotOn(frame1,Slice(sample,"physics"),ProjWData(sample,combData), LineColor(kGreen)) ;
//////  //gx1.plotOn(frame1,Slice(sample,"physics"),ProjWData(sample,combData), LineColor(kRed)) ;
  RooPlot* frame2_Stop = ak08Pruned_1_mass.frame(Bins(nbins),Title("Single Top Pass/UnMatch sample")) ;

  TCanvas* c_Stop = new TCanvas("c_Stop","c_Stop",1);//800,400) ;
  //c_Stop->Divide(2) ;
  /*c_Stop->cd(1) ;*/ gPad->SetLeftMargin(0.15) ; frame1_Stop->GetYaxis()->SetTitleOffset(1.4) ; frame1_Stop->Draw() ;
//  c_Stop->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_Stop->GetYaxis()->SetTitleOffset(1.4) ; frame2_Stop->Draw() ;
   //return; 
   TCanvas* c_Stop_fail = new TCanvas("c_Stop_fail","c_Stop_unmatch",1);//,800,400) ;
   gPad->SetLeftMargin(0.15) ; frame2_Stop->GetYaxis()->SetTitleOffset(1.4) ; frame2_Stop->Draw() ;


  std::cout<<"+++++++++++++++++++"<<std::endl;
  std::cout<<"      DiBoson      "<<std::endl;
  std::cout<<"+++++++++++++++++++"<<std::endl;
  TFile *file_WW;
  file_WW=TFile::Open("../ttbar_output_new/WW_20160703_141826/rootfile_ww__20160703_141826_0_reduced_skim.root");//testMuonEff/firstTest_reduced_skim.root");
  double nEvt_WW= ((TH1D *)(file_WW->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_WW = 118.7;
  TFile *file_WW1;
  file_WW1=TFile::Open("../ttbar_output_new/WW_20160703_141826/rootfile_ww__20160703_141826_0_reduced_skim_redu.root");
  TTree *tree_WW=(TTree *)file_WW1->Get("mio");
  
  RooDataSet ds_WW("ds_WW","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_WW));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_WW = (RooDataSet*) ds_WW.reduce(pass.c_str()) ;
  dsC_WW->SetName("dsC_WW");
  dsC_WW->Print();
  double weight_WW=xSec_WW*Lumi/nEvt_WW;
  if(!doWeight) weight_WW=1;
  RooConstVar w_WW ("w_WW", "w_WW", weight_WW);
  RooFormulaVar wFunc_WW("wFunc_WW","event weight","w_WW",RooArgList(w_WW,ak08Pruned_1_mass)) ;
  RooRealVar* W_WW = (RooRealVar*) dsC_WW->addColumn(wFunc_WW) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_WW.getValV()<<std::endl;
  RooDataSet wdata_WW(dsC_WW->GetName(),dsC_WW->GetTitle(),dsC_WW,*dsC_WW->get(),0,W_WW->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
  

  TFile *file_WZ;
  file_WZ=TFile::Open("../ttbar_output_new/WZ_20160703_141826/rootfile_wz__20160703_141826_0_reduced_skim.root");//testMuonEff/firstTest_reduced_skim.root");
  double nEvt_WZ= ((TH1D *)(file_WZ->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_WZ = 16.5;
  TFile *file_WZ1;
  file_WZ1=TFile::Open("../ttbar_output_new/WZ_20160703_141826/rootfile_wz__20160703_141826_0_reduced_skim_redu.root");
  TTree *tree_WZ=(TTree *)file_WZ1->Get("mio");
 
  RooDataSet ds_WZ("ds_WZ","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_WZ));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_WZ = (RooDataSet*) ds_WZ.reduce(pass.c_str()) ; 
  dsC_WZ->SetName("dsC_WZ");
  dsC_WZ->Print();
  double weight_WZ=xSec_WZ*Lumi/nEvt_WZ;
  if(!doWeight) weight_WZ=1;
  RooConstVar w_WZ ("w_WZ", "w_WZ", weight_WZ);
  RooFormulaVar wFunc_WZ("wFunc_WZ","event weight","w_WZ",RooArgList(w_WZ,ak08Pruned_1_mass)) ;
  RooRealVar* W_WZ = (RooRealVar*) dsC_WZ->addColumn(wFunc_WZ) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_WZ.getValV()<<std::endl;
  RooDataSet wdata_WZ(dsC_WZ->GetName(),dsC_WZ->GetTitle(),dsC_WZ,*dsC_WZ->get(),0,W_WZ->GetName()) ; 
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;


  TFile *file_ZZ;
  file_ZZ=TFile::Open("../ttbar_output_new/ZZ_20160703_141827/rootfile_zz__20160703_141827_0_reduced_skim.root");//testMuonEff/firstTest_re    duced_skim.r        oot");
  double nEvt_ZZ= ((TH1D *)(file_ZZ->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_ZZ = 47.13;
  TFile *file_ZZ1;
  file_ZZ1=TFile::Open("../ttbar_output_new/ZZ_20160703_141827/rootfile_zz__20160703_141827_0_reduced_skim_redu.root");
  TTree *tree_ZZ=(TTree *)file_ZZ1->Get("mio");
  
  RooDataSet ds_ZZ("ds_ZZ","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_ZZ));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_ZZ = (RooDataSet*) ds_ZZ.reduce(pass.c_str()) ; 
  dsC_ZZ->SetName("dsC_ZZ");
  dsC_ZZ->Print();
  double weight_ZZ=xSec_ZZ*Lumi/nEvt_ZZ;
  if(!doWeight) weight_ZZ=1;
  RooConstVar w_ZZ ("w_ZZ", "w_ZZ", weight_ZZ);
  RooFormulaVar wFunc_ZZ("wFunc_ZZ","event weight","w_ZZ",RooArgList(w_ZZ,ak08Pruned_1_mass)) ;
  RooRealVar* W_ZZ = (RooRealVar*) dsC_ZZ->addColumn(wFunc_ZZ) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_ZZ.getValV()<<std::endl;
  RooDataSet wdata_ZZ(dsC_ZZ->GetName(),dsC_ZZ->GetTitle(),dsC_ZZ,*dsC_ZZ->get(),0,W_ZZ->GetName()) ; 
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;

  
  
  

  std::cout<<"+++++++++++++++++++"<<std::endl;
  std::cout<<"      W+Jets       "<<std::endl;
  std::cout<<"+++++++++++++++++++"<<std::endl;
  TFile *file_W_100_200;
  file_W_100_200=TFile::Open("../ttbar_output_new/WJets_100_200_20160703_141822/rootfile_wjets100_200__20160703_141822_0_reduced_skim.root");
  double nEvt_W_100_200= ((TH1D *)(file_W_100_200->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_100_200 = 1629.87;
  TFile *file_W_100_200_1;
  file_W_100_200_1=TFile::Open("../ttbar_output_new/WJets_100_200_20160703_141822/rootfile_wjets100_200__20160703_141822_0_reduced_skim_redu.root");
  TTree *tree_W_100_200=(TTree *)file_W_100_200_1->Get("mio");
  
  RooDataSet ds_W_100_200("ds_W_100_200","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_W_100_200));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_W_100_200 = (RooDataSet*) ds_W_100_200.reduce(pass.c_str()) ;
  dsC_W_100_200->SetName("dsC_W_100_200");
  dsC_W_100_200->Print();
  double weight_W_100_200=xSec_W_100_200*Lumi/nEvt_W_100_200;
  if(!doWeight) weight_W_100_200=1;
  RooConstVar w_W_100_200 ("w_W_100_200", "w_W_100_200", weight_W_100_200);
  RooFormulaVar wFunc_W_100_200("wFunc_W_100_200","event weight","w_W_100_200",RooArgList(w_W_100_200,ak08Pruned_1_mass)) ;
  RooRealVar* W_W_100_200 = (RooRealVar*) dsC_W_100_200->addColumn(wFunc_W_100_200) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_W_100_200.getValV()<<std::endl;
  RooDataSet wdata_W_100_200(dsC_W_100_200->GetName(),dsC_W_100_200->GetTitle(),dsC_W_100_200,*dsC_W_100_200->get(),0,W_W_100_200->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;


  TFile *file_W_200_400;
  file_W_200_400=TFile::Open("../ttbar_output_new/WJets_200_400_20160703_141824/rootfile_wjets200_400__20160703_141824_0_reduced_skim.root");
  double nEvt_W_200_400= ((TH1D *)(file_W_200_400->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_200_400 = 435.6;
  TFile *file_W_200_400_1;
  file_W_200_400_1=TFile::Open("../ttbar_output_new/WJets_200_400_20160703_141824/rootfile_wjets200_400__20160703_141824_0_reduced_skim_redu.root");
  TTree *tree_W_200_400=(TTree *)file_W_200_400_1->Get("mio");

  RooDataSet ds_W_200_400("ds_W_200_400","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_W_200_400));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_W_200_400 = (RooDataSet*) ds_W_200_400.reduce(pass.c_str()) ; 
  dsC_W_200_400->SetName("dsC_W_200_400");
  dsC_W_200_400->Print();
  double weight_W_200_400=xSec_W_200_400*Lumi/nEvt_W_200_400;
  if(!doWeight) weight_W_200_400=1;
  RooConstVar w_W_200_400 ("w_W_200_400", "w_W_200_400", weight_W_200_400);
  RooFormulaVar wFunc_W_200_400("wFunc_W_200_400","event weight","w_W_200_400",RooArgList(w_W_200_400,ak08Pruned_1_mass)) ;
  RooRealVar* W_W_200_400 = (RooRealVar*) dsC_W_200_400->addColumn(wFunc_W_200_400) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_W_200_400.getValV()<<std::endl;
  RooDataSet wdata_W_200_400(dsC_W_200_400->GetName(),dsC_W_200_400->GetTitle(),dsC_W_200_400,*dsC_W_200_400->get(),0,W_W_200_400->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;


  TFile *file_W_400_600;
  file_W_400_600=TFile::Open("../ttbar_output_new/WJets_400_600_20160703_141825/rootfile_wjets400_600__20160703_141825_0_reduced_skim.root");
  double nEvt_W_400_600= ((TH1D *)(file_W_400_600->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_400_600 = 56.17;
  TFile *file_W_400_600_1;
  file_W_400_600_1=TFile::Open("../ttbar_output_new/WJets_400_600_20160703_141825/rootfile_wjets400_600__20160703_141825_0_reduced_skim_redu.root");
  TTree *tree_W_400_600=(TTree *)file_W_400_600_1->Get("mio");
  
  RooDataSet ds_W_400_600("ds_W_400_600","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_W_400_600));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_W_400_600 = (RooDataSet*) ds_W_400_600.reduce(pass.c_str()) ;
  dsC_W_400_600->SetName("dsC_W_400_600");
  dsC_W_400_600->Print();
  double weight_W_400_600=xSec_W_400_600*Lumi/nEvt_W_400_600;
  if(!doWeight) weight_W_400_600=1;
  RooConstVar w_W_400_600 ("w_W_400_600", "w_W_400_600", weight_W_400_600);
  RooFormulaVar wFunc_W_400_600("wFunc_W_400_600","event weight","w_W_400_600",RooArgList(w_W_400_600,ak08Pruned_1_mass)) ;
  RooRealVar* W_W_400_600 = (RooRealVar*) dsC_W_400_600->addColumn(wFunc_W_400_600) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_W_400_600.getValV()<<std::endl;
  RooDataSet wdata_W_400_600(dsC_W_400_600->GetName(),dsC_W_400_600->GetTitle(),dsC_W_400_600,*dsC_W_400_600->get(),0,W_W_400_600->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;


  TFile *file_W_600_800;
  file_W_600_800=TFile::Open("../ttbar_output_new/WJets_600_800_20160703_141825/rootfile_wjets600_800__20160703_141825_0_reduced_skim.root");
  double nEvt_W_600_800= ((TH1D *)(file_W_600_800->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_600_800 = 14.61;
  TFile *file_W_600_800_1;
  file_W_600_800_1=TFile::Open("../ttbar_output_new/WJets_600_800_20160703_141825/rootfile_wjets600_800__20160703_141825_0_reduced_skim_redu.root");
  TTree *tree_W_600_800=(TTree *)file_W_600_800_1->Get("mio");
  
  RooDataSet ds_W_600_800("ds_W_600_800","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_W_600_800));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_W_600_800 = (RooDataSet*) ds_W_600_800.reduce(pass.c_str()) ;
  dsC_W_600_800->SetName("dsC_W_600_800");
  dsC_W_600_800->Print();
  double weight_W_600_800=xSec_W_600_800*Lumi/nEvt_W_600_800;
  if(!doWeight) weight_W_600_800=1;
  RooConstVar w_W_600_800 ("w_W_600_800", "w_W_600_800", weight_W_600_800);
  RooFormulaVar wFunc_W_600_800("wFunc_W_600_800","event weight","w_W_600_800",RooArgList(w_W_600_800,ak08Pruned_1_mass)) ;
  RooRealVar* W_W_600_800 = (RooRealVar*) dsC_W_600_800->addColumn(wFunc_W_600_800) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_W_600_800.getValV()<<std::endl;
  RooDataSet wdata_W_600_800(dsC_W_600_800->GetName(),dsC_W_600_800->GetTitle(),dsC_W_600_800,*dsC_W_600_800->get(),0,W_W_600_800->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;



  TFile *file_W_800_1200;
  file_W_800_1200=TFile::Open("../ttbar_output_new/WJets_800_1200_20160703_141825/rootfile_wjets800_1200__20160703_141825_0_reduced_skim.root");
  double nEvt_W_800_1200= ((TH1D *)(file_W_800_1200->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_800_1200 = 6.36;
  TFile *file_W_800_1200_1;
  file_W_800_1200_1=TFile::Open("../ttbar_output_new/WJets_800_1200_20160703_141825/rootfile_wjets800_1200__20160703_141825_0_reduced_skim_redu.root");
  TTree *tree_W_800_1200=(TTree *)file_W_800_1200_1->Get("mio");
  
  RooDataSet ds_W_800_1200("ds_W_800_1200","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_W_800_1200));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_W_800_1200 = (RooDataSet*) ds_W_800_1200.reduce(pass.c_str()) ;
  dsC_W_800_1200->SetName("dsC_W_800_1200");
  dsC_W_800_1200->Print();
  double weight_W_800_1200=xSec_W_800_1200*Lumi/nEvt_W_800_1200;
  if(!doWeight) weight_W_800_1200=1;
  RooConstVar w_W_800_1200 ("w_W_800_1200", "w_W_800_1200", weight_W_800_1200);
  RooFormulaVar wFunc_W_800_1200("wFunc_W_800_1200","event weight","w_W_800_1200",RooArgList(w_W_800_1200,ak08Pruned_1_mass)) ;
  RooRealVar* W_W_800_1200 = (RooRealVar*) dsC_W_800_1200->addColumn(wFunc_W_800_1200) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_W_800_1200.getValV()<<std::endl;
  RooDataSet wdata_W_800_1200(dsC_W_800_1200->GetName(),dsC_W_800_1200->GetTitle(),dsC_W_800_1200,*dsC_W_800_1200->get(),0,W_W_800_1200->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;


  TFile *file_W_1200_2500;
  file_W_1200_2500=TFile::Open("../ttbar_output_new/WJets_1200_2500_20160703_141826/rootfile_wjets1200_2500__20160703_141826_0_reduced_skim.root");
  double nEvt_W_1200_2500= ((TH1D *)(file_W_1200_2500->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_1200_2500 = 1.61;
  TFile *file_W_1200_2500_1;
  file_W_1200_2500_1=TFile::Open("../ttbar_output_new/WJets_1200_2500_20160703_141826/rootfile_wjets1200_2500__20160703_141826_0_reduced_skim_redu.root");
  TTree *tree_W_1200_2500=(TTree *)file_W_1200_2500_1->Get("mio");
  
  RooDataSet ds_W_1200_2500("ds_W_1200_2500","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_W_1200_2500));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_W_1200_2500 = (RooDataSet*) ds_W_1200_2500.reduce(pass.c_str()) ;
  dsC_W_1200_2500->SetName("dsC_W_1200_2500");
  dsC_W_1200_2500->Print();
  double weight_W_1200_2500=xSec_W_1200_2500*Lumi/nEvt_W_1200_2500;
  if(!doWeight) weight_W_1200_2500=1;
  RooConstVar w_W_1200_2500 ("w_W_1200_2500", "w_W_1200_2500", weight_W_1200_2500);
  RooFormulaVar wFunc_W_1200_2500("wFunc_W_1200_2500","event weight","w_W_1200_2500",RooArgList(w_W_1200_2500,ak08Pruned_1_mass)) ;
  RooRealVar* W_W_1200_2500 = (RooRealVar*) dsC_W_1200_2500->addColumn(wFunc_W_1200_2500) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_W_1200_2500.getValV()<<std::endl;
  RooDataSet wdata_W_1200_2500(dsC_W_1200_2500->GetName(),dsC_W_1200_2500->GetTitle(),dsC_W_1200_2500,*dsC_W_1200_2500->get(),0,W_W_1200_2500->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;


  TFile *file_W_2500_Inf;
  file_W_2500_Inf=TFile::Open("../ttbar_output_new/WJets_2500_Inf_20160703_141826/rootfile_wjets2500_Inf__20160703_141826_0_reduced_skim.root");
  double nEvt_W_2500_Inf= ((TH1D *)(file_W_2500_Inf->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_2500_Inf = 37;
  TFile *file_W_2500_Inf_1;
  file_W_2500_Inf_1=TFile::Open("../ttbar_output_new/WJets_2500_Inf_20160703_141826/rootfile_wjets2500_Inf__20160703_141826_0_reduced_skim_redu.root");
  TTree *tree_W_2500_Inf=(TTree *)file_W_2500_Inf_1->Get("mio");
  
  RooDataSet ds_W_2500_Inf("ds_W_2500_Inf","ds_ttbar",RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21,genW_genBquark1_DR,genW_genBquark2_DR,ak08Ungroomed_WGen_DR),Import(*tree_W_2500_Inf));//,Cut(Form("%s", cut_matched_pass.c_str())));//, (xSec_ttbar*Lumi/nEvt_ttbar)) ;
  
  RooDataSet* dsC_W_2500_Inf = (RooDataSet*) ds_W_2500_Inf.reduce(pass.c_str()) ;
  dsC_W_2500_Inf->SetName("dsC_W_2500_Inf");
  dsC_W_2500_Inf->Print();
  double weight_W_2500_Inf=xSec_W_2500_Inf*Lumi/nEvt_W_2500_Inf;
  if(!doWeight) weight_W_2500_Inf=1;
  RooConstVar w_W_2500_Inf ("w_W_2500_Inf", "w_W_2500_Inf", weight_W_2500_Inf);
  RooFormulaVar wFunc_W_2500_Inf("wFunc_W_2500_Inf","event weight","w_W_2500_Inf",RooArgList(w_W_2500_Inf,ak08Pruned_1_mass)) ;
  RooRealVar* W_W_2500_Inf = (RooRealVar*) dsC_W_2500_Inf->addColumn(wFunc_W_2500_Inf) ;
  std::cout<<"+++++++++++++++++++"<<wFunc_W_2500_Inf.getValV()<<std::endl;
  RooDataSet wdata_W_2500_Inf(dsC_W_2500_Inf->GetName(),dsC_W_2500_Inf->GetTitle(),dsC_W_2500_Inf,*dsC_W_2500_Inf->get(),0,W_W_2500_Inf->GetName()) ;
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;


  TCanvas* c_Wjet = new TCanvas("c_Wjet","c_Wjet",1);//,800,400) ;

  RooRealVar a_pass_bkg("a_pass_bkg","a_pass1", 1,-100,100);
  RooRealVar a1_pass_bkg("a1_pass_bkg","a1_pass1", -0.1,-100,100);
  RooRealVar a2_pass_bkg("a2_pass_bkg","a2_pass1", -0.1,-100,100);
  /*RooRealVar a_pass1("a_pass1","a_Stop", 1,-1,1);
  RooRealVar a1_pass1("a1_pass1","a1_Stop", 0.1,-1,1);
  RooRealVar a2_pass1("a2_pass1","a2_Stop", -0.1,-1,1);
  */RooChebychev cheby_pass_bkg("cheby_pass_bkg","cheby_pass",ak08Pruned_1_mass,RooArgSet(a_pass_bkg,a1_pass_bkg));//, a2_pass_bkg)) ;
  RooRealVar a_bkg("a_bkg","a_bkg", 0.1,-1,1);
  RooRealVar a1_bkg("a1_bkg","a1_bkg", 0.1,-1,1);
  RooRealVar a2_bkg("a2_bkg","a2_bkg", 0.1,-1,1);
  //RooRealVar ap_bkg("ap_bkg","ap_bkg", 0.5, -10,10);
  RooRealVar mean_bkg("mean_bkg","mean_bkg", 80, 60,95);
  RooRealVar sigma_bkg("sigma_bkg","sigma_bkg", 5, 0,15);
  RooGaussian gx_bkg("gx_bkg","gx_bkg",ak08Pruned_1_mass,mean_bkg,sigma_bkg) ;
  RooRealVar mean_bkg1("mean_bkg1","mean_bkg1", 80, 60,95);
  RooRealVar sigma_bkg1("sigma_bkg1","sigma_bkg1", 30, 0,100);
  RooGaussian gx_bkg1("gx_bkg1","gx_bkg1",ak08Pruned_1_mass,mean_bkg1,sigma_bkg1) ;
  RooRealVar k_bkg("k_bkg","k_bkg", 2, 0,5);
  RooAddPdf model_bkg("model_bkg", "model_bkg", RooArgList(gx_bkg,gx_bkg1), RooArgList(k_bkg));
  RooChebychev p2("p2","p2",ak08Pruned_1_mass,RooArgSet(a_bkg,a1_bkg,a2_bkg));//,0) ;
  //RooGenericPdf passBkg("passBkg", "passBkg", "(1+(TMath::Erf((ak08Pruned_1_mass-a_bkg)/b_bkg)))*0.5*exp(-(mean_bkg-ak08Pruned_1_mass)*(mean_bkg-ak08Pruned_1_mass)/(2*sigma_bkg*sigma_bkg))", RooArgSet(ak08Pruned_1_mass, a_bkg, b_bkg, mean_bkg, sigma_bkg));
  RooRealVar a_bkg_exp("a_bkg_exp","a_bkg_exp", .1,-1,1);
  RooExponential exp_bkg("exp_bkg","exp_bkg", ak08Pruned_1_mass, a_bkg_exp);

  RooArgusBG argus("argus","argus", a_pass_bkg,a1_pass_bkg, a2_bkg);
//////  cheby_pass_bkg/*exp_bkg*/.fitTo(*dh_bkgPass) ;
  
  
  RooDCBShape dcb_unmatch1("dcb_unmatch1", "double crystal ball", ak08Pruned_1_mass, mean_ttbar_unmatch, sigma_ttbar_unmatch, dCBCutL_ttbar_unmatch, dCBPowL_ttbar_unmatch, dCBCutR_ttbar_unmatch,dCBPowR_ttbar_unmatch);
  //dcb_unmatch1.fitTo(*dh_bkgPass) ;

  wdata_ttbar_unmatch.append(wdata_WW);
  wdata_ttbar_unmatch.append(wdata_WZ);
  wdata_ttbar_unmatch.append(wdata_ZZ);
  wdata_ttbar_unmatch.append(wdata_W_100_200);
  wdata_ttbar_unmatch.append(wdata_W_200_400);
  wdata_ttbar_unmatch.append(wdata_W_400_600);
  wdata_ttbar_unmatch.append(wdata_W_600_800);
  wdata_ttbar_unmatch.append(wdata_W_800_1200);
  wdata_ttbar_unmatch.append(wdata_W_1200_2500);
  wdata_ttbar_unmatch.append(wdata_W_2500_Inf);
  wdata_ttbar_unmatch.append(wdata_Stop1_unmatch);
  wdata_ttbar_unmatch.append(wdata_Stop2_unmatch);
  wdata_ttbar_unmatch.append(wdata_Stop3_unmatch);
  //wdata_ttbar_unmatch.append(wdata_Stop4_unmatch);
  //wdata_ttbar_unmatch.append(wdata_Stop5_unmatch);
  RooPlot* frame1_bkg = ak08Pruned_1_mass.frame(Bins(nbins),Title("BackGround Pass sample")) ;
  cheby_pass_bkg.fitTo(wdata_ttbar_unmatch, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2));
  cheby_pass_bkg.fitTo(wdata_ttbar_unmatch, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2));
  cheby_pass_bkg.fitTo(wdata_ttbar_unmatch, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2));
  //dh_bkgPass->plotOn(frame1_bkg) ;
  wdata_ttbar_unmatch.plotOn(frame1_bkg, DataError(RooAbsData::SumW2));

  cheby_pass_bkg.plotOn(frame1_bkg);
  cheby_pass_bkg.paramOn(frame1_bkg);
  //RooPlot* frame2_ttbar = ak08Pruned_1_mass.frame(Bins(nbins),Title("TopTop Pass/UnMatch sample")) ;
  //dh_ttbarUnMatch->plotOn(frame2_ttbar) ;
  TCanvas* c_bkg = new TCanvas("c_bkg","c_bkg",1);//,800,400) ;
  //c_ttbar->Divide(2) ;
  /*c_ttbar->cd(1) ;*/ gPad->SetLeftMargin(0.15) ; frame1_bkg->GetYaxis()->SetTitleOffset(1.4) ; frame1_bkg->Draw() ;
  //c_ttbar->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_ttbar->GetYaxis()->SetTitleOffset(1.4) ; frame2_ttbar->Draw() ;

  std::cout<<"+++++++++++++++++++++++++++++++++++"<<std::endl;
  //std::cout<<"ttbar Integral (match): "<<hh_ttbarMatch->Integral()<<" "<<matched_ttbar<<std::endl;
  //return;


  RooRealVar mean_pass("mean_pass","mean_pass", 80, 60,95);
  RooRealVar sigma_pass("sigma_pass","sigma_pass", 1, 0,20);
  RooGaussian gx_pass("gx_pass","gx_pass",ak08Pruned_1_mass,mean_pass,sigma_pass) ;
  
  RooRealVar mean_pass1("mean_pass1","mean_pass1", 120, 70,130);
  RooRealVar sigma_pass1("sigma_pass1","sigma_pass1", 10, 0,40);
  RooGaussian gx_pass1("gx_pass1","gx_pass1",ak08Pruned_1_mass,mean_pass1,sigma_pass1) ;

  RooRealVar a_pass("a_pass","a_pass", 0.1,-1,1);
  RooRealVar a1_pass("a1_pass","a1_pass", 0.1,-1,1);
  RooRealVar a2_pass("a2_pass","a2_pass", 0.1,-1,1);
  RooPolynomial p2_pass("p2_pass","p2_pass",ak08Pruned_1_mass,RooArgList(a_pass,a1_pass,a2_pass),0) ;

  RooRealVar Nsig("Nsig", "Nsig", 1200,10,10000);//(Nttbar_true.getValV()+NStop_true.getValV())+1000);
  RooRealVar Nbkg("Nbkg", "Nbkg", 100,10,100000);
  RooExtendPdf gx_pass_norm("gx_pass_norm", "gx_pass_norm", gx_pass, Nsig);
  RooExtendPdf p2_pass_norm("p2_pass_norm","p2_pass_norm", p2_pass, Nbkg);
  //RooAddPdf gauss2_pass("gauss2_pass", "gauss2_pass",RooArgList(gx_pass,gx_pass1) );

  RooRealVar a_ttbar1x("a_ttbar1x","a_ttbar", 2., 0.2, 50.);
  RooRealVar a1_ttbar1x("a1_ttbar1x","a1_ttbar", 2., 0.2, 50.);//0.1,-1,1);
  RooRealVar a_ttbarx("a_ttbarx","a_ttbar", 1., 0.1, 50.);
  RooRealVar a1_ttbarx("a1_ttbarx","a1_ttbar", 1., 0.1, 50.);//0.1,-1,1);
  RooRealVar a2_ttbarx("a2_ttbarx","a2_ttbar", -.1,-10,10);
  RooRealVar mean_ttbarx("mean_ttbarx","mean",75,74,81) ;
  RooRealVar sigma_ttbarx("sigma_ttbarx","sigma",4,0,10) ;
  RooDCBShape dcbx("dcbx", "double crystal ball", ak08Pruned_1_mass, mean_ttbarx, sigma_ttbarx, a1_ttbarx, a_ttbarx, a1_ttbar1x, a_ttbar1x);
  

  RooRealVar a_pass1("a_pass1","a_pass1", 1,-100,100);
  RooRealVar a1_pass1("a1_pass1","a1_pass1", 0.1,-100,100);
  RooRealVar a2_pass1("a2_pass1","a2_pass1", -0.1,-100,100);
  /*RooRealVar a_pass1("a_pass1","a_Stop", 1,-1,1);
  RooRealVar a1_pass1("a1_pass1","a1_Stop", 0.1,-1,1);
  RooRealVar a2_pass1("a2_pass1","a2_Stop", -0.1,-1,1);
  */RooChebychev cheby_pass("cheby_pass","cheby_pass",ak08Pruned_1_mass,RooArgSet(a_pass1,a1_pass1, a2_pass1, a_pass)) ;
  RooAddPdf modelPass("modelPass", "modelPass", RooArgList(/*gx_pass*/dcbx,cheby_pass), RooArgList( Nsig, Nbkg));
  //RooProdPdf modelPass("modelPass", "modelPass", RooArgList(gx_pass_norm,p2_pass_norm));
///////  modelPass.fitTo(*dh_totalPass);
  wdata_ttbar_match->append(wdata_ttbar_unmatch);
  std::cout<<"FITTING"<<std::endl;
  modelPass.fitTo(*wdata_ttbar_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2)) ;
  modelPass.fitTo(*wdata_ttbar_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2)) ;
  modelPass.fitTo(*wdata_ttbar_match, RooFit::SumW2Error(kTRUE), RooFit::Strategy(2)) ;
  std::cout<<"After fit"<<std::endl;
  RooPlot* frame1_pass = ak08Pruned_1_mass.frame(Bins(nbins),Title("Pass sample")) ;
  //gx_pass.plotOn(frame1_pass, LineColor(kRed));
  std::cout<<"PLOTTING!!!"<<std::endl;
  wdata_ttbar_match->plotOn(frame1_pass,DataError(RooAbsData::SumW2));
  modelPass.plotOn(frame1_pass, Components(cheby_pass), LineColor(kGreen), Normalization(1.0,RooAbsReal::RelativeExpected));
  modelPass.plotOn(frame1_pass, Components(dcbx), LineColor(kRed), Normalization(1.0,RooAbsReal::RelativeExpected));
  //cheby_pass.plotOn(frame1_pass, LineColor(kGreen));
  modelPass.plotOn(frame1_pass, Normalization(1.0,RooAbsReal::RelativeExpected));
  modelPass.paramOn(frame1_pass);

  TCanvas* c_pass = new TCanvas("c_pass","c_pass",800,400) ;
  gPad->SetLeftMargin(0.15) ; frame1_pass->GetYaxis()->SetTitleOffset(1.4) ; frame1_pass->Draw() ;

  return;
}
