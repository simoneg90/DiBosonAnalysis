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
#include <TF1.h>
#include <TMath.h>
using namespace RooFit ;

#define nbins 22
#define xmin 45
#define xMax 130
#define Lumi 2.2*1000 //[pb-1]

void fit_Signal()
{


  //Take variable from tree
  
  std::string cut_matched_pass = "((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1&&lepton_goodNumber==1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&btag>0&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Ungroomed_1_tau21<0.6&&ak08Ungroomed_WGen_DR<.1";
  std::string cut_unmatched_pass = "((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1&&lepton_goodNumber==1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&btag>0&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Ungroomed_1_tau21<0.6&&ak08Ungroomed_WGen_DR>.3";//&&lepton_WGen_DR>1";
  std::string cut_pass = "((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1&&lepton_goodNumber==1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&btag>0&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Ungroomed_1_tau21<0.6";//"((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1&&lepton_goodNumber==1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&btag>0&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Ungroomed_1_tau21<0.6";
  std::string cut_fail = "((abs(lepton_pdgID)==11&&((abs(lepton_eta)>0&&abs(lepton_eta)<1.442)||(abs(lepton_eta)>1.56&&abs(lepton_eta)<2.5))&&metType1>80&&lepton_pt>120&&ak08Ungroomed_lepton_DR>1&&lepton_goodNumber==1)||(abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&lepton_pt>40))&&btag>0&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Ungroomed_1_tau21>0.6";
  std::string cut = "abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&btag>0&&lepton_pt>40&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Pruned_1_mass>40&&ak08Pruned_1_mass<150&&ak08Ungroomed_1_tau21<0.6&&ak08Ungroomed_WGen_DR<0.1";
  std::string cutCtl = "abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&btag>0&&lepton_pt>40&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Pruned_1_mass>40&&ak08Pruned_1_mass<150&&ak08Ungroomed_1_tau21>0.6";
  std::string variable="ak08Pruned_1_mass";

  std::cout<<"String 1: "<<cut_matched_pass.c_str()<<std::endl;
  std::cout<<"String 2: "<<cut.c_str()<<std::endl;
 
  std::cout<<"+++++++++++++++++++"<<std::endl;
  std::cout<<"      TopTop       "<<std::endl;
  std::cout<<"+++++++++++++++++++"<<std::endl;


  TH1* hh_ttbarMatch= new TH1D ("hh_ttbarMatch", "hh_ttbarMatch", nbins, xmin, xMax); //histo for HP&&matched cuts
  TH1* hh_ttbarUnMatch= new TH1D ("hh_ttbarUnMatch", "hh_ttbarUnMatch", nbins, xmin, xMax); //histo for HP&&unmatched cuts
  TH1* hh_ttbarPass= new TH1D ("hh_ttbarPass", "hh_ttbarPass", nbins, xmin, xMax);
  TH1* hh_ttbarMatch_fail= new TH1D ("hh_ttbarMatch_fail", "hh_ttbarMatch_fail", nbins, xmin, xMax);//histo for LP&&matched cuts
  TH1* hh_ttbarUnMatch_fail= new TH1D ("hh_ttbarUnMatch_fail", "hh_ttbarUnMatch_fail", nbins, xmin, xMax);//histo for LP&&unmatched cuts
  TH1* hh_ttbarFail= new TH1D ("hh_ttbarFail", "hh_ttbarFail", nbins, xmin, xMax);
  TFile *file_ttbar;
  file_ttbar=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/TTBar_Inclusive.root");
  double nEvt_ttbar= ((TH1D *)(file_ttbar->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_ttbar = 831.76;
  TTree* tree_ttbar;
  TDirectory * dir_ttbar = (TDirectory*)file_ttbar->Get("rootTupleTree");
  dir_ttbar->GetObject("tree", tree_ttbar);
  TH1* hh_ttbar1= new TH1D ("hh_ttbar1", "hh_ttbar1", nbins, xmin, xMax);
  TH1* hh_ttbar2= new TH1D ("hh_ttbar2", "hh_ttbar2", nbins, xmin, xMax);
  TH1* hh_ttbar3= new TH1D ("hh_ttbar3", "hh_ttbar3", nbins, xmin, xMax);
  TH1* hh_ttbar4= new TH1D ("hh_ttbar4", "hh_ttbar4", nbins, xmin, xMax);
  TH1* hh_ttbar5= new TH1D ("hh_ttbar5", "hh_ttbar5", nbins, xmin, xMax);
  TH1* hh_ttbar6= new TH1D ("hh_ttbar6", "hh_ttbar6", nbins, xmin, xMax);
  tree_ttbar->Project("hh_ttbar1", Form("%s", variable.c_str()), Form("%s", cut_matched_pass.c_str()));
  tree_ttbar->Project("hh_ttbar2", Form("%s", variable.c_str()), Form("%s", cut_unmatched_pass.c_str()));
  tree_ttbar->Project("hh_ttbar3", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  hh_ttbar1->Scale(xSec_ttbar*Lumi/nEvt_ttbar);
  hh_ttbar2->Scale(xSec_ttbar*Lumi/nEvt_ttbar);
  hh_ttbar3->Scale(xSec_ttbar*Lumi/nEvt_ttbar);

  //tree_ttbar->Project("hh_Stop1Ctl", Form("%s", variable.c_str()), Form("%s", cutCtl.c_str()));
  hh_ttbarMatch->Add(hh_ttbar1);
  hh_ttbarUnMatch->Add(hh_ttbar2);
  hh_ttbarPass->Add(hh_ttbar3);


//  //PUT here normalization parameters!!!
//
  std::cout<<hh_ttbarMatch->GetEntries()<<std::endl;
  std::cout<<hh_ttbarMatch->Integral()<<std::endl;
  RooRealVar ak08Pruned_1_mass_ttbar("ak08Pruned_1_mass_ttbar","ak08Pruned_1_mass_ttbar",xmin,xMax);
  //Create RooDataHist for the fit
  RooDataHist *dh_ttbarMatch= new RooDataHist("dh_ttbarMatch","dh_ttbar",ak08Pruned_1_mass_ttbar,Import(*hh_ttbarMatch)) ;
  RooDataHist *dh_ttbarUnMatch= new RooDataHist("dh_ttbarUnMatch","dh",ak08Pruned_1_mass_ttbar,Import(*hh_ttbarUnMatch)) ;

  RooRealVar Nttbar_true("Nttbar_true","Nttbar_true", 1000,0,10000);
  RooRealVar mean_ttbar("mean_ttbar","mean",80,60,100) ;
  RooRealVar sigma_ttbar("sigma_ttbar","sigma",5,0,10) ;
  RooGaussian gx_ttbar("gx_ttbar","gx",ak08Pruned_1_mass_ttbar,mean_ttbar,sigma_ttbar) ;
  RooExtendPdf gx_ttbar_norm("gx_ttbar_norm", "gx_ttbar_norm", gx_ttbar, Nttbar_true);
  RooRealVar NttbarBKG_true("NttbarBKG_true","NttbarBKG_true", 1000,0,10000);
  RooRealVar mean_ttbar1("mean_ttbar1","mean",80,60,100) ;
  RooRealVar sigma_ttbar1("sigma_ttbar1","sigma",50,0,150) ;
  RooGaussian gx_ttbar1("gx_ttbar1","gx",ak08Pruned_1_mass_ttbar,mean_ttbar1,sigma_ttbar1) ;
  RooExtendPdf gx_ttbar1_norm("gx_ttbar1_norm", "gx_ttbar1_norm", gx_ttbar1, NttbarBKG_true);
  RooRealVar k("k", "k", 100, 0,10000); 

  RooRealVar a_ttbar("a_ttbar","a_ttbar", 1,-100,100);
  RooRealVar a1_ttbar("a1_ttbar","a1_ttbar", 0.1,-1,1);
  RooRealVar a2_ttbar("a2_ttbar","a2_ttbar", -0.1,-1,1);
  RooPolynomial p2_ttbar("p2_ttbar","p2",ak08Pruned_1_mass_ttbar,RooArgList(a_ttbar,a1_ttbar, a2_ttbar),0) ;
  RooChebychev cheby_ttbar("cheby_ttbar","cheby_ttbar",ak08Pruned_1_mass_ttbar,RooArgSet(a_ttbar,a1_ttbar, a2_ttbar)) ;

  RooAddPdf model_ttbar("model_ttbar","model",RooArgList(gx_ttbar, cheby_ttbar/*gx_ttbar1*/), RooArgList(Nttbar_true,k));//, px_ctl)) ;  
  model_ttbar.fitTo(*dh_ttbarMatch) ;
 
  RooPlot* frame1_ttbar = ak08Pruned_1_mass_ttbar.frame(Bins(nbins),Title("TopTop Pass/Match sample")) ;
  dh_ttbarMatch->plotOn(frame1_ttbar) ;
  //gx_ttbar.plotOn(frame1_ttbar, LineColor(kRed)) ;
  //p2_ttbar.plotOn(frame1_ttbar, LineColor(kGreen));
  //gx_ttbar1.plotOn(frame1_ttbar, LineColor(kGreen));
  //cheby_ttbar.plotOn(frame1_ttbar, LineColor(kGreen));
  //model_ttbar.plotOn(frame1_ttbar);//, LineStyle(8));
  //model_ttbar.plotOn(frame1_ttbar, Components(gx_ttbar), LineColor(kRed), Normalization(1.0,RooAbsReal::RelativeExpected));
  //model_ttbar.plotOn(frame1_ttbar, Components(cheby_ttbar), LineColor(kGreen), Normalization(1.0,RooAbsReal::RelativeExpected));
  //model_ttbar.paramOn(frame1_ttbar);
  RooPlot* frame2_ttbar = ak08Pruned_1_mass_ttbar.frame(Bins(nbins),Title("TopTop Pass/UnMatch sample")) ;
  dh_ttbarUnMatch->plotOn(frame2_ttbar) ;
  TCanvas* c_ttbar = new TCanvas("c_ttbar","c_ttbar",800,400) ;
  //c_ttbar->Divide(2) ;
  /*c_ttbar->cd(1) ;*/ gPad->SetLeftMargin(0.15) ; frame1_ttbar->GetYaxis()->SetTitleOffset(1.4) ; frame1_ttbar->Draw() ;
  //c_ttbar->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_ttbar->GetYaxis()->SetTitleOffset(1.4) ; frame2_ttbar->Draw() ;
  TCanvas* c_ttbar_fail = new TCanvas("c_ttbar_fail","c_ttbar",800,400) ;
  gPad->SetLeftMargin(0.15) ; frame2_ttbar->GetYaxis()->SetTitleOffset(1.4) ; frame2_ttbar->Draw() ;



  //return;


  std::cout<<"+++++++++++++++++++"<<std::endl;
  std::cout<<"    Single Top     "<<std::endl;
  std::cout<<"+++++++++++++++++++"<<std::endl;
  TH1* hh_Stop_match= new TH1D ("hh_Stop_match", "hh_Stop", nbins, xmin, xMax); //histo for HP&&match cuts
  TH1* hh_Stop_unmatch= new TH1D ("hh_Stop_unmatch", "hh_Stop_unmatch", nbins, xmin, xMax); //histo for HP&unmatch cuts
  TH1* hh_Stop_pass= new TH1D ("hh_Stop_pass", "hh_Stop_pass", nbins, xmin, xMax);
  TFile *file_Stop1;
  file_Stop1=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_s-channel_4f.root");//testMuonEff/firstTest_reduced_skim.root");
  double nEvt_Stop1= ((TH1D *)(file_Stop1->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_Stop1 = 47.13;
  TTree* tree_Stop1;
  TDirectory * dir_Stop1 = (TDirectory*)file_Stop1->Get("rootTupleTree");
  dir_Stop1->GetObject("tree", tree_Stop1);
  TH1* hh_Stop1= new TH1D ("hh_Stop1", "hh_Stop1", nbins, xmin, xMax);
  TH1* hh_Stop1Ctl= new TH1D ("hh_Stop1Ctl", "hh_Stop1Ctl", nbins, xmin, xMax);
  TH1* hh_Stop1pass= new TH1D ("hh_Stop1pass", "hh_Stop1pass", nbins, xmin, xMax);
  tree_Stop1->Project("hh_Stop1", Form("%s", variable.c_str()), Form("%s", cut_matched_pass.c_str()));//"abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metType1>40&&btag>0&&lepton_pt>40&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Pruned_1_mass>40&&ak08Pruned_1_mass<150&&ak08Ungroomed_1_tau21<0.6");
  tree_Stop1->Project("hh_Stop1Ctl", Form("%s", variable.c_str()), Form("%s", cut_unmatched_pass.c_str()));
  tree_Stop1->Project("hh_Stop1pass", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  hh_Stop1->Scale(xSec_Stop1*Lumi/nEvt_Stop1);
  hh_Stop1Ctl->Scale(xSec_Stop1*Lumi/nEvt_Stop1);
  hh_Stop1Ctl->Scale(xSec_Stop1*Lumi/nEvt_Stop1);
  hh_Stop_match->Add(hh_Stop1);
  hh_Stop_unmatch->Add(hh_Stop1Ctl);
  //hh_Stop1pass->Scale(1000);
  hh_Stop_pass->Add(hh_Stop1pass);
  std::cout<<"Scale Factor: "<<xSec_Stop1*Lumi/nEvt_Stop1<<std::endl;
  TFile *file_Stop2; 
  file_Stop2=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_tW-channel_top.root");//testMuonEff/firstTest_reduced_skim.root");
  double nEvt_Stop2= ((TH1D *)(file_Stop2->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_Stop2 = 35.6;
  TTree* tree_Stop2;
  TDirectory * dir_Stop2 = (TDirectory*)file_Stop2->Get("rootTupleTree");
  dir_Stop2->GetObject("tree", tree_Stop2);
  TH1* hh_Stop2= new TH1D ("hh_Stop2", "hh_Stop2", nbins, xmin, xMax);
  TH1* hh_Stop2Ctl= new TH1D ("hh_Stop2Ctl", "hh_Stop2Ctl", nbins, xmin, xMax);
  TH1* hh_Stop2pass= new TH1D ("hh_Stop2pass", "hh_Stop2pass", nbins, xmin, xMax);
  tree_Stop2->Project("hh_Stop2", Form("%s", variable.c_str()), Form("%s", cut_matched_pass.c_str()));
  tree_Stop2->Project("hh_Stop2Ctl", Form("%s", variable.c_str()), Form("%s", cut_unmatched_pass.c_str()));
  tree_Stop2->Project("hh_Stop2pass", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  hh_Stop2->Scale(xSec_Stop2*Lumi/nEvt_Stop2);
  hh_Stop2Ctl->Scale(xSec_Stop2*Lumi/nEvt_Stop2);
  hh_Stop2pass->Scale(xSec_Stop2*Lumi/nEvt_Stop2);
  hh_Stop_match->Add(hh_Stop2);
  hh_Stop_unmatch->Add(hh_Stop2Ctl);
  hh_Stop_pass->Add(hh_Stop2pass);
  std::cout<<"Scale Factor: "<<xSec_Stop2*Lumi/nEvt_Stop2<<std::endl;
  TFile *file_Stop3; 
  file_Stop3=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_tW-channel_antitop.root");//testMuonEff/firstTest_reduced_skim.root");
  double nEvt_Stop3= ((TH1D *)(file_Stop3->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_Stop3 = 35.6;
  TTree* tree_Stop3;
  TDirectory * dir_Stop3 = (TDirectory*)file_Stop3->Get("rootTupleTree");
  dir_Stop3->GetObject("tree", tree_Stop3);
  TH1* hh_Stop3= new TH1D ("hh_Stop3", "hh_Stop3", nbins, xmin, xMax);
  TH1* hh_Stop3Ctl= new TH1D ("hh_Stop3Ctl", "hh_Stop3Ctl", nbins, xmin, xMax);
  TH1* hh_Stop3pass= new TH1D ("hh_Stop3pass", "hh_Stop3pass", nbins, xmin, xMax);
  tree_Stop3->Project("hh_Stop3", Form("%s", variable.c_str()), Form("%s", cut_matched_pass.c_str()));
  tree_Stop3->Project("hh_Stop3Ctl", Form("%s", variable.c_str()), Form("%s", cut_unmatched_pass.c_str()));
  tree_Stop3->Project("hh_Stop3pass", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  hh_Stop3->Scale(xSec_Stop3*Lumi/nEvt_Stop3);
  hh_Stop3Ctl->Scale(xSec_Stop3*Lumi/nEvt_Stop3);
  hh_Stop3pass->Scale(xSec_Stop3*Lumi/nEvt_Stop3);
  hh_Stop_match->Add(hh_Stop3);
  hh_Stop_unmatch->Add(hh_Stop3Ctl);
  hh_Stop_pass->Add(hh_Stop3pass);
  std::cout<<"Scale Factor: "<<xSec_Stop3*Lumi/nEvt_Stop3<<std::endl;
  TFile *file_Stop4; 
  file_Stop4=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_t-channel_antitop.root");//testMuonEff/firstTest_reduced_skim.root");
  double nEvt_Stop4= ((TH1D *)(file_Stop4->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_Stop4 = 26;
  TTree* tree_Stop4;
  TDirectory * dir_Stop4 = (TDirectory*)file_Stop4->Get("rootTupleTree");
  dir_Stop4->GetObject("tree", tree_Stop4);
  TH1* hh_Stop4= new TH1D ("hh_Stop4", "hh_Stop4", nbins, xmin, xMax);
  TH1* hh_Stop4Ctl= new TH1D ("hh_Stop4Ctl", "hh_Stop4Ctl", nbins, xmin, xMax);
  TH1* hh_Stop4pass= new TH1D ("hh_Stop4pass", "hh_Stop4pass", nbins, xmin, xMax);
  tree_Stop4->Project("hh_Stop4", Form("%s", variable.c_str()), Form("%s", cut_matched_pass.c_str()));
  tree_Stop4->Project("hh_Stop4Ctl", Form("%s", variable.c_str()), Form("%s", cut_unmatched_pass.c_str()));
  tree_Stop4->Project("hh_Stop4pass", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  hh_Stop4->Scale(xSec_Stop4*Lumi/nEvt_Stop4);
  hh_Stop4Ctl->Scale(xSec_Stop4*Lumi/nEvt_Stop4);
  hh_Stop4pass->Scale(xSec_Stop4*Lumi/nEvt_Stop4);
  hh_Stop_match->Add(hh_Stop4);
  hh_Stop_unmatch->Add(hh_Stop4Ctl);
  hh_Stop_pass->Add(hh_Stop4pass);
  std::cout<<"Scale Factor: "<<xSec_Stop4*Lumi/nEvt_Stop4<<std::endl;
  TFile *file_Stop5; 
  file_Stop5=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ST_t-channel_top.root");//testMuonEff/firstTest_reduced_skim.root");
  double nEvt_Stop5= ((TH1D *)(file_Stop5->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_Stop5 = 43.8;
  TTree* tree_Stop5;
  TDirectory * dir_Stop5 = (TDirectory*)file_Stop5->Get("rootTupleTree");
  dir_Stop5->GetObject("tree", tree_Stop5);
  TH1* hh_Stop5= new TH1D ("hh_Stop5", "hh_Stop5", nbins, xmin, xMax);
  TH1* hh_Stop5Ctl= new TH1D ("hh_Stop5Ctl", "hh_Stop5Ctl", nbins, xmin, xMax);
  TH1* hh_Stop5pass= new TH1D ("hh_Stop5pass", "hh_Stop5pass", nbins, xmin, xMax);
  tree_Stop5->Project("hh_Stop5", Form("%s", variable.c_str()), Form("%s", cut_matched_pass.c_str()));
  tree_Stop5->Project("hh_Stop5Ctl", Form("%s", variable.c_str()), Form("%s", cut_unmatched_pass.c_str()));
  tree_Stop5->Project("hh_Stop5pass", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  hh_Stop5->Scale(xSec_Stop5*Lumi/nEvt_Stop5);
  hh_Stop5Ctl->Scale(xSec_Stop5*Lumi/nEvt_Stop5);
  hh_Stop5pass->Scale(xSec_Stop5*Lumi/nEvt_Stop5);
  hh_Stop_match->Add(hh_Stop5);
  hh_Stop_unmatch->Add(hh_Stop5Ctl);
  hh_Stop_pass->Add(hh_Stop5pass);
  std::cout<<"Scale Factor: "<<xSec_Stop5*Lumi/nEvt_Stop5<<std::endl;
  std::cout<<"Integral: "<<hh_Stop_match->Integral()<<std::endl;
  //Create RooRealVar to Fit (the Pruned Mass)
  RooRealVar ak08Pruned_1_mass_Stop("ak08Pruned_1_mass_Stop","ak08Pruned_1_mass",xmin,xMax);
  //Create RooDataHist for the fit
  RooDataHist *dh_StopMatch= new RooDataHist("dh_StopMatch","dh",ak08Pruned_1_mass_Stop,Import(*hh_Stop_match)) ;
  RooDataHist *dh_StopUnMatch= new RooDataHist("dh_StopUnMatch","dh",ak08Pruned_1_mass_Stop,Import(*hh_Stop_unmatch)) ;

  RooRealVar mean_Stop("mean_Stop","mean",80,60,100) ;
  RooRealVar sigma_Stop("sigma_Stop","sigma",2,0,8) ;
  RooGaussian gx_Stop("gx_Stop","gx",ak08Pruned_1_mass_Stop,mean_Stop,sigma_Stop) ;
  RooRealVar NStop_true("NStop_true","NStop_true", 100,0,10000);
  RooExtendPdf gx_Stop_norm("gx_Stop_norm", "gx_Stop_norm", gx_Stop, NStop_true);
  RooRealVar mean_Stop1("mean_Stop1","mean",80,60,100) ;
  RooRealVar sigma_Stop1("sigma_Stop1","sigma",50,0,100) ;
  RooGaussian gx_Stop1("gx_Stop1","gx",ak08Pruned_1_mass_Stop,mean_Stop1,sigma_Stop1) ;
  RooRealVar k_Stop("k_Stop", "k", 100, 0,10000);
  RooRealVar a_Stop("a_Stop","a_Stop", 1,-100,100);
  RooRealVar a1_Stop("a1_Stop","a1_Stop", 0.1,-1,1);
  RooRealVar a2_Stop("a2_Stop","a2_Stop", -0.1,-1,1);
  //RooPolynomial p2_ttbar("p2_ttbar","p2",ak08Pruned_1_mass_ttbar,RooArgList(a_ttbar,a1_ttbar, a2_ttbar),0) ;
  RooChebychev cheby_Stop("cheby_Stop","cheby_Stop",ak08Pruned_1_mass_Stop,RooArgSet(a_Stop,a1_Stop, a2_Stop)) ;

  RooAddPdf model_Stop("model_Stop","model",RooArgList(gx_Stop, cheby_Stop/*gx_Stop1*/), RooArgList(NStop_true,k_Stop));//, px_ctl)) ;
  model_Stop.fitTo(*dh_StopMatch) ;
  //gx_Stop_norm.fitTo(*dh_StopMatch) ;

  //Plotting
  RooPlot* frame1_Stop = ak08Pruned_1_mass_Stop.frame(Bins(nbins),Title("Single Top Pass/Match sample")) ;
  dh_StopMatch->plotOn(frame1_Stop) ;
  //gx_Stop.plotOn(frame1_Stop, LineColor(kRed)) ;
  //cheby_Stop.plotOn(frame1_Stop, LineColor(kGreen)) ;
  //gx_Stop1.plotOn(frame1_Stop, LineColor(kGreen)) ;
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
  RooPlot* frame2_Stop = ak08Pruned_1_mass_Stop.frame(Bins(nbins),Title("Single Top Pass/UnMatch sample")) ;
////  dh_StopCtl->plotOn(frame2_Stop) ;
////////  combData.plotOn(frame2,Cut("sample==sample::control")) ;
  dh_StopUnMatch->plotOn(frame2_Stop) ;
////  simPdf_Stop.plotOn(frame2_Stop,Slice(sampleStop,"fail"),ProjWData(sampleStop,combDataStop)) ;
//////  //px_ctl.plotOn(frame2,Slice(sample,"control"),ProjWData(sample,combData), LineColor(kGreen)) ;
//////  //gx_ctl.plotOn(frame2,Slice(sample,"control"),ProjWData(sample,combData), LineColor(kRed)) ;
//////
  TCanvas* c_Stop = new TCanvas("c_Stop","c_Stop",800,400) ;
  //c_Stop->Divide(2) ;
  /*c_Stop->cd(1) ;*/ gPad->SetLeftMargin(0.15) ; frame1_Stop->GetYaxis()->SetTitleOffset(1.4) ; frame1_Stop->Draw() ;
//  c_Stop->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_Stop->GetYaxis()->SetTitleOffset(1.4) ; frame2_Stop->Draw() ;
  
   TCanvas* c_Stop_fail = new TCanvas("c_Stop_fail","c_Stop",800,400) ;
   gPad->SetLeftMargin(0.15) ; frame2_Stop->GetYaxis()->SetTitleOffset(1.4) ; frame2_Stop->Draw() ;


  std::cout<<"+++++++++++++++++++"<<std::endl;
  std::cout<<"      DiBoson      "<<std::endl;
  std::cout<<"+++++++++++++++++++"<<std::endl;
  TH1* hh_DiBoson_pass= new TH1D ("hh_DiBoson_pass", "hh_Diboson", nbins, xmin, xMax); //histo for HP&&match cuts
  TH1* hh_DiBoson_fail= new TH1D ("hh_DiBoson_fail", "hh_Diboson_unmatch", nbins, xmin, xMax); //histo for HP&&unmatch cuts
  TFile *file_WW;
  file_WW=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/WW.root");//testMuonEff/firstTest_re    duced_skim.root");
  double nEvt_WW= ((TH1D *)(file_WW->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_WW = 118.7;
  TTree* tree_WW;
  TDirectory * dir_WW = (TDirectory*)file_WW->Get("rootTupleTree");
  dir_WW->GetObject("tree", tree_WW);
  TH1* hh_WW= new TH1D ("hh_WW", "hh_WW", nbins, xmin, xMax);
  TH1* hh_WWCtl= new TH1D ("hh_WWCtl", "hh_WWCtl", nbins, xmin, xMax);
  tree_WW->Project("hh_WW", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));//"abs(lepton_pdgID)==13&&muonRelIso03<0.1&&metTy    pe1>40&&btag>0&&lepton_pt>40&&ak08Ungroomed_1_pt>200&&WType1_pt>200&&ak08Pruned_1_mass>40&&ak08Pruned_1_mass<150&&ak08Ungroomed_1_tau21<0.6");
  tree_WW->Project("hh_WWCtl", Form("%s", variable.c_str()), Form("%s", cut_fail.c_str()));
  hh_WW->Scale(xSec_WW*Lumi/nEvt_WW);
  hh_WWCtl->Scale(xSec_WW*Lumi/nEvt_WW);
  hh_DiBoson_pass->Add(hh_WW);
  hh_DiBoson_fail->Add(hh_WWCtl);
  std::cout<<"Scale Factor: "<<xSec_WW*Lumi/nEvt_WW<<std::endl;

  TFile *file_WZ;
  file_WZ=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/WZ.root");//testMuonEff/firstTest_re    duced_skim.r    oot");
  double nEvt_WZ= ((TH1D *)(file_WZ->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_WZ = 16.5;
  TTree* tree_WZ;
  TDirectory * dir_WZ = (TDirectory*)file_WZ->Get("rootTupleTree");
  dir_WZ->GetObject("tree", tree_WZ);
  TH1* hh_WZ= new TH1D ("hh_WZ", "hh_WZ", nbins, xmin, xMax);
  TH1* hh_WZCtl= new TH1D ("hh_WZCtl", "hh_WZCtl", nbins, xmin, xMax);
  tree_WZ->Project("hh_WZ", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  tree_WZ->Project("hh_WZCtl", Form("%s", variable.c_str()), Form("%s", cut_fail.c_str()));
  hh_WZ->Scale(xSec_WZ*Lumi/nEvt_WZ);
  hh_WZCtl->Scale(xSec_WZ*Lumi/nEvt_WZ);
  hh_DiBoson_pass->Add(hh_WZ);
  hh_DiBoson_fail->Add(hh_WZCtl);
  std::cout<<"Scale Factor: "<<xSec_WZ*Lumi/nEvt_WZ<<std::endl;

  TFile *file_ZZ;
  file_ZZ=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/ZZ.root");//testMuonEff/firstTest_re    duced_skim.r        oot");
  double nEvt_ZZ= ((TH1D *)(file_ZZ->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_ZZ = 47.13;
  TTree* tree_ZZ;
  TDirectory * dir_ZZ = (TDirectory*)file_ZZ->Get("rootTupleTree");
  dir_ZZ->GetObject("tree", tree_ZZ);
  TH1* hh_ZZ= new TH1D ("hh_ZZ", "hh_ZZ", nbins, xmin, xMax);
  TH1* hh_ZZCtl= new TH1D ("hh_ZZCtl", "hh_ZZCtl", nbins, xmin, xMax);
  tree_ZZ->Project("hh_ZZ", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  tree_ZZ->Project("hh_ZZCtl", Form("%s", variable.c_str()), Form("%s", cut_fail.c_str()));
  hh_ZZ->Scale(xSec_ZZ*Lumi/nEvt_ZZ);
  hh_ZZCtl->Scale(xSec_ZZ*Lumi/nEvt_ZZ);
  hh_DiBoson_pass->Add(hh_ZZ);
  hh_DiBoson_fail->Add(hh_ZZCtl);
  std::cout<<"Scale Factor: "<<xSec_ZZ*Lumi/nEvt_ZZ<<std::endl;
  TCanvas* c_Diboson = new TCanvas("c_Diboson","c_Diboson",800,400) ;
  hh_DiBoson_pass->Draw();



  std::cout<<"+++++++++++++++++++"<<std::endl;
  std::cout<<"      W+Jets       "<<std::endl;
  std::cout<<"+++++++++++++++++++"<<std::endl;
  TH1* hh_Wjet_pass= new TH1D ("hh_Wjet_pass", "hh_Wjet", nbins, xmin, xMax); //histo for HP&&match cuts
  TH1* hh_Wjet_fail= new TH1D ("hh_Wjet_fail", "hh_Wjet_unmatch", nbins, xmin, xMax); //histo for HP&&unmatch cuts
  TFile *file_W_100_200;
  file_W_100_200=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/WJets_100_200.root");
  double nEvt_W_100_200= ((TH1D *)(file_W_100_200->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_100_200 = 1629.87;
  TTree* tree_W_100_200;
  TDirectory * dir_W_100_200 = (TDirectory*)file_W_100_200->Get("rootTupleTree");
  dir_W_100_200->GetObject("tree", tree_W_100_200);
  TH1* hh_W_100_200= new TH1D ("hh_W_100_200", "hh_W_100_200", nbins, xmin, xMax);
  TH1* hh_W_100_200Ctl= new TH1D ("hh_W_100_200Ctl", "hh_W_100_200Ctl", nbins, xmin, xMax);
  tree_W_100_200->Project("hh_W_100_200", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  tree_W_100_200->Project("hh_W_100_200Ctl", Form("%s", variable.c_str()), Form("%s", cut_fail.c_str()));
  hh_W_100_200->Scale(xSec_W_100_200*Lumi/nEvt_W_100_200);
  hh_W_100_200Ctl->Scale(xSec_W_100_200*Lumi/nEvt_W_100_200);
  hh_Wjet_pass->Add(hh_W_100_200);
  hh_Wjet_fail->Add(hh_W_100_200Ctl);
  std::cout<<"Scale Factor: "<<xSec_W_100_200*Lumi/nEvt_W_100_200<<std::endl;

  TFile *file_W_200_400;
  file_W_200_400=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/WJets_200_400.root");
  double nEvt_W_200_400= ((TH1D *)(file_W_200_400->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_200_400 = 435.6;
  TTree* tree_W_200_400;
  TDirectory * dir_W_200_400 = (TDirectory*)file_W_200_400->Get("rootTupleTree");
  dir_W_200_400->GetObject("tree", tree_W_200_400);
  TH1* hh_W_200_400= new TH1D ("hh_W_200_400", "hh_W_200_400", nbins, xmin, xMax);
  TH1* hh_W_200_400Ctl= new TH1D ("hh_W_200_400Ctl", "hh_W_200_400Ctl", nbins, xmin, xMax);
  tree_W_200_400->Project("hh_W_200_400", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  tree_W_200_400->Project("hh_W_200_400Ctl", Form("%s", variable.c_str()), Form("%s", cut_fail.c_str()));
  hh_W_200_400->Scale(xSec_Stop1*Lumi/nEvt_Stop1);
  hh_W_200_400->Scale(xSec_Stop1*Lumi/nEvt_Stop1);
  hh_Wjet_pass->Add(hh_W_200_400);
  hh_Wjet_fail->Add(hh_W_200_400Ctl);
  std::cout<<"Scale Factor: "<<xSec_W_200_400*Lumi/nEvt_W_200_400<<std::endl;

  TFile *file_W_400_600;
  file_W_400_600=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/WJets_400_600.root");
  double nEvt_W_400_600= ((TH1D *)(file_W_400_600->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_400_600 = 56.17;
  TTree* tree_W_400_600;
  TDirectory * dir_W_400_600 = (TDirectory*)file_W_400_600->Get("rootTupleTree");
  dir_W_400_600->GetObject("tree", tree_W_400_600);
  TH1* hh_W_400_600= new TH1D ("hh_W_400_600", "hh_W_400_600", nbins, xmin, xMax);
  TH1* hh_W_400_600Ctl= new TH1D ("hh_W_400_600Ctl", "hh_W_400_600Ctl", nbins, xmin, xMax);
  tree_W_400_600->Project("hh_W_400_600", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  tree_W_400_600->Project("hh_W_400_600Ctl", Form("%s", variable.c_str()), Form("%s", cut_fail.c_str()));
  hh_W_400_600->Scale(xSec_W_400_600*Lumi/nEvt_W_400_600);
  hh_W_400_600Ctl->Scale(xSec_W_400_600*Lumi/nEvt_W_400_600);
  hh_Wjet_pass->Add(hh_W_400_600);
  hh_Wjet_fail->Add(hh_W_400_600Ctl);
  std::cout<<"Scale Factor: "<<xSec_W_400_600*Lumi/nEvt_W_400_600<<std::endl;


  TFile *file_W_600_800;
  file_W_600_800=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/WJets_600_800.root");
  double nEvt_W_600_800= ((TH1D *)(file_W_600_800->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_600_800 = 14.61;
  TTree* tree_W_600_800;
  TDirectory * dir_W_600_800 = (TDirectory*)file_W_600_800->Get("rootTupleTree");
  dir_W_600_800->GetObject("tree", tree_W_600_800);
  TH1* hh_W_600_800= new TH1D ("hh_W_600_800", "hh_W_600_800", nbins, xmin, xMax);
  TH1* hh_W_600_800Ctl= new TH1D ("hh_W_600_800Ctl", "hh_W_600_800Ctl", nbins, xmin, xMax);
  tree_W_600_800->Project("hh_W_600_800", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  tree_W_600_800->Project("hh_W_600_800Ctl", Form("%s", variable.c_str()), Form("%s", cut_fail.c_str()));
  hh_W_600_800->Scale(xSec_W_600_800*Lumi/nEvt_W_600_800);
  hh_W_600_800Ctl->Scale(xSec_W_600_800*Lumi/nEvt_W_600_800);
  hh_Wjet_pass->Add(hh_W_600_800);
  hh_Wjet_fail->Add(hh_W_600_800Ctl);
  std::cout<<"Scale Factor: "<<xSec_W_600_800*Lumi/nEvt_W_600_800<<std::endl;

  TFile *file_W_800_1200;
  file_W_800_1200=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/WJets_800_1200.root");
  double nEvt_W_800_1200= ((TH1D *)(file_W_800_1200->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_800_1200 = 6.36;
  TTree* tree_W_800_1200;
  TDirectory * dir_W_800_1200 = (TDirectory*)file_W_800_1200->Get("rootTupleTree");
  dir_W_800_1200->GetObject("tree", tree_W_800_1200);
  TH1* hh_W_800_1200= new TH1D ("hh_W_800_1200", "hh_W_800_1200", nbins, xmin, xMax);
  TH1* hh_W_800_1200Ctl= new TH1D ("hh_W_800_1200Ctl", "hh_W_800_1200Ctl", nbins, xmin, xMax);
  tree_W_800_1200->Project("hh_W_800_1200", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  tree_W_800_1200->Project("hh_W_800_1200Ctl", Form("%s", variable.c_str()), Form("%s", cut_fail.c_str()));
  hh_W_800_1200->Scale(xSec_W_800_1200*Lumi/nEvt_W_800_1200); 
  hh_W_800_1200Ctl->Scale(xSec_W_800_1200*Lumi/nEvt_W_800_1200);
  hh_Wjet_pass->Add(hh_W_800_1200);
  hh_Wjet_fail->Add(hh_W_800_1200Ctl);
  std::cout<<"Scale Factor: "<<xSec_W_800_1200*Lumi/nEvt_W_800_1200<<std::endl;


  TFile *file_W_1200_2500;
  file_W_1200_2500=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/WJets_1200_2500.root");
  double nEvt_W_1200_2500= ((TH1D *)(file_W_1200_2500->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_1200_2500 = 1.61;
  TTree* tree_W_1200_2500;
  TDirectory * dir_W_1200_2500 = (TDirectory*)file_W_1200_2500->Get("rootTupleTree");
  dir_W_1200_2500->GetObject("tree", tree_W_1200_2500);
  TH1* hh_W_1200_2500= new TH1D ("hh_W_1200_2500", "hh_W_1200_2500", nbins, xmin, xMax);
  TH1* hh_W_1200_2500Ctl= new TH1D ("hh_W_1200_2500Ctl", "hh_W_1200_2500Ctl", nbins, xmin, xMax);
  tree_W_1200_2500->Project("hh_W_1200_2500", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  tree_W_1200_2500->Project("hh_W_1200_2500Ctl", Form("%s", variable.c_str()), Form("%s", cut_fail.c_str()));
  hh_W_1200_2500->Scale(xSec_W_1200_2500*Lumi/nEvt_W_1200_2500); 
  hh_W_1200_2500Ctl->Scale(xSec_W_1200_2500*Lumi/nEvt_W_1200_2500);
  hh_Wjet_pass->Add(hh_W_1200_2500);
  hh_Wjet_fail->Add(hh_W_1200_2500Ctl);
  std::cout<<"Scale Factor: "<<xSec_W_1200_2500*Lumi/nEvt_W_1200_2500<<std::endl;

  TFile *file_W_2500_Inf;
  file_W_2500_Inf=TFile::Open("/cmshome/gellisim/Diboson/CMSSW_7_4_15/src/DiBosonAnalysis/rootFolder9/WJets_2500_Inf.root");
  double nEvt_W_2500_Inf= ((TH1D *)(file_W_2500_Inf->Get("DijetFilter/EventCount/EventCounter")))->GetBinContent(1);
  double xSec_W_2500_Inf = 37;
  TTree* tree_W_2500_Inf;
  TDirectory * dir_W_2500_Inf = (TDirectory*)file_W_2500_Inf->Get("rootTupleTree");
  dir_W_2500_Inf->GetObject("tree", tree_W_2500_Inf);
  TH1* hh_W_2500_Inf= new TH1D ("hh_W_2500_Inf", "hh_W_2500_Inf", nbins, xmin, xMax);
  TH1* hh_W_2500_InfCtl= new TH1D ("hh_W_2500_InfCtl", "hh_W_2500_InfCtl", nbins, xmin, xMax);
  tree_W_2500_Inf->Project("hh_W_2500_Inf", Form("%s", variable.c_str()), Form("%s", cut_pass.c_str()));
  tree_W_2500_Inf->Project("hh_W_2500_InfCtl", Form("%s", variable.c_str()), Form("%s", cut_fail.c_str()));
  hh_W_2500_Inf->Scale(xSec_W_2500_Inf*Lumi/nEvt_W_2500_Inf); 
  hh_W_2500_InfCtl->Scale(xSec_W_2500_Inf*Lumi/nEvt_W_2500_Inf);
  hh_Wjet_pass->Add(hh_W_2500_Inf);
  hh_Wjet_fail->Add(hh_W_2500_InfCtl);
  std::cout<<"Scale Factor: "<<xSec_W_2500_Inf*Lumi/nEvt_W_2500_Inf<<std::endl;

  TCanvas* c_Wjet = new TCanvas("c_Wjet","c_Wjet",800,400) ;
  hh_Wjet_pass->Draw();


  TH1* hh_bkgPass= new TH1D ("hh_bkgPass", "hh_bkgPass", nbins, xmin, xMax); //histo for HP&&match cuts
  hh_bkgPass->Add(hh_Wjet_pass);
  hh_bkgPass->Add(hh_DiBoson_pass);
  hh_bkgPass->Add(hh_ttbarUnMatch);
  hh_bkgPass->Add(hh_Stop_unmatch);
  //TCanvas* c_bkg = new TCanvas("c_bkg", "c_bkg",800,400);
  //hh_bkgPass->Draw();

  RooRealVar ak08Pruned_1_mass_bkg("ak08Pruned_1_mass_bkg","ak08Pruned_1_mass_bkg",xmin,xMax) ;
  RooDataHist *dh_bkgPass= new RooDataHist("dh_bkgPass","dh_bkgPass",ak08Pruned_1_mass_bkg,Import(*hh_bkgPass)) ;
  RooRealVar a_bkg("a_bkg","a_bkg", 0.1,-1,1);
  RooRealVar a1_bkg("a1_bkg","a1_bkg", 0.1,-1,1);
  RooRealVar a2_bkg("a2_bkg","a2_bkg", 0.1,-1,1);
  //RooRealVar ap_bkg("ap_bkg","ap_bkg", 0.5, -10,10);
  RooRealVar mean_bkg("mean_bkg","mean_bkg", 80, 60,95);
  RooRealVar sigma_bkg("sigma_bkg","sigma_bkg", 5, 0,15);
  RooGaussian gx_bkg("gx_bkg","gx_bkg",ak08Pruned_1_mass_bkg,mean_bkg,sigma_bkg) ;
  RooRealVar mean_bkg1("mean_bkg1","mean_bkg1", 80, 60,95);
  RooRealVar sigma_bkg1("sigma_bkg1","sigma_bkg1", 30, 0,100);
  RooGaussian gx_bkg1("gx_bkg1","gx_bkg1",ak08Pruned_1_mass_bkg,mean_bkg1,sigma_bkg1) ;
  RooRealVar k_bkg("k_bkg","k_bkg", 2, 0,5);
  RooAddPdf model_bkg("model_bkg", "model_bkg", RooArgList(gx_bkg,gx_bkg1), RooArgList(k_bkg));
  RooChebychev p2("p2","p2",ak08Pruned_1_mass_bkg,RooArgSet(a_bkg,a1_bkg,a2_bkg));//,0) ;
  //RooGenericPdf passBkg("passBkg", "passBkg", "(1+(TMath::Erf((ak08Pruned_1_mass_bkg-a_bkg)/b_bkg)))*0.5*exp(-(mean_bkg-ak08Pruned_1_mass_bkg)*(mean_bkg-ak08Pruned_1_mass_bkg)/(2*sigma_bkg*sigma_bkg))", RooArgSet(ak08Pruned_1_mass_bkg, a_bkg, b_bkg, mean_bkg, sigma_bkg));
  p2.fitTo(*dh_bkgPass) ;


  RooPlot* frame1_bkg = ak08Pruned_1_mass_bkg.frame(Bins(nbins),Title("BackGround Pass sample")) ;
  dh_bkgPass->plotOn(frame1_bkg) ;
  /*model_bkg*/p2.plotOn(frame1_bkg);
  p2.paramOn(frame1_bkg);
  //RooPlot* frame2_ttbar = ak08Pruned_1_mass_ttbar.frame(Bins(nbins),Title("TopTop Pass/UnMatch sample")) ;
  //dh_ttbarUnMatch->plotOn(frame2_ttbar) ;
  TCanvas* c_bkg = new TCanvas("c_bkg","c_bkg",800,400) ;
  //c_ttbar->Divide(2) ;
  /*c_ttbar->cd(1) ;*/ gPad->SetLeftMargin(0.15) ; frame1_bkg->GetYaxis()->SetTitleOffset(1.4) ; frame1_bkg->Draw() ;
  //c_ttbar->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_ttbar->GetYaxis()->SetTitleOffset(1.4) ; frame2_ttbar->Draw() ;


  TH1* hh_Pass= new TH1D ("hh_Pass", "hh_Pass", nbins, xmin, xMax);
  hh_Pass->Add(hh_ttbarPass);
  hh_Pass->Add(hh_Stop_pass);
  hh_Pass->Add(hh_DiBoson_pass);
  hh_Pass->Add(hh_Wjet_pass);
  RooRealVar ak08Pruned_1_mass_pass("ak08Pruned_1_mass_pass","ak08Pruned_1_mass_pass",xmin,xMax) ;
  RooDataHist *dh_totalPass = new RooDataHist("dh_totalPass", "dh_totalPass", ak08Pruned_1_mass_pass, Import(*hh_Pass));

  RooRealVar mean_pass("mean_pass","mean_pass", 80, 60,95);
  RooRealVar sigma_pass("sigma_pass","sigma_pass", 10, 0,20);
  RooGaussian gx_pass("gx_pass","gx_pass",ak08Pruned_1_mass_pass,mean_pass,sigma_pass) ;
  
  RooRealVar mean_pass1("mean_pass1","mean_pass1", 120, 70,130);
  RooRealVar sigma_pass1("sigma_pass1","sigma_pass1", 10, 0,40);
  RooGaussian gx_pass1("gx_pass1","gx_pass1",ak08Pruned_1_mass_pass,mean_pass1,sigma_pass1) ;

  RooRealVar a_pass("a_pass","a_pass", 0.1,-1,1);
  RooRealVar a1_pass("a1_pass","a1_pass", 0.1,-1,1);
  RooRealVar a2_pass("a2_pass","a2_pass", 0.1,-1,1);
  RooPolynomial p2_pass("p2_pass","p2_pass",ak08Pruned_1_mass_pass,RooArgList(a_pass,a1_pass,a2_pass),0) ;

  RooRealVar Nsig("Nsig", "Nsig", 1200,10,10000);//(Nttbar_true.getValV()+NStop_true.getValV())+1000);
  RooRealVar Nbkg("Nbkg", "Nbkg", 100,10,100000);
  RooExtendPdf gx_pass_norm("gx_pass_norm", "gx_pass_norm", gx_pass, Nsig);
  RooExtendPdf p2_pass_norm("p2_pass_norm","p2_pass_norm", p2_pass, Nbkg);
  //RooAddPdf gauss2_pass("gauss2_pass", "gauss2_pass",RooArgList(gx_pass,gx_pass1) );

  RooRealVar a_pass1("a_pass1","a_pass1", 1,-100,100);
  RooRealVar a1_pass1("a1_pass1","a1_pass1", 0.1,-1,1);
  RooRealVar a2_pass1("a2_pass1","a2_pass1", -0.1,-1,1);
  /*RooRealVar a_pass1("a_pass1","a_Stop", 1,-1,1);
  RooRealVar a1_pass1("a1_pass1","a1_Stop", 0.1,-1,1);
  RooRealVar a2_pass1("a2_pass1","a2_Stop", -0.1,-1,1);
  */RooChebychev cheby_pass("cheby_pass","cheby_pass",ak08Pruned_1_mass_pass,RooArgSet(a_pass1,a1_pass1, a2_pass1)) ;
  RooAddPdf modelPass("modelPass", "modelPass", RooArgList(gx_pass,cheby_pass), RooArgList( Nsig, Nbkg));
  //RooProdPdf modelPass("modelPass", "modelPass", RooArgList(gx_pass_norm,p2_pass_norm));
  modelPass.fitTo(*dh_totalPass);
  
  RooPlot* frame1_pass = ak08Pruned_1_mass_pass.frame(Bins(nbins),Title("Pass sample")) ;
  dh_totalPass->plotOn(frame1_pass);
  //gx_pass.plotOn(frame1_pass, LineColor(kRed));
  modelPass.plotOn(frame1_pass, Components(cheby_pass), LineColor(kGreen), Normalization(1.0,RooAbsReal::RelativeExpected));
  modelPass.plotOn(frame1_pass, Components(gx_pass), LineColor(kRed), Normalization(1.0,RooAbsReal::RelativeExpected));
  //cheby_pass.plotOn(frame1_pass, LineColor(kGreen));
  modelPass.plotOn(frame1_pass, Normalization(1.0,RooAbsReal::RelativeExpected));
  modelPass.paramOn(frame1_pass);

  TCanvas* c_pass = new TCanvas("c_pass","c_pass",800,400) ;
  gPad->SetLeftMargin(0.15) ; frame1_pass->GetYaxis()->SetTitleOffset(1.4) ; frame1_pass->Draw() ;

  std::cout<<"+++++++++++++++++++++++++++"<<std::endl;
  std::cout<<"  Starting Toy Generation  "<<std::endl;
  std::cout<<"+++++++++++++++++++++++++++"<<std::endl;
  //RooDataSet protoData("protoData", "protoData", ak08Pruned_1_mass_pass);
  RooHistPdf histpdf1("histpdf1","histpdf1",ak08Pruned_1_mass_pass,*dh_totalPass,0) ;
  RooMCStudy* mcstudy = new RooMCStudy(/*histpdf1,*/ modelPass,ak08Pruned_1_mass_pass,Binned(kTRUE),Silence(),Extended(), FitOptions(Save(kTRUE),PrintEvalErrors(0)));
  int Ngen=10000;
  mcstudy->generateAndFit(Ngen, /*Nsig.getValV()+Nbkg.getValV()*/0, kTRUE);

  bool kTrue= 1;
  std::cout<<"plotting first pull"<<std::endl;
  RooPlot* frame1_MCStudy = mcstudy->plotPull(Nsig, Bins(40), FitGauss(kTrue));
  RooPlot* frame3_MCStudy = mcstudy->plotPull(Nbkg, Bins(40), FitGauss(kTrue));
  RooPlot* frame4_MCStudy = mcstudy->plotPull(mean_pass, Bins(40), FitGauss(kTrue));
  RooPlot* frame5_MCStudy = mcstudy->plotPull(sigma_pass, Bins(40), FitGauss(kTrue));
  RooPlot* frame6_MCStudy = mcstudy->plotPull(a_pass1, Bins(40), FitGauss(kTrue));
  RooPlot* frame7_MCStudy = mcstudy->plotPull(a1_pass1, Bins(40), FitGauss(kTrue));
  RooPlot* frame8_MCStudy = mcstudy->plotPull(a2_pass1, Bins(40), FitGauss(kTrue));
  TCanvas* c_MCStudy1 = new TCanvas("c_MCStudy1","c_toy",800,400) ;
  gPad->SetLeftMargin(0.15) ; frame1_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame1_MCStudy->Draw() ;
  TCanvas* c_MCStudy3 = new TCanvas("c_MCStudy3","c_toy",800,400) ;
  gPad->SetLeftMargin(0.15) ; frame3_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame3_MCStudy->Draw() ;
  TCanvas* c_MCStudy4 = new TCanvas("c_MCStudy4","c_toy",800,400) ;
  gPad->SetLeftMargin(0.15) ; frame4_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame4_MCStudy->Draw() ;
  TCanvas* c_MCStudy5 = new TCanvas("c_MCStudy5","c_toy",800,400) ;
  gPad->SetLeftMargin(0.15) ; frame5_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame5_MCStudy->Draw() ;
  TCanvas* c_MCStudy6 = new TCanvas("c_MCStudy6","c_toy",800,400) ;
  gPad->SetLeftMargin(0.15) ; frame6_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame6_MCStudy->Draw() ;
  TCanvas* c_MCStudy7 = new TCanvas("c_MCStudy7","c_toy",800,400) ;
  gPad->SetLeftMargin(0.15) ; frame7_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame7_MCStudy->Draw() ;
  TCanvas* c_MCStudy8 = new TCanvas("c_MCStudy8","c_toy",800,400) ;
  gPad->SetLeftMargin(0.15) ; frame8_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame8_MCStudy->Draw() ;

  std::cout<<"plotting second pull"<<std::endl;
  RooPlot* frame2_MCStudy = ak08Pruned_1_mass_pass.frame(Bins(nbins),Title("Pass sample")) ;//mcstudy->plotPull(Nbkg, Bins(40), FitGauss(kTrue));
  TCanvas* c_toy = new TCanvas("c_toy","c_toy",800,400) ;
  //mcstudy->fitParDataSet().get(1)->writeToFile("../counter.txt");//Print("v");
  //c_toy->Divide(2) ;
  //c_toy->cd(1) ;*/ gPad->SetLeftMargin(0.15) ; frame1_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame1_MCStudy->Draw() ;
  //c_toy->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame2_MCStudy->Draw() ;
  //RooGaussian gx_test("gx_test","gx_pass1",ak08Pruned_1_mass_pass,(RooRealVar &)mean_fitresult,(RooRealVar &)sigma_fitresult) ;
  //gx_test.plotOn(frame2_MCStudy);
  TH1D* pull= new TH1D("pull","pull", 400, -30,30);
  //pull->SetStats(0);
  std::cout<<"Nttbar True: "<<Nttbar_true.getValV()<</*" Nttbar BKG: "<<NttbarBKG_true<<*/" NStop true: "<<NStop_true.getValV()<<std::endl;

  int testCounter=0;


  //(mcstudy->fitResult(0)->createHessePdf(RooArgSet(mcstudy->fitParams(0)))).plotOn(frame2_MCStudy);
  //mcstudy->fitParDataSet().plotOn(frame2_MCStudy, LineColor(kBlue));
  for(int i=0; i<Ngen;++i){
    //mcstudy->genData(i)->plotOn(frame2_MCStudy, LineColor(kBlue), DrawOption("C"));
    std::cout<<mcstudy->genData(i)->numEntries()<<std::endl;
    std::cout<<mcstudy->genData(i)->sumEntries()<<std::endl;
    //mcstudy->fitResult(i)->
    RooRealVar* par1_fitresult = (RooRealVar*) mcstudy->fitResult(i)->floatParsFinal().find("Nsig"); 
    std::cout<<"Variable value..." << par1_fitresult->getValV()<<std::endl;
    std::cout<<"Variable error..." << par1_fitresult->getError()<<std::endl;
    modelPass.fitTo((RooAbsData &)*(mcstudy->genData(i)));
    pull->Fill((Nsig.getValV()-(Nttbar_true.getValV()/*hh_ttbarMatch->GetEntries()*/+NStop_true.getValV()))/Nsig.getError());
    if(((par1_fitresult->getValV()-(Nttbar_true.getValV()+NStop_true.getValV()))/par1_fitresult->getError())>-15. && ((par1_fitresult->getValV()-(Nttbar_true.getValV()+NStop_true.getValV()))/par1_fitresult->getError())<14.8){
      testCounter=i;
      std::cout<<"testCounter: "<<testCounter<<std::endl;
    }
    std::cout<<"Pull: "<<(par1_fitresult->getValV()-(Nttbar_true.getValV()+NStop_true.getValV()))/par1_fitresult->getError()<<std::endl;
  }//end loop for pull
  mcstudy->genData(testCounter)->plotOn(frame2_MCStudy);
  //const RooFitResult* r = (mcstudy->fitResult(0));
  //RooFitResult* r = modelPass.fitTo((RooAbsData &)*(mcstudy->genData(0)), Save());//modelPass.fitTo((mcstudy->genData(0)));//,Save()) ;
  //RooAbsPdf* parabPdf = r->createHessePdf(RooArgSet(Nsig, mean_pass, sigma_pass, Nbkg, a_pass, a1_pass, a2_pass)) ;
  modelPass.fitTo((RooAbsData &)*(mcstudy->genData(testCounter)));
  std::cout<<"After pdf!!!!"<<std::endl;
  std::cout<<"testCounter: "<<testCounter<<std::endl;
  modelPass.plotOn(frame2_MCStudy);
  /*TF1 *myTestF= new TF1("myTestF", "0.023*(([0]*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])))+[3]*([4]+[5]*x+[6]*x*x))", xmin, xMax);
  RooAbsReal* rfa1 = bindFunction(myTestF,ak08Pruned_1_mass_pass, RooArgList(Nsig, mean_pass, sigma_pass, Nbkg, a_pass, a1_pass, a2_pass)) ;
  rfa1->Print();
  myTestF->SetParameter(0, ((RooRealVar*) mcstudy->fitResult(0)->floatParsFinal().find("Nsig"))->getValV());
  myTestF->SetParameter(1, ((RooRealVar*) mcstudy->fitResult(0)->floatParsFinal().find("mean_pass"))->getValV());
  myTestF->SetParameter(2, ((RooRealVar*) mcstudy->fitResult(0)->floatParsFinal().find("sigma_pass"))->getValV());
  myTestF->SetParameter(3, ((RooRealVar*) mcstudy->fitResult(0)->floatParsFinal().find("Nbkg"))->getValV());
  myTestF->SetParameter(4, ((RooRealVar*) mcstudy->fitResult(0)->floatParsFinal().find("a_pass"))->getValV());
  myTestF->SetParameter(5, ((RooRealVar*) mcstudy->fitResult(0)->floatParsFinal().find("a1_pass"))->getValV());
  myTestF->SetParameter(6, ((RooRealVar*) mcstudy->fitResult(0)->floatParsFinal().find("a2_pass"))->getValV());
  rfa1->Print();
*/
  //rfa1->plotOn(frame2_MCStudy);
  /*c_toy->cd(2) ;*/ gPad->SetLeftMargin(0.15) ; frame2_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame2_MCStudy->Draw() ;
  
  TCanvas* c_pull = new TCanvas("c_pull", "c_pull", 800,400) ;
  //RooPlot* frame1_toy = ak08Pruned_1_mass_pass.frame(Bins(nbins),Title("Pass sample")) ;
  TH1* hh_cor_a0_s1f = mcstudy->fitParDataSet().createHistogram("hh",mean_pass,YVar(sigma_pass)) ;
  //myTestF->Draw();
  pull->Draw();
  std::cout<<"Nttbar True: "<<Nttbar_true.getValV()<</*" Nttbar BKG: "<<NttbarBKG_true<<*/" NStop true: "<<NStop_true.getValV()<<std::endl;


  return;
}
