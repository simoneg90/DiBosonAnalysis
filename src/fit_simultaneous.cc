//USAGE: ./fit_simultaneous fit_rootList.txt fit_cutList.txt 2.2 test 1



#include "doPlots.h"
#include "fit_functions.cc"
#include "setTDRStyle.h"
#include "utility.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "time.h"
#include "RooRealVar.h"
#include "RooExtendPdf.h"
#include "RooHist.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
//#include "RooDCBShape.h"
//#include "RooDCBShape.cxx"
#include "RooGlobalFunc.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooSimPdfBuilder.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "RooWorkspace.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooMCStudy.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooRealProxy.h"
#include "TMath.h"
#include "RooMath.h"
//#include "LinkDef.h"
//#include "otherPdf.h"

#include "TROOT.h"
#include "TSystem.h"

#include <fstream>

//#define doTOY 1

using namespace RooFit;

#define MAX_bkg 8 //also the n of colors
#define MAX_sgn 8
#define mass_min 40
#define mass_max 130
#define toyGen 1000


int main(int argc, char* argv[]){

  Color_t COLOR[] = {kAzure-1, kRed, kCyan, kGreen+2, kOrange+10, kSpring-9, kCyan+0, kBlue+0};

  frame("Fit environment - VTagging");
  std::cout<<"Preliminary infos: in his program we use a dedicated dictionary to associate the right BKG"<<std::endl;
  std::cout<<"1: ttbar"<<std::endl;
  std::cout<<"2: SingleTop"<<std::endl;
  std::cout<<"3: WJets"<<std::endl;
  std::cout<<"4: Diboson"<<std::endl;

  frame("OCCHIO!");
  RooRealVar rrv_x("rrv_x","rrv_x", -100.,100.);
  RooRealVar rrv_c_ErfExp      /*= new RooRealVar*/("rrv_c_ErfExp","rrv_c_ErfExp",-0.026,-0.05, 0.05);
  RooRealVar rrv_offset_ErfExp /*= new RooRealVar*/("rrv_offset_ErfExp","rrv_offset_ErfExp",41.,0.,100);
  RooRealVar rrv_width_ErfExp  /*= new RooRealVar*/("rrv_width_ErfExp","rrv_width_ErfExp",30.,1.,100.);
  
  RooErfExpPdf model_pdf ("model_pdf","model_pdf",rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

  //=============================
  //   FIT FUNCTION DECLARATION
  //=============================
  
  //ttbar matched events
  RooRealVar ak08Pruned_1_mass_ttbar("ak08Pruned_1_mass_ttbar","ak08Pruned_1_mass_ttbar",mass_min,mass_max);
  RooRealVar Nttbar_true("Nttbar_true","Nttbar_true", 1000,0,10000);
  RooRealVar mean_ttbar("mean_ttbar","mean",80,60,100) ;
  RooRealVar sigma_ttbar("sigma_ttbar","sigma",5,0,10) ;
  RooGaussian gx_ttbar("gx_ttbar","gx_ttbar",ak08Pruned_1_mass_ttbar,mean_ttbar,sigma_ttbar) ;
  RooExtendPdf gx_ttbar_norm("gx_ttbar_norm", "gx_ttbar_norm", gx_ttbar, Nttbar_true);
  RooRealVar NttbarBKG_true("NttbarBKG_true","NttbarBKG_true", 1000,0,10000);
  RooRealVar mean_ttbar1("mean_ttbar1","mean",80,60,100) ;
  RooRealVar sigma_ttbar1("sigma_ttbar1","sigma",50,0,150) ;
  RooGaussian gx_ttbar1("gx_ttbar1","gx_ttbar1",ak08Pruned_1_mass_ttbar,mean_ttbar1,sigma_ttbar1) ;
  RooExtendPdf gx_ttbar1_norm("gx_ttbar1_norm", "gx_ttbar1_norm", gx_ttbar1, NttbarBKG_true);
  RooRealVar k("k", "k", 0.1, 0,1);
  
  RooAddPdf model_ttbar("model_ttbar","model",RooArgList(gx_ttbar, gx_ttbar1), RooArgList(Nttbar_true,k));//, px_ctl)) ; 

  //single top events
  RooRealVar ak08Pruned_1_mass_Stop("ak08Pruned_1_mass_Stop","ak08Pruned_1_mass",mass_min,mass_max);
  RooRealVar mean_Stop("mean_Stop","mean",80,60,100) ;
  RooRealVar sigma_Stop("sigma_Stop","sigma",2,0,8) ;
  RooGaussian gx_Stop("gx_Stop","gx",ak08Pruned_1_mass_Stop,mean_Stop,sigma_Stop) ;
  RooRealVar NStop_true("NStop_true","NStop_true", 100,0,10000);
  RooExtendPdf gx_Stop_norm("gx_Stop_norm", "gx_Stop_norm", gx_Stop, NStop_true);
  RooRealVar mean_Stop1("mean_Stop1","mean",80,60,100) ;
  RooRealVar sigma_Stop1("sigma_Stop1","sigma",50,0,100) ;
  RooGaussian gx_Stop1("gx_Stop1","gx",ak08Pruned_1_mass_Stop,mean_Stop1,sigma_Stop1) ;
  RooRealVar k_Stop("k_Stop", "k", 0.1, 0,1);
  RooAddPdf model_Stop("model_Stop","model",RooArgList(gx_Stop, gx_Stop1), RooArgList(NStop_true,k_Stop));//, px_ctl)) ;


  //total function for passed events [MC]
  RooRealVar ak08Pruned_1_mass_pass_MC("ak08Pruned_1_mass_pass_MC","ak08Pruned_1_mass_pass_MC",mass_min,mass_max) ;
  
  RooRealVar mean_pass_MC("mean_pass_MC","mean_pass_MC", 80, 60,95);
  RooRealVar sigma_pass_MC("sigma_pass_MC","sigma_pass_MC", 10, 0,20);
  RooGaussian gx_pass_MC("gx_pass_MC","gx_pass_MC",ak08Pruned_1_mass_pass_MC,mean_pass_MC,sigma_pass_MC) ;
  
  RooRealVar a_pass_MC("a_pass_MC","a_pass_MC", 0.1,-1,1);
  RooRealVar a1_pass_MC("a1_pass_MC","a1_pass_MC", 0.1,-1,1);
  RooRealVar a2_pass_MC("a2_pass_MC","a2_pass_MC", 0.1,-1,1);
  RooPolynomial p2_pass_MC("p2_pass_MC","p2_pass_MC",ak08Pruned_1_mass_pass_MC,RooArgList(a_pass_MC,a1_pass_MC,a2_pass_MC),0) ;
  
/*  RooRealVar a_pass_MC("a_pass_MC","a_pass_MC", 1,-100,100);
  RooRealVar a1_pass_MC("a1_pass_MC","a1_pass_MC", 0.1,-1,1);
  RooRealVar a2_pass_MC("a2_pass_MC","a2_pass_MC", -0.1,-1,1);
  RooChebychev cheby_pass_MC("cheby_pass_MC","cheby_pass_MC",ak08Pruned_1_mass_pass_MC,RooArgSet(a_pass_MC,a1_pass_MC, a2_pass_MC)) ;
  */
  RooRealVar Nsig_MC("Nsig_MC", "Nsig_MC", 1200,10,10000);
  RooRealVar Nbkg_MC("Nbkg_MC", "Nbkg_MC", 100,10,100000);
  RooAddPdf modelPass_MC("modelPass_MC", "modelPass_MC", RooArgList(gx_pass_MC,/*cheby_pass_MC*/p2_pass_MC), RooArgList(Nsig_MC, Nbkg_MC));

  //total function for passed events [data]
  RooRealVar ak08Pruned_1_mass_pass_data("ak08Pruned_1_mass_pass_data","ak08Pruned_1_mass_pass_data",mass_min,mass_max) ;
  
  RooRealVar mean_pass_data("mean_pass_data","mean_pass_data", 80, 60,95);
  RooRealVar sigma_pass_data("sigma_pass_data","sigma_pass_data", 10, 0,20);
  RooGaussian gx_pass_data("gx_pass_data","gx_pass_data",ak08Pruned_1_mass_pass_data,mean_pass_data,sigma_pass_data) ;
  
  RooRealVar a_pass_data("a_pass_data","a_pass_data", 0.1,-1,1);
  RooRealVar a1_pass_data("a1_pass_data","a1_pass_data", 0.1,-1,1);
  RooRealVar a2_pass_data("a2_pass_data","a2_pass_data", 0.1,-1,1);
  RooPolynomial p2_pass_data("p2_pass_data","p2_pass_data",ak08Pruned_1_mass_pass_data,RooArgList(a_pass_data,a1_pass_data,a2_pass_data),0) ;
  
  RooRealVar Nsig_data("Nsig_data", "Nsig_data", 1200,10,10000);
  RooRealVar Nbkg_data("Nbkg_data", "Nbkg_data", 100,10,100000);
  RooAddPdf modelPass_data("modelPass_data", "modelPass_data", RooArgList(gx_pass_data,p2_pass_data), RooArgList(Nsig_data, Nbkg_data));


  //total function for failed events [MC]
  RooRealVar ak08Pruned_1_mass_fail_MC("ak08Pruned_1_mass_fail_MC","ak08Pruned_1_mass_fail_MC",mass_min,mass_max) ;
  
  RooRealVar mean_fail_MC("mean_fail_MC","mean_fail_MC", 80, 78,105);
  RooRealVar sigma_fail_MC("sigma_fail_MC","sigma_fail_MC", 10, 0,20);
  RooGaussian gx_fail_MC("gx_fail_MC","gx_fail_MC",ak08Pruned_1_mass_fail_MC,mean_fail_MC,sigma_fail_MC) ;
  
  RooRealVar a_fail_MC("a_fail_MC","a_fail_MC", 0.1,-1,1);
  RooRealVar a1_fail_MC("a1_fail_MC","a1_fail_MC", 0.1,-1,1);
  RooRealVar a2_fail_MC("a2_fail_MC","a2_fail_MC", 0.1,-1,1);
  RooRealVar a3_fail_MC("a3_fail_MC","a3_fail_MC", 0.1,-1,1);
  RooPolynomial p2_fail_MC("p2_fail_MC","p2_fail_MC",ak08Pruned_1_mass_fail_MC,RooArgList(a_fail_MC,a1_fail_MC,a2_fail_MC),0) ;
  RooRealVar k_fail_MC("k_fail_MC","k_fail_MC", 0.7,0,1);
  RooExponential exp_fail_MC("exp_fail_MC", "exp_fail_MC", ak08Pruned_1_mass_fail_MC, a3_fail_MC);
  RooAddPdf testModel_fail_MC("testModel_fail_MC", "", RooArgList(exp_fail_MC,p2_fail_MC), RooArgList(k_fail_MC));
  
  RooChebychev cheby_fail_MC("cheby_fail_MC","cheby_fail_MC",ak08Pruned_1_mass_fail_MC,RooArgSet(a_fail_MC,a1_fail_MC, a2_fail_MC)) ;
  
  RooRealVar Nsig_fail_MC("Nsig_fail_MC", "Nsig_fail_MC", 1200,10,10000);
  RooRealVar Nbkg_fail_MC("Nbkg_fail_MC", "Nbkg_fail_MC", 100,10,100000);
  RooExtendPdf gx_pass_norm_fail_MC("gx_pass_norm_fail_MC", "gx_pass_norm_fail_MC", gx_fail_MC, Nsig_fail_MC);
  RooAddPdf model_fail_MC("model_fail_MC", "model_fail_MC", RooArgList(gx_fail_MC,cheby_fail_MC), RooArgList(Nsig_fail_MC, Nbkg_fail_MC));

  //total function for failed events [data]
  RooRealVar ak08Pruned_1_mass_fail_data("ak08Pruned_1_mass_fail_data","ak08Pruned_1_mass_fail_data",mass_min,mass_max) ;
  
  RooRealVar mean_fail_data("mean_fail_data","mean_fail_data", 80, 78,105);
  RooRealVar sigma_fail_data("sigma_fail_data","sigma_fail_data", 10, 0,20);
  RooGaussian gx_fail_data("gx_fail_data","gx_fail_data",ak08Pruned_1_mass_fail_data,mean_fail_data,sigma_fail_data) ;
  
  RooRealVar a_fail_data("a_fail_data","a_fail_data", 0.1,-1,1);
  RooRealVar a1_fail_data("a1_fail_data","a1_fail_data", 0.1,-1,1);
  RooRealVar a2_fail_data("a2_fail_data","a2_fail_data", 0.1,-1,1);
  RooRealVar a3_fail_data("a3_fail_data","a3_fail_data", 0.1,-1,1);
  RooPolynomial p2_fail_data("p2_fail_data","p2_fail_data",ak08Pruned_1_mass_fail_data,RooArgList(a_fail_data,a1_fail_data,a2_fail_data),0) ;
  RooRealVar k_fail_data("k_fail_data","k_fail_data", 0.7,0,1);
  RooExponential exp_fail_data("exp_fail_data", "exp_fail_data", ak08Pruned_1_mass_fail_data, a3_fail_data);
  RooAddPdf testModel_fail_data("testModel_fail_data", "", RooArgList(exp_fail_data,p2_fail_data), RooArgList(k_fail_data));
  
  RooChebychev cheby_fail_data("cheby_fail_data","cheby_fail_data",ak08Pruned_1_mass_fail_data,RooArgSet(a_fail_data,a1_fail_data, a2_fail_data)) ;
  
  RooRealVar Nsig_fail_data("Nsig_fail_data", "Nsig_fail_data", 1200,10,10000);
  RooRealVar Nbkg_fail_data("Nbkg_fail_data", "Nbkg_fail_data", 100,10,100000);
  RooExtendPdf gx_pass_norm_fail_data("gx_pass_norm_fail_data", "gx_pass_norm_fail_data", gx_fail_data, Nsig_fail_data);
  RooAddPdf model_fail_data("model_fail_data", "model_fail_data", RooArgList(gx_fail_data,cheby_fail_data), RooArgList(Nsig_fail_data, Nbkg_fail_data));

  //simultaneous pdf declaration [MC]
  RooRealVar efficiency ("efficiency", "efficiency", 0,1.);
  //RooRealVar efficiency_ctl ("efficiency_ctl", "efficiency", 0,1.);
  //RooRealVar Nsig_sim_MC_ctl("Nsig_sim_MC_ctl", "Nsig_sim_MC", 1200,10,10000);

  //+++ Pass part +++
  RooRealVar ak08Pruned_1_mass("ak08Pruned_1_mass","ak08Pruned_1_mass",mass_min,mass_max) ;
  RooRealVar ak08Ungroomed_1_tau21("ak08Ungroomed_1_tau21","ak08Ungroomed_1_tau21", 0,1);
  
  /*RooRealVar mean_pass_sim_MC("mean_pass_sim_MC","mean_pass_sim_MC", 80, 60,95);
  RooRealVar sigma_pass_sim_MC("sigma_pass_sim_MC","sigma_pass_sim_MC", 10, 0,20);
  RooGaussian gx_pass_sim_MC("gx_pass_sim_MC","gx_pass_sim_MC",ak08Pruned_1_mass,mean_pass_sim_MC,sigma_pass_sim_MC) ;*/


  //constant part

  RooRealVar const_efficiency ("const_efficiency", "efficiency",  .91, 0 ,1.);
  RooRealVar const_mean_pass_sim_MC("const_mean_pass_sim_MC", "Double CB mean", 81., 78., 87.);//, 75.,85.);//78.7);//78., 75.,85.);
  RooRealVar const_sigma_pass_sim_MC("const_sigma_pass_sim_MC", "Double CB Width", 7.,5.,8.5);//,5.,10.);//8.4);//8.,5.,10.);
  RooConstVar const_dCBCutL_pass_sim_MC("const_dCBCutL_pass_sim_MC", "Double CB Cut left", 1.77);//3.7);//1.85);//***** tau21_04 1.73);//***** tau21_06 1.5);
  RooConstVar const_dCBCutR_pass_sim_MC("const_dCBCutR_pass_sim_MC", "Double CB Cut right", 2.);//4.1);//***** tau21_04  2.);//***** tau21_06 1.9);//1.9);
  RooConstVar const_dCBPowerL_pass_sim_MC("const_dCBPowerL_pass_sim_MC", "Double CB Power left", 1.16);//2.3);//***** tau21_04 1.14);//***** tau21_06 1.16);
  RooConstVar const_dCBPowerR_pass_sim_MC("const_dCBPowerR_pass_sim_MC", "Double CB Power right", 6.09);//3.9);//5.4);//***** tau21_04 2.26);//***** tau21_06 2.2);// 4.7);
  RooDCBShape const_dcb_pass_sim_MC("const_dcb_pass_sim_MC", "double crystal ball", ak08Pruned_1_mass, const_mean_pass_sim_MC, const_sigma_pass_sim_MC, const_dCBCutL_pass_sim_MC, const_dCBCutR_pass_sim_MC, const_dCBPowerL_pass_sim_MC, const_dCBPowerR_pass_sim_MC);
  RooRealVar const_a_pass_sim_MC("const_a_pass_sim_MC","a_pass_sim_MC",0.47,0.3,0.5);//.01, .001,.1);//.06,.001,.2);//, 0.001,.1);//.1,.01,1.);//***** tau21_04  .48);//***** tau21_06 .1);//.12);//1,-100,100);
  RooRealVar const_a1_pass_sim_MC("const_a1_pass_sim_MC","a1_pass_sim_MC",-.54,-.6,-.4);//-.4,-1.,-.1);//,-1.,-.1);//-.5,-1.,-.01);//***** tau21_04 -.61);//***** tau21_06 -.5);//-.32);//0.1,-1,1);
  RooConstVar const_a2_pass_sim_MC("const_a2_pass_sim_MC","a2_pass_sim_MC",-.1);//,-.1,.1);//***** tau21_04 -.19);//***** tau21_06  .07);//.07);//0.1,-1,1);
  //RooPolynomial p2_pass_sim_MC("p2_pass_sim_MC","p2_pass_sim_MC",ak08Pruned_1_mass,RooArgList(a_pass_sim_MC,a1_pass_sim_MC,a2_pass_sim_MC),0) ;
  RooChebychev const_cheby_pass_sim_MC("const_cheby_pass_sim_MC","cheby_pass_sim_MC",ak08Pruned_1_mass,RooArgSet(const_a_pass_sim_MC,const_a1_pass_sim_MC, const_a2_pass_sim_MC)) ;
  RooGaussian gx_fail_sim_MC("gx_fail_sim_MC","gx_fail_sim_MC",ak08Pruned_1_mass,const_mean_pass_sim_MC, const_sigma_pass_sim_MC);
  RooRealVar const_Nsig_sim_MC("const_Nsig_sim_MC", "Nsig_sim_MC", 1000.,0.,100000000.);
  RooFormulaVar const_k_pass_sim("const_k_pass_sim", "pass norm", "(const_Nsig_sim_MC*const_efficiency)", RooArgList(const_Nsig_sim_MC, const_efficiency));
  RooRealVar const_Nbkg_pass_sim_MC("const_Nbkg_pass_sim_MC", "Nbkg_pass_sim_MC",1000.,0.,100000000.);
  RooConstVar a_const("a_const","a_const", 10.);
  RooPolynomial p2_pass_sim_MC("p2_pass_sim_MC","p2_pass_sim_MC",ak08Pruned_1_mass,RooArgList(a_const),0) ;
  RooAddPdf const_modelPass_sim_MC("const_modelPass_sim_MC", "modelPass_sim_MC", RooArgList(const_dcb_pass_sim_MC,/*p2_pass_sim_MC*/const_cheby_pass_sim_MC), RooArgList(const_k_pass_sim, const_Nbkg_pass_sim_MC));


  //+++ Fail part +++
  RooConstVar const_sim_sigmaRatio("const_sim_sigmaRatio", "ratio data sigma", 1.51);
  RooConstVar const_sim_meanRatio("const_sim_meanRatio", "ratio mean data", .96);//.96);
  RooFormulaVar const_mean_fail_sim_MC("const_mean_fail_sim_MC", "Double CB mean", "(const_sim_meanRatio*const_mean_pass_sim_MC)", RooArgList(const_sim_meanRatio,const_mean_pass_sim_MC));
  //RooRealVar const_mean_fail_sim_MC("const_mean_fail_sim_MC", "Double CB mean", 76,75,85);
  //RooFormulaVar const_sigma_fail_sim_MC("const_sigma_fail_sim_MC", "Double CB Width","(const_sim_sigmaRatio*const_sigma_pass_sim_MC", RooArgList(const_sim_sigmaRatio,const_sigma_pass_sim_MC));

  //RooConstVar const_mean_fail_sim_MC("const_mean_fail_sim_MC", "Double CB mean", 75.);//***** tau21_04 76.53);//***** tau21_06 75.5);
  RooRealVar const_sigma_fail_sim_MC("const_sigma_fail_sim_MC", "Double CB Width", 9.6,8.,11.);//5.7, 5., 15.);//10.,5,15); //working points 10, 0, 40 for tau21<.4
  RooConstVar const_dCBCutL_fail_sim_MC("const_dCBCutL_fail_sim_MC", "Double CB Cut left",.8);//.48);//***** tau21_04 .7);//***** tau21_06 .8);//.42);// 0.18); //for tau21<0.4 1., 0.1, 50. works
  RooConstVar const_dCBCutR_fail_sim_MC("const_dCBCutR_fail_sim_MC", "Double CB Cut right",.82);//1.);//***** tau21_04 1.4);//***** tau21_06 1.5);//.2);//1.5);
  RooConstVar const_dCBPowerL_fail_sim_MC("const_dCBPowerL_fail_sim_MC", "Double CB Power left",1.);// 1.);//***** tau21_04 1.5);//***** tau21_06 2.5 .72);//.82); //.72);
  RooConstVar const_dCBPowerR_fail_sim_MC("const_dCBPowerR_fail_sim_MC", "Double CB Power right",2.24);//***** tau21_04  2.);//***** tau21_06 2.3);//2.8);//2.3); //working points 2., -0.2, 50. for tau21<0.6 using 2 DCB for LP and HP when same mu and sigma while 2., 0.2, 50. for tau21<0.4
  RooDCBShape const_dcb_fail_sim_MC("const_dcb_fail_sim_MC", "double crystal ball", ak08Pruned_1_mass, const_mean_fail_sim_MC, const_sigma_fail_sim_MC, const_dCBCutL_fail_sim_MC, const_dCBCutR_fail_sim_MC, const_dCBPowerL_fail_sim_MC, const_dCBPowerR_fail_sim_MC);
  RooRealVar const_a_fail_sim_MC("const_a_fail_sim_MC","a_fail_sim_MC", -.7,-1.,0.);//-.69, -1.,0.);//-1.06,-2.,-.5);//, -3.,-.1);//-1.1,-2.,0.);//***** tau21_04 -.66);//-.7);//***** tau21_06 -1.1);
  RooRealVar const_a1_fail_sim_MC("const_a1_fail_sim_MC","a1_fail_sim_MC", .02,0.,.1);//.31,0.,.6);//,.1,2.);//.34,-1.,.5);//***** tau21_04 -.01);//.04);//***** tau21_06 .34);
  RooConstVar const_a2_fail_sim_MC("const_a2_fail_sim_MC","a2_fail_sim_MC", .067);//-.05,-1.,.0);//***** tau21_04 .07);//-.013);//***** tau21_06 -.05);
  RooChebychev const_cheby_fail_sim_MC("const_cheby_fail_sim_MC","cheby_fail_sim_MC",ak08Pruned_1_mass,RooArgSet(const_a_fail_sim_MC,const_a1_fail_sim_MC, const_a2_fail_sim_MC)) ;
  
  RooRealVar const_Nbkg_fail_sim_MC("const_Nbkg_fail_sim_MC", "Nbkg_fail_sim_MC", 1000.,0.,100000000.);
  RooFormulaVar const_k_fail_sim("const_k_fail_sim", "fail norm", "(const_Nsig_sim_MC*(1-const_efficiency))", RooArgList(const_Nsig_sim_MC, const_efficiency));
  RooAddPdf const_model_fail_sim_MC("const_model_fail_sim_MC", "model_fail_sim_MC", RooArgList(const_dcb_fail_sim_MC,const_cheby_fail_sim_MC), RooArgList(const_k_fail_sim, const_Nbkg_fail_sim_MC));

  //end of constant part
  //+++++++++++++++++++++++++++++++++++



  //RooRealVar ak08Pruned_1_mass_data("ak08Pruned_1_mass_data","ak08Pruned_1_mass_data",mass_min,mass_max) ;
  //constant part data
  RooRealVar const_data_efficiency ("const_data_efficiency", "efficiency",  .91, 0 ,1.);
  RooRealVar const_data_mean_pass_sim_MC("const_data_mean_pass_sim_MC", "Double CB mean", 81., 78., 87.);//78., 75.,85.);
  RooRealVar const_data_sigma_pass_sim_MC("const_data_sigma_pass_sim_MC", "Double CB Width", 7.,5.,8.5);
  RooConstVar const_data_dCBCutL_pass_sim_MC("const_data_dCBCutL_pass_sim_MC", "Double CB Cut left", 1.77); 
  RooConstVar const_data_dCBCutR_pass_sim_MC("const_data_dCBCutR_pass_sim_MC", "Double CB Cut right", 2.);
  RooConstVar const_data_dCBPowerL_pass_sim_MC("const_data_dCBPowerL_pass_sim_MC", "Double CB Power left", 1.16);
  RooConstVar const_data_dCBPowerR_pass_sim_MC("const_data_dCBPowerR_pass_sim_MC", "Double CB Power right", 6.09);
  RooDCBShape const_data_dcb_pass_sim_MC("const_data_dcb_pass_sim_MC", "double crystal ball", ak08Pruned_1_mass, const_data_mean_pass_sim_MC, const_data_sigma_pass_sim_MC, const_data_dCBCutL_pass_sim_MC, const_data_dCBCutR_pass_sim_MC, const_data_dCBPowerL_pass_sim_MC, const_data_dCBPowerR_pass_sim_MC);
  RooRealVar const_data_a_pass_sim_MC("const_data_a_pass_sim_MC","a_pass_sim_MC",.47, .3,.5);
  RooRealVar const_data_a1_pass_sim_MC("const_data_a1_pass_sim_MC","a1_pass_sim_MC",-.54,-.6,-.4);
  RooConstVar const_data_a2_pass_sim_MC("const_data_a2_pass_sim_MC","a2_pass_sim_MC",-.1);//.07,-0.1,.1);//***** tau21_04 -.19);//***** tau21_06  .07);//.07);//0.1,-1,1);
  RooChebychev const_data_cheby_pass_sim_MC("const_data_cheby_pass_pass_MC","cheby_pass_sim_MC",ak08Pruned_1_mass,RooArgSet(const_data_a_pass_sim_MC,const_data_a1_pass_sim_MC, const_data_a2_pass_sim_MC)) ;
  RooRealVar const_data_Nsig_sim_MC("const_data_Nsig_sim_MC", "Nsig_sim_MC", 300.,100,1000000000000);
  RooFormulaVar const_data_k_pass_sim("const_data_k_pass_sim", "pass norm", "(const_data_Nsig_sim_MC*const_data_efficiency)", RooArgList(const_data_Nsig_sim_MC, const_data_efficiency));
  RooRealVar const_data_Nbkg_pass_sim_MC("const_data_Nbkg_pass_sim_MC", "Nbkg_pass_sim_MC", 575., 400.,1000000000000.);
  RooAddPdf const_data_modelPass_sim_MC("const_data_modelPass_sim_MC", "modelPass_sim_MC", RooArgList(const_data_dcb_pass_sim_MC,const_data_cheby_pass_sim_MC), RooArgList(const_data_k_pass_sim, const_data_Nbkg_pass_sim_MC));


  //+++ Fail part +++
  //RooConstVar const_data_sigmaRatio("const_data_sigmaRatio", "ratio data sigma", 1.25);
  RooConstVar const_data_meanRatio("const_data_meanRatio", "ratio mean data", .96);
  RooFormulaVar const_data_mean_fail_sim_MC("const_data_mean_fail_sim_MC", "Double CB mean", "(const_data_meanRatio*const_data_mean_pass_sim_MC)", RooArgList(const_data_meanRatio,const_data_mean_pass_sim_MC));
  //RooFormulaVar const_data_sigma_fail_sim_MC("const_data_sigma_fail_sim_MC", "Double CB Width","(const_data_sigmaRatio*const_data_sigma_pass_sim_MC", RooArgList(const_data_sigmaRatio,const_data_sigma_pass_sim_MC));
  //RooConstVar const_data_mean_fail_sim_MC("const_data_mean_fail_sim_MC", "Double CB mean", 75.5);//***** tau21_04 76.53);//***** tau21_06 75.5);
  RooRealVar const_data_sigma_fail_sim_MC("const_data_sigma_fail_sim_MC", "Double CB Width", 9.6,8.,11.);
  RooConstVar const_data_dCBCutL_fail_sim_MC("const_data_dCBCutL_fail_sim_MC", "Double CB Cut left",.8);
  RooConstVar const_data_dCBCutR_fail_sim_MC("const_data_dCBCutR_fail_sim_MC", "Double CB Cut right", .82);
  RooConstVar const_data_dCBPowerL_fail_sim_MC("const_data_dCBPowerL_fail_sim_MC", "Double CB Power left", 1.);//***** tau21_04 1.5);//***** tau21_06 .72);//.82); //.72);
  RooConstVar const_data_dCBPowerR_fail_sim_MC("const_data_dCBPowerR_fail_sim_MC", "Double CB Power right",2.24);
  RooDCBShape const_data_dcb_fail_sim_MC("const_data_dcb_fail_sim_MC", "double crystal ball", ak08Pruned_1_mass, const_data_mean_fail_sim_MC, const_data_sigma_fail_sim_MC, const_data_dCBCutL_fail_sim_MC, const_data_dCBCutR_fail_sim_MC, const_data_dCBPowerL_fail_sim_MC, const_data_dCBPowerR_fail_sim_MC);
  RooRealVar const_data_a_fail_sim_MC("const_data_a_fail_sim_MC","a_fail_sim_MC", -.7,-2.,0.);
  RooRealVar const_data_a1_fail_sim_MC("const_data_a1_fail_sim_MC","a1_fail_sim_MC", .02,0.,.1);
  RooConstVar const_data_a2_fail_sim_MC("const_data_a2_fail_sim_MC","a2_fail_sim_MC", .067);//-.05,-1.,0.);//***** tau21_04 .07);//-.013);//***** tau21_06 -.05);
  RooChebychev const_data_cheby_fail_sim_MC("const_data_cheby_fail_sim_MC","cheby_fail_sim_MC",ak08Pruned_1_mass,RooArgSet(const_data_a_fail_sim_MC,const_data_a1_fail_sim_MC, const_data_a2_fail_sim_MC)) ;
  
  RooRealVar const_data_Nbkg_fail_sim_MC("const_data_Nbkg_fail_sim_MC", "Nbkg_fail_sim_MC", 292., 100., 1000000000.);
  RooFormulaVar const_data_k_fail_sim("const_data_k_fail_sim", "fail norm", "(const_data_Nsig_sim_MC*(1-const_data_efficiency))", RooArgList(const_data_Nsig_sim_MC, const_data_efficiency));
  RooAddPdf const_data_model_fail_sim_MC("const_data_model_fail_sim_MC", "model_fail_sim_MC", RooArgList(const_data_dcb_fail_sim_MC,const_data_cheby_fail_sim_MC), RooArgList(const_data_k_fail_sim, const_data_Nbkg_fail_sim_MC));

  //end of constant part
  RooRealVar mean_pass_sim_MC("mean_pass_sim_MC", "Double CB mean", 80, 60, 95);
  RooRealVar sigma_pass_sim_MC("sigma_pass_sim_MC", "Double CB Width", 10, 0, 20);
  //Used before putting constraints!
  //RooRealVar dCBCutL_pass_sim_MC("dCBCutL_pass_sim_MC", "Double CB Cut left", 1., 0.1, 50.);
  //RooRealVar dCBCutR_pass_sim_MC("dCBCutR_pass_sim_MC", "Double CB Cut right", 1., 0.1, 50.);
  //RooRealVar dCBPowerL_pass_sim_MC("dCBPowerL_pass_sim_MC", "Double CB Power left", 2., 0.2, 50.);
  //RooRealVar dCBPowerR_pass_sim_MC("dCBPowerR_pass_sim_MC", "Double CB Power right", 2., 0.2, 50.);
  //new! with constraints
  RooRealVar dCBCutL_pass_sim_MC("dCBCutL_pass_sim_MC", "Double CB Cut left",  1.4, 1.3, 1.5);
  RooRealVar dCBCutR_pass_sim_MC("dCBCutR_pass_sim_MC", "Double CB Cut right",  2.5, 2.4, 2.6);
  RooRealVar dCBPowerL_pass_sim_MC("dCBPowerL_pass_sim_MC", "Double CB Power left", .8, 0.7, .9);
  RooRealVar dCBPowerR_pass_sim_MC("dCBPowerR_pass_sim_MC", "Double CB Power right",  9., 8., 10.);
  RooDCBShape dcb_pass_sim_MC("dcb_pass_sim_MC", "double crystal ball", ak08Pruned_1_mass, mean_pass_sim_MC, sigma_pass_sim_MC, dCBCutL_pass_sim_MC, dCBCutR_pass_sim_MC, dCBPowerL_pass_sim_MC, dCBPowerR_pass_sim_MC);
  RooRealVar a_pass_sim_MC("a_pass_sim_MC","a_pass_sim_MC", .1,-10,10);//1,-100,100);
  RooRealVar a1_pass_sim_MC("a1_pass_sim_MC","a1_pass_sim_MC", .1,-1,1);//0.1,-1,1);
  RooRealVar a2_pass_sim_MC("a2_pass_sim_MC","a2_pass_sim_MC", .1,-1,1);//0.1,-1,1);
  //RooPolynomial p2_pass_sim_MC("p2_pass_sim_MC","p2_pass_sim_MC",ak08Pruned_1_mass,RooArgList(a_pass_sim_MC,a1_pass_sim_MC,a2_pass_sim_MC),0) ;
  RooChebychev cheby_pass_sim_MC("cheby_pass_pass_MC","cheby_pass_sim_MC",ak08Pruned_1_mass,RooArgSet(a_pass_sim_MC,a1_pass_sim_MC, a2_pass_sim_MC)) ;
  
  RooRealVar Nsig_sim_MC("Nsig_sim_MC", "Nsig_sim_MC", 1200,10,10000);
  RooFormulaVar k_pass_sim("k_pass_sim", "pass norm", "(Nsig_sim_MC*efficiency)", RooArgList(Nsig_sim_MC, efficiency));
  RooRealVar Nbkg_pass_sim_MC("Nbkg_pass_sim_MC", "Nbkg_pass_sim_MC", 100,10,100000);
  RooAddPdf modelPass_sim_MC("modelPass_sim_MC", "modelPass_sim_MC", RooArgList(dcb_pass_sim_MC/*gx_pass_sim_MC*/,cheby_pass_sim_MC), RooArgList(k_pass_sim, Nbkg_pass_sim_MC));


  //+++ Fail part +++
  /*RooRealVar mean_fail_sim_MC("mean_fail_sim_MC","mean_fail_sim_MC", 80, 78,105);
  RooRealVar sigma_fail_sim_MC("sigma_fail_sim_MC","sigma_fail_sim_MC", 10, 0,20);
  RooGaussian gx_fail_sim_MC("gx_fail_sim_MC","gx_fail_sim_MC",ak08Pruned_1_mass,mean_fail_sim_MC,sigma_fail_sim_MC) ;*/

  RooRealVar mean_fail_sim_MC("mean_fail_sim_MC", "Double CB mean", 80, 70, 95);
  RooRealVar sigma_fail_sim_MC("sigma_fail_sim_MC", "Double CB Width", 10, 0, 20); //working points 10, 0, 40 for tau21<.4
  RooRealVar dCBCutL_fail_sim_MC("dCBCutL_fail_sim_MC", "Double CB Cut left",  1.4, 1.1, 1.7); //for tau21<0.4 1., 0.1, 50. works
  RooRealVar dCBCutR_fail_sim_MC("dCBCutR_fail_sim_MC", "Double CB Cut right",  2.5, 2.2, 2.8);
  RooRealVar dCBPowerL_fail_sim_MC("dCBPowerL_fail_sim_MC", "Double CB Power left",  .8, 0.5, 1.1);
  RooRealVar dCBPowerR_fail_sim_MC("dCBPowerR_fail_sim_MC", "Double CB Power right", 9., 6., 12.); //working points 2., -0.2, 50. for tau21<0.6 using 2 DCB for LP and HP when same mu and sigma while 2., 0.2, 50. for tau21<0.4
  RooDCBShape dcb_fail_sim_MC("dcb_fail_sim_MC", "double crystal ball", ak08Pruned_1_mass, mean_pass_sim_MC, sigma_pass_sim_MC, dCBCutL_fail_sim_MC, dCBCutR_fail_sim_MC, dCBPowerL_fail_sim_MC, dCBPowerR_fail_sim_MC);
  RooRealVar a_fail_sim_MC("a_fail_sim_MC","a_fail_sim_MC", 0.1,-10,10);
  RooRealVar a1_fail_sim_MC("a1_fail_sim_MC","a1_fail_sim_MC", 0.1,-1,1);
  RooRealVar a2_fail_sim_MC("a2_fail_sim_MC","a2_fail_sim_MC", 0.1,-1,1);
  RooChebychev cheby_fail_sim_MC("cheby_fail_sim_MC","cheby_fail_sim_MC",ak08Pruned_1_mass,RooArgSet(a_fail_sim_MC,a1_fail_sim_MC, a2_fail_sim_MC)) ;
  
  RooRealVar Nbkg_fail_sim_MC("Nbkg_fail_sim_MC", "Nbkg_fail_sim_MC", 100,10,100000);
  RooFormulaVar k_fail_sim("k_fail_sim", "fail norm", "(Nsig_sim_MC*(1-efficiency))", RooArgList(Nsig_sim_MC, efficiency));
  RooAddPdf model_fail_sim_MC("model_fail_sim_MC", "model_fail_sim_MC", RooArgList(dcb_fail_sim_MC,cheby_fail_sim_MC), RooArgList(k_fail_sim, Nbkg_fail_sim_MC));

  RooCategory sample_MC("sample_MC","sample_MC") ;
  sample_MC.defineType("passed") ;
  sample_MC.defineType("failed") ;

  RooCategory sample_data("sample_data","sample_data") ;
  sample_data.defineType("passed") ;
  sample_data.defineType("failed") ;

  
  //=== end of function declaration ===

  TH1D *histoD_pass[MAX_NUMBER];//histos for project High Purity
  TH1D *histoD_fail[MAX_NUMBER];//histos for project Low Purity
  TH1D *bkgCounts[MAX_NUMBER];//histos to count events to rescale
  TH1D *histoBkg_pass[MAX_bkg];    //1 histo for each background HP
  TH1D *histoBkg_fail[MAX_bkg];    //1 histo for each background LP
  TH1D *dataHisto_pass; //data histo HP
  TH1D *dataHisto_fail; //data histo LP
  TH1D *allBkgHisto_pass;
  TH1D *allBkgHisto_fail;
  TH1D *allBkgHisto_pass_toy;
  TH1D *allBkgHisto_fail_toy;
  TH1D *signalHisto[MAX_sgn]; //signal histos (if we want to superimpose them)

  ////////////DataSet Variables declaration
  TFile *file_reduced[MAX_NUMBER];
  TTree *tree_reduced[MAX_NUMBER];
  RooDataSet *completeDataset[MAX_NUMBER]; //dataset from ttree
  RooDataSet *passDataset[MAX_NUMBER]; //only passed events
  RooDataSet *failDataset[MAX_NUMBER]; //only failed events
  RooDataSet *passWeightDataset[MAX_NUMBER]; //only passed events with weight
  RooDataSet *failWeightDataset[MAX_NUMBER]; //only failed events with weight
  RooFormulaVar *weightFunc[MAX_NUMBER]; //needed to take into account weights
  RooRealVar ds_weight ("ds_weight", "ds_weight", 1.); //real weight used by weightFunc
  RooRealVar *final_Wpass[MAX_NUMBER];    //the column added to the dataset
  RooRealVar *final_Wfail[MAX_NUMBER];
  RooDataSet *dataPass;
  RooDataSet *dataFail;
  RooDataSet *MC_Pass;
  RooDataSet *MC_Fail;
  /////////////////////////////////////////

  //TH2D *histo2D[MAX_NUMBER];
  //TH2F *histo2F[MAX_NUMBER];
  TTree *tree[MAX_NUMBER];
  TCanvas *c[MAX_NUMBER];
  //TLegend *leg[MAX_NUMBER];
  TFile *file[MAX_NUMBER];
  TLegend *leg = new TLegend(0.65,0.65,0.9,0.9);
  TStopwatch timer;
  timer.Start(true);

  if(argc<=2){
    breakLine();
    std::cout<<"Attention! Not enough inputs!"<<std::endl;
    std::cout<<"Digit './fit_simultaneous --help' for informations!"<<std::endl;
    breakLine();
    exit(-1);
  }
   
  if(strcmp(argv[1],"--help")==0){
    breakLine();
    std::cout<<"This program superimposes or Stacks the same histogram from different files"<<std::endl;
    std::cout<<"Program usage:"<<std::endl;
    std::cout<<"./fit_simultaneous rootList cutList lumi[fb-1] output weight"<<std::endl;
    //argv[1]=inputList 
    breakLine();
    exit(-1);
  }
  
  if(argc>6){
    breakLine();
    std::cout<<"ATTENTION! TOO MANY ARGUMENTS ADDED!"<<std::endl;
    std::cout<<"Exiting program"<<std::endl;
    breakLine();
    exit(-1);
  }

  std::ifstream inputList; //list of rootFiles used for the analysis
  inputList.open(argv[1]);
  if(inputList==NULL){
    std::cout<<"ERROR! File: "<<argv[1]<<" doesn't exist!"<<std::endl;
    exit(-1);
  }
  std::cout<<"Opened root List: " <<argv[1]<<std::endl;

  std::ifstream cutList; //list of cuts to be applied: 1. match HP; 2. unmatch HP; 3. match LP; 4. unmatch LP; 5. Pass; 6. Fail
  cutList.open(argv[2]);
  if(cutList==NULL){
    std::cout<<"ERROR! File: "<<argv[2]<<" doesn't exist!"<<std::endl;
    exit(-1);
  }
  std::cout<<"Opened cut list: "<<argv[2]<<std::endl;

  //=== variable definition ===
  std::string prefix=argv[4];
  double scaleFactor=1.;
  double lumi=(atof(argv[3]))*1000;
  double rescale=(atoi(argv[5]));
  if(rescale==0){
    frame("No rescale applied");
  }else if(rescale==1){
    frame("Samples will be rescaled for Lumi and xsec");
  }

  //=== assigning the cut to the right string ===
  std::string cut_string, match_HP, unmatch_HP, match_LP, unmatch_LP, pass, fail;
  cut_string=match_HP=unmatch_HP=match_LP=unmatch_LP=pass=fail="";
  int cut_counter=0;
  while(cutList>>cut_string){

    if(cut_counter==0) match_HP=cut_string;
    if(cut_counter==1) unmatch_HP=cut_string;
    if(cut_counter==2) match_LP=cut_string;
    if(cut_counter==3) unmatch_LP=cut_string;
    if(cut_counter==4) pass=cut_string;
    if(cut_counter==5) fail=cut_string;
    std::cout<<cut_string.c_str()<<std::endl;
    ++cut_counter;
  }//end loop over cut file 'cutList'
  
  if((strcmp(match_HP.c_str(), "")==0)||(strcmp(unmatch_HP.c_str(), "")==0)||(strcmp(match_LP.c_str(), "")==0)||(strcmp(unmatch_LP.c_str(), "")==0)||(strcmp(pass.c_str(), "")==0)||(strcmp(fail.c_str(), "")==0)){
    breakLine();
    std::cout<<"ATTENTION! one or more cut strings are not filled!"<<std::endl;
    std::cout<<"Please control that your cutList file " <<argv[2]<<" has 6 rows"<<std::endl;
    breakLine();
    exit(-1);
  }

  frame("List of cuts");
  
  std::cout<<"Match High Purity: "<<match_HP.c_str()<<std::endl;
  std::cout<<"UnMatch High Purity: "<<unmatch_HP.c_str()<<std::endl;
  std::cout<<"Match Low Purity: "<<match_LP.c_str()<<std::endl;
  std::cout<<"UnMatch Low Purity: "<<unmatch_LP.c_str()<<std::endl;
  std::cout<<"Pass: "<<pass.c_str()<<std::endl;
  std::cout<<"Fail: "<<fail.c_str()<<std::endl;

  std::string bkg_nameTMP, bkg_file, bkg_name;
  //double scaleFactor=1.;
  int bins, min, max;
  int file_counter=0;
  int bkg_counter=-1;
  bool isData=0;
  double counts=0;
  std::string variable="ak08Pruned_1_mass";
  bkg_name="";
  std::string redu_suff= "_redu.root";
  std::string pass_unbin="ak08Ungroomed_1_tau21<.45";
  std::string fail_unbin="ak08Ungroomed_1_tau21>.45";
  while(inputList>>bkg_nameTMP>>bkg_file>>scaleFactor>>bins>>min>>max){
   std::cout<<"BKG: "<<bkg_nameTMP.c_str()<<" File: "<<bkg_file.c_str()<<" xSection: "<<scaleFactor<<std::endl;

   file[file_counter]=TFile::Open(bkg_file.c_str());
   if(file[file_counter]==0){
     breakLine();
     std::cout<<"ATTENTION! "<<bkg_file.c_str()<<" doesn't exists!"<<std::endl;
     breakLine();
     exit(-1);
   }
   bkgCounts[file_counter]=((TH1D *)(file[file_counter]->Get("DijetFilter/EventCount/EventCounter")));
   std::cout<<"Total Events: "<<bkgCounts[file_counter]->GetBinContent(1)<<std::endl;
   counts=bkgCounts[file_counter]->GetBinContent(1);
  // TDirectory * dir = (TDirectory*)file[file_counter]->Get("rootTupleTree");
  // dir->GetObject("tree", tree[file_counter]);
   std::cout<<"Creating Histo"<<std::endl;
   histoD_pass[file_counter] = new TH1D(Form("histo_%d_pass", file_counter), Form("histo_%d_pass", file_counter),bins,min,max);
   histoD_fail[file_counter] = new TH1D(Form("histo_%d_fail", file_counter), Form("histo_%d_fail", file_counter),bins,min,max);
  // std::cout<<"Projecting"<<std::endl;
  // tree[file_counter]->Project(Form("histo_%d_pass", file_counter),Form("%s", variable.c_str()),pass.c_str());
  // tree[file_counter]->Project(Form("histo_%d_fail", file_counter),Form("%s", variable.c_str()),fail.c_str());
    
   frame("Unbinned part");
   file_reduced[file_counter]=TFile::Open(Form("%s%s",(bkg_file.substr(0,bkg_file.size()-5)).c_str(),redu_suff.c_str()));
   if(file_reduced[file_counter]==0){
     breakLine();
     std::cout<<"ATTENTION! "<<(bkg_file.substr(0,bkg_file.size()-5)).c_str()<<redu_suff.c_str()<<" doesn't exists!"<<std::endl;
     breakLine();
     exit(-1);
   }
   tree_reduced[file_counter]=(TTree *)file_reduced[file_counter]->Get("mio");
   completeDataset[file_counter] = new RooDataSet(Form("ds_%d",file_counter),Form("ds_%d",file_counter),RooArgSet(ak08Pruned_1_mass, ak08Ungroomed_1_tau21),Import(*tree_reduced[file_counter]));

   //reducing datasets to fail and pass components
   passDataset[file_counter] = (RooDataSet*) completeDataset[file_counter]->reduce(pass_unbin.c_str()) ;
   failDataset[file_counter] = (RooDataSet*) completeDataset[file_counter]->reduce(fail_unbin.c_str()) ;
   passDataset[file_counter]->SetName(Form("passed_%d",file_counter));
   passDataset[file_counter]->Print();
   failDataset[file_counter]->SetName(Form("failed_%d",file_counter));
   failDataset[file_counter]->Print();

   ds_weight.setVal(scaleFactor*lumi/counts);
   if(strcmp(bkg_nameTMP.c_str(), "data")==0) ds_weight.setVal(1.);
   std::cout<<"UNBINNED WEIGHT: "<<ds_weight.getValV()<<std::endl;
   weightFunc[file_counter]= new RooFormulaVar(Form("weightFunc_%d",file_counter),"event weight","ds_weight",RooArgList(ds_weight,ak08Pruned_1_mass)) ; 
    final_Wpass[file_counter]= (RooRealVar*) passDataset[file_counter]->addColumn(((RooAbsArg &)*weightFunc[file_counter]));
    final_Wfail[file_counter]= (RooRealVar*) failDataset[file_counter]->addColumn(((RooAbsArg &)*weightFunc[file_counter]));
    passWeightDataset[file_counter]= new RooDataSet(passDataset[file_counter]->GetName(),passDataset[file_counter]->GetTitle(),passDataset[file_counter],*passDataset[file_counter]->get(),0,final_Wpass[file_counter]->GetName()) ;
    failWeightDataset[file_counter]= new RooDataSet(failDataset[file_counter]->GetName(),failDataset[file_counter]->GetTitle(),failDataset[file_counter],*failDataset[file_counter]->get(),0,final_Wfail[file_counter]->GetName()) ;    



   if(strcmp(bkg_nameTMP.c_str(), "data")==0){
     std::cout<<"DATA!"<<std::endl;
     if(isData==0) {
       allBkgHisto_pass= new TH1D("allBkgHisto_pass","allBkgHisto_pass", bins,min,max);
       allBkgHisto_fail= new TH1D("allBkgHisto_fail","allBkgHisto_fail", bins,min,max);
       allBkgHisto_pass->Sumw2();
       allBkgHisto_fail->Sumw2();
       dataHisto_pass= new TH1D("dataHisto_pass","dataHisto_pass", bins,min,max);
       dataHisto_fail= new TH1D("dataHisto_fail","dataHisto_fail", bins,min,max);
       dataPass = new RooDataSet (passDataset[file_counter]->GetName(),passDataset[file_counter]->GetTitle(),passDataset[file_counter],*passDataset[file_counter]->get(),0,final_Wpass[file_counter]->GetName());
       dataFail = new RooDataSet (failDataset[file_counter]->GetName(),failDataset[file_counter]->GetTitle(),failDataset[file_counter],*failDataset[file_counter]->get(),0,final_Wfail[file_counter]->GetName()) ;
       std::cout<<"Data counts: "<<counts<<std::endl;
       //dataHisto->GetXaxis()->SetTitle(xTitle.c_str());
       //dataHisto_pass->GetYaxis()->SetTitle(yTitle.c_str());
       dataHisto_pass->SetMarkerStyle(8);
       dataHisto_pass->SetMarkerColor(1);
       dataHisto_pass->SetLineColor(1);
       //dataHisto_fail->GetYaxis()->SetTitle(yTitle.c_str());
       dataHisto_fail->SetMarkerStyle(8);
       dataHisto_fail->SetMarkerColor(1);
       dataHisto_fail->SetLineColor(1);
       //leg->AddEntry(dataHisto, "Data", "pe");
      }else{
        dataPass->append(*passWeightDataset[file_counter]);
        dataFail->append(*failWeightDataset[file_counter]);
      }
      isData=1;
      dataHisto_pass->Add(histoD_pass[file_counter]);
      dataHisto_fail->Add(histoD_fail[file_counter]);

    }else{
      if(bkg_counter<0){
        MC_Pass = new RooDataSet(passDataset[file_counter]->GetName(),passDataset[file_counter]->GetTitle(),passDataset[file_counter],*passDataset[file_counter]->get(),0,final_Wpass[file_counter]->GetName());
        MC_Fail= new RooDataSet(failDataset[file_counter]->GetName(),failDataset[file_counter]->GetTitle(),failDataset[file_counter],*failDataset[file_counter]->get(),0,final_Wfail[file_counter]->GetName()) ;
      }else{
        MC_Pass->append(*passWeightDataset[file_counter]);
        MC_Fail->append(*failWeightDataset[file_counter]);
      }
      if(rescale==1){
        histoD_pass[file_counter]->Scale(scaleFactor*lumi/counts);
        histoD_fail[file_counter]->Scale(scaleFactor*lumi/counts);
        std::cout<<" " <<scaleFactor<<" "<<lumi<<" "<<counts<<std::endl;
        std::cout<<"Total scale factor: "<<scaleFactor*lumi/counts<<std::endl;
      }//end if rescale

      if(strcmp(bkg_name.c_str(),bkg_nameTMP.c_str())==0){//if the bkg is always the same
        std::cout<<"Background: "<<bkg_name.c_str()<<std::endl;
        histoBkg_pass[bkg_counter]->Add(histoD_pass[file_counter]);
        histoBkg_fail[bkg_counter]->Add(histoD_fail[file_counter]);
      }//end for the same bkg
      else{

        bkg_name=bkg_nameTMP;
        //if(bkg_counter!=-1) bkgStack->Add(histoBkg[bkg_counter]);
        ++bkg_counter;
        histoBkg_pass[bkg_counter]= new TH1D(Form("histo_%s_pass", bkg_name.c_str()), Form("histo_%s_pass", bkg_name.c_str()),bins,min,max);
        histoBkg_fail[bkg_counter]= new TH1D(Form("histo_%s_fail", bkg_name.c_str()), Form("histo_%s_fail", bkg_name.c_str()),bins,min,max);
        histoBkg_pass[bkg_counter]->Sumw2();
        histoBkg_fail[bkg_counter]->Sumw2();
        histoBkg_pass[bkg_counter]->SetStats(0);
        histoBkg_fail[bkg_counter]->SetStats(0);
        //histoBkg[bkg_counter]->SetLineColor(COLOR[bkg_counter]);
        histoBkg_pass[bkg_counter]->SetFillColor(COLOR[bkg_counter]);
        histoBkg_pass[bkg_counter]->SetLineWidth(1);
        histoBkg_pass[bkg_counter]->SetLineColor(1);
        histoBkg_pass[bkg_counter]->Add(histoD_pass[file_counter]);
        histoBkg_fail[bkg_counter]->SetFillColor(COLOR[bkg_counter]);
        histoBkg_fail[bkg_counter]->SetLineWidth(1);
        histoBkg_fail[bkg_counter]->SetLineColor(1);
        histoBkg_fail[bkg_counter]->Add(histoD_fail[file_counter]);

        //leg->AddEntry(histoBkg[bkg_counter], bkg_nameTMP.c_str(), "f");
        std::cout<<"New Background: "<<bkg_name.c_str()<<std::endl;



      }//end if bkg already known

      ++file_counter;

    }//end if data or NotData


  }//end loop over rootList file

  std::cout<<"Background counter: "<<bkg_counter+1<<std::endl;
  double integralBKG_all_pass=0;
  double integralBKG_all_fail=0;

  for(int i=0; i<=bkg_counter; ++i){
    //c[i]=new TCanvas(Form("c%d",i), "Grafico1", 200, 10, 600, 400);
    //histoBkg[i]->Draw("hist");
    //c[i]->SaveAs(Form("c%d.png",i));
    std::cout<<"Integral [pass]: "<<histoBkg_pass[i]->Integral()<<std::endl;
    std::cout<<"Integral [fail]: "<<histoBkg_fail[i]->Integral()<<std::endl;
    integralBKG_all_pass+=histoBkg_pass[i]->Integral();
    integralBKG_all_fail+=histoBkg_fail[i]->Integral();
    std::cout<<"All [pass]: "<<integralBKG_all_pass<<" "<<histoBkg_pass[i]->Integral()<<std::endl;
    std::cout<<"All [fail]: "<<integralBKG_all_fail<<" "<<histoBkg_fail[i]->Integral()<<std::endl;
    std::cout<<"Entries [pass]: "<<histoBkg_pass[i]->GetEntries()<<std::endl;
    std::cout<<"Entries [fail]: "<<histoBkg_fail[i]->GetEntries()<<std::endl;
    if(isData==1){
      allBkgHisto_pass->Add(histoBkg_pass[i]);
      allBkgHisto_fail->Add(histoBkg_fail[i]);
    }
    //integralBKG=allBkgHisto->Integral();
    std::cout<<"adding back"<<std::endl;
    //histoBkg[i]->Write();
    //bkgStack->Add(histoBkg[i]);
  }

  TCanvas *c_passed_MC = new TCanvas("c_passed_MC", "Grafico1", 200, 10, 600, 400);
  allBkgHisto_pass_toy= new TH1D("allBkgHisto_pass_toy","generated with allBkgHisto_pass", bins,min,max);
  allBkgHisto_pass_toy->FillRandom(allBkgHisto_pass,(int)allBkgHisto_pass->Integral());
  RooDataHist *dh_totalPass = new RooDataHist("dh_totalPass", "dh_totalPass", ak08Pruned_1_mass_pass_MC, Import(*allBkgHisto_pass));
  RooPlot* frame1_pass = ak08Pruned_1_mass_pass_MC.frame(Bins(22),Title("Pass sample")) ;
  modelPass_MC.fitTo(*dh_totalPass, SumW2Error(kTRUE));
  dh_totalPass->plotOn(frame1_pass);
  modelPass_MC.plotOn(frame1_pass, Components(gx_pass_MC), LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));
  modelPass_MC.plotOn(frame1_pass, Components(p2_pass_MC), LineColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));
  modelPass_MC.plotOn(frame1_pass);
  modelPass_MC.paramOn(frame1_pass);
  /*histoBkg_pass[0]->Draw("");
  allBkgHisto_pass->Draw("SAME");
  */
  gPad->SetLeftMargin(0.15) ; frame1_pass->GetYaxis()->SetTitleOffset(1.4) ; frame1_pass->Draw() ;
  
  c_passed_MC->SaveAs(Form("%s_provaRoofit.png", prefix.c_str()));
  
  TCanvas *c_fail_MC = new TCanvas("c_fail_MC", "Grafico1", 200, 10, 600, 400);
  allBkgHisto_fail_toy= new TH1D("allBkgHisto_fail_toy","generated with allBkgHisto_fail", bins,min,max);
  allBkgHisto_fail_toy->FillRandom(allBkgHisto_fail,(int)allBkgHisto_fail->Integral());
  RooDataHist *dh_totalFail = new RooDataHist("dh_totalFail", "dh_totalFail", ak08Pruned_1_mass_fail_MC, Import(*allBkgHisto_fail_toy));
  RooPlot* frame1_fail = ak08Pruned_1_mass_fail_MC.frame(Bins(22),Title("Fail sample")) ;
  model_fail_MC.fitTo(*dh_totalFail, SumW2Error(kTRUE));
  dh_totalFail->plotOn(frame1_fail);
  model_fail_MC.plotOn(frame1_fail, Components(gx_fail_MC), LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));
  model_fail_MC.plotOn(frame1_fail, Components(cheby_fail_MC), LineColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));
  model_fail_MC.plotOn(frame1_fail);
  //model_fail_MC.paramOn(frame1_fail);
  gPad->SetLeftMargin(0.15) ; frame1_fail->GetYaxis()->SetTitleOffset(1.4) ; frame1_fail->Draw() ;
  
  c_fail_MC->SaveAs(Form("%s_provaRoofit_1.png",prefix.c_str()));


  //Simultaneous fit
  frame("Simultaneous fit");

  RooDataHist *dh_totalPass_sim = new RooDataHist("dh_totalPass_sim", "dh_totalPass_sim", ak08Pruned_1_mass, Import(*allBkgHisto_pass));
  RooDataHist *dh_totalFail_sim = new RooDataHist("dh_totalFail_sim", "dh_totalFail_sim", ak08Pruned_1_mass, Import(*allBkgHisto_fail));
  RooDataHist combData_MC("combData_MC","combined data for MC",ak08Pruned_1_mass,Index(sample_MC),Import("passed",*dh_totalPass_sim), Import("failed", *dh_totalFail_sim));
  RooDataSet combData_MC_unbinned("combData_MC_unbinned","combData_MC_unbinned", ak08Pruned_1_mass,Index(sample_MC),Import("passed",*MC_Pass), Import("failed",*MC_Fail)); 
  RooSimultaneous simPdf_MC("simPdf_MC","simultaneous pdf for MC",sample_MC) ;
  simPdf_MC.addPdf(const_modelPass_sim_MC,"passed") ;
  simPdf_MC.addPdf(const_model_fail_sim_MC, "failed");
  //simPdf_MC.addParameters(efficiency, mean_pass_sim_MC, sigma_pass_sim_MC,a_pass_sim_MC,a1_pass_sim_MC,a2_pass_sim_MC,Nsig_sim_MC,Nbkg_pass_sim_MC,mean_fail_sim_MC,sigma_fail_sim_MC,a_fail_sim_MC,a1_fail_sim_MC,a2_fail_sim_MC,Nbkg_fail_sim_MC);
  //gx_fail_sim_MC.plotOn(frame2_sim, LineColor(kRed));

  simPdf_MC.fitTo(combData_MC_unbinned, SumW2Error(kTRUE));//, RooFit::Strategy(2));
  TCanvas *c_sim_MC_pass = new TCanvas("c_sim_MC_pass", "Simultaneous PDF Pass",1);// 200, 10, 600, 400);
  //c_sim_MC->Divide(2) ;
  RooPlot* frame1_sim = ak08Pruned_1_mass.frame(Bins(bins),Title("Pass sample")) ;
  RooPlot* frame2_sim = ak08Pruned_1_mass.frame(Bins(bins),Title("Fail sample")) ;
  //TF1* myFunc=(TF1 *)const_modelPass_sim_MC->asTF(ak08Pruned_1_mass,efficiency, RooArgList(const_mean_pass_sim_MC, const_sigma_pass_sim_MC,const_a_pass_sim_MC,const_a1_pass_sim_MC,const_a2_pass_sim_MC,const_Nsig_sim_MC,const_Nbkg_pass_sim_MC,const_mean_fail_sim_MC,const_sigma_fail_sim_MC,const_a_fail_sim_MC,const_a1_fail_sim_MC,const_a2_fail_sim_MC,const_Nbkg_fail_sim_MC));
  //dh_totalPass_sim->plotOn(frame1_sim, DataError(RooAbsData::SumW2), LineColor(2));//, LineStyle(1),LineWidth(2) /* DrawOption("l") , XErrorSize(0)*/); const_modelPass_sim_MC.plotOn(frame1_sim);
  MC_Pass->plotOn(frame1_sim, DataError(RooAbsData::SumW2));
  const_modelPass_sim_MC.plotOn(frame1_sim);
  const_modelPass_sim_MC.plotOn(frame1_sim,Components(const_dcb_pass_sim_MC),LineColor(kRed));//, Normalization(1.0,RooAbsReal::RelativeExpected));
  const_modelPass_sim_MC.plotOn(frame1_sim,Components(const_cheby_pass_sim_MC),LineColor(kGreen));//, Normalization(1.0,RooAbsReal::RelativeExpected));
 
  //simPdf_MC.plotOn(frame1_sim,Slice(sample_MC,"passed"),ProjWData(sample_MC,combData_MC)) ;
  //simPdf_MC.paramOn(frame1_sim,Slice(sample_MC,"passed"),ProjWData(sample_MC,combData_MC)) ;
  
  /*gx_fail_sim_MC.plotOn(frame2_sim, LineColor(kRed));
  cheby_fail_sim_MC.plotOn(frame2_sim, LineColor(kGreen));*/
  //dh_totalFail_sim->plotOn(frame2_sim, LineColor(2), DataError(RooAbsData::SumW2));
  MC_Fail->plotOn(frame2_sim, DataError(RooAbsData::SumW2));
  const_model_fail_sim_MC.plotOn(frame2_sim);
  const_model_fail_sim_MC.plotOn(frame2_sim,Components(const_dcb_fail_sim_MC),LineColor(kRed));//, Normalization(1.0,RooAbsReal::RelativeExpected) );
  const_model_fail_sim_MC.plotOn(frame2_sim,Components(const_cheby_fail_sim_MC),LineColor(kGreen));//, Normalization(1.0,RooAbsReal::RelativeExpected));
  /*c_sim_MC->cd(1) ;*/ gPad->SetLeftMargin(0.15) ; frame1_sim->GetYaxis()->SetTitleOffset(1.4) ; frame1_sim->Draw() ;
  //c_sim_MC->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_sim->GetYaxis()->SetTitleOffset(1.4) ; frame2_sim->Draw() ;
  c_sim_MC_pass->SaveAs(Form("%s_provaRoofitSimPass.png",prefix.c_str()));
  c_sim_MC_pass->SaveAs(Form("%s_provaRoofitSimPass.pdf",prefix.c_str()));
  
  TCanvas *c_sim_MC_fail = new TCanvas("c_sim_MC_fail", "Simultaneous PDF Fail",1);// 200, 10, 600, 400);
  //TH1 * mio = dh_totalFail_sim->createHistogram("mio",ak08Pruned_1_mass);
  gPad->SetLeftMargin(0.15) ; frame2_sim->GetYaxis()->SetTitleOffset(1.4) ; frame2_sim->Draw() ;
  c_sim_MC_fail->SaveAs(Form("%s_provaRoofitSimFail.png",prefix.c_str()));
  c_sim_MC_fail->SaveAs(Form("%s_provaRoofitSimFail.pdf",prefix.c_str()));
  RooRealVar t("t", "t", 2, 0,10);
  std::cout<<t.getValV()<<std::endl;


  //Data
  

  RooDataHist *dh_totalPass_data = new RooDataHist("dh_totalPass_data", "dh_totalPass_data", ak08Pruned_1_mass, Import(*dataHisto_pass));
  RooDataHist *dh_totalFail_data = new RooDataHist("dh_totalFail_data", "dh_totalFail_data", ak08Pruned_1_mass, Import(*dataHisto_fail));
  RooDataHist combData_data("combData_data","combined data for data",ak08Pruned_1_mass,Index(sample_data),Import("passed",*dh_totalPass_data), Import("failed", *dh_totalFail_data));
 
  RooDataSet combData_data_unbin("combData_data_unbin","combData_data_unbin",ak08Pruned_1_mass,Index(sample_data),Import("passed",*dataPass), Import("failed", *dataFail));
  RooSimultaneous simPdf_data("simPdf_data","simultaneous pdf for data",sample_data) ;
  simPdf_data.addPdf(const_data_modelPass_sim_MC,"passed") ;
  simPdf_data.addPdf(const_data_model_fail_sim_MC, "failed");
/////////  simPdf_data.fitTo(combData_data_unbin);// DATA are not weighted!, SumW2Error(kTRUE));
  simPdf_data.fitTo(combData_data_unbin);
  TCanvas *c_sim_data_pass = new TCanvas("c_sim_data_pass", "Simultaneous PDF Pass",1);// 200, 10, 600, 400);
  //c_sim_MC->Divide(2) ;
  RooPlot* frame1_data = ak08Pruned_1_mass.frame(Bins(bins),Title("Pass sample")) ;
  RooPlot* frame2_data = ak08Pruned_1_mass.frame(Bins(bins),Title("Fail sample")) ;
  //dh_totalPass_data->plotOn(frame1_data);//,Cut("sample_MC==sample_MC::passed"));
  dataPass->plotOn(frame1_data);
  const_data_modelPass_sim_MC.plotOn(frame1_data);
  const_data_modelPass_sim_MC.plotOn(frame1_data,Components(const_data_dcb_pass_sim_MC),LineColor(kRed), Normalization(1.0,RooAbsReal::RelativeExpected));
  const_data_modelPass_sim_MC.plotOn(frame1_data,Components(const_data_cheby_pass_sim_MC),LineColor(kGreen), Normalization(1.0,RooAbsReal::RelativeExpected));
  //dh_totalFail_data->plotOn(frame2_data, DataError(RooAbsData::SumW2));
  dataFail->plotOn(frame2_data);
  /*gx_fail_sim_MC.plotOn(frame2_sim, LineColor(kRed));
  cheby_fail_sim_MC.plotOn(frame2_sim, LineColor(kGreen));*/
  const_data_model_fail_sim_MC.plotOn(frame2_data);
  const_data_model_fail_sim_MC.plotOn(frame2_data,Components(const_data_dcb_fail_sim_MC),LineColor(kRed));//, Normalization(1.0,RooAbsReal::RelativeExpected) );
  const_data_model_fail_sim_MC.plotOn(frame2_data,Components(const_data_cheby_fail_sim_MC),LineColor(kGreen));//, Normalization(1.0,RooAbsReal::RelativeExpected));
  /*c_sim_MC->cd(1) ;*/ gPad->SetLeftMargin(0.15) ; frame1_data->GetYaxis()->SetTitleOffset(1.4) ; frame1_data->Draw() ;
  //c_sim_MC->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_sim->GetYaxis()->SetTitleOffset(1.4) ; frame2_sim->Draw() ;
  c_sim_data_pass->SaveAs(Form("%s_provaRoofitDataPass.png",prefix.c_str()));
  c_sim_data_pass->SaveAs(Form("%s_provaRoofitDataPass.pdf",prefix.c_str()));
  
  TCanvas *c_sim_data_fail = new TCanvas("c_sim_data_fail", "Simultaneous PDF Fail",1);// 200, 10, 600, 400);
  gPad->SetLeftMargin(0.15) ; frame2_data->GetYaxis()->SetTitleOffset(1.4) ; frame2_data->Draw() ;
  c_sim_data_fail->SaveAs(Form("%s_provaRoofitDataFail.png",prefix.c_str()));
  c_sim_data_fail->SaveAs(Form("%s_provaRoofitDataFail.pdf",prefix.c_str()));

#ifdef doTOY
  //definition only for TOY fit
  RooRealVar const_efficiency_toy ("const_efficiency_toy", "efficiency",  .91, 0 ,1.);
  RooRealVar const_mean_pass_toy("const_mean_pass_toy", "Double CB mean", 78., 75.,85.);
  RooRealVar const_sigma_pass_toy("const_sigma_pass_toy", "Double CB Width", 7.,6.,8.); //it was 8.,5.,10.
  RooConstVar const_dCBCutL_pass_toy("const_dCBCutL_pass_toy", "Double CB Cut left", 1.77);//***** tau21_04 1.73);//***** tau21_06 1.5);
  RooConstVar const_dCBCutR_pass_toy("const_dCBCutR_pass_toy", "Double CB Cut right", 2.);//***** tau21_04  2.);//***** tau21_06 1.9);//1.9);
  RooConstVar const_dCBPowerL_pass_toy("const_dCBPowerL_pass_toy", "Double CB Power left", 1.16);//***** tau21_04 1.14);//***** tau21_06 1.16);
  RooConstVar const_dCBPowerR_pass_toy("const_dCBPowerR_pass_toy", "Double CB Power right", 6.09);//***** tau21_04 2.26);//***** tau21_06 2.2);// 4.7);
  RooDCBShape const_dcb_pass_toy("const_dcb_pass_toy", "double crystal ball", ak08Pruned_1_mass, const_mean_pass_toy, const_sigma_pass_toy, const_dCBCutL_pass_toy, const_dCBCutR_pass_toy, const_dCBPowerL_pass_toy, const_dCBPowerR_pass_toy);
  RooRealVar const_a_pass_toy("const_a_pass_toy","a_pass_sim_MC",.47,.3,.5);//, 0.01,100.);//***** tau21_04  .48);//***** tau21_06 .1);//.12);//1,-100,100);
  RooRealVar const_a1_pass_toy("const_a1_pass_toy","a1_pass_sim_MC",-.54, -.6,-.4);//,-100,-.1);//***** tau21_04 -.61);//***** tau21_06 -.5);//-.32);//0.1,-1,1);
  RooConstVar const_a2_pass_toy("const_a2_pass_toy","a2_pass_sim_MC",-.1);//***** tau21_04 -.19);//***** tau21_06  .07);//.07);//0.1,-1,1);
  RooChebychev const_cheby_pass_toy("const_cheby_pass_toy","cheby_pass_sim_MC",ak08Pruned_1_mass,RooArgSet(const_a_pass_toy,const_a1_pass_toy ,const_a2_pass_toy)) ;
  RooRealVar const_Nsig_toy("const_Nsig_toy", "Nsig_sim_MC", 1000.,0.,10000000.);
  RooFormulaVar const_k_pass_toy("const_k_pass_toy", "pass norm", "(const_Nsig_toy*const_efficiency_toy)", RooArgList(const_Nsig_toy, const_efficiency_toy));
  RooRealVar const_Nbkg_pass_toy("const_Nbkg_pass_toy", "Nbkg_pass_sim_MC", 1000., 0.,100000000.);
  RooAddPdf const_modelPass_toy("const_modelPass_toy", "modelPass_sim_MC", RooArgList(const_dcb_pass_toy,const_cheby_pass_toy), RooArgList(const_k_pass_toy, const_Nbkg_pass_toy));


  //+++ Fail part +++
  RooConstVar const_sigmaRatio("const_sigmaRatio", "ratio MC sigma", 1.25);
  //RooRealVar const_mean_fail_toy("const_mean_fail_toy", "Double CB mean", 75.4, 75.,85.);//***** tau21_04 76.53);//***** tau21_06 75.5);
  RooConstVar const_meanRatio("const_meanRatio", "ratio mean MC", .96);
  RooFormulaVar const_mean_fail_toy("const_mean_fail_toy", "Double CB mean", "(const_meanRatio*const_mean_pass_toy)", RooArgList(const_meanRatio, const_mean_pass_toy));
  RooRealVar const_sigma_fail_toy("const_sigma_fail_toy", "Double CB Width", 9.6,8,11); //working points 10, 0, 40 for tau21<.4
  //RooFormulaVar const_sigma_fail_toy("const_sigma_fail_toy", "Double CB Width","(const_sigmaRatio*const_sigma_pass_toy)", RooArgList(const_sigmaRatio, const_sigma_pass_toy));
  RooConstVar const_dCBCutL_fail_toy("const_dCBCutL_fail_toy", "Double CB Cut left",.8);//***** tau21_04 .7);//***** tau21_06 .8);//.42);// 0.18); //for tau21<0.4 1., 0.1, 50. works
  RooConstVar const_dCBCutR_fail_toy("const_dCBCutR_fail_toy", "Double CB Cut right", .82);//***** tau21_04 1.4);//***** tau21_06 1.5);//.2);//1.5);
  RooConstVar const_dCBPowerL_fail_toy("const_dCBPowerL_fail_toy", "Double CB Power left", 1.);//***** tau21_04 1.5);//***** tau21_06 2.5 .72);//.82); //.72);
  RooConstVar const_dCBPowerR_fail_toy("const_dCBPowerR_fail_toy", "Double CB Power right",2.24);//***** tau21_04  2.);//***** tau21_06 2.3);//2.8);//2.3); //working points 2., -0.2, 50. for tau21<0.6 using 2 DCB for LP and HP when same mu and sigma while 2., 0.2, 50. for tau21<0.4
  RooDCBShape const_dcb_fail_toy("const_dcb_fail_toy", "double crystal ball", ak08Pruned_1_mass, const_mean_fail_toy, const_sigma_fail_toy, const_dCBCutL_fail_toy, const_dCBCutR_fail_toy, const_dCBPowerL_fail_toy, const_dCBPowerR_fail_toy);
  RooRealVar const_a_fail_toy("const_a_fail_toy","a_fail_sim_MC", -.7,-2.,0.);//,-100,-.1);//***** tau21_04 -.66);//-.7);//***** tau21_06 -1.1);
  RooRealVar const_a1_fail_toy("const_a1_fail_toy","a1_fail_sim_MC", .02,0.,.1);//,.1,100);//***** tau21_04 -.01);//.04);//***** tau21_06 .34);
  RooConstVar const_a2_fail_toy("const_a2_fail_toy","a2_fail_sim_MC", .067);//***** tau21_04 .07);//-.013);//***** tau21_06 -.05);
  RooChebychev const_cheby_fail_toy("const_cheby_fail_toy","cheby_fail_sim_MC",ak08Pruned_1_mass,RooArgSet(const_a_fail_toy,const_a1_fail_toy, const_a2_fail_toy)) ;
  
  RooRealVar const_Nbkg_fail_toy("const_Nbkg_fail_toy", "Nbkg_fail_sim_MC", 1000., 0., 10000000.);
  RooFormulaVar const_k_fail_toy("const_k_fail_toy", "fail norm", "(const_Nsig_toy*(1-const_efficiency_toy))", RooArgList(const_Nsig_toy, const_efficiency_toy));
  RooAddPdf const_model_fail_toy("const_model_fail_toy", "model_fail_sim_MC", RooArgList(const_dcb_fail_toy,const_cheby_fail_toy), RooArgList(const_k_fail_toy, const_Nbkg_fail_toy));

  //end of constant part
  TTree* toyTree = new TTree("toyTree","data from my toy MC");

  //Toy generation part
  frame("Generating toys!!!");

  RooMCStudy* mcstudy_pass = new RooMCStudy(const_modelPass_sim_MC,ak08Pruned_1_mass,Binned(kTRUE),Silence(),Extended(), FitOptions(Save(kTRUE),PrintEvalErrors(0)));
  RooMCStudy* mcstudy_fail = new RooMCStudy(const_model_fail_sim_MC,ak08Pruned_1_mass,Binned(kTRUE),Silence(),Extended(), FitOptions(Save(kTRUE),PrintEvalErrors(0)));
  mcstudy_pass->generate(toyGen,MC_Pass->sumEntries()/*950*0.91+1600*/, kTRUE);//(const_Nbkg_pass_sim_MC.getValV()+(const_Nsig_sim_MC.getValV()*const_efficiency.getValV()))/10, kTRUE);
  mcstudy_fail->generate(toyGen,MC_Fail->sumEntries()/*950*0.09+563*/,kTRUE);//(const_Nbkg_fail_sim_MC.getValV()+(const_Nsig_sim_MC.getValV()*(1-const_efficiency.getValV())))/10, kTRUE);


  RooPlot* frame2_MCStudy = ak08Pruned_1_mass.frame(Bins(bins),Title("Pass sample"));
  RooPlot* frame3_MCStudy = ak08Pruned_1_mass.frame(Bins(bins),Title("Pass sample"));
  //(RooAbsData &)*(mcstudy->genData(i))
  double eff= const_efficiency.getValV();
  double mean=const_mean_pass_sim_MC.getValV();
  double sigma=const_sigma_pass_sim_MC.getValV();
  double a0=const_a_pass_sim_MC.getValV();
  double a1=const_a1_pass_sim_MC.getValV();
  double a2=const_a2_pass_sim_MC.getValV();
  double sigma_f=const_sigma_fail_sim_MC.getValV();
  double a0_f=const_a_fail_sim_MC.getValV();
  double a1_f=const_a1_fail_sim_MC.getValV();
  double a2_f=const_a2_fail_sim_MC.getValV();
  double Nbkg_f=const_Nbkg_fail_sim_MC.getValV();
  double Nbkg=const_Nbkg_pass_sim_MC.getValV();
  double Nsig=const_Nsig_sim_MC.getValV();
  //variables for tree
  double eff_toy, eff_toy_err;
  double mean_toy, mean_toy_err;
  double mean1_toy, mean1_toy_err;
  double sigma_toy, sigma_toy_err;
  double a0_toy, a0_toy_err;
  double a1_toy, a1_toy_err;
  double a2_toy, a2_toy_err;
  double sigma_f_toy, sigma_f_toy_err;
  double a0_f_toy, a0_f_toy_err;
  double a1_f_toy, a1_f_toy_err;
  double a2_f_toy, a2_f_toy_err;
  double Nbkg_f_toy, Nbkg_f_toy_err;
  double Nbkg_toy, Nbkg_toy_err;
  double Nsig_toy, Nsig_toy_err;
  toyTree->Branch("eff_toy", &eff_toy, "eff_toy/D");
  toyTree->Branch("eff_toy_err",&eff_toy_err,"eff_toy_err/D");
  toyTree->Branch("mean_toy", &mean_toy,"mean_toy/D");
  toyTree->Branch("mean_toy_err",&mean_toy_err,"mean_toy_err/D");
  toyTree->Branch("mean1_toy", &mean1_toy,"mean1_toy/D");
  toyTree->Branch("mean1_toy_err",&mean1_toy_err,"mean1_toy_err/D");
  toyTree->Branch("sigma_toy", &sigma_toy,"sigma_toy/D");
  toyTree->Branch("sigma_toy_err", &sigma_toy_err,"sigma_toy_err/D");
  toyTree->Branch("a0_toy", &a0_toy,"a0_toy/D");
  toyTree->Branch("a0_toy_err",&a0_toy_err,"a0_toy_err/D");
  toyTree->Branch("a1_toy",&a1_toy,"a1_toy/D");
  toyTree->Branch("a1_toy_err",&a1_toy_err,"a1_toy_err/D");
  //toyTree->Branch("a2_toy",&a2_toy, "a2_toy/D");
  //toyTree->Branch("a2_toy_err", &a2_toy_err,"a2_toy_err/D");
  toyTree->Branch("sigma_f_toy",&sigma_f_toy,"sigma_f_toy/D");
  toyTree->Branch("sigma_f_toy_err",&sigma_f_toy_err,"sigma_f_toy_err/D");
  toyTree->Branch("a0_f_toy",&a0_f_toy,"a0_f_toy/D");
  toyTree->Branch("a0_f_toy_err",&a0_f_toy_err,"a0_f_toy_err/D");
  toyTree->Branch("a1_f_toy",&a1_f_toy,"a1_f_toy/D");
  toyTree->Branch("a1_f_toy_err",&a1_f_toy_err,"a1_f_toy_err/D");
  //toyTree->Branch("a2_f_toy",&a2_f_toy,"a2_f_toy/D");
  //toyTree->Branch("a2_f_toy_err",&a2_f_toy_err,"a2_f_toy_err/D");
  toyTree->Branch("Nbkg_f_toy",&Nbkg_f_toy,"Nbkg_f_toy/D");
  toyTree->Branch("Nbkg_f_toy_err",&Nbkg_f_toy_err,"Nbkg_f_toy_err/D");
  toyTree->Branch("Nbkg_toy",&Nbkg_toy,"Nbkg_toy/D");
  toyTree->Branch("Nbkg_toy_err",&Nbkg_toy_err,"Nbkg_toy_err/D");
  toyTree->Branch("Nsig_toy",&Nsig_toy,"Nsig_toy/D");
  toyTree->Branch("Nsig_toy_err",&Nsig_toy_err,"Nsig_toy_err/D");
  TH1 *allBkgHisto_pass_test= new TH1D("allBkgHisto_fail_pass","allBkgHisto_pass", bins,min,max);
  TH1 *allBkgHisto_fail_test= new TH1D("allBkgHisto_fail_test","allBkgHisto_fail", bins,min,max);

  TH1F* pull_eff= new TH1F("pull_eff","pull_eff", 40,-10,10);
  TH1F* pull_mean= new TH1F("pull_mean","pull_mean", 40,-10,10);
  TH1F* pull_sigma= new TH1F("pull_sigma","pull_sigma", 40,-10,10);
  TH1F* pull_a0= new TH1F("pull_a0","pull_a0", 40,-10,10);
  TH1F* pull_a1= new TH1F("pull_a1","pull_a1", 40,-10,10);
  TH1F* pull_a2= new TH1F("pull_a2","pull_a2", 40,-10,10);
  TH1F* pull_sigma_f= new TH1F("pull_sigma_f","pull_sigma_f", 40,-10,10);
  TH1F* pull_a0_f= new TH1F("pull_a0_f","pull_a0_f", 40,-10,10);
  TH1F* pull_a1_f= new TH1F("pull_a1_f","pull_a1_f", 40,-10,10);
  TH1F* pull_a2_f= new TH1F("pull_a2_f","pull_a2_f", 40,-10,10);
  TH1F* pull_Nbkg_f= new TH1F("pull_Nbkg_f","pull_Nbkg_f", 40,-10,10);
  TH1F* pull_Nbkg= new TH1F("pull_Nbkg","pull_Nbkg", 40,-10,10);
  TH1F* pull_Nsig= new TH1F("pull_Nsig","pull_Nsig", 40,-10,10);
  RooDataHist *dh_totalPass_toy;
  RooDataHist *dh_totalFail_toy;
  RooDataSet *ds_toy_pass;
  RooDataSet *ds_toy_fail;
  RooSimultaneous simPdf_toy("simPdf_toy","simultaneous pdf for MC",sample_MC) ;
  simPdf_toy.addPdf(const_modelPass_toy,"passed") ;
  simPdf_toy.addPdf(const_model_fail_toy, "failed");
  int index=0;

  for(int i=0; i<toyGen; ++i){ 

    if((i%100)==0) std::cout<<"++++ Toy: "<<i<<std::endl;
    const_modelPass_toy.Clear();
    const_model_fail_toy.Clear();
    
    allBkgHisto_pass_test=((RooAbsData &)*(mcstudy_pass->genData(i))).createHistogram("ak08Pruned_1_mass",bins);
    allBkgHisto_fail_test=((RooAbsData &)*(mcstudy_fail->genData(i))).createHistogram("ak08Pruned_1_mass",bins);

    dh_totalPass_toy = new RooDataHist("dh_totalPass_toy", "dh_totalPass_data_test", ak08Pruned_1_mass, Import(*allBkgHisto_pass_test));
    dh_totalFail_toy = new RooDataHist("dh_totalPass_toy", "dh_totalPass_data_test", ak08Pruned_1_mass, Import(*allBkgHisto_fail_test));
    
//    ds_toy_pass =const_modelPass_sim_MC.generate(RooArgSet(ak08Pruned_1_mass),*MC_Pass, 1000);//(RooDataSet &)(mcstudy_pass->genData(i)) ;//new RooDataSet("ds_toy_pass","ds_toy_pass", ak08Pruned_1_mass,*(mcstudy_pass->genData(i)));
//    ds_toy_fail =const_model_fail_sim_MC.generate(RooArgSet(ak08Pruned_1_mass),*MC_Fail, 1000);//new RooDataSet("ds_toy_fail","ds_toy_fail", ak08Pruned_1_mass,Import((RooDataSet &)*mcstudy_fail->genData(i))); 
//    RooDataSet combData_MC_toy_unbin("combData_MC_toy_unbin","combData_MC_toy_unbin",ak08Pruned_1_mass,Index(sample_MC),Import("passed",*(ds_toy_pass)), Import("failed", *(ds_toy_fail)));
    RooDataHist combData_MC_toy("combData_MC_toy","combined data for MC",ak08Pruned_1_mass,Index(sample_MC),Import("passed",*(dh_totalPass_toy)), Import("failed",*(dh_totalFail_toy)));
    RooMsgService::instance().setSilentMode(true);
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    
    simPdf_toy.fitTo(combData_MC_toy, PrintEvalErrors(-1), SumW2Error(kTRUE)); 
    std::cout<<"++++++++++++++++++++ "<<const_efficiency_toy.getValV()<<std::endl;
    eff_toy=const_efficiency_toy.getValV(); eff_toy_err=const_efficiency_toy.getError();
    mean_toy=const_mean_pass_toy.getValV(); mean_toy_err=const_mean_pass_toy.getError();
    //mean1_toy=const_mean_fail_toy.getValV(); mean1_toy_err=const_mean_fail_toy.getError();
    sigma_toy=const_sigma_pass_toy.getValV(); sigma_toy_err=const_sigma_pass_toy.getError();
    a0_toy=const_a_pass_toy.getValV(); a0_toy_err=const_a_pass_toy.getError();
    a1_toy=const_a1_pass_toy.getValV(); a1_toy_err=const_a1_pass_toy.getError();
    //a2_toy=const_a2_pass_toy.getValV(); a2_toy_err=const_a2_pass_toy.getError();
    sigma_f_toy=const_sigma_fail_toy.getValV(); sigma_f_toy_err=const_sigma_fail_toy.getError();
    a0_f_toy=const_a_fail_toy.getValV(); a0_f_toy_err=const_a_fail_toy.getError();
    a1_f_toy=const_a1_fail_toy.getValV(); a1_f_toy_err=const_a1_fail_toy.getError();
    //a2_f_toy=const_a2_fail_toy.getValV(); a2_f_toy_err=const_a2_fail_toy.getError();
    Nbkg_f_toy=const_Nbkg_fail_toy.getValV(); Nbkg_f_toy_err=const_Nbkg_fail_toy.getError();
    Nbkg_toy=const_Nbkg_pass_toy.getValV(); Nbkg_toy_err=const_Nbkg_pass_toy.getError();
    Nsig_toy=const_Nsig_sim_MC.getValV(); Nsig_toy_err=const_Nsig_sim_MC.getError();
    //std::cout<<(const_efficiency.getValV()-eff)/const_efficiency.getError()<<" "<<(const_mean_pass_sim_MC.getValV()-mean)/const_mean_pass_sim_MC.getError()<<" "<<(const_sigma_pass_sim_MC.getValV()-sigma)/const_sigma_pass_sim_MC.getError()<<" "<<(const_a_pass_sim_MC.getValV()-a0)/const_a_pass_sim_MC.getError()<<" "<<(const_a1_pass_sim_MC.getValV()-a1)/const_a1_pass_sim_MC.getError()<<" "<<(const_a2_pass_sim_MC.getValV()-a2)/const_a2_pass_sim_MC.getError()<<" "<<(const_sigma_fail_sim_MC.getValV()-sigma_f)/const_sigma_fail_sim_MC.getError()<<" "<<std::endl; 
    pull_eff->Fill((const_efficiency_toy.getValV()-eff)/const_efficiency_toy.getError());
    std::cout<<eff<<" "<<eff_toy<<" "<<eff_toy_err<<std::endl;
    std::cout<<"Pull: "<<(const_efficiency_toy.getValV()-eff)/const_efficiency_toy.getError()<<std::endl;
    pull_mean->Fill((const_mean_pass_toy.getValV()-mean)/const_mean_pass_toy.getError());
    pull_sigma->Fill((const_sigma_pass_toy.getValV()-sigma)/const_sigma_pass_toy.getError());
    pull_a0->Fill((const_a_pass_toy.getValV()-a0)/const_a_pass_toy.getError());
    pull_a1->Fill((const_a1_pass_toy.getValV()-a1)/const_a1_pass_toy.getError());
    //pull_a2->Fill((const_a2_pass_toy.getValV()-a2)/const_a2_pass_toy.getError());
    pull_sigma_f->Fill((const_sigma_fail_toy.getValV()-sigma_f)/const_sigma_fail_toy.getError());
    pull_a0_f->Fill((const_a_fail_toy.getValV()-a0_f)/const_a_fail_toy.getError());
    pull_a1_f->Fill((const_a1_fail_toy.getValV()-a1_f)/const_a1_fail_toy.getError());
    ////pull_a2_f->Fill((const_a2_fail_toy.getValV()-a2_f)/const_a2_fail_toy.getError());
    pull_Nbkg_f->Fill((const_Nbkg_fail_toy.getValV()-Nbkg_f)/const_Nbkg_fail_toy.getError());
    pull_Nbkg->Fill((const_Nbkg_pass_toy.getValV()-Nbkg)/const_Nbkg_pass_toy.getError());
    pull_Nsig->Fill((const_Nsig_toy.getValV()-Nsig)/const_Nsig_toy.getError());
    //if(const_sigma_pass_toimPdf_MClV()<9) index=i;
    toyTree->Fill();

    dh_totalPass_toy->Clear();
    dh_totalFail_toy->Clear();
    combData_MC_toy.Clear();
    allBkgHisto_pass_test->Clear();
    allBkgHisto_fail_test->Clear();
    delete allBkgHisto_pass_test;
    delete allBkgHisto_fail_test;
  }//end loop over toys 

  

  dh_totalPass_toy->plotOn(frame2_MCStudy);
  //(mcstudy_pass->genData(0))->plotOn(frame2_MCStudy);
  const_modelPass_toy.plotOn(frame2_MCStudy);
  const_modelPass_toy.plotOn(frame2_MCStudy,Components(const_dcb_pass_toy),LineColor(kRed));//, Normalization(1.0,RooAbsReal::RelativeExpected));
  const_modelPass_toy.plotOn(frame2_MCStudy,Components(const_cheby_pass_toy),LineColor(kGreen));//, Normalization(1.0,RooAbsReal::RelativeExpected));

  dh_totalFail_toy->plotOn(frame3_MCStudy);
  const_model_fail_toy.plotOn(frame3_MCStudy);
  const_model_fail_toy.plotOn(frame3_MCStudy,Components(const_dcb_fail_toy),LineColor(kRed));//, Normalization(1.0,RooAbsReal::RelativeExpected) );
  const_model_fail_toy.plotOn(frame3_MCStudy,Components(const_cheby_fail_toy),LineColor(kGreen));//, Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas* c_toy = new TCanvas("c_toy","c_toy",1) ;
  
  //allBkgHisto_pass_test->Draw();
  gPad->SetLeftMargin(0.15) ; frame2_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame2_MCStudy->Draw();
  c_toy->SaveAs(Form("%s_toyMC.pdf",prefix.c_str()));
  TCanvas* c_toy1 = new TCanvas("c_toy1","c_toy",1) ;
  gPad->SetLeftMargin(0.15) ; frame2_MCStudy->GetYaxis()->SetTitleOffset(1.4) ; frame3_MCStudy->Draw();
  c_toy1->SaveAs(Form("%s_toyMC1.pdf",prefix.c_str())); 

  TFile *outFile = new TFile(Form("%s.root",prefix.c_str()), "RECREATE");

  TCanvas* c_MCStudy1 = new TCanvas("c_MCStudy1","c_toy",1) ;
  pull_eff->Draw("hist");
  pull_eff->Write();
  TCanvas* c_MCStudy3 = new TCanvas("c_MCStudy3","c_toy",1) ;
  pull_mean->Draw("hist");
  pull_mean->Write();
  TCanvas* c_MCStudy4 = new TCanvas("c_MCStudy4","c_toy",1) ;
  pull_sigma->Draw("hist");
  pull_sigma->Write();
  TCanvas* c_MCStudy5 = new TCanvas("c_MCStudy5","c_toy",1) ;
  pull_a0->Draw("hist");
  pull_a0->Write();
  TCanvas* c_MCStudy6 = new TCanvas("c_MCStudy6","c_toy",1) ;
  pull_a1->Draw("hist");
  pull_a1->Write();
  TCanvas* c_MCStudy7 = new TCanvas("c_MCStudy7","c_toy",1) ;
  pull_a2->Draw("hist");
  pull_a2->Write();
  TCanvas* c_MCStudy8 = new TCanvas("c_MCStudy8","c_toy",1) ;
  pull_sigma_f->Draw("hist");
  pull_sigma_f->Write();
  TCanvas* c_MCStudy9 = new TCanvas("c_MCStudy9","c_toy",1) ;
  pull_a0_f->Draw("hist");
  pull_a0_f->Write();
  TCanvas* c_MCStudy10 = new TCanvas("c_MCStudy10","c_toy",1) ;
  pull_a1_f->Draw("hist");
  pull_a1_f->Write();
  TCanvas* c_MCStudy11 = new TCanvas("c_MCStudy11","c_toy",1) ;
  pull_a2_f->Draw("hist");
  pull_a2_f->Write();
  TCanvas* c_MCStudy12 = new TCanvas("c_MCStudy12","c_toy",1) ;
  pull_Nbkg_f->Draw("hist");
  pull_Nbkg_f->Write();
  TCanvas* c_MCStudy13 = new TCanvas("c_MCStudy13","c_toy",1) ;
  pull_Nbkg->Draw("hist");
  pull_Nbkg->Write();
  TCanvas* c_MCStudy14 = new TCanvas("c_MCStudy14","c_toy",1) ;
  pull_Nsig->Draw("hist");
  pull_Nsig->Write();

  toyTree->Write();
  outFile->Close();
  c_MCStudy1->SaveAs(Form("%s_mcstudy1.pdf",prefix.c_str()));
  c_MCStudy3->SaveAs(Form("%s_mcstudy3.pdf",prefix.c_str()));
  c_MCStudy4->SaveAs(Form("%s_mcstudy4.pdf",prefix.c_str()));
  c_MCStudy5->SaveAs(Form("%s_mcstudy5.pdf",prefix.c_str()));
  c_MCStudy6->SaveAs(Form("%s_mcstudy6.pdf",prefix.c_str()));
  c_MCStudy7->SaveAs(Form("%s_mcstudy7.pdf",prefix.c_str()));
  c_MCStudy8->SaveAs(Form("%s_mcstudy8.pdf",prefix.c_str()));
  c_MCStudy9->SaveAs(Form("%s_mcstudy9.pdf",prefix.c_str()));
  c_MCStudy10->SaveAs(Form("%s_mcstudy10.pdf",prefix.c_str()));
  c_MCStudy11->SaveAs(Form("%s_mcstudy11.pdf",prefix.c_str()));
  c_MCStudy12->SaveAs(Form("%s_mcstudy12.pdf",prefix.c_str()));
  c_MCStudy13->SaveAs(Form("%s_mcstudy13.pdf",prefix.c_str()));
  c_MCStudy14->SaveAs(Form("%s_mcstudy14.pdf",prefix.c_str()));
  
  std::cout<<"Generated: "<<(Nsig*eff)+Nbkg<<" for passed toy"<<std::endl;
  std::cout<<"Generated: "<<(Nsig*(1-eff))+Nbkg_f<<" for passed toy"<<std::endl;

#endif //doTOY

  std::cout<<"Pass: "<<MC_Pass->numEntries()<<" "<<MC_Pass->sumEntries()<<std::endl;
  std::cout<<"Fail: "<<MC_Fail->numEntries()<<" "<<MC_Fail->sumEntries()<<std::endl;


  std::cout<<""<<std::endl;
  std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << timer.RealTime() << " s" <<std::endl;
  std::cout<<""<<std::endl;
  
  return 0;

}
