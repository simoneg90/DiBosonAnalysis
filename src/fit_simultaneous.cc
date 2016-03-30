#include "doPlots.h"
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
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"


#include <fstream>


using namespace RooFit;

#define MAX_bkg 8 //also the n of colors
#define MAX_sgn 8
#define mass_min 45
#define mass_max 130



int main(int argc, char* argv[]){

  Color_t COLOR[] = {kAzure-1, kRed, kCyan, kGreen+2, kOrange+10, kSpring-9, kCyan+0, kBlue+0};

  frame("Fit environment - VTagging");
  std::cout<<"Preliminary infos: in his program we use a dedicated dictionary to associate the right BKG"<<std::endl;
  std::cout<<"1: ttbar"<<std::endl;
  std::cout<<"2: SingleTop"<<std::endl;
  std::cout<<"3: WJets"<<std::endl;
  std::cout<<"4: Diboson"<<std::endl;

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
  
  RooRealVar Nsig_MC("Nsig_MC", "Nsig_MC", 1200,10,10000);
  RooRealVar Nbkg_MC("Nbkg_MC", "Nbkg_MC", 100,10,100000);
  RooAddPdf modelPass_MC("modelPass_MC", "modelPass_MC", RooArgList(gx_pass_MC,p2_pass_MC), RooArgList(Nsig_MC, Nbkg_MC));



  TH1D *histoD_pass[MAX_NUMBER];//histos for project High Purity
  TH1D *histoD_fail[MAX_NUMBER];//histos for project Low Purity
  TH1D *bkgCounts[MAX_NUMBER];//histos to count events to rescale
  TH1D *histoBkg_pass[MAX_bkg];    //1 histo for each background HP
  TH1D *histoBkg_fail[MAX_bkg];    //1 histo for each background LP
  TH1D *dataHisto_pass; //data histo HP
  TH1D *dataHisto_fail; //data histo LP
  TH1D *allBkgHisto;
  TH1D *signalHisto[MAX_sgn]; //signal histos (if we want to superimpose them)
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
  
  while(inputList>>bkg_nameTMP>>bkg_file>>scaleFactor>>bins>>min>>max){
   std::cout<<"BKG: "<<bkg_nameTMP.c_str()<<" File: "<<bkg_file.c_str()<<" scaleFactor: "<<scaleFactor<<std::endl;

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
   TDirectory * dir = (TDirectory*)file[file_counter]->Get("rootTupleTree");
   dir->GetObject("tree", tree[file_counter]);
   std::cout<<"Creating Histo"<<std::endl;
   histoD_pass[file_counter] = new TH1D(Form("histo_%d_pass", file_counter), Form("histo_%d_pass", file_counter),bins,min,max);
   histoD_fail[file_counter] = new TH1D(Form("histo_%d_fail", file_counter), Form("histo_%d_fail", file_counter),bins,min,max);
   std::cout<<"Projecting"<<std::endl;
   tree[file_counter]->Project(Form("histo_%d_pass", file_counter),Form("%s", variable.c_str()),pass.c_str());
   tree[file_counter]->Project(Form("histo_%d_fail", file_counter),Form("%s", variable.c_str()),fail.c_str());

   if(strcmp(bkg_nameTMP.c_str(), "data")==0){
     std::cout<<"DATA!"<<std::endl;
     if(isData==0) {
       allBkgHisto= new TH1D("allBkgHisto","allBkgHisto", bins,min,max);
       dataHisto_pass= new TH1D("dataHisto_pass","dataHisto_pass", bins,min,max);
       dataHisto_fail= new TH1D("dataHisto_fail","dataHisto_fail", bins,min,max);
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
      }
      isData=1;
      dataHisto_pass->Add(histoD_pass[file_counter]);
      dataHisto_fail->Add(histoD_fail[file_counter]);

    }else{
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





    }//end if data or NotData


  }//end loop over rootList file

  RooRealVar t("t", "t", 2, 0,10);
  std::cout<<t.getValV()<<std::endl;




  return 0;


}
