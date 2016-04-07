#include "doPlots.h"
#include "setTDRStyle.h"
#include "utility.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "time.h"
#include "CMS_lumi.h"
//#include "CMS_lumi.C"

#include <fstream>

#define MAX_bkg 8 //also the n of colors
#define MAX_sgn 8
/*
TPaveText* get_labelCMStop(bool wide) const {

  //TLatex *latex = new TLatex();
  //latex->SetNDC();
  TPaveText* label_cmstop = new TPaveText(0.10, 0.94, 0.96, 0.98, "brNDC");
  label_cmstop->SetTextSize(0.045);
  label_cmstop->SetFillColor(0);
  label_cmstop->SetTextFont(42);

  label_cmstop->SetTextAlign(31); // align right
  //latex->DrawLatex(wide ? 0.98 : 0.95, 0.96, "#sqrt{s} = 7 TeV");
  label_cmstop->AddText("#sqrt{s} = 13 TeV");
  std::string leftText;
  if (dataFiles_.size() == 0) {
    leftText = "CMS Simulation 2015";
  } else {
    if (isCMSArticle_) {
      leftText = "CMS 2015";
    } else {
      leftText = "CMS Preliminary 2015";
    }
  }

  if (lumi_ > 0.) {
    label_cmstop->SetTextAlign(11); // align left
    std::string lumiText = this->get_lumiText();
    if (dataFiles_.size() == 0) {
      lumiText = "L = " + lumiText;
    }
    //label_cmstop->DrawLatex(wide ? 0.06 : 0.15, 0.96, Form("%s, %s", leftText.c_str(), lumiText.c_str()));
    label_cmstop->AddText(Form("%s, %s", leftText.c_str(), lumiText.c_str()));
  } else {
    label_cmstop->SetTextAlign(11); // align left
    //label_cmstop->DrawLatex(wide ? 0.06 : 0.15, 0.96, Form("%s", leftText.c_str()));
    label_cmstop->AddText(Form("%s", leftText.c_str()));
  }

  return label_cmstop;

} // cmsPrel



TPaveText* get_labelCMS(int legendQuadrant, bool hasRatio) const {

  if (legendQuadrant != 0 && legendQuadrant != 1 && legendQuadrant != 2 && legendQuadrant != 3) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for CMS label. Using 2." << std::endl;
    legendQuadrant = 2;
  }

  float x1, y1, x2, y2;
  if (legendQuadrant == 1) {
    x1 = 0.63;
    y1 = 0.86;
    x2 = 0.8;
    y2 = 0.92;
  } else if (legendQuadrant == 2) {
    x1 = 0.10;
    y1 = 0.86;
    x2 = (isCMSArticle_) ? 0.39 : 0.42;
    y2 = 0.92;
  } else if (legendQuadrant == 3) {
    x1 = 0.25;
    y1 = 0.2;
    x2 = 0.42;
    y2 = 0.24;
  } else if (legendQuadrant == 0) {
    x1 = hasRatio ? 0.25 : 0.30;
    y1 = 0.963;
    x2 = 0.65;
    y2 = 0.985;
  }


  TPaveText* cmslabel = new TPaveText(x1, y1, x2, y2, "brNDC");
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  if (legendQuadrant == 0) {
    //cmslabel->SetTextAlign(11);
  }
  cmslabel->SetTextFont(42);
  std::string label_CMS_text = this->get_CMSText();
  if (legendQuadrant != 0) {
    cmslabel->AddText(label_CMS_text.c_str());
  } else {
    std::string leftText;
    if (dataFiles_.size() == 0) {
      leftText = "CMS Simulation 2015";
    } else {
      if (isCMSArticle_) {
	leftText = "CMS 2015";
      } else {
	leftText = "CMS Preliminary 2015";
      }
    }
    if (lumi_ > 0.) {
      //cmslabel->SetTextAlign(11); // align left
      std::string lumiText = this->get_lumiText();
      cmslabel->AddText(Form("%s, %s", leftText.c_str(), lumiText.c_str()));
    } else {
      //cmslabel->SetTextAlign(11); // align left
      cmslabel->AddText(Form("%s", leftText.c_str()));
    }
  }

  return cmslabel;

}

TPaveText* get_labelSqrt(int legendQuadrant) const {

  if (legendQuadrant != 0 && legendQuadrant != 1 && legendQuadrant != 2 && legendQuadrant != 3) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for Sqrt label. Using 2." << std::endl;
    legendQuadrant = 2;
  }


  float x1, y1, x2, y2;
  if (legendQuadrant == 1) {
    x1 = 0.63;
    y1 = 0.82;
    x2 = 0.8;
    y2 = 0.86;
  } else if (legendQuadrant == 2) {
    x1 = (isCMSArticle_) ? 0.22 : 0.25;
    y1 = 0.82;
    x2 = (isCMSArticle_) ? 0.39 : 0.42;
    y2 = 0.86;
  } else if (legendQuadrant == 3) {
    x1 = 0.25;
    y1 = 0.16;
    x2 = 0.42;
    y2 = 0.2;
  } else if (legendQuadrant == 0) {
    x1 = 0.7;
    y1 = 0.953;
    x2 = 0.96;
    y2 = 0.975;
  }


  TPaveText* label_sqrt = new TPaveText(x1, y1, x2, y2, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  std::string label_sqrt_text = this->get_sqrtText();
  if (legendQuadrant != 0) {
    label_sqrt->AddText(label_sqrt_text.c_str());
  } else {
    label_sqrt->SetTextAlign(31); // align right
    label_sqrt->AddText("#sqrt{s} = 13 TeV");
  }

  return label_sqrt;

}


TPaveText* get_labelAlgo(int legendQuadrant) const {


  float x1, y1, x2, y2;
  if (legendQuadrant == 1) {
    x1 = 0.77;
    y1 = 0.88;
    x2 = 0.82;
    y2 = 0.92;
  } else if (legendQuadrant == 2) {
    x1 = 0.30;
    y1 = 0.86;
    x2 = 0.35;
    y2 = 0.92;
  } else if (legendQuadrant == 3) {
    x1 = 0.3;
    //y1 = 0.15;
    y1 = 0.18;
    x2 = 0.35;
    //y2 = 0.2;
    y2 = 0.21;
  } else if (legendQuadrant == 4) {
    x1 = 0.75;
    y1 = 0.18;
    x2 = 0.8;
    y2 = 0.21;
  } else {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for Algo label. Using 3." << std::endl;
    x1 = 0.27;
    y1 = 0.15;
    x2 = 0.32;
    y2 = 0.2;
  }

  Float_t labelTextSize = 0.035;
  std::string jetAlgoName = (recoType_ != "" && jetAlgo_ != "") ? get_algoName() : "";
  TPaveText* label_algo = new TPaveText(x1, y1, x2, y2, "brNDC");
  label_algo->SetTextFont(42);
  label_algo->SetFillColor(kWhite);
  label_algo->SetTextSize(labelTextSize);
  label_algo->AddText(jetAlgoName.c_str());
  //label_algo->SetTextAlign(11);

  return label_algo;

}

*/
int main(int argc, char* argv[]){

  frame("Plotter env");
  
  setTDRStyle(); 
  //Color_t COLOR[] = {kRed+1, kGreen+3, kAzure-1, kViolet-3, kOrange+10, kSpring-9, kCyan+0, kBlue+0};
  Color_t COLOR[] = {kGreen+2,kAzure-1, kCyan, kSpring-9,kOrange+9,kGreen+2,kBlue+0, kOrange+8, kCyan+0};

  TH1D *histoD[MAX_NUMBER];//histos for project
  TH1D *bkgCounts[MAX_NUMBER];//histos to count events to rescale
  TH1D *histoBkg[MAX_bkg];    //1 histo for each background
  TH1F *histoF[MAX_NUMBER];//histos for project
  TH1D *dataHisto; //data histo
  TH1D *allBkgHisto; 
  TH1D *signalHisto[MAX_sgn]; //signal histos
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
    std::cout<<"Digit './doPlots --help' for informations!"<<std::endl;
    breakLine();
    exit(-1);
  }

  if(strcmp(argv[1],"--help")==0){
     breakLine();
     std::cout<<"This program superimposes or Stacks the same histogram from different files"<<std::endl;
     std::cout<<"Program usage:"<<std::endl;
     std::cout<<"./doPlots inputList lumi[fb-1] output eff_file"<<std::endl;
     //argv[1]=inputList 
     breakLine();
     exit(-1);
  }

  if(argc>5){
     breakLine();
     std::cout<<"ATTENTION! TOO MANY ARGUMENTS ADDED!"<<std::endl;
     std::cout<<"Exiting program"<<std::endl;
     breakLine();
     exit(-1);
  }
  
  std::ifstream inputList;
  inputList.open(argv[1]);
  if(inputList==NULL){
    std::cout<<"ERROR! File: "<<argv[1]<<" doesn't exist!"<<std::endl;
    exit(-1);
  }
  std::cout<<"Opened: " <<argv[1]<<std::endl;
  std::cout<<"OutPutFile: "<< argv[3]<<std::endl;
  std::string outString = argv[3];
  
  std::cout<<"Will be writing .pdf, .png, .root: "<<outString.c_str()<<std::endl;
  std::string filePath;
  std::string legendName[MAX_NUMBER];
  //Input file variables
  std::string bkg_name, bkg_nameTMP, bkg_file, variable, cut;
  double scaleFactor=1.;
  double lumi=(atof(argv[2]))*1000;
  bkg_name="";
  std::string xTitle="";
  std::string yTitle="";
  int file_counter;
  file_counter=0;
  double bins, min, max;
  int bkg_counter=-1;
  double integralData, integralBKG;
  integralData=integralBKG=0;
  double counts=0;
  double integralBKG_all=0;
  int sgn_counter=0;
  bool isData=0;
  THStack *bkgStack = new THStack("bkgStack","");
  while(inputList>>bkg_nameTMP>>bkg_file>>scaleFactor>>variable>>cut>>bins>>min>>max>>xTitle>>yTitle){
    std::cout<<"BKG: "<<bkg_nameTMP.c_str()<<" File: "<<bkg_file.c_str()<<" scaleFactor: "<<scaleFactor<<std::endl;
    if(strcmp(cut.c_str(), "-")==0) cut="";
    std::cout<<"CUT: "<<cut.c_str()<<std::endl;
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
    histoD[file_counter] = new TH1D(Form("histo_%d", file_counter), Form("histo_%d", file_counter),bins,min,max);
    std::cout<<"Projecting"<<std::endl;
    tree[file_counter]->Project(Form("histo_%d", file_counter),Form("%s", variable.c_str()),cut.c_str());
    //std::cout<<"Scaling"<<std::endl;
    //histoD[file_counter]->Scale(scaleFactor);//rescale for the crossSection and the number of events

    if(strcmp(bkg_nameTMP.c_str(), "data")==0){
      std::cout<<"DATA!"<<std::endl;
      if(isData==0) {
        allBkgHisto= new TH1D("allBkgHisto","allBkgHisto", bins,min,max);
        dataHisto= new TH1D("dataHisto","dataHisto", bins,min,max);
        std::cout<<"Data counts: "<<counts<<std::endl;
        //dataHisto->GetXaxis()->SetTitle(xTitle.c_str());
        dataHisto->GetYaxis()->SetTitle(yTitle.c_str());
        dataHisto->SetMarkerStyle(8);
        dataHisto->SetMarkerColor(1);
        dataHisto->SetLineColor(1);
        leg->AddEntry(dataHisto, "Data", "pe");
      }
      isData=1;
      dataHisto->Add(histoD[file_counter]);
    }else if(strcmp(bkg_nameTMP.c_str(), "signal")==0){
      signalHisto[sgn_counter]= new TH1D(Form("signalHisto_%d",sgn_counter),Form("signalHisto_%d",sgn_counter), bins,min,max);
      leg->AddEntry(signalHisto[sgn_counter],"Signal", "l");
      signalHisto[sgn_counter]->Add(histoD[file_counter]);
      ++sgn_counter;

    }else{
      histoD[file_counter]->Scale(scaleFactor*lumi/counts);
      std::cout<<" " <<scaleFactor<<" "<<lumi<<" "<<counts<<std::endl;
      std::cout<<"Total scale factor: "<<scaleFactor*lumi/counts<<std::endl;
      if(strcmp(bkg_name.c_str(),bkg_nameTMP.c_str())==0){//if the bkg is always the same
        std::cout<<"Background: "<<bkg_name.c_str()<<std::endl;
        histoBkg[bkg_counter]->Add(histoD[file_counter]);
      }//end for the same bkg
      else{
        bkg_name=bkg_nameTMP;
        //if(bkg_counter!=-1) bkgStack->Add(histoBkg[bkg_counter]);
        ++bkg_counter;
        histoBkg[bkg_counter]= new TH1D(Form("histo_%s", bkg_name.c_str()), Form("histo_%s", bkg_name.c_str()),bins,min,max);
        histoBkg[bkg_counter]->SetStats(0);
        //histoBkg[bkg_counter]->SetLineColor(COLOR[bkg_counter]);
        histoBkg[bkg_counter]->SetFillColor(COLOR[bkg_counter]);
        histoBkg[bkg_counter]->SetLineWidth(1);
        histoBkg[bkg_counter]->SetLineColor(1);
        histoBkg[bkg_counter]->Add(histoD[file_counter]);
        leg->AddEntry(histoBkg[bkg_counter], bkg_nameTMP.c_str(), "f");
        std::cout<<"New Background: "<<bkg_name.c_str()<<std::endl;
      }

    }
    //file[file_counter]=TFile::Close("R");
    ++file_counter;
   
  }//end while over inputFile
  std::cout<<"counter del bkg "<<bkg_counter<<std::endl; 
  //TCanvas *c[10];
  TFile *outFile = new TFile(Form("%s.root",outString.c_str()), "RECREATE");
  for(int i=0; i<=bkg_counter; ++i){
    c[i]=new TCanvas(Form("c%d",i), "Grafico1",1);// 200, 10, 600, 400);
    histoBkg[i]->Draw("hist");
    c[i]->SaveAs(Form("c%d.png",i));
    std::cout<<"Integral: "<<histoBkg[i]->Integral()<<std::endl;
    integralBKG_all+=histoBkg[i]->Integral();
    std::cout<<"All: "<<integralBKG_all<<" "<<histoBkg[i]->Integral()<<std::endl;
    std::cout<<"Entries: "<<histoBkg[i]->GetEntries()<<std::endl;
    if(isData==1) allBkgHisto->Add(histoBkg[i]);
    //integralBKG=allBkgHisto->Integral();
    std::cout<<"adding back"<<std::endl;
    histoBkg[i]->Write();
    bkgStack->Add(histoBkg[i]);
  }
  
  //setTDRStyle();
  TCanvas *c_histo = new TCanvas("c_histo", "Grafico1", 1);//200, 10, 600, 400);
  std::cout<<"Defined canvas"<<std::endl; 

  std::ofstream effList;
  std::string effFile = argv[4];
  std::cout<<"Writing: "<<effFile.c_str()<<std::endl;
  effList.open(Form("%s.txt",effFile.c_str()), std::fstream::app);
  
  if(isData==1){
    TPad *pad2=new TPad("pad2", "bottom pad", 0.,0.05,1,0.3);
    pad2->SetBottomMargin(0.2);
    pad2->SetTopMargin(0);
    pad2->Draw();
    TPad *pad1=new TPad("pad1", "top pad", 0,0.3,1,1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    dataHisto->SetStats(0);
    dataHisto->Write();
    std::cout<<"DataHisto Entries: "<<dataHisto->GetEntries()<<std::endl;   
    std::cout<<"DataHisto integral: "<<dataHisto->Integral()<<std::endl;
    std::cout<<"BKGHisto integral: "<<integralBKG_all<<std::endl;
    dataHisto->GetXaxis()->SetLabelSize(0);
    dataHisto->GetXaxis()->SetLabelOffset(999);
    dataHisto->GetYaxis()->SetLabelFont(43);
    dataHisto->GetYaxis()->SetLabelSize(14);
    //dataHisto->Draw("EP");
    dataHisto->SetTitle("");
    std::cout<<"Drawing data"<<std::endl;
    allBkgHisto->Write();
    bkgStack->Write();
    //bkgStack->Draw("");
    std::cout<<"Drawing bkg"<<std::endl;
    dataHisto->Draw("EP");
    bkgStack->Draw("SAME");
    dataHisto->Draw("EPSAME");
    std::cout<<"Redrawing data"<<std::endl;
    TGraphErrors *gr = new TGraphErrors(0);
    integralData=dataHisto->Integral();
    integralBKG=allBkgHisto->Integral();
    double error, ratio;
    frame("Graph values");
    for(int w=1; w<bins; ++w){
      if((dataHisto->GetBinContent(w)==0&&allBkgHisto->GetBinContent(w)==0) || (allBkgHisto->GetBinContent(w)==0)) continue;
      gr->SetPoint(w, dataHisto->GetBinCenter(w),(dataHisto->GetBinContent(w))/(allBkgHisto->GetBinContent(w)));
      ratio= (dataHisto->GetBinContent(w))/(allBkgHisto->GetBinContent(w));
      error= (dataHisto->GetBinContent(w)*sqrt(allBkgHisto->GetBinContent(w)) + allBkgHisto->GetBinContent(w)*sqrt(dataHisto->GetBinContent(w)))/(allBkgHisto->GetBinContent(w)*allBkgHisto->GetBinContent(w));
      std::cout<<"VALUE: "<<ratio<<" ERROR: "<<error<<std::endl;
      gr->SetPointError(w, dataHisto->GetBinWidth(w)/2,error);

    }
    gr->GetHistogram()->SetMaximum(2);//1.5);
    gr->GetHistogram()->SetMinimum(0.1);
    gr->GetXaxis()->SetLimits(min, max);
    //TGraphAsymmErrors *gr= new TGraphAsymmErrors();
    //gr->Divide(dataHisto,allBkgHisto, "pois");
    pad2->cd();
    gStyle->SetTextSize(14);
    gROOT->ForceStyle();
    pad2->SetGrid();
    gr->Draw("ZAP");
    gr->GetXaxis()->SetLabelFont(43);
    gr->GetXaxis()->SetLabelSize(14);
    gr->GetYaxis()->SetLabelFont(43);
    gr->GetYaxis()->SetLabelSize(14);
    allBkgHisto->GetXaxis()->SetLabelFont(43);
    allBkgHisto->GetXaxis()->SetLabelSize(14);
    allBkgHisto->GetYaxis()->SetLabelFont(43);
    allBkgHisto->GetYaxis()->SetLabelSize(14);
    //allBkgHisto->Draw("hist");
    pad1->cd();
    effList<<integralData<<" ";
    for(int ii=0; ii<=bkg_counter; ++ii){
      effList<<histoBkg[ii]->Integral()<<" ";
    }
    effList<<integralBKG_all<<std::endl;
  }else{
    std::cout<<"entered if only BKG"<<std::endl;
  /*  bkgStack->GetXaxis()->SetTitle(xTitle.c_str());
    bkgStack->GetYaxis()->SetTitle(yTitle.c_str());
    */bkgStack->Draw();
    std::cout<<"Drawing BKG"<<std::endl;
  }
  
  for(int i=0; i<sgn_counter;++i){
    std::cout<<"Drawing signal n. "<<sgn_counter<<std::endl;
    signalHisto[i]->SetStats(0);
    signalHisto[i]->Draw("histSAME");

  }
  leg->Draw("SAME");

  /*Float_t cmsTextSize = 0.043;
  TPaveText* label_cms = get_labelCMS(1);
  label_cms->SetTextSize(cmsTextSize);
  TPaveText* label_sqrt = get_labelSqrt(1);
  TPaveText* label_algo = get_labelAlgo(3);

  label_cms->Draw("SAME");
  label_sqrt->Draw("SAME");
  label_algo->Draw("SAME");*/
  //std::cout<<"Writing outTree"<<std::endl;
  //TTree *T2 = tree[file_counter]->CopyTree(0);
  //std::cout<<"Drawing legend"<<std::endl;
  //T2->Write();
  c_histo->SaveAs(Form("%s.png", outString.c_str()));
  c_histo->SaveAs(Form("%s.pdf", outString.c_str()));
  //c_histo->SaveAs("histograms.png");
  //c_histo->SaveAs("histograms.pdf");

  //gROOT->ProcessLine(".L src/CMS_lumi.C");//
  //gROOT->LoadMacro("src/CMS_lumi.C");
  //TPad *pad1 = new TPad("pad1", "pad1",0,0.15,1,1) ;
  //CMS_lumi(c2, 3,10);
  //gROOT->ProcessLine("CMS_lumi(c2,1,0)");
  //c2->Update();
  //c2->RedrawAxis();
  //c2->GetFrame()->Draw();
  //c2->SaveAs("histograms_CMS.png");

  //std::string outName = "outfile.root";
  

  std::cout<<integralData<<" "<<integralBKG_all<<std::endl;
  std::cout<<"Wrote: "<<outString.c_str()<<std::endl; 
  c_histo->Write();
  std::cout<<"Closing output file"<<std::endl;
  outFile->Close();

  std::cout<<""<<std::endl;
  std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << timer.RealTime() << " s" <<std::endl;
  std::cout<<""<<std::endl;
}
