#include "doPlots.h"
#include "setTDRStyle.h"
#include "utility.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "TLatex.h"
#include "time.h"
//#include "CMS_lumi.h"
//#include "CMS_lumi.C"

#include <fstream>

#define MAX_bkg 8 //also the n of colors
#define MAX_sgn 8


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




int main(int argc, char* argv[]){

  TH1::SetDefaultSumw2();
  frame("Plotter env");
  
  //setTDRStyle(); 
  //Color_t COLOR[] = {kRed+1, kGreen+3, kAzure-1, kViolet-3, kOrange+10, kSpring-9, kCyan+0, kBlue+0};
  //Color_t COLOR[] = {kGreen+2,kAzure-1, kCyan, kSpring-9,kOrange+9,kGreen+2,kBlue+0, kOrange+8, kCyan+0};
  //Color_t COLOR[] = {kAzure-1, kRed, kCyan, kSpring+7,kGreen+2, kViolet, kOrange+10, kSpring-9, kCyan+0, kBlue+0};
  //Color_t COLOR[] = {kAzure-1, kViolet, kCyan, kSpring+7,kGreen+2, kOrange+1, kSpring+9, kRed+1, kBlue+0}; //ZGamma MIO
  Color_t COLOR[] = {kOrange-3,kGreen-7,kMagenta-7,kAzure-3,kBlue,kViolet+2 ,kGreen-3, kRed-3};
  //Color_t COLOR[] = {kAzure+8,kAzure-2, kTeal+9, kSpring+7, kPink+7, kPink+6,kOrange-3,kOrange-2,kOrange-1, kMagenta, kCyan};

  //gStyle->SetStripDecimals(kTRUE);
  //gStyle->SetTickLength(0.03, "XYZ");
  //gStyle->SetNdivisions(510, "XYZ");
  //gStyle->SetPadTickY(.1);
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
  TLegend *leg = new TLegend(0.65,0.5,0.8,0.85);
  leg->SetFillStyle(3004);
  leg->SetTextSize(.04);
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
     std::cout<<"./doPlots inputList lumi[fb-1] output eff_file doScale superimpose"<<std::endl;
     //argv[1]=inputList 
     breakLine();
     exit(-1);
  }

  if(argc>8){
     breakLine();
     std::cout<<"ATTENTION! TOO MANY ARGUMENTS ADDED!"<<std::endl;
     std::cout<<"Exiting program"<<std::endl;
     breakLine();
     exit(-1);
  }
  
  std::ifstream inputList;
  inputList.open(argv[1]);
  if(!inputList){//==NULL){
    std::cout<<"ERROR! File: "<<argv[1]<<" doesn't exist!"<<std::endl;
    exit(-1);
  }
  std::cout<<"Opened: " <<argv[1]<<std::endl;
  std::cout<<"OutPutFile: "<< argv[3]<<std::endl;
  std::string outString = argv[3];
  int doScale = atoi(argv[5]);
  int superimpose = atoi(argv[6]);
  int doLogY = atoi(argv[7]);
  
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
  std::cout<<"Before while"<<std::endl;
  int line_counter=0;
  //When working at Rome T2
  //while(inputList>>bkg_nameTMP>>bkg_file>>scaleFactor>>variable>>cut>>bins>>min>>max>>xTitle>>yTitle){
  //When working at lxplus
  while(true){
    
    inputList>>bkg_nameTMP>>bkg_file>>scaleFactor>>variable>>cut>>bins>>min>>max>>xTitle>>yTitle;
    ++line_counter;
    if(inputList.eof()) break;
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
    //exit(1);
    bkgCounts[file_counter]=((TH1D *)(file[file_counter]->Get("DijetFilter/EventCount/EventCounter")));
    std::cout<<"Total Events: "<<bkgCounts[file_counter]->GetBinContent(1)<<std::endl;
    counts=bkgCounts[file_counter]->GetBinContent(1);
    TDirectory * dir = (TDirectory*)file[file_counter]->Get("rootTupleTree");
    dir->GetObject("tree", tree[file_counter]);
    std::cout<<"Creating Histo"<<std::endl;
    histoD[file_counter] = new TH1D(Form("histo_%d", file_counter), Form("histo_%d", file_counter),bins,min,max);
    //histoD[file_counter]->Sumw2();
    std::cout<<"Projecting"<<std::endl;
    //tree[file_counter]->Project(Form("histo_%d", file_counter),Form("%s", variable.c_str()),cut.c_str());
    tree[file_counter]->Draw(Form("%s>>histo_%d", variable.c_str(), file_counter),Form("photonCorrection*(%s)",cut.c_str()));
    //std::cout<<"Scaling"<<std::endl;
    //histoD[file_counter]->Scale(scaleFactor);//rescale for the crossSection and the number of events

    if(strcmp(bkg_nameTMP.c_str(), "data")==0){
      std::cout<<"DATA!"<<std::endl;
      if(isData==0) {
        allBkgHisto= new TH1D("allBkgHisto","allBkgHisto", bins,min,max);
        //allBkgHisto->Sumw2();
        dataHisto= new TH1D("dataHisto","dataHisto", bins,min,max);
        //dataHisto->Sumw2();
        std::cout<<"Data counts: "<<counts<<std::endl;
        //dataHisto->GetXaxis()->SetTitle(xTitle.c_str());
        dataHisto->GetYaxis()->SetTitle(yTitle.c_str());
        //dataHisto->GetYaxis()->SetRangeUser(-1000.,dataHisto->GetMaximum()*1.3);
        dataHisto->SetMarkerStyle(8);
        dataHisto->SetMarkerColor(1);
        dataHisto->SetLineColor(1);
        leg->AddEntry(dataHisto, "Data", "pe");
      }
      isData=1;
      dataHisto->Add(histoD[file_counter]);
    }else if(strncmp(bkg_nameTMP.c_str(), "signal", 6)==0){
      signalHisto[sgn_counter]= new TH1D(Form("signalHisto_%d",sgn_counter),Form("signalHisto_%d",sgn_counter), bins,min,max);
      //signalHisto[sgn_counter]->Sumw2();
      leg->AddEntry(signalHisto[sgn_counter],bkg_nameTMP.c_str()/*"Signal"*/, "l");
      if(doScale==0){
        scaleFactor=1; 
        lumi=1; 
        counts=1;
      }
      histoD[file_counter]->Scale(scaleFactor*lumi/counts);
      signalHisto[sgn_counter]->Add(histoD[file_counter]);
      ++sgn_counter;

    }else{
      if(doScale==0) {
        scaleFactor=1; 
        lumi=1; 
        counts=1;
      }
      histoD[file_counter]->Scale(scaleFactor*lumi/counts);
      std::cout<<" " <<scaleFactor<<" "<<lumi<<" "<<counts<<std::endl;
      std::cout<<"Total scale factor: "<<scaleFactor*lumi/counts<<std::endl;
      if(strcmp(bkg_name.c_str(),bkg_nameTMP.c_str())==0){//if the bkg is always the same
        std::cout<<"Background: "<<bkg_name.c_str()<<std::endl;
        histoBkg[bkg_counter]->Add(histoD[file_counter]);
        for(int m=0; m<histoD[file_counter]->GetNbinsX();++m){
          std::cout<<"BIN: "<<m<<" "<<histoD[file_counter]->GetBinContent(m)<<std::endl;
          //std::cout<<"BIN1: "<<m<<" "<<histoBkg[file_counter]->GetBinContent(m)<<std::endl;
        }
      }//end for the same bkg
      else{
        bkg_name=bkg_nameTMP;
        //if(bkg_counter!=-1) bkgStack->Add(histoBkg[bkg_counter]);
        ++bkg_counter;
        histoBkg[bkg_counter]= new TH1D(Form("histo_%s", bkg_name.c_str()), Form("histo_%s", bkg_name.c_str()),bins,min,max);
        //histoBkg[bkg_counter]->Sumw2();
        histoBkg[bkg_counter]->SetStats(0);
        //histoBkg[bkg_counter]->SetLineColor(COLOR[bkg_counter]);
        histoBkg[bkg_counter]->SetFillColor(COLOR[bkg_counter]);
        histoBkg[bkg_counter]->SetLineWidth(1);
        histoBkg[bkg_counter]->SetLineColor(1);
        histoBkg[bkg_counter]->Add(histoD[file_counter]);
        leg->AddEntry(histoBkg[bkg_counter], bkg_nameTMP.c_str(), "f");
        std::cout<<"New Background: "<<bkg_name.c_str()<<std::endl;
        std::cout<<"Color: "<<COLOR[bkg_counter]<<std::endl;
      }

    }
    //file[file_counter]=TFile::Close("R");
    ++file_counter;
   
  }//end while over inputFile
  std::cout<<"counter del bkg "<<bkg_counter<<std::endl; 
  //TCanvas *c[10];
  std::string outFileName = "";
  if(outString.size()<100) {
    outFileName = outString;
  }else{
    outFileName = "long_variable";
  }

  TFile *outFile = new TFile(Form("%s.root",outFileName.c_str()), "RECREATE");
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
    if(doLogY==1) {
      c_histo->SetLogy();
    }
    TPad *pad2=new TPad("pad2", "bottom pad", 0.,0.05,1,0.25);
    pad2->SetBottomMargin(0.2);
    pad2->SetTopMargin(0);
    //pad2->SetNdivisions(10);
    pad2->Draw();
    TPad *pad1=new TPad("pad1", "top pad", 0,0.3,1,1);
    pad1->SetBottomMargin(0);
    //if(doLogY==1) pad1->SetLogy();
    pad1->Draw();
    pad1->cd();
    if(doLogY==1) gPad-> SetLogy();
    float H = pad1->GetWh();
    float W = pad1->GetWw();
    float l = pad1->GetLeftMargin();
    float t = pad1->GetTopMargin();
    float r = pad1->GetRightMargin();
    float b = pad1->GetBottomMargin();
    TString lumiText = Form("#it{L}=%.01f fb^{-1} (2016) (13 TeV)", lumi/1000);
    TString cmsLogo = "#bf{CMS}  #it{Preliminary}";
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
    if(doLogY==1){
      dataHisto->GetYaxis()->SetRangeUser(1.,dataHisto->GetMaximum()*10);
    }else{
      dataHisto->GetYaxis()->SetRangeUser(0.,dataHisto->GetMaximum()*1.3);
    }

    //TPaveText *ptText = new TPaveText(.5,.5,.6,.6);
    //ptText->AddText("stocazzo");
    //ptText->Draw();

    dataHisto->Draw("EP");
    bkgStack->Draw("histoSAME");
    dataHisto->Draw("EPSAME");
    //gPad->Modified();
    //ptText->Draw();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);    
    latex.SetTextFont(42);
    latex.SetTextAlign(31); 
    float lumiTextOffset = 0.2;
    float lumiTextSize = 0.04;
    latex.SetTextSize(lumiTextSize);    
    latex.DrawLatex(.9,.92,lumiText);//1-r,1-t+lumiTextOffset*t,lumiText);
    latex.DrawLatex(.25,.92,cmsLogo);
    latex.Draw();
    gPad->Modified();
    gPad->Update();
    pad1->Update();
    std::cout<<"Redrawing data"<<std::endl;
    TGraphErrors *gr = new TGraphErrors(0);
    TH1D *fake_plotHisto = new TH1D("fake_plotHisto","fake_plotHisto", bins,min,max);
    integralData=dataHisto->Integral();
    integralBKG=allBkgHisto->Integral();
    double error, ratio;
    frame("Graph values");
    for(int w=1; w<bins; ++w){
      if((dataHisto->GetBinContent(w)==0&&allBkgHisto->GetBinContent(w)==0) || (allBkgHisto->GetBinContent(w)==0)) continue;
      gr->SetPoint(w, dataHisto->GetBinCenter(w),(dataHisto->GetBinContent(w))/(allBkgHisto->GetBinContent(w)));
      ratio= (dataHisto->GetBinContent(w))/(allBkgHisto->GetBinContent(w));
      error= (dataHisto->GetBinContent(w)*sqrt(allBkgHisto->GetBinContent(w)) + allBkgHisto->GetBinContent(w)*sqrt(dataHisto->GetBinContent(w)))/(allBkgHisto->GetBinContent(w)*allBkgHisto->GetBinContent(w));
      std::cout<<"XAXIS: "<<dataHisto->GetBinCenter(w)<<" VALUE: "<<ratio<<" ERROR: "<<error<<std::endl;
      gr->SetPointError(w, dataHisto->GetBinWidth(w)/2,error);

    }
    gr->GetHistogram()->SetMaximum(2.);//1.5);
    gr->GetHistogram()->SetMinimum(0.1);
    gr->GetXaxis()->SetLimits(min, max);
    gr->GetYaxis()->SetTitleOffset(.2);
    gr->GetYaxis()->SetTitleFont(42);
    gr->GetYaxis()->SetTitleSize(.15);
    gr->GetXaxis()->SetTitleFont(42);
    gr->GetXaxis()->SetTitleSize(.15);
    gr->GetXaxis()->SetTitle(xTitle.c_str());
    gr->GetYaxis()->SetTitle("DATA/MC");
    //TGraphAsymmErrors *gr= new TGraphAsymmErrors();
    //gr->Divide(dataHisto,allBkgHisto, "pois");
    pad2->cd();
    gStyle->SetTextSize(14);
    gROOT->ForceStyle();
    pad2->SetGrid();
    fake_plotHisto->GetYaxis()->SetRangeUser(0.,2.);
    //fake_plotHisto->SetNdivisions(400);
    fake_plotHisto->SetTitle("");
    fake_plotHisto->SetStats(0);
    fake_plotHisto->SetLineColor(0);
   // fake_plotHisto->Draw("hist");
    gr->Draw("ZAP");//SAME");
    gr->GetXaxis()->SetLabelFont(43);
    gr->GetXaxis()->SetLabelSize(14);
    gr->GetYaxis()->SetLabelFont(43);
    gr->GetYaxis()->SetLabelSize(14);
    gr->GetYaxis()->SetNdivisions(5);
    gr->SetMarkerStyle(8);
    gr->SetTitle("");
    //Chebyshev * cheb = new Chebyshev(4,min,max);
    //TF1 * f1 = new TF1("f1",cheb,min,1230,5,"Chebyshev");
    //TFile* gJet_correction = TFile::Open("/cmshome/gellisim/CMSSW_VGamma/src/DiBosonAnalysis/uncertainties_EWK_24bins.root", "r");
    //if (!gJet_correction ){
    //   frame("GammaJets corrections file not found! ERROR!");
    //   exit(-1);
    //}
    //TH1F* gCorrNominal = (TH1F*)gJet_correction->Get("GJets_1j_NLO/nominal_G");
    //TH1F* gCorrInv = (TH1F*)gJet_correction->Get("GJets_LO/inv_pt_G");
    //TGraphAsymmErrors *grCorr = new TGraphAsymmErrors(0);
    //grCorr->Divide(gCorrNominal,gCorrInv, "pois");
    //for (int i = 0; i <=5; ++i) f1->SetParameter(i,1);
    //grCorr->Fit(f1, "R");
    //if(strcmp(variable.c_str(),"photon_pt")==0){
    //  //f1->Draw("SAME");
    //  //grCorr->SetMarkerColor(2);
    //  //grCorr->Draw("PSAME");
    //}
    std::cout<<"Drawing chebyshev"<<std::endl;
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
    if (superimpose==1){
      for(int i=0; i<=bkg_counter; ++i){
        histoBkg[i]->SetFillColor(0);
        histoBkg[i]->SetLineColor(COLOR[i]);
        histoBkg[i]->SetLineWidth(2.5);
        histoBkg[i]->Scale(1/(histoBkg[i]->Integral()));
        if(i==0) {
          histoBkg[i]->GetXaxis()->SetTitle(xTitle.c_str());
          histoBkg[i]->GetYaxis()->SetTitle(yTitle.c_str());
          histoBkg[i]->Draw("hist");
        }
        if(i>0) histoBkg[i]->Draw("histSAME");
      }
    }else{
      bkgStack->GetXaxis()->SetTitle(xTitle.c_str());
      bkgStack->GetYaxis()->SetTitle(yTitle.c_str());
      bkgStack->Draw("h");
    }
    TLatex latex;
    TString lumiText = Form("#it{L}=%.01f fb^{-1} (2016) (13 TeV)", lumi/1000);//"#bf{CMS Preliminary} #it{L}=12.9 fb^{-1}";
    TString cmsLogo = "#bf{CMS}  #it{Preliminary}";
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    float lumiTextOffset = 0.2;
    float lumiTextSize = 0.04; 
    latex.SetTextSize(lumiTextSize);
    latex.DrawLatex(.9,.92,lumiText);//1-r,1-t+lumiTextOffset*t,lumiText);
    latex.DrawLatex(.35/*.25*/,.92,cmsLogo);
    latex.Draw();
    gPad->Modified();
    gPad->Update();

    std::cout<<"Drawing BKG"<<std::endl;
  }
  
  for(int i=0; i<sgn_counter;++i){
    std::cout<<"Drawing signal n. "<<sgn_counter<<std::endl;
    std::cout<<"++++++++++++++++++++++ "<<(signalHisto[i]->Integral())<<" "<<integralData<<std::endl;
    //signalHisto[i]->Scale((signalHisto[i]->Integral())/integralData);
    std::cout<<"++++++++++++++++++++++ "<<(signalHisto[i]->Integral())<<std::endl;
    signalHisto[i]->SetStats(0);
    signalHisto[i]->SetLineColor(1);
    signalHisto[i]->SetLineWidth(2);
    signalHisto[i]->SetLineStyle(i+1);
    signalHisto[i]->Draw("histSAME");
    signalHisto[i]->Write();

  }
  leg->SetBorderSize(0);
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
  if(doLogY==1){
    std::cout<<"LOGY scale"<<std::endl;
    if(isData==0){
      c_histo->SetLogy();
    }
  }
  gPad->RedrawAxis();
  c_histo->SaveAs(Form("%s.png", outFileName.c_str()));
  c_histo->SaveAs(Form("%s.pdf", outFileName.c_str()));
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
  leg->Write();
  std::cout<<"Closing output file"<<std::endl;
  outFile->Close();

  std::cout<<""<<std::endl;
  std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << timer.RealTime() << " s" <<std::endl;
  std::cout<<""<<std::endl;
}
