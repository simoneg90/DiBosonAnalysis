//#include "utility.cc"

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include<math.h>
#include "fit_functions.cc"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TPaveLabel.h"
#include "TText.h"
#include "TLegend.h"
#include "TAxis.h"
#include <iomanip>
#include <fstream>
#include "RooGlobalFunc.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooSimPdfBuilder.h"
#include "RooGenericPdf.h"
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
#include "RooChi2Var.h"
#include "TMath.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooRealProxy.h"


#include "TStopwatch.h"

#define NUMBER 1000
#define maxIterations 200
using namespace RooFit;

int main(int argc, char* argv[]){

  //Usage: ./dijetFunctionFitter tau21Quantiles/ak08_photon_mass_All_tau090.root allBkgHisto 3000 13. range_min range_max
  TH1::SetDefaultSumw2();
  TStopwatch time;
  time.Start(true);

  breakLine();
  std::cout<<"This program uses the name given in INPUT to localize the file and the histogram to analyze"<<std::endl;
  breakLine();

  std::cout<<"File: "<<argv[1]<<std::endl;
  std::cout<<"Histogram: "<<argv[2]<<std::endl;
  TFile *file= TFile::Open(argv[1]);
  TH1D *histo= ((TH1D *)(file->Get(argv[2])));
  TH1D *histo2 = (TH1D *)histo->Clone();

  for(int i = 0; i < histo->GetNbinsX(); ++i){
//        if(histo->GetBinContent(i) <1 ){
//           histo->SetBinError(i,0);
//           std::cout<<"<1"<<std::endl;
//        }
 //       std::cout<<"----- "<<histo->GetBinContent(i)<<std::endl;
        if(histo->GetBinContent(i) == 0){
   //       std::cout<<"negativo!"<<std::endl;
           histo->SetBinError(i, histo->Integral());//9999999999);
        }
  }
  

  std::cout << "\n ************** ITERATIVE DCB FIT *************\n " << std::endl;

  double max_range = atof(argv[3]);
  double p1_value = atof(argv[4]);
  RooRealVar y("y","y",600,max_range);//4500);
  RooDataHist dh("dh","dh",y,Import(*histo));
  RooRealVar bg_p0("bg_p0", "bg_p0", 3000, 100, 100000.);
  RooRealVar bg_p1("bg_p1", "bg_p1", -13., -16.,-12.);//-13., -1000., 1000.);
  //bg_p1.setVal(p1_value);
  RooRealVar bg_p2("bg_p2", "bg_p2", -1.4, -2.,-1.);//-1000, 1000.);
  RooRealVar bg_p3("bg_p3", "bg_p3", .2,0.01,10.);
  //RooGenericPdf bg_pure = RooGenericPdf((std::string("bg_pure_")+blah).c_str(),"(pow(@0/13000,@1+@2*log(@0/13000)))",RooArgList(x,bg_p1,bg_p2));
  //RooGenericPdf bg = RooGenericPdf("bg","(pow((1-(@0/13000)),@3)*pow(@0/13000,@1+@2*log(@0/13000)))",RooArgList(x,bg_p1,bg_p2,bg_p3));
  RooGenericPdf bg = RooGenericPdf("bg","bg_p0*(pow((1-(y/13000)),bg_p3)*pow(y/13000,(bg_p1+bg_p2*log(y/13000))))",RooArgList(y,bg_p1,bg_p2,bg_p3,bg_p0));
  //RooGenericPdf bg = RooGenericPdf("bg","@3*(pow(1-@0/13000,@4))*pow(@0/13000,@1+@2*log(@0/13000))",RooArgList(x,bg_p1,bg_p2,bg_p0,bg_p3));
  //RooGenericPdf bg = RooGenericPdf("bg","@1*exp(1-@0)", RooArgList(x, bg_p3));
  //TF1* f1 = new TF1( "gaussian", "gaus",  histo->GetMean() - sigma * histo->GetRMS(),  histo->GetMean() + sigma * histo->GetRMS() );

  //double x;
  TF1* f1 = new TF1( "f1", "(TMath::Power((1-(x/13000)),[2])*TMath::Power(x/13000,[0]+[1]*TMath::Log(x/13000)))", 600,2800);
  f1->SetParameters(-13,-1.4, .2);
  histo->Fit("f1","R");
  f1->SetParameters(f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2));
  histo->Fit("f1","R");

  std::cout<<"++++++++ "<<f1->Integral(600,1000)<<std::endl;
  std::cout<<"++++++++ "<<f1->Integral(660,820)<<std::endl;

  //bg.fitTo(dh, /*RooFit::Strategy(2),*/ SumW2Error(kTRUE));

  ////dh.add(x, 1.,1.);
  //RooChi2Var chi2("chi2","chi2",bg,dh,DataError(RooAbsData::SumW2)) ;
  //RooMinuit m(chi2) ;
  //m.migrad() ;
  //m.hesse() ;

  //// Plot chi^2 fit result on frame as well
  //RooFitResult* r_chi2_wgt = m.save() ;
//  bg.plotOn(frame,LineStyle(kDashed),LineColor(kRed)) ;

////  bg.chi2FitTo(dh, RooFit::Strategy(2), SumW2Error(kTRUE));
//  bg.chi2FitTo(dh, RooFit::Strategy(2), SumW2Error(kTRUE));
//  bg.fitTo(dh, /*RooFit::Strategy(2),*/ SumW2Error(kTRUE));
  TCanvas *c = new TCanvas("c", "Double Crystal Ball Fit",1);
  RooPlot* frame1 = y.frame(Bins(histo->GetNbinsX()),Title("dijet"));//, Range(500,5000)) ;
  dh.plotOn(frame1, DataError(RooAbsData::SumW2),LineColor(0), MarkerColor(0));
  bg.plotOn(frame1, LineColor(2), LineWidth(2));//,Normalization(1.0,RooAbsReal::RelativeExpected) );
  bg.paramOn(frame1);
  //x.setRange("my_range", 600, 5000);
  //bg.plotOn(frame1, LineColor(3), LineWidth(2));
  //TPad *pad2=new TPad("pad2", "bottom pad", 0.,0.05,1,0.25);
  //pad2->SetBottomMargin(0.2);
  //pad2->SetTopMargin(0);
  //pad2->Draw();
  //TPad *pad1=new TPad("pad1", "top pad", 0,0.3,1,1);
  //pad1->SetBottomMargin(0);
  //pad1->Draw();
  //pad1->cd();
  c->SetLogy();

  std::cout<<"chi square! "<< f1->GetChisquare()<<std::endl;//frame1->chiSquare()<<std::endl;
  std::cout<<"NDF: "<<f1->GetNDF()<<std::endl;

  TPaveText *pt = new TPaveText(.05,.1,.95,.8);//.7,.75,.9,.85);
  //TText *t1 = pt->AddText(Form("#Chi^{2} = %f", frame1->chiSquare()),"brNDC");
  TPaveLabel *t1 = new TPaveLabel(0.7,0.6,0.9,0.68, Form("#chi^{2} = %f", f1->GetChisquare()/f1->GetNDF()/*frame1->chiSquare()*/),"brNDC");//"bg","dh")),"brNDC"); 
  frame1->addObject(t1);

  //TPaveLabel *t1 = new TPaveLabel(0.7,0.6,0.9,0.68, Form("#chi^{2} = %f", frame1->chiSquare())); 
  pt->Draw();
  histo->SetStats(0);
  histo->SetLineColor(1);
  histo->SetMarkerColor(1);
  histo->SetMarkerStyle(8);
  histo->GetXaxis()->SetRangeUser(600,5000);
  histo->Draw();
  //frame1->Draw("SAME");
  f1->Draw("SAME");
  f1->SetRange(600,4000);
  f1->SetLineStyle(2);
  f1->SetLineColor(3);
  f1->Draw("SAME");
  histo->Draw("E1same");
  t1->Draw("SAME");

  pt->Draw("SAME");

  //pad2->cd();
  TH1D *histoRatio = (TH1D *) histo->Clone();
  TGraphErrors *grRatio = new TGraphErrors(0);

  for(int w=1; w<histo->GetNbinsX(); ++w){
//    x.setRange(Form("my_range%d",w), histo->GetBinCenter(w)-.5*histo->GetBinWidth(w), histo->GetBinCenter(w)+.5*histo->GetBinWidth(w)); 
//    RooAbsReal* fun = bg.createIntegral(x, RooFit::NormSet(x), RooFit::Range(Form("my_range%d",w)));
//    std::cout<<"---------------------- "<<fun->getVal()*bg.getNorm()<<std::endl;
    //histoRatio->SetBinContent(w/*, histo->GetBinCenter(w)*/,(fun->getVal()*histo->Integral()-histo->GetBinContent(w)/histo->GetBinError(w)));//(histo->GetBinContent(w)-fun->getVal())/sqrt(fun->getVal()));
  //  grRatio->SetPoint(w,histo->GetBinCenter(w),(histo->GetBinContent(w)-bg.Eval(w)*histo->Integral()/sqrt(bg.Eval(w)*histo->Integral())));
     //ratio= (dataHisto->GetBinContent(w))/(allBkgHisto->GetBinContent(w));
     //error= (dataHisto->GetBinContent(w)*sqrt(allBkgHisto->GetBinContent(w)) + allBkgHisto->GetBinContent(w)*sqrt(dataHisto->GetBinContent(w)))/(allBkgHisto->GetBinContent(w)*allBkgHisto->GetBinContent(w));
     //std::cout<<"XAXIS: "<<dataHisto->GetBinCenter(w)<<" VALUE: "<<ratio<<" ERROR: "<<error<<std::endl;
     //gr->SetPointError(w, dataHisto->GetBinWidth(w)/2,error);

  }
  //histoRatio->Draw("histo");

  //RooPlot* frame2 = y.frame(Bins(histo->GetNbinsX()),Title(""));
  //frame2->addObject(frame1->pullHist()) ; 
  //frame2->SetMinimum(-5.);
  //frame2->SetMaximum(5.);
  //frame2->getAttText()->SetTextSize(.15) ;
  //frame2->Draw();
  //grRatio->SetFillColor(2);
  //grRatio->Draw("ac");

  c->SaveAs("dijet_test.png");

  std::cout<<""<<std::endl;
  std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << time.RealTime() << " s" <<std::endl;
  std::cout<<""<<std::endl;

  //RooRealVar x1 ("x1");
  //RooAbsPdf* pdf = w->pdf("pdf");
  Double_t x0 = atof(argv[5]);
  Double_t x_max = atof(argv[6]);

  y.setRange("my_range", x0, x_max); // create range to integrate over
  RooAbsReal* i = bg.createIntegral(y, RooFit::NormSet(y), RooFit::Range("my_range"));
  // alternatively: RooAbsReal* i = pdf->createIntegral(*x, *x, "my_range");
  // // in my example, I am also using RooArgSet(*x) as a normalisation set; you might want to use a different one; be sure to check out the documentation for RooAbsReal::createIntegral
  std::cout << "Integral value: " << i->getVal() << std::endl;
  std::cout<< "bg_p0 value " <<bg_p0.getValV()<<std::endl;
  std::cout<< "histo integral: "<<histo->Integral(11,19)<<std::endl;
  std::cout<< "histo total integral: "<<histo->Integral(0,150)<<std::endl;
  std::cout<< "n bins " <<histo->GetNbinsX()<<std::endl;
  std::cout<< "entries histo: "<<histo->GetEntries()<<std::endl;
  std::cout<< "sumOfWeights: "<<histo->GetSumOfWeights()<<std::endl;
  //std::cout<<" bg norm "<<bg.getNormIntegral()<<std::endl;

  std::ofstream outList;
  outList.open("integrallist.txt", std::fstream::app);
  outList<<f1->Integral(x0, x_max)<<std::endl;
  return 0;

  TFile *file1 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau020.root");
  TFile *file2 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau025.root");
  TFile *file3 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau030.root");
  TFile *file4 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau035.root");
  TFile *file5 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau040.root");
  TFile *file6 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau045.root");
  TFile *file7 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau050.root");
  TFile *file8 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau055.root");
  TFile *file9 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau060.root");
  TFile *file10 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau065.root");
  TFile *file11 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau070.root");
  TFile *file12 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau075.root");
  TFile *file13 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau080.root");
  TFile *file14 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau085.root");
  TFile *file15 = TFile::Open ("tau21Quantiles/ak08_photon_mass_All_tau090.root");

  TH1D *h1 = (TH1D *)file1->Get("allBkgHisto");
  TH1D *h2 = (TH1D *)file2->Get("allBkgHisto");
  TH1D *h3 = (TH1D *)file3->Get("allBkgHisto");
  TH1D *h4 = (TH1D *)file4->Get("allBkgHisto");
  TH1D *h5 = (TH1D *)file5->Get("allBkgHisto");
  TH1D *h6 = (TH1D *)file6->Get("allBkgHisto");
  TH1D *h7 = (TH1D *)file7->Get("allBkgHisto");
  TH1D *h8 = (TH1D *)file8->Get("allBkgHisto");
  TH1D *h9 = (TH1D *)file9->Get("allBkgHisto");
  TH1D *h10 = (TH1D *)file10->Get("allBkgHisto");
  TH1D *h11 = (TH1D *)file11->Get("allBkgHisto");
  TH1D *h12 = (TH1D *)file12->Get("allBkgHisto");
  TH1D *h13 = (TH1D *)file13->Get("allBkgHisto");
  TH1D *h14 = (TH1D *)file14->Get("allBkgHisto");
  TH1D *h15 = (TH1D *)file15->Get("allBkgHisto");

  //h1->SetFillColor(1);
  //h2->SetFillColor(2);
  //h3->SetFillColor(3);
  //h4->SetFillColor(4);
  //h5->SetFillColor(5);
  //h6->SetFillColor(6);
  //h7->SetFillColor(7);
  //h8->SetFillColor(8);
  //h9->SetFillColor(9);
  //h10->SetFillColor(10);
  //h11->SetFillColor(11);
  //h12->SetFillColor(12);
  //h13->SetFillColor(13);
  //h14->SetFillColor(14);
  //h15->SetFillColor(15);
  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(3);
  h4->SetLineColor(4);
  h5->SetLineColor(5);
  h6->SetLineColor(6);
  h7->SetLineColor(7);
  h8->SetLineColor(8);
  h9->SetLineColor(9);
  h10->SetLineColor(10);
  h11->SetLineColor(11);
  h12->SetLineColor(12);
  h13->SetLineColor(13);
  h14->SetLineColor(14);
  h15->SetLineColor(15);
  h1->Scale(1/h1->Integral());
  h2->Scale(1/h2->Integral());
  h3->Scale(1/h3->Integral());
  h4->Scale(1/h4->Integral());
  h5->Scale(1/h5->Integral());
  h6->Scale(1/h6->Integral());
  h7->Scale(1/h7->Integral());
  h8->Scale(1/h8->Integral());
  h9->Scale(1/h9->Integral());
  h10->Scale(1/h10->Integral());
  h11->Scale(1/h11->Integral());
  h12->Scale(1/h12->Integral());
  h13->Scale(1/h13->Integral());
  h14->Scale(1/h14->Integral());
  h15->Scale(1/h15->Integral());
  TCanvas *c0 = new TCanvas("c0","c0",1);
  c0->SetLogy();

  TGraphAsymmErrors *gr1 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr2 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr3 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr4 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr5 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr6 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr7 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr8 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr9 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr10 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr11 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr12 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr13 = new TGraphAsymmErrors();
  TGraphAsymmErrors *gr14 = new TGraphAsymmErrors();

  for(int i=0; i<15; ++i){

    if(h1->GetBinContent(i)>0) gr1->SetPoint(i, i*.05+.2,h2->GetBinContent(i)/h1->GetBinContent(i));
    if(h2->GetBinContent(i)>0) gr2->SetPoint(i, i*.05+.2,h3->GetBinContent(i)/h2->GetBinContent(i));
    if(h3->GetBinContent(i)>0) gr3->SetPoint(i, i*.05+.2,h4->GetBinContent(i)/h3->GetBinContent(i));
    if(h4->GetBinContent(i)>0) gr4->SetPoint(i, i*.05+.2,h5->GetBinContent(i)/h4->GetBinContent(i));
    if(h5->GetBinContent(i)>0) gr5->SetPoint(i, i*.05+.2,h6->GetBinContent(i)/h5->GetBinContent(i));
    if(h6->GetBinContent(i)>0) gr6->SetPoint(i, i*.05+.2,h7->GetBinContent(i)/h6->GetBinContent(i));
    if(h7->GetBinContent(i)>0) gr6->SetPoint(i, i*.05+.2,h8->GetBinContent(i)/h7->GetBinContent(i));
    if(h8->GetBinContent(i)>0) gr7->SetPoint(i, i*.05+.2,h9->GetBinContent(i)/h8->GetBinContent(i));
    if(h9->GetBinContent(i)>0) gr8->SetPoint(i, i*.05+.2,h10->GetBinContent(i)/h9->GetBinContent(i));
    if(h10->GetBinContent(i)>0) gr9->SetPoint(i, i*.05+.2,h11->GetBinContent(i)/h10->GetBinContent(i));
    if(h11->GetBinContent(i)>0) gr10->SetPoint(i, i*.05+.2,h12->GetBinContent(i)/h11->GetBinContent(i));
    if(h12->GetBinContent(i)>0) gr11->SetPoint(i, i*.05+.2,h13->GetBinContent(i)/h12->GetBinContent(i));
    if(h13->GetBinContent(i)>0) gr12->SetPoint(i, i*.05+.2,h14->GetBinContent(i)/h13->GetBinContent(i));
    if(h14->GetBinContent(i)>0) gr13->SetPoint(i, i*.05+.2,h15->GetBinContent(i)/h14->GetBinContent(i));
    std::cout<<"i: "<<i*.05+.2<< " "<<h2->GetBinContent(i)/h1->GetBinContent(i)<<std::endl;

  }

//  gr1->Divide(h2,h1,"pois");
//  gr2->Divide(h3,h2,"pois");
//  gr3->Divide(h4,h3,"pois");
//  gr4->Divide(h5,h4,"pois");
//  gr5->Divide(h6,h5,"pois");
//  gr6->Divide(h7,h6,"pois");
//  gr7->Divide(h8,h7,"pois");
//  gr8->Divide(h9,h8,"pois");
//  gr9->Divide(h10,h9, "pois");
//  gr10->Divide(h11,h10,"pois");
//  gr11->Divide(h12,h11,"pois");
//  gr12->Divide(h13,h12,"pois");
//  gr13->Divide(h14,h13,"pois");
//  gr14->Divide(h15,h14,"pois");

  gr1->SetMarkerStyle(8);
  gr2->SetMarkerStyle(8);
  gr3->SetMarkerStyle(8);
  gr4->SetMarkerStyle(8);
  gr5->SetMarkerStyle(8);
  gr6->SetMarkerStyle(8);
  gr7->SetMarkerStyle(8);
  gr8->SetMarkerStyle(8);
  gr9->SetMarkerStyle(8);
  gr10->SetMarkerStyle(8);
  gr11->SetMarkerStyle(8);
  gr12->SetMarkerStyle(8);
  gr13->SetMarkerStyle(8);
  gr14->SetMarkerStyle(8);
  gr1->SetMarkerColor(1);
  gr2->SetMarkerColor(2);
  gr3->SetMarkerColor(3);
  gr4->SetMarkerColor(4);
  gr5->SetMarkerColor(5);
  gr6->SetMarkerColor(6);
  gr7->SetMarkerColor(7);
  gr8->SetMarkerColor(8);
  gr9->SetMarkerColor(9);
  gr10->SetMarkerColor(10);
  gr11->SetMarkerColor(11);
  gr12->SetMarkerColor(12);
  gr13->SetMarkerColor(13);
  gr14->SetMarkerColor(14);

//  h15->Draw("p");
//  h1->Draw("pSAME");
//  h13->Draw("pSAME");
//  h12->Draw("pSAME");
//  h11->Draw("pSAME");
//  h10->Draw("pSAME");
//  h9->Draw("pSAME");
//  h8->Draw("pSAME");
//  h7->Draw("pSAME");
//  h6->Draw("pSAME");
//  h5->Draw("pSAME");
//  h4->Draw("pSAME");
//  h3->Draw("pSAME");
//  h2->Draw("pSAME");
//  h14->Draw("pSAME");
//
//  c0->SaveAs("superimposeTau21.png");


  return 0;

}
