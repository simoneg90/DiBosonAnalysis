//#include "utility.cc"

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include "fit_functions.cc"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include <iomanip>
#include <fstream>
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


#include "TStopwatch.h"

#define NUMBER 1000
#define maxIterations 200
using namespace RooFit;

int main(int argc, char* argv[]){

  TStopwatch time;
  time.Start(true);

  breakLine();
  std::cout<<"This program uses the name given in INPUT to localize the file and the histogram to analyze"<<std::endl;
  breakLine();

  std::cout<<"File: "<<argv[1]<<std::endl;
  std::cout<<"Histogram: "<<argv[2]<<std::endl;
  TFile *file= TFile::Open(argv[1]);
  TH1D *histo= ((TH1D *)(file->Get(argv[2])));
  RooRealVar x("x","x",0,99999);
  RooDataHist dh("dh","dh",x,Import(*histo));

  std::cout << "\n ************** ITERATIVE DCB FIT *************\n " << std::endl;


  RooRealVar mean("mean", "Double CB mean", 501., 400., 800.);
  RooRealVar sigma("sigma", "Double CB Width", 11.,10.,100.);
  RooRealVar dCBCutL("dCBCutL", "Double CB Cut left", 1.,.1,50.);
  RooRealVar dCBCutR("dCBCutR", "Double CB Cut right", 1.,.1,50.);
  RooRealVar dCBPowerL("dCBPowerL", "Double CB Power left", 0.,-.001,.001);
  RooRealVar dCBPowerR("dCBPowerR", "Double CB Power right", 1.,.01,50.);
  RooDCBShape dcb("dcb", "double crystal ball", x, mean, sigma, dCBCutL, dCBCutR, dCBPowerL, dCBPowerR);
  //TF1* f1 = new TF1( "gaussian", "gaus",  histo->GetMean() - sigma * histo->GetRMS(),  histo->GetMean() + sigma * histo->GetRMS() );

  dcb.fitTo(dh, /*RooFit::Strategy(2),*/ SumW2Error(kTRUE));
  dcb.fitTo(dh, /*RooFit::Strategy(2),*/ SumW2Error(kTRUE));
  dcb.fitTo(dh, /*RooFit::Strategy(2),*/ SumW2Error(kTRUE));
  TCanvas *c = new TCanvas("c", "Double Crystal Ball Fit",1);
  RooPlot* frame1 = x.frame(Bins(histo->GetNbinsX()),Title("double crystal ball")) ;
  dh.plotOn(frame1, LineColor(1), LineWidth(2));
  dcb.plotOn(frame1, LineColor(2), LineWidth(2));//,Normalization(1.0,RooAbsReal::RelativeExpected) );

  std::cout<<"chi square! "<<frame1->chiSquare()<<std::endl;
  frame1->Draw();
  c->SaveAs("dcb_test.png");

  std::cout<<""<<std::endl;
  std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << time.RealTime() << " s" <<std::endl;
  std::cout<<""<<std::endl;

  //RooRealVar x1 ("x1");
  //RooAbsPdf* pdf = w->pdf("pdf");
  Double_t x0 = atof(argv[3]);
  Double_t x_max = atof(argv[4]);

  x.setRange("my_range", x0, x_max); // create range to integrate over
  RooAbsReal* i = dcb.createIntegral(x, RooFit::NormSet(x), RooFit::Range("my_range"));
  // alternatively: RooAbsReal* i = pdf->createIntegral(*x, *x, "my_range");
  // // in my example, I am also using RooArgSet(*x) as a normalisation set; you might want to use a different one; be sure to check out the documentation for RooAbsReal::createIntegral
  std::cout << "Integral value: " << i->getVal() << std::endl;

  std::ofstream outList;
  outList.open("integrallist.txt", std::fstream::app);
  outList<<i->getVal()<<std::endl;

  return 0;

}
