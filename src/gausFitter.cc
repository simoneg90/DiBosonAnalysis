//#include "utility.cc"

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
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

#include "TStopwatch.h"

#define NUMBER 1000
#define maxIterations 200

int main(int argc, char* argv[]){

    TStopwatch time;
    time.Start(true);

    std::cout<<"This program uses the name given in INPUT to localize the file and the histogram to analyze"<<std::endl;

    std::cout<<"File: "<<argv[1]<<std::endl;
    std::cout<<"Histogram: "<<argv[2]<<std::endl;
    TFile *file= TFile::Open(argv[1]);
    TH1D *histo= ((TH1D *)(file->Get(argv[2])));
    std::cout << "\n ************** ITERATIVE GAUSSIAN FIT *************\n " << std::endl;

    double rangeMin = atof(argv[3]);
    double rangeMax = atof(argv[4]);
    double sigma;
    double integral = histo->Integral();
    double thisNorm = -999.;
    double thisMean = -999.;
    double thisSigma = -999.;
    int iter=1;
    sigma=2.;

    TF1* f1 = new TF1( "gaussian", "gaus", rangeMin, rangeMax);//histo->GetMean() - sigma * histo->GetRMS(),  histo->GetMean() + sigma * histo->GetRMS() );
    std::cout<<"Range: "<<histo->GetMean() - sigma * histo->GetRMS()<<" - "<<histo->GetMean() + sigma * histo->GetRMS()<<std::endl;
    f1->SetParameters(histo->GetMaximum(), histo->GetMean(), histo->GetRMS());
    f1->SetParLimits( 2, 0.2*histo->GetRMS(), 3*histo->GetRMS() );
    //histo->Scale(1/histo->Integral());


    histo->Fit( f1, "R+" );

    while ( (fabs((f1->GetParameter(1) - thisMean ) / f1->GetParameter(1) ) > 0.001) && (fabs((f1->GetParameter(2) - thisSigma ) / f1->GetParameter(2) ) > 0.001) && (iter < maxIterations) ) {
       thisNorm = f1->GetParameter(0);
       thisMean = f1->GetParameter(1);
       thisSigma = f1->GetParameter(2);
       f1 = new TF1( "gaussian", "gaus", thisMean - sigma*thisSigma, thisMean + sigma*thisSigma );
       f1->SetParameters(thisNorm, thisMean, thisSigma);
       histo->Fit( f1, "R+" );
       ++iter;

   }

   std::cout << "\n ************** EXIT AFTER " << iter << " ITERATIONS *************\n " << std::endl;

   std::cout<<f1->Integral(rangeMin,rangeMax)<<"++++++++++++"<<std::endl;

   std::ofstream outList;
   outList.open("integrallistSignal.txt", std::fstream::app);
   outList<<histo->Integral(/*histo->FindBin(rangeMin), histo->FindBin(rangeMax)*/)<<std::endl;//f1->Integral(rangeMin, rangeMax)<<std::endl;
   std::cout<<"---- "<<histo->FindBin(rangeMin)<<std::endl;

   std::cout<<""<<std::endl;
   std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << time.RealTime() << " s" <<std::endl;
   std::cout<<""<<std::endl;

   TCanvas *c1 = new TCanvas ("c1","c1",1);

   histo->Draw("h");
   f1->Draw("SAME");

   c1->SaveAs(Form("gaus_%s.png",argv[5]));



   return 0;

}
