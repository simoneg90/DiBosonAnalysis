#include "histoPlotter.h"
#include "setTDRStyle.h"
#include "utility.h"

#include <fstream>

int main(int argc, char* argv[]){

  frame("Removed the scale part");
  frame("Histo Plotter env");

  THStack *hs= new THStack ("hs","");
  frame("1");
  TH1D *histoD[MAX_NUMBER];
  TH1F *histoF[MAX_NUMBER];
  TH2D *histo2D[MAX_NUMBER];
  TH2F *histo2F[MAX_NUMBER];
  TCanvas *c[MAX_NUMBER];
  TLegend *leg[MAX_NUMBER];
  TFile *file[MAX_NUMBER];
  std::string xAxis, yAxis, histoName;
  int thstack;
  TStopwatch time;
  time.Start(true);

  histoName="";
  xAxis="";
  yAxis="";
  histoName+=argv[2];
  xAxis += argv[3];
  yAxis += argv[4];
  thstack = 0;
  std::ifstream inputList;
  inputList.open(argv[1]);
  std::string filePath;
  std::string legendName[MAX_NUMBER];
  bool isData;
  int evts[MAX_NUMBER]; //fileCounter;
  float xsec[MAX_NUMBER];
  int maxFiles;
  //fileCounter=0;
  for(int j=0; j<MAX_NUMBER;++j){
    if(!inputList.good()) break;
      maxFiles = j;
      inputList>>filePath>>legendName[j];
      std::cout<<"File: "<<filePath.c_str()<<" legend: "<<legendName[j]<<std::endl;
      file[j]=TFile::Open(filePath.c_str());
    if(file[j]==0){
      breakLine();
      std::cout<<"ATTENTION! "<<argv[j+1]<<" doesn't exists!"<<std::endl;
      breakLine();
      exit(-1);
    }//end if (file exists)
    //++fileCounter;
  }//end loop over files=>from list

  //== Variables declaration ==

  int i=0;
  int histo_type=0; // =1: TH1D; =-1: TH1F; =2: TH2D; =-2 TH2F;
  int key_costraint=0;
  int histoD_counter,histoF_counter,histo2D_counter,histo2F_counter;
  int canvas_counter=0;
  double histo_max=0;
  int histo_index=0;
  std::string word;
  std::string title;
  histoD_counter=histoF_counter=histo2D_counter=histo2F_counter=0;

  gStyle->SetLegendBorderSize(0);
  TLegend *legSingle=new TLegend(0.67,0.67,0.87,0.87);
  TCanvas *cSingle = new TCanvas("cSingle","cSingle",200,10,600,400);
  setTDRStyle();
  float max=0;
  int max_index;

  frame("starting loop over files");
  int myVar; 
  for(int j=0;j<maxFiles;++j){//Loop over files

          frame("reading histo");         
          histoF[j]=((TH1F *)(file[j]->Get(histoName.c_str())));
          frame("read histo");
	  histoF[j]->SetStats(0);
          histoF[j]->Scale(1/(histoF[j]->Integral()));
          histoF[j]->GetXaxis()->SetTitle("x");
          if(j==0){
             histoF[j]->SetLineColor(kRed + 1);
             histoF[j]->SetMarkerColor(kRed + 1);
            // histoF[j]->SetFillColor(kRed + 1);
          }else if(j==1){
             histoF[j]->SetLineColor(kGreen + 1);
             histoF[j]->SetMarkerColor(kGreen + 1);
            // histoF[j]->SetFillColor(kGreen + 1);
          }else if(j==2){
            histoF[j]->SetLineColor(kBlue + 1);
            histoF[j]->SetMarkerColor(kBlue + 1);
            //histoF[j]->SetFillColor(kBlue + 1);
          }else if(j==3){
            histoF[j]->SetLineColor(kViolet);
            histoF[j]->SetMarkerColor(kViolet);
            //histoF[j]->SetFillColor(kYellow + 1);
          }
          //TH1F_config(histoF[j],"",xAxis.c_str(),yAxis.c_str(),j+2);
          //histoF[j]->SetFillColor(kRed +1);

          if(max<histoF[j]->GetMaximum()){
              max=histoF[j]->GetMaximum();
              histoF[j]->Draw("histo");
              max_index=j;
          }

          legSingle->AddEntry(histoF[j],legendName[j].c_str(),"l");
          hs->Add(histoF[j]);
          //if(thstack==1) hs->Add(histoF[j]);//added                     
	  //if(j==0 && thstack==0) histoF[j]->Draw("histo");
	  //if(j>0 && thstack==0) histoF[j]->Draw("histoSAME");
   
  }//end while inFile.good()
  for(int j=0; j<maxFiles;++j){
     
    //if(thstack==1) hs->Add(histoF[j]);//added                     
    if(j!=max_index)  histoF[j]->Draw("histoSAME");
  }

  frame("ended loop over files");

  //if(thstack==0)hs->Draw("histoNOSTACK");
  //if(thstack==0)hs->Draw("histo");
  legSingle->Draw("SAME");
  //hs->GetXaxis()->SetTitle(xAxis.c_str());
  //hs->GetYaxis()->SetTitle(yAxis.c_str());
  //cSingle->Modified();
  cSingle->SaveAs(Form("/afs/cern.ch/user/s/sgelli/www/DiBoson/%s.pdf",histoName.c_str()));
  cSingle->SaveAs(Form("/afs/cern.ch/user/s/sgelli/www/DiBoson/%s.png",histoName.c_str()));



  std::cout<<""<<std::endl;
  std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << time.RealTime() << " s" <<std::endl;
  std::cout<<""<<std::endl;
}
