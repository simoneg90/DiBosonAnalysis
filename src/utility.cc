#include "utility.h"

#include <stdlib.h>
#include <iostream>
#include <string.h>
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TCanvas.h"

//This function prints a +++ line
void breakLine(){

  std::cout<<""<<std::endl;
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<<""<<std::endl;

}

//This function prints a string inside a + box
void frame(std::string stringa){

  std::cout<<""<<std::endl;
  for (int i=0; i<(stringa.size()+6); ++i){
    std::cout<<"+";
  }
  std::cout<<""<<std::endl;
  std::cout<<"+  "<<stringa.c_str()<<"  +"<<std::endl;
  for(int i=0; i<(stringa.size()+6); ++i){
    std::cout<<"+";
  }
  std::cout<<"\n"<<std::endl;
  

}

//This function 'labels' a TH1D
void TH1D_config(TH1D *histo, const char* title, const char* title_x, const char* title_y, int color){

  std::string x="";
  x+=title_x;

  if(strcmp(x.c_str(),"Eta")==0)x="#eta";
  if(strcmp(x.c_str(),"Phi")==0)x="#Phi";
  if(strcmp(x.c_str(),"Pt")==0)x="Pt [GeV]";
  histo->SetTitle(title);
  histo->GetXaxis()->SetTitle(x.c_str());
  histo->GetYaxis()->SetTitle(title_y);
  if(color==1){
    histo->SetLineColor(kRed + 1);
    histo->SetMarkerColor(kRed + 1);
    histo->SetFillColor(kRed + 1);
  }else if(color==2){
    histo->SetLineColor(kGreen + 1);
    histo->SetMarkerColor(kGreen + 1);
    histo->SetFillColor(kGreen + 1);
  }

}

//This function 'labels' a TH1F
void TH1F_config(TH1F *histo, const char* title, const char* title_x, const char* title_y, int color){
 
  std::string x="";

  if(strcmp(x.c_str(),"Eta")==0)x="#eta";
  if(strcmp(x.c_str(),"Phi")==0)x="#Phi";
  if(strcmp(x.c_str(),"Pt")==0)x="Pt [GeV]";
  histo->SetTitle(title);
  histo->GetXaxis()->SetTitle(x.c_str());
  histo->GetYaxis()->SetTitle(title_y);
  histo->SetFillColor(color);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  
}

//This function 'labels' a TH2D
void TH2D_config(TH2D *histo, const char* title, const char* title_x, const char* title_y, int color){
 
  histo->SetTitle(title);
  histo->GetXaxis()->SetTitle(title_x);
  histo->GetYaxis()->SetTitle(title_y);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);

}

//This function 'labels' a TH2F
void TH2F_config(TH2F *histo, const char* title, const char* title_x, const char* title_y, int color){
    
  histo->SetTitle(title);
  histo->GetXaxis()->SetTitle(title_x);
  histo->GetYaxis()->SetTitle(title_y);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);

}

//This function analyses whole the histograms contained in a .root file and writes the output in a .txt/.dat file
//Still need to import the script...
void LS_config(std::string rootFile, std::string outfile){

  std::cout<<"Removing old config file: "<<outfile.c_str()<<std::endl;
  system(Form("rm -rf %s",outfile.c_str()));
  std::cout<<"Root File analyzed: "<<rootFile.c_str()<<std::endl;
  system(Form("./scripts/ls_script.sh %s | tee %s",rootFile.c_str(),outfile.c_str()));

}
