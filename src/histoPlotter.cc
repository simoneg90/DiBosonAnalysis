#include "histoPlotter.h"
#include "setTDRStyle.h"
#include "utility.h"

#include <fstream>

int main(int argc, char* argv[]){

  frame("Removed the scale part");
  frame("Histo Plotter env");


  TH1D *histoD[MAX_NUMBER];
  TH1F *histoF[MAX_NUMBER];
  TH2D *histo2D[MAX_NUMBER];
  TH2F *histo2F[MAX_NUMBER];
  TCanvas *c[MAX_NUMBER];
  TLegend *leg[MAX_NUMBER];
  TFile *file[MAX_NUMBER];

  TStopwatch time;
  time.Start(true);

  if(argc<=1){
    breakLine();
    std::cout<<"Attention! Not enough inputs!"<<std::endl;
    std::cout<<"Digit './histoPlotter --help' for informations!"<<std::endl;
    breakLine();
    exit(-1);
  }


  if(strcmp(argv[1],"--help")==0){
    breakLine();
    std::cout<<"This program superimposes or Stacks the same histogram from different files"<<std::endl;
    std::cout<<"Program usage:"<<std::endl;
    std::cout<<"./histoPlotter inputList"<<std::endl;
    //argv[1]=inputList 
    breakLine();
    exit(-1);

  }

  //if(argc>MAX_NUMBER+1){
  if(argc>6){
    breakLine();
    std::cout<<"ATTENTION! TOO MANY ARGUMENTS ADDED!"<<std::endl;
    std::cout<<"Exiting program"<<std::endl;
    breakLine();
    exit(-1);
  }
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
  //while(inputList.good()){
    if(!inputList.good()) break;
    maxFiles = j;
    inputList>>filePath>>isData>>evts[j]>>xsec[j]>>legendName[j];
    std::cout<<"File: "<<filePath.c_str()<<" IsData: "<<isData<<" evt: "<<evts[j]<<" xsec: "<<xsec[j]<<" legend: "<<legendName[j]<<std::endl;
    if(j==0){
        breakLine();
        std::cout<<"Launching histo_config script!"<<std::endl;
        LS_config(filePath.c_str(),"config_file.txt");
        breakLine();
    }
    file[j]=TFile::Open(filePath.c_str());
    if(file[j]==0){
      breakLine();
      std::cout<<"ATTENTION! "<<argv[j+1]<<" doesn't exists!"<<std::endl;
      breakLine();
      exit(-1);
    }//end if (file exists)
    //++fileCounter;
  }//end loop over files=>from list

  std::ifstream inFile;
  inFile.open("config_file.txt");

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

  setTDRStyle();
  //setTDRStyle.tdrStyle->SetMarkerStyle(19);
  std::string x_axis="";
  for (int counter=0; counter<MAX_NUMBER; ++counter){
    if(!inFile.good()) break;
  //while(inFile.good()){
    inFile>>word;

    //== Segmentation Fault protection ==
    if(histoD_counter>(MAX_NUMBER-argc+1) || histoF_counter>(MAX_NUMBER-argc+1) || histo2D_counter>(MAX_NUMBER-argc+1) || histo2F_counter>(MAX_NUMBER-argc+1)){
      breakLine();
      std::cout<<"ATTENTION!!! TOO MANY ARGUMENTS INSERTED!!! PLEASE CHANGE VARIABLE 'MAX_NUMBER'"<<std::endl;
      std::cout<<"# TH1D = "<<histoD_counter<<" # TH1F = "<<histoF_counter<<" # TH2D = "<<histo2D_counter<<" # TH2F = "<<histo2F_counter<<std::endl;
      std::cout<<"# Canvas = "<<canvas_counter<<" # File = "<<argc-1<<std::endl;
      breakLine();
           
      std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << time.RealTime() << " s" <<std::endl;
        
      exit(-1);
    }//end if for Segmentation Fault protection
        
        
    if(histo_type!=0){//reading histo name
      std::cout<<"Histo name: "<<(word.substr(0,word.size()-2)).c_str()<<std::endl;
      c[canvas_counter]=new TCanvas(Form("c%d",canvas_counter),Form("c%d",canvas_counter), 200,10,600,400);
      leg[canvas_counter]=new TLegend(0.7,0.7,0.9,0.9);
      for(int y=0; word[y]!='_'&& word[y]!='.';++y){
	x_axis+=word[y];
      }

      for(int j=0;j<maxFiles;++j){//Loop over files

	//if(histo_type==1)histo[histo_counter]=((TH1D *)(file1->Get((word.substr(0,word.size()-2)).c_str())));
                   
	if(histo_type==1){
	  histoD[histoD_counter]=((TH1D *)(file[j]->Get((word.substr(0,word.size()-2)).c_str())));
	  if((((TH1D *)(file[j]->Get((word.substr(0,word.size()-2)).c_str())))->GetMaximum())>histo_max){
	    histo_max=((TH1D *)(file[j]->Get((word.substr(0,word.size()-2)).c_str())))->GetMaximum(); 
	    histo_index=j;
	  }
	  histoD[histoD_counter]->SetStats(0);
	  TH1D_config(histoD[histoD_counter],"",x_axis.c_str(),"# Events",j+1);
	  //leg[canvas_counter]->AddEntry(histoD[histoD_counter],Form("file%d",j),"lp");
	  title=argv[j+1];
	  leg[canvas_counter]->AddEntry(histoD[histoD_counter],(title.substr(0,title.size()-5)).c_str(),"lp");
          //leg[canvas_counter]->AddEntry(histoD[histoD_counter-1],legendName[j].c_str(),"lp");

	  ++histoD_counter;
	}//end if TH1D
	if(histo_type==-1){
	  histoF[histoF_counter]=((TH1F *)(file[j]->Get((word.substr(0,word.size()-2)).c_str())));
          //histoF[histoF_counter]=((TH1F *)(file[j]->Get("bGen_pt")));
	  if((((TH1F *)(file[j]->Get((word.substr(0,word.size()-2)).c_str())))->GetMaximum())>histo_max){
	    histo_max=((TH1F *)(file[j]->Get((word.substr(0,word.size()-2)).c_str())))->GetMaximum();
	    histo_index=j;
                                                                               
	  }

	  histoF[histoF_counter]->SetStats(0);
	  TH1F_config(histoF[histoF_counter],"",x_axis.c_str(),"# Events",j+1);
	  //title=argv[j+1];
	  //leg[canvas_counter]->AddEntry(histoF[histoF_counter],(title.substr(0,title.size()-5)).c_str(),"lp");
          leg[canvas_counter]->AddEntry(histoF[histoF_counter-1],legendName[j].c_str(),"lp");
	  if(j==0) histoF[histoF_counter]->Draw();
	  if(j>0) histoF[histoF_counter]->Draw("SAME");

	  ++histoF_counter;
	}//end if TH1F
	if(histo_type==2){
	  histo2D[histo2D_counter]=((TH2D *)(file[j]->Get((word.substr(0,word.size()-2)).c_str())));
	  histo2D[histo2D_counter]->SetStats(0);
	  TH2D_config(histo2D[histo2D_counter],"Titolo","x-axis","y-axis",j+1);
	  //title=argv[j+1];
	  leg[canvas_counter]->AddEntry(histo2D[histo2D_counter],(title.substr(0,title.size()-5)).c_str(),"lp");
          leg[canvas_counter]->AddEntry(histo2D[histo2D_counter],legendName[counter].c_str(),"lp");

	  if(j==0) histo2D[histo2D_counter]->Draw();
	  if(j>0) histo2D[histo2D_counter]->Draw("SAME");
	  ++histo2D_counter;  
	}//end if TH2D

	if(histo_type==-2){
	  histo2F[histo2F_counter]=((TH2F *)(file[j]->Get((word.substr(0,word.size()-2)).c_str())));
	  histo2F[histo2F_counter]->SetStats(0);
	  TH2F_config(histo2F[histo2F_counter],"","x-axis","y-axis",j+1);
	  //leg[canvas_counter]->AddEntry(histo2F[histo2F_counter],Form("file%d",j),"lp");
          leg[canvas_counter]->AddEntry(histo2F[histo2F_counter],legendName[histo2D_counter].c_str(),"lp");

	  if(j==0) histo2F[histo2F_counter]->Draw();
	  if(j>0) histo2F[histo2F_counter]->Draw("SAME");
	  ++histo2F_counter;
	}//end if TH2F
      }//end for file_number

      double integral, integral0;
      //== Draw INFO ==
      //== Scaling to the histo with more events ==

      if(histo_type==1){
         integral0=histoD[histoD_counter-argc+1+histo_index]->Integral();
         histoD[histoD_counter-argc+1+histo_index]->Scale(1/integral0);
         histoD[histoD_counter-argc+1+histo_index]->Draw();

         for(int j=0; j<(argc-1);++j){
             integral=histoD[histoD_counter-argc+1+j]->Integral();
             if(j!=histo_index){
                histoD[histoD_counter-argc+1+j]->Scale(1/integral);
                //histoD[histoD_counter-argc+1+j]->Scale(integral0/integral);
                histoD[histoD_counter-argc+1+j]->Draw("SAME");
             }
         }
      }
                                                                            
  /*    if(histo_type==-1){
          histoF[histoF_counter-argc+1+histo_index]->Draw();
          integral0=histoF[histoF_counter-argc+1+histo_index]->Integral();
                             
          for(int j=0; j<(argc-1);++j){
              integral=histoF[histoF_counter-argc+1+j]->Integral();
              if(j!=histo_index){
                 histoF[histoF_counter-argc+1+j]->Scale(integral0/integral);
                 histoF[histoF_counter-argc+1+j]->Draw("SAME");
              }
          }
      }
*/
      //== End Draw INFO==                                                                                                   

      x_axis="";
      leg[canvas_counter]->Draw("SAME");
      c[canvas_counter]->SaveAs(Form("grafici/%s.pdf",(word.substr(0,word.size()-2)).c_str()));
      ++canvas_counter;
      histo_type=0;
                                  
    }//end if histo name

    if(key_costraint==1){//We are reading the histo type
        if(strcmp(word.c_str(),"TH1D")==0) histo_type=1;
        if(strcmp(word.c_str(),"TH1F")==0) histo_type=-1;
        if(strcmp(word.c_str(),"TH2D")==0) histo_type=2;
        if(strcmp(word.c_str(),"TH2F")==0) histo_type=-2;

        std::cout<<"Histo type: "<<histo_type<<" "<<word.c_str()<<std::endl;

        key_costraint=0;

    }//end if for histo type

    if(strcmp(word.c_str(),"KEY:")==0){
        key_costraint=1;
        ++i;
    }//end if (KEY: costraint)

    histo_max=0;
      
  }//end while inFile.good()

  breakLine();
  std::cout<<"# TH1D = "<<histoD_counter<<" # TH1F = "<<histoF_counter<<" # TH2D = "<<histo2D_counter<<" # TH2F = "<<histo2F_counter<<std::endl;
  std::cout<<"# Canvas = "<<canvas_counter<<" # File = "<<argc-1<<std::endl;
  breakLine();
  
  std::cout<<""<<std::endl;
  std::cout << "------ TIME ELAPSED DURING ANALYSIS  ----- " << time.RealTime() << " s" <<std::endl;
  std::cout<<""<<std::endl;
}
