#!/bin/bash
###usage script.sh fileToReduce
cutString="((abs(lepton_pdgID)==11&&met>80&&lepton_pt>120)||(abs(lepton_pdgID)==13&&met>40&&lepton_pt>40))&&(btag_medium>0)&&ak08Ungroomed_1_pt>200&&WType1_pt>200"
echo $PWD
echo "This script will reduce a rootFile using this string: "${cutString}
rootFile=${1}
rootFile=${rootFile:0:(${#rootFile}-5)}
echo ${rootFile}
root -l <<EOF
TFile *_file0 = TFile::Open("${rootFile}.root")
TDirectory * dir = (TDirectory*)_file0->Get("rootTupleTree");
TTree *oldtree = (TTree *)dir->Get("tree");
TFile *newfile = new TFile("${rootFile}_redu.root", "recreate")
TTree *COPY1= oldtree->CopyTree("${cutString}")
COPY1->Write("mio")
newfile->Close()
_file0->Close()


EOF
