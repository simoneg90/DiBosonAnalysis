#!/bin/bash
source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh
export X509_USER_PROXY=$PWD/X509_USER_PROXY
COPYCMD=xrdcp
EOSPREFIX=root://eoscms//eos/cms/
EOSDIR=/store/group/dpg_ecal/alca_ecalcalib/EcalTiming/
AFSDIR=/afs/cern.ch/work/s/sgelli/private/CMSSW_Heppy/src/CMGTools/MonoXAnalysis/cfg/MCTrees/
FOLDER=/store/group/phys_jetmet/sgelli/MCsamples_V/
USER=$(whoami)
EOSCMD=/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select
echo $USER
query=`$EOSCMD ls /store/group/phys_jetmet/sgelli/MCsamples_V/`
#$FOLDER`
#eos find -f /eos/atlas/user/t/test/ find files
#$EOSCMD ls /store/
echo $query
for folder in $query
do 
  files=`$EOSCMD ls $FOLDER/$folder`
  #echo "File List +++++++++++++++++++++"
  #echo $files
  #echo "+++++++++++++++++++++++++++++++"
  for rootFile in $files
  do
    echo "File to analyze "${FOLDER}"/"${folder}"/"${rootFile}
    newfile=${rootFile:35:10}
    echo "----------------"${#rootFile}"-"${newfile:0:(${#newfile}-5)}
    if [ "ls $AFSDIR"/"$folder"_Chunk"${newfile:0:(${#newfile}-5)}"/skimAnalyzerCount/SkimReport.txt"" ]; then
      echo "exist"
      number=`awk '{print $3}' $AFSDIR"/"$folder"_Chunk"${newfile:0:(${#newfile}-5)}"/skimAnalyzerCount/SkimReport.txt"`
      counter=0
      for number1 in ${number}
      do
        if [ ${counter} = "2" ]; then  #beware. look for data!
          echo "-------"${number1}
          evtProcessed=${number1}
        fi
        echo "counter " $counter 
        counter=$((counter+1))
      done #end loop while reading skimReport

    else
      echo "Warning! File do NOT exists!"
      echo ${FOLDER}"/"${folder}"/"${rootFile} >> problematic.txt
    fi
    echo "Copying file in local "${FOLDER}"/"${folder}"/"${rootFile}
    $COPYCMD root://xrootd.unl.edu/${FOLDER}/${folder}/${rootFile} .
    root -l <<EOF
    TFile *_file0 = TFile::Open("${rootFile}","update")
    TH1F* Count= new TH1F("Count","Count", 1, 0, 2)
    Count->SetBinContent(1,${evtProcessed})
    Count->Write()
    _file0->Close()
EOF
    echo "Removing old file from eos"
    $EOSCMD rm ${FOLDER}/${folder}/${rootFile}
    echo "Copying new file to eos "${FOLDER}"/"${folder}"/"${rootFile}
    $COPYCMD ${rootFile} root://eoscms.cern.ch//eos/cms/${FOLDER}/${folder}/${rootFile}
    echo "Removing local file"
    rm -r ${rootFile}
#exit
  done #looping over files in the folder
done
echo "bella!"
