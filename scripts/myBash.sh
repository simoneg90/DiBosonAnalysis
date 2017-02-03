#!/bin/bash -l
export PATH=$PATH:/opt/exp_soft/local/software/root-5.28/bin
source /cvmfs/cms.cern.ch/cmsset_default.sh

folder=$PWD
name=$1
outFolder=$2
fileList=$3
cutFile=$4
echo ${folder}
#cd /cmshome/gellisim/MCP/CMSSW_6_2_0/
#eval `scramv1 runtime -sh`
#cd ..
#cd BTF_TB_SW1
#python reco_script.py


echo "export PATH=\$PATH:/opt/exp_soft/local/software/root-5.28/bin"  > ${name}.sh
echo "source /cvmfs/cms.cern.ch/cmsset_default.sh" >> ${name}.sh
echo "cd "${folder} >> ${name}.sh
echo "eval \`scramv1 runtime -sh\`" >> ${name}.sh
echo "./main "${folder}"/"${fileList}" "${folder}"/"${cutFile}" tree "${folder}"/"${outFolder}"/outFile "${folder}"/"${outFolder}"/cutEfficiencyList" >> ${name}.sh 
chmod +x ${name}.sh
echo "bsub -q cmsan source ${folder}"/"${name}.sh -o ${folder}"/"${name}.log"
bsub -q cmslong source ${folder}/${name}.sh -o ${folder}/${name}.log
#./${name}.sh
