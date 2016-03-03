#!/bin/bash
#USAGE ./scripts/makeList.sh /pnfsPATH/0000/ listName

rm -rf ${2}.txt
touch ${2}.txt

echo ${2}".txt"

for file in $1*.root
do
  echo "dcap://cmsrm-se01.roma1.infn.it/"${file} >> ${2}.txt
  echo ${1}"/"${file}

done

echo "Finished"
