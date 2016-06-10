#!/bin/bash
#USAGE ./scripts/makeList.sh /pnfsPATH/0000/ listName

rm -rf ${2}.txt
touch ${2}.txt

echo ${2}".txt"

for file in $1*.root
do
  echo "root://xrootd.unl.edu/"${file} >> ${2}.txt
  echo ${1}"/"${file}

done

echo "Finished"
