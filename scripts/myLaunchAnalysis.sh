#!/bin/bash

#USAGE ./scripts/myLaunchAnalysis ttbar_analysis/ ttbar_output/
#where ttbar_analysis/ contains all the folder with the list.txt of the /pnfs/ files inside

echo "REMEMBER TO DO CMSENV!"

for file in $(ls ${1})
do
  echo ${1}$file"/"
  echo "OUTPUT" ${2}"/"$file
  python scripts/submit_batch_T2_split.py -q 8nh -i  ${1}$file"/" -o ${2}"/"$file

done

echo "Finished"
