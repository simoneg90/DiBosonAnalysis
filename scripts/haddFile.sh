#!/bin/bash
#USAGE ./scripts/haddFile.sh ttbar_output/ 
#where ttbar_output/ contains all the folders with the .root files
echo "REMEMBER TO DO CMSENV!"

for file in $(ls ${1})
do
  fileName=${file:0:(${#file}-16)}
  echo "OUTPUT " ${fileName}".root"
  rm -rf ${1}${file}/${fileName}.root
  echo "Removed " ${1}${file}"/"${fileName}".root"
  hadd ${1}${file}/${fileName}.root ${1}/${file}/*reduced*
  echo $file
done

echo "Finished"
