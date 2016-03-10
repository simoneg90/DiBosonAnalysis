#!/bin/bash
#USAGE ./scripts/haddFile.sh ttbar_output/ fileFolder/ 
#where ttbar_output/ contains all the folders with the .root files
echo "REMEMBER TO DO CMSENV!"

echo "Creating folder " ${2}
mkdir ${2}

for file in $(ls ${1})
do
  fileName=${file:0:(${#file}-16)}
  echo "OUTPUT " ${fileName}".root"
  rm -rf ${2}/${fileName}.root
  echo "Removed " ${1}${file}"/"${fileName}".root"
  hadd ${2}/${fileName}.root ${1}/${file}/*reduced*
  echo $file
done

echo "Finished"
