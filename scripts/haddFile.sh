#!/bin/bash
#USAGE ./scripts/haddFile.sh ttbar_output/ fileFolder/ user[for t3 folder]
#where ttbar_output/ contains all the folders with the .root files
echo "REMEMBER TO DO CMSENV!"

echo "Creating folder " ${2}
mkdir -v -p ${2}

for file in $(ls ${1})
do
  fileName=${file:0:(${#file}-16)}
  echo "OUTPUT " ${fileName}".root"
  rm -rf ${2}/${fileName}.root
  echo "Removed " ${2}"/"${fileName}".root"
  hadd ${2}/${fileName}.root ${1}/${file}/*reduced*
  echo $file
done

timestamp="$(date +"%m_%d_%Y_%H_%M_%S")"
echo ${timestamp}
#Only when working in Rome's T2
#cp -r ${2} /t3/users/${3}/rootFolder_${timestamp}

echo "Finished"
