#!/bin/bash

#USAGE ./scripts/scpToCern.sh user whatToCopy whereIfYouWant
#$2 can be omitted!

userLetter="$(echo $1 | head -c 1)"

echo "+++++++++++++++++"
echo "Hi there!"
echo "I'm going to copy ${2} here /afs/cern.ch/user/${userLetter}/${1}/www/${3}"
echo "Remember to add \'\' before and after what you want to copy"
echo "Are you sure???"
echo "+++++++++++++++++"
echo "[Enter to confirm, Ctrl+C to exit]"
read

scp -r ${2} ${1}@lxplus.cern.ch:/afs/cern.ch/user/${userLetter}/${1}/www/${3}
scp ${1}@lxplus.cern.ch:/afs/cern.ch/user/${userLetter}/${1}/www/index.php ${1}@lxplus.cern.ch:/afs/cern.ch/user/${userLetter}/${1}/www/${3}

