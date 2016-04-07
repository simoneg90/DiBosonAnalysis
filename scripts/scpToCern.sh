#!/bin/bash

#USAGE ./scripts/scpToCern.sh whatToCopy whereIfYouWant
#$2 can be omitted!

echo "+++++++++++++++++"
echo "Hi there!"
echo "I'm going to copy ${1} here /afs/cern.ch/user/s/sgelli/www/${2}"
echo "Remember to add \'\' before and after what you want to copy"
echo "Are you sure???"
echo "+++++++++++++++++"
echo "[Enter to confirm, Ctrl+C to exit]"
read
scp -r ${1} sgelli@lxplus.cern.ch:/afs/cern.ch/user/s/sgelli/www/${2}

