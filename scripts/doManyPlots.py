#! /usr/bin/env python


import os
import sys
import optparse
import datetime
import time

start = time.time()

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-i", "--input", action='store', dest="rootList", help="write report to FILE", metavar="FILE", default="files.txt")
parser.add_option("-l", "--list", action='store', dest="cutList", help="write report to FILE", metavar="FILE", default="list.txt")
parser.add_option("-o", "--outFile", action='store', dest="outFile", help="write report to FILE", metavar="FILE", default="list_tmp.txt")
parser.add_option("-L", "--Lumi", action='store', dest="Lumi", help="write report to FILE", metavar="FILE", default="list.txt")
parser.add_option("-s", "--suffix", action='store', dest="suffix", help="write report to FILE", metavar="FILE", default="")
parser.add_option("-P", "--prefix", action='store', dest="prefix",help="write report to FILE", metavar="FILE", default="outFolder")
parser.add_option("-e", "--efficiency", action='store', dest="efficiency",help="write report to FILE", metavar="FILE", default="effFile")
parser.add_option("-d", "--doScale", action='store', dest="doScale", help="write report to FILE", metavar="FILE", default="1")
parser.add_option("-S", "--superimpose",action='store', dest="superimpose", help="write report to FILE", metavar="FILE", default="0")
parser.add_option("-Y", "--logY",action='store', dest="logY", help="write report to FILE", metavar="FILE", default="0")
(opt, args) = parser.parse_args()

rootList=opt.rootList
cutList=opt.cutList
outFile=opt.outFile
Lumi=opt.Lumi
doScale=opt.doScale
suffix=opt.suffix
prefix=opt.prefix
superimpose=opt.superimpose
logY=opt.logY
efficiency=opt.efficiency
print "Suffix ", suffix
print "OutFile ", outFile
os.system("mkdir -p "+prefix)
i=0
print rootList
#for subdir, dirs,files in os.walk(rootFolder):
#  for f in files:
#    print i," ", f
#    i+=1
#    #with open(cutList) as 
#for files in open(rootList):
#    myList = [L for L in files.split()]
#    print myList[0], " ", myList[1]
#    for line in open(cutList):#source:
#      cols = [x for x in line.split()]
#      print cols[0], " ", cols[1]
#      os.system("echo ") 

for line in open(cutList):#source:
    print "removing tmp file"
    os.system("rm -rf "+outFile)
    cols = [x for x in line.split()]
    print cols[0], " ", cols[1], " ",cols[2]
    for files in open(rootList):
      myList = [L for L in files.split()]
      print myList[0], " ", myList[1]
      os.system("echo "+myList[0]+" "+myList[1]+" "+myList[2]+" "+cols[0]+" "+cols[1]+" "+cols[2]+" "+cols[3]+" "+cols[4]+" "+cols[5]+" "+cols[6]+">>"+ outFile)
      #if  myList[1].find('SinglePhoton')!=-1:
      #  print "SinglePhoton dataset"
      #  trigger='\&\&\(HLT_BIT_HLT_Photon165_HE10_v==1\|\|HLT_BIT_HLT_Photon175_v==1\)'
      #  print trigger
      #  #cols[1]+="&&HLT_BIT_HLT_Photon165_HE10_v==1&&HLT_BIT_HLT_Photon175_v==1"
      #  os.system("echo "+myList[0]+" "+myList[1]+" "+myList[2]+" "+cols[0]+" "+cols[1]+trigger+" "+cols[2]+" "+cols[3]+" "+cols[4]+" "+cols[5]+" "+cols[6]+">>"+ outFile)
      #elif myList[1].find('JetHT')!=-1:
      #  print "JetHT dataset"
      #  trigger='\&\&\(\(HLT_BIT_HLT_Photon165_HE10_v==0\|\|HLT_BIT_HLT_Photon175_v==0\)\&\&HLT_BIT_HLT_PFHT800_v==1\)'
      #  print trigger
      #  os.system("echo "+myList[0]+" "+myList[1]+" "+myList[2]+" "+cols[0]+" "+cols[1]+trigger+" "+cols[2]+" "+cols[3]+" "+cols[4]+" "+cols[5]+" "+cols[6]+">>"+ outFile)
      #else:
      #  print cols[0], " ", cols[1], " ",cols[2], " ",cols[3], " ",cols[4], " ",cols[5], " ",cols[6]
      #  os.system("echo "+myList[0]+" "+myList[1]+" "+myList[2]+" "+cols[0]+" "+cols[1]+" "+cols[2]+" "+cols[3]+" "+cols[4]+" "+cols[5]+" "+cols[6]+">>"+ outFile)

    print "File recorded"
    os.system("./doPlots "+outFile+" "+Lumi+" "+prefix+"/"+cols[0]+suffix+" "+efficiency+" "+doScale+" "+superimpose+" "+logY)




end = time.time()
print ">TIME ELAPSED: ",(end - start)
