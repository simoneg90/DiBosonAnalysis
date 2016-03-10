#! /usr/bin/env python


import os
import sys
import optparse
import datetime


from optparse import OptionParser

parser = OptionParser()

parser.add_option("-i", "--input", action='store', dest="rootList", help="write report to FILE", metavar="FILE", default="files.txt")
parser.add_option("-l", "--list", action='store', dest="cutList", help="write report to FILE", metavar="FILE", default="list.txt")
parser.add_option("-o", "--outFile", action='store', dest="outFile", help="write report to FILE", metavar="FILE", default="list_tmp.txt")
parser.add_option("-L", "--Lumi", action='store', dest="Lumi", help="write report to FILE", metavar="FILE", default="list.txt")
parser.add_option("-s", "--suffix", action='store', dest="suffix", help="write report to FILE", metavar="FILE", default="")
parser.add_option("-P", "--prefix", action='store', dest="prefix",help="write report to FILE", metavar="FILE", default="outFolder")
(opt, args) = parser.parse_args()

rootList=opt.rootList
cutList=opt.cutList
outFile=opt.outFile
Lumi=opt.Lumi
suffix=opt.suffix
prefix=opt.prefix
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
      print cols[0], " ", cols[1], " ",cols[2], " ",cols[3], " ",cols[4], " ",cols[5], " ",cols[6]
      os.system("echo "+myList[0]+" "+myList[1]+" "+myList[2]+" "+cols[0]+" "+cols[1]+" "+cols[2]+" "+cols[3]+" "+cols[4]+" "+cols[5]+" "+cols[6]+">>"+ outFile)

    print "File recorded"
    os.system("./doPlots "+outFile+" "+Lumi+" "+prefix+"/"+cols[0]+suffix)

