#! /usr/bin/env python

import os
import sys
import optparse
import datetime

usage = "usage: To be run from DijetRootTreeAnalyzer/ :  python scripts/submit_batch_T2.py -i directory_containing_lists_to_run -o output_directory"

parser = optparse.OptionParser("submitAllGJetsID.py")
parser.add_option('-q', '--queue',       action='store',     dest='queue',       
    help='run in batch in queue specified as option (default -q cmslong)', 
    default='cmsan',
    metavar="QUEUE")

parser.add_option("-i", "--input", dest="input",
    help="directory containing the lists of samples to be analyzed",
    )

parser.add_option("-o", "--output", dest="output",
    help="the directory OUTDIR contains the output of the program",
    metavar="OUTDIR")

parser.add_option("-m", "--match", dest="match",
    help="run only the samples containing this string in the name",
    default="")

parser.add_option("--tag", dest="tag",
    help="useful tag to identify the version",
    default="")

parser.add_option("--split", dest="filesperjob", type=int,
    help="files to analyze per job ",
    default=2)

parser.add_option('-I', '--interactive',      
    action='store_true',
    dest='interactive',      
    help='run the jobs interactively, 2 jobs at a time',
    default=False)

parser.add_option("-c", "--cut", dest="cutfile",
    help="cutfile",                  
    default="config/cutFile_ttbarSelection.txt",#cutFile_mainDiBosonSelection.txt ",
    )

(opt, args) = parser.parse_args()
################################################

###
current_time = datetime.datetime.now()
simpletimeMarker = "_%04d%02d%02d_%02d%02d%02d" % (current_time.year,current_time.month,current_time.day,current_time.hour,current_time.minute,current_time.second) 
timeMarker = "mycutFile_%04d%02d%02d_%02d%02d%02d__" % (current_time.year,current_time.month,current_time.day,current_time.hour,current_time.minute,current_time.second) 
cutfileName = timeMarker+os.path.split(opt.cutfile)[1]
print cutfileName
###
os.system("mkdir -p "+opt.output)#+simpletimeMarker)
#os.system("rm -rf batch")
os.system("mkdir -p batch")
pwd = os.environ['PWD']
###
os.system("cp "+opt.cutfile+" batch/"+cutfileName)
###
newTag = opt.tag#+simpletimeMarker

if not opt.match:
  print opt.input
  os.system("ls "+opt.input+" | grep txt > config/lists_to_run.txt")
  os.system("cat config/lists_to_run.txt")
else:
  os.system("ls "+opt.input+" | grep "+opt.match+"  > config/lists_to_run.txt")

ins = open("config/lists_to_run.txt", "r") 


#lists of lists (in each position there is a list of commads)
inputlists = []
njobs_list = []
#commands = []
#hadd_cmd = []
#filenames_skim = []

##create a directory to store logfiles containing the "tag" in the name
os.system("mkdir batch/"+newTag)
os.system("mkdir "+opt.input+"/splitted")

##loop over lists (one for datasets) to create splitted lists
for line in  ins:
  lists_dataset = []
  sample = os.path.splitext(line)[0]
  print ("process %s" % sample)
  line = line.rstrip('\n')
  sample = sample.rstrip('\n')
  # name of final rootfile merged
  #filenames_skim.append(opt.output+simpletimeMarker+"/rootfile_"+sample+"_"+newTag+"_"+str(jj)+"_reduced_skim.root")
  # list of commands for current dataset
  #commands_dataset = []
  
  jf=0  ##counter for num of files
  jj=0  ##counter for num of jobs
  #open list
  list = open(opt.input+"/"+line,"r") 
  ## remove splitted lists if they already exist (necessary beacuse we append to txt file)
  os.system("rm "+opt.input+"/splitted/"+sample+"_"+newTag+"*.txt")
  for file in list:
    #print "file:%i  filesperjob:%i  job:%i op.modulo:%i  list %s " % (jf, opt.filesperjob,jj,(jf+1 % opt.filesperjob), opt.input+"/"+line)
    #print file
    modulo = int(jf+1) % int(opt.filesperjob)
    #print "modulo = %i" % modulo
    splittedlist = open(opt.input+"/splitted/"+sample+"_"+newTag+"_"+str(jj)+".txt","a+")
    splittedlist.write(file)
    if ( modulo == 0 ):
      lists_dataset.append(opt.input+"/splitted/"+sample+"_"+newTag+"_"+str(jj)+".txt")
      print "job "+str(jj)+"   appending "+opt.input+"/splitted/"+sample+"_"+newTag+"_"+str(jj)+".txt"
      jj += 1 #increment counter of jobs  
    jf += 1   #increment counter of files         
  print "job "+str(jj)+"   appending "+opt.input+"/splitted/"+sample+"_"+newTag+"_"+str(jj)+".txt"
  lists_dataset.append(opt.input+"/splitted/"+sample+"_"+newTag+"_"+str(jj)+".txt")  
  njobs_list.append(jj)
  inputlists.append(lists_dataset)

print ""
print njobs_list
print ""
print inputlists
print "inputlists size = "+str(len(inputlists))

i_f = 0
ins = open("config/lists_to_run.txt", "r") 
for line in  ins:
  lists_dataset = []
  sample = os.path.splitext(line)[0]
  line = line.rstrip('\n')
  sample = sample.rstrip('\n')
  splittedlist = inputlists[i_f]
  print splittedlist
  for jj in range(0,njobs_list[i_f]+1):
    print sample+"  job "+str(jj)
    #command = "./main "+splittedlist[jj]+" config/cutFile_mainDijetSelection.txt dijets/events "+opt.output+simpletimeMarker+"/rootfile_"+sample+"_"+newTag+"_"+str(jj)+" "+opt.output+simpletimeMarker+"/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)
    #command = "./main "+splittedlist[jj]+" config/cutFile_mainDijetSelection.txt dijets/events /tmp/rootfile_"+sample+"_"+newTag+"_"+str(jj)+" /tmp/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)
    command = "./main "+splittedlist[jj]+" batch/"+cutfileName+" tree /tmp/rootfile_"+sample+"_"+newTag+"_"+str(jj)+" /tmp/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)
    print "submit "+command
    print ""
    
    logfile = "batch/"+newTag+"/logfile_"+sample+"_"+newTag+"_"+str(jj)+".log"
    outputname = "batch/"+newTag+"/submit_"+sample+"_"+newTag+"_"+str(jj)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    #outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc481\n')
    outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc491\n')
    outputfile.write('cd '+pwd+' \n')
    outputfile.write('eval `scramv1 runtime -sh`\n')
    outputfile.write(command+"\n")
    outputfile.write("dccp /tmp/rootfile_"+sample+"_"+newTag+"_"+str(jj)+"_reduced_skim.root "+opt.output+"\n")#simpletimeMarker+"\n")
    outputfile.write("dccp /tmp/rootfile_"+sample+"_"+newTag+"_"+str(jj)+".root "+opt.output+"\n")#simpletimeMarker+"\n")
    outputfile.write("dccp /tmp/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)+".dat "+opt.output+"\n")#simpletimeMarker+"\n")
    outputfile.write("rm /tmp/rootfile_"+sample+"_"+newTag+"_"+str(jj)+"_reduced_skim.root\n")
    outputfile.write("rm /tmp/rootfile_"+sample+"_"+newTag+"_"+str(jj)+".root\n")
    outputfile.write("rm /tmp/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)+".dat\n")
    
    print outputname 
    if opt.interactive==False:
      os.system("bsub -q "+opt.queue+" -o "+logfile+" source "+pwd+"/"+outputname)
    else:
      print logfile
      if imc==0: os.system(command+" >&! "+logfile+"&")
      else: os.system(command+" >&! "+logfile)
  i_f += 1
      
