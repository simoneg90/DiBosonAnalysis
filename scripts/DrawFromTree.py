#!usr/bin/python

import sys, os, subprocess, string, re
from ROOT import *
from array import array
import CMS_lumi
import optparse
from setTDRStyle import setTDRStyle

CMS_lumi.extraText = "Simulation Preliminary"
CMS_lumi.lumi_sqrtS = "1 fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0

gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadBottomMargin(0.15)
gROOT.ForceStyle()

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')


#######################################################

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--var",action="store",type="string",dest="var",default='ptHat')
parser.add_option("--xmin",action="store",type="float",dest="xmin",default=1)
parser.add_option("--xmax",action="store",type="float",dest="xmax",default=1)
parser.add_option("--xtitle",action="store",type="string",dest="xtitle",default='')
parser.add_option("--bins",action="store",type="int",dest="bins",default=11111111111)
parser.add_option("--rebin",action="store",type="int",dest="rebin",default=1)
parser.add_option("--logy",action="store_true",default=False,dest="logy")
parser.add_option("--outputDir",action="store",type="string",default="./",dest="outputDir")
parser.add_option("--inputList",action="store",type="string",default="list.txt",dest="inputList")
parser.add_option("--lumi",action="store",type="float",default="1000.",dest="lumi")

(options, args) = parser.parse_args()

var = options.var
xmin = options.xmin
xmax = options.xmax
bins = options.bins
xtitle = options.xtitle
rebin = options.rebin
logy = options.logy
outputDir = options.outputDir
inputList = options.inputList
lumi = options.lumi

#fileNames = ['QCD_Pt-301to470','QCD_Pt-470to600','QCD_Pt-600to800', 'QCD_Pt-800to1000', 'QCD_Pt-1000to1400', 'QCD_Pt-1400to1800', 'QCD_Pt-1800to2400', 'QCD_Pt-2400to3200', 'QCD_Pt-3200']
#xsections = [7475 ,587., 167, 28.25, 8.195, 0.7346, 0.102, 0.00644, 0.000163]
#colorF    = [ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9,ROOT.kBlue-9]
#colorL    = [ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack,ROOT.kBlack]
hist_allCuts      = []

LUMI      = lumi
#PATH      = inputDir

#---- read the list -----------------
lines = [line.strip() for line in open(inputList)]

#---- split sample name and xsec
fileNames = []
xsecs = []
ii = 0
for line in lines:
  parts = line.split()
  fileNames.append(parts[0])
  xsecs.append(parts[1])
  print ("dataset : %s    xsec : %s" % (fileNames[ii], xsecs[ii]))


  ii+=1

#---- open the files --------------------
i_f = 0
for f in fileNames:
  inf = TFile.Open(f)
  print inf.GetName()
  
  Nev = inf.Get('DijetFilter/EventCount/EventCounter').GetBinContent(1)
  print ('processed events: %s' % Nev)
  wt = 1.0
  #if i_f < 3:

  h_allCuts = TH1F("h_allCuts", "", bins, xmin, xmax)
  tree = inf.Get('rootTupleTree/tree')
  tree.Project(h_allCuts.GetName(), var,'deltaETAjj < 1.3')
  Npassed = h_allCuts.GetEntries()
  eff = float(Npassed)/Nev
  print('eff : %f' % eff)
  print('(not using efficiency in the weight)')
  wt = LUMI*float(xsecs[i_f])/Nev
  print('weight : %f' % wt)
  h_allCuts.Scale(wt)
  h_allCuts.SetDirectory(0)
  h_allCuts.SetFillColor(9)
  h_allCuts.SetLineColor(9)
  h_allCuts.SetMarkerColor(9)
  hist_allCuts.append(h_allCuts)
   
  i_f += 1

h_sig = hist_allCuts[9].Clone()
h_sig.SetLineStyle(2)
h_sig.SetLineWidth(2)
h_sig.SetLineColor(kRed)

#kFactor = NDAT/NQCD
kFactor = 1.3
print ("kFactor = %f" % kFactor)

                                      
#for i in range(0,len(fileNames)) :
for i in range(0,9) :
  hist_allCuts[i].Scale(kFactor)
 
NQCD_allCuts = hist_allCuts[0].Integral()

#for i in range(0,len(fileNames)) :
for i in range(0,9) :
  NQCD_allCuts += hist_allCuts[i].Integral()
    
hist_allCutsQCD = hist_allCuts[0].Clone('hist_allCutsQCD')

#for i in range(1,len(fileNames)):
for i in range(1,9):
  hist_allCutsQCD.Add(hist_allCuts[i])

hsQCD_allCuts = THStack('QCD_allCuts','QCD_allCuts')

#for i in range(0,len(fileNames)) :
for i in range(0,9) :
  hsQCD_allCuts.Add(hist_allCuts[i])


print ("---- After scaling signal to bkg (if not plotting mjj) -----")
print ("bkg integral all cuts = %f" % NQCD_allCuts)


#----- Drawing  and save on file -----------------------

os.system("mkdir -p "+options.outputDir)

outFile = TFile(outputDir+"histo_signal_bkg_"+var+"_fromTree.root", "recreate")
outFile.cd()
hist_allCutsQCD.Rebin(rebin)
hist_allCutsQCD.Write()

can_allCuts = TCanvas('can_allCuts_'+var,'can_allCuts_'+var,900,600)

leg = TLegend(0.6, 0.7, 0.85, 0.85)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.AddEntry(hist_allCutsQCD, "QCD", "l")
leg.AddEntry(h_sig, "q* m = 5 TeV", "l")

can_allCuts.cd()
max = hist_allCutsQCD.GetBinContent(hist_allCutsQCD.GetMaximumBin())

if logy:
  gPad.SetLogy()
  hist_allCutsQCD.SetMaximum(1000.*max)  

#hist_allCutsQCD.Reset()
hist_allCutsQCD.GetXaxis().SetRangeUser(xmin,xmax)
hist_allCutsQCD.GetXaxis().SetTitle(xtitle)
#maximumBin = array('f',  [hist_allCutsQCD.GetBinContent(hist_allCutsQCD.GetMaximumBin()), hist_allCutsSig_1000.GetBinContent(hist_allCutsSig_1000.GetMaximumBin()), hist_allCutsSig_5000.GetBinContent(hist_allCutsSig_5000.GetMaximumBin()), hist_allCutsSig_8000.GetBinContent(hist_allCutsSig_8000.GetMaximumBin())])
#max = TMath.MaxElement(4, maximumBin)
#hist_allCutsQCD.SetMaximum(1.2*max)
hist_allCutsQCD.SetMinimum(0.0001)
hist_allCutsQCD.Draw("hist")
h_sig.Draw("same")
leg.Draw()

#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(can_allCuts, iPeriod, iPos)

gPad.RedrawAxis()


can_allCuts.Write()   
can_allCuts.SaveAs(outputDir+var+'_allCuts.png')
can_allCuts.SaveAs(outputDir+var+'_allCuts.svg')
can_allCuts.Close()

outFile.Close()

##----- keep the GUI alive ------------
#if __name__ == '__main__':
#  rep = ''
#  while not rep in ['q','Q']:
#    rep = raw_input('enter "q" to quit: ')
#    if 1 < len(rep):
#      rep = rep[0]
