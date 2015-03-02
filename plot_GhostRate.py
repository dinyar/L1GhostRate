#!/usr/bin/python

from ROOT import *
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../L1AnalysisHelpers"))
from CreateHistograms import *

gROOT.Reset()
gROOT.SetBatch(kTRUE);

efficiencyList = []
# Entries: Label for histogram (Will be used for filename and title) | binning | parameters used for project functions
efficiencyList.append(["mu_recoEta", binningDict["eta"], "Eta1_reco", cutDict["diMu-gmtPt1"], cutDict["recoPt1"], [0, 0.05]])
efficiencyList.append(["mu_recoPhi", binningDict["phi"], "Phi1_reco", cutDict["diMu-gmtPt1"], cutDict["recoPt1"], [0, 0.05]])
efficiencyList.append(["mu_recoPt", binningDict["pt"], "pT1_reco", cutDict["diMu-gmtPt1"], cutDict["recoPt1"], [0, 0.05]])
# TODO: Add invariant mass calculation.

rateList = []
# Entries: Label for histogram (Will be used for filename and title) | binning | parameters used for project functions
rateList.append(["mu_recoPt", binningDict["ptFine"], "pT1_reco", cutDict["recoPt1"]])
rateList.append(["mu_recoPt", binningDict["ptFine"], "pT1_reco", cutDict["diMu-gmtPt1"]])
rateList.append(["mu_recoEta", binningDict["etaFine"], "Eta1_reco", cutDict["recoPt1"]])
rateList.append(["mu_recoEta", binningDict["etaFine"], "Eta1_reco", cutDict["diMu-gmtPt1"]])
rateList.append(["mu_recoPhi", binningDict["phiFine"], "Phi1_reco", cutDict["recoPt1"]])
rateList.append(["mu_recoPhi", binningDict["phiFine"], "Phi1_reco", cutDict["diMu-gmtPt1"]])

f = TFile.Open("SingleMuNtuple.root")

ntuple = f.Get("ntuple")

for varList in efficiencyList:
    generateGhostPercHist(varList)

for varList in rateList:
    generateRateHist(varList)
