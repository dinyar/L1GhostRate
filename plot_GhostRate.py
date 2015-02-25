#!/usr/bin/python

from ROOT import *
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../L1AnalysisHelpers"))
from CreateHistograms import *

gROOT.Reset()
gROOT.SetBatch(kTRUE);

efficiencyList = []
# Entries: Label for histogram (Will be used for filename and title) | binning | parameters used for project functions
efficiencyList.append(["mu_recoEta", 50, -2.5, 2.5, "Eta_reco", cutDict["diMu-gmtPt1"], cutDict["recoPt1"]])
# TODO: Add invariant mass calculation.

rateList = []
# Entries: Label for histogram (Will be used for filename and title) | binning | parameters used for project functions
rateList.append(["mu_recoPt", 25, 0, 50, "pT_reco", cutDict["recoPt1"]])
rateList.append(["mu_recoPt", 25, 0, 50, "pT_reco", cutDict["gmtPt1"]])
rateList.append(["mu_recoPt", 25, 0, 50, "pT_reco", cutDict["diMu-gmtPt1"]])
rateList.append(["mu_recoEta", 25, -2.5, 2.5, "Eta_reco", cutDict["recoPt1"]])
rateList.append(["mu_recoEta", 25, -2.5, 2.5, "Eta_reco", cutDict["gmtPt1"]])
rateList.append(["mu_recoEta", 25, -2.5, 2.5, "Eta_reco", cutDict["diMu-gmtPt1"]])
rateList.append(["mu_recoPhi", 25, -3.2, 3.2, "Phi_reco", cutDict["recoPt1"]])
rateList.append(["mu_recoPhi", 25, -3.2, 3.2, "Phi_reco", cutDict["gmtPt1"]])
rateList.append(["mu_recoPhi", 25, -3.2, 3.2, "Phi_reco", cutDict["diMu-gmtPt1"]])

f = TFile.Open("SingleMuNtuple.root")

ntuple = f.Get("ntuple")

for varList in efficiencyList:
    generateEfficiencyHist(varList)

for varList in rateList:
    generateRateHist(varList)
