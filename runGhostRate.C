{
    gROOT->ProcessLine(".x ../L1TriggerDPG/L1Ntuples/macros/initL1Analysis.C");
    gROOT->ProcessLine(".L GMTGhostRateNtupleizer.C+");
    gROOT->ProcessLine("GMTGhostRateNtupleizer macro = GMTGhostRateNtupleizer(\"file_list\");");
    gROOT->ProcessLine("macro.run(-1);");
}
