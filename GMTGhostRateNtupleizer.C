#include "L1Ntuple.h"
#include "hist.C"
#include "Style.C"
#include <sstream>
#include "TNtuple.h"
#include "TFile.h"

// --------------------------------------------------------------------
//                       GMTGhostRateNtupleizer macro definition
// --------------------------------------------------------------------

#define PHI 0
#define ETA 1
#define PT 2
#define PTETA 3

#define DT 0
#define RPCb 1
#define CSC 2
#define RPCf 3
#define GMT 4
#define RECOPASS 5
enum muSysEnum {eDTRPC, eCSCRPC, eDT, eCSC, eRPCb, eRPCf, none};
int Nsys = 5;
const TString msystem[6]={"DT","RPCb","CSC","RPCf","GMT","Reco"};
Double_t pig = acos(-1.);

class GMTGhostRateNtupleizer : public L1Ntuple
{
  public :

  //constructor
  GMTGhostRateNtupleizer(std::string filename) : L1Ntuple(filename) {}
  GMTGhostRateNtupleizer() {}
  ~GMTGhostRateNtupleizer() {}

  //main function macro : arguments can be adpated to your need
  void run(Long64_t nevents);

  private :

  //your private methods can be declared here
  TNtuple *ntuple;
  void toggleBranches();
  double dphi(int iRecoMu, int iL1Mu, int iL1Sys); // calculate delta phi between iRecoMu muon and iL1Mu trigcand of type iL1Sys
  double deta(int iRecoMu, int iL1Mu, int iL1Sys); // calculate delta eta between iRecoMu muon and iL1Mu trigcand of type iL1Sys
  void fillNtuple(int recoMu, int gmtMu1, int gmtMu2, std::pair<bool, bool> muMatch, std::vector<std::string> contDict, Float_t ntupleValues[]);
  double bestL1match(int iRecoMu, int &iL1Muint, int iL1Sys, float ptcut, int exclMu); // finds best match between iRecoMu muon and trig cands of type iL1Sys
  std::pair<bool, bool> matchDiMuons(int iRecoMu1, int iRecoMu2, int &L1Mu1, int &L1Mu2, int iL1Sys, float ptcut, float dRmax); // Find best match for two reco muons
  muSysEnum whichSubsystem(int mu);
  bool glmucuts(int imu);
};


// --------------------------------------------------------------------
//                             run function
// --------------------------------------------------------------------
void GMTGhostRateNtupleizer::run(Long64_t nevents)
{

  toggleBranches();

  // Create ntuple
  std::vector<std::string> varList;
  varList.push_back("Eta");
  varList.push_back("Phi");
  varList.push_back("pT");
  std::ostringstream ntupleContStream;
  std::vector<std::string> contDict;
  ntupleContStream << "N_reco:N_GMT:Qual1_GMT:Qual2_GMT:SubsysID1_GMT:SubsysID2_GMT";
  contDict.push_back("N_reco");
  contDict.push_back("N_GMT");
  contDict.push_back("Qual1_GMT");
  contDict.push_back("Qual2_GMT");
  contDict.push_back("SubsysID1_GMT");
  contDict.push_back("SubsysID2_GMT");
  for(auto name:varList) {
    ntupleContStream << ":" << name << "1_reco:" << name << "1_GMT";
    contDict.push_back(name + "1_reco");
    contDict.push_back(name + "1_GMT");
    ntupleContStream << ":" << name << "2_GMT";
    contDict.push_back(name + "2_GMT");
  }
  std::string fname("SingleMuNtuple.root");
  TFile *out = new TFile(fname.c_str(), "RECREATE");
  out->cd();
  std::string ntupleContent(ntupleContStream.str());
  ntuple = new TNtuple("ntuple", "ntupledump", ntupleContent.c_str());

  //load TDR style
  setTDRStyle();

  //number of events to process
  if (nevents==-1 || nevents>GetEntries()) {
    nevents=GetEntries();
  }
  std::cout << nevents << " to process ..." << std::endl;

  // DEBUG:
  int fails = 0;
  int skips = 0;
  int eventsRun = 0;
  int filledTuples = 0;

  //loop over the events
  for (Long64_t i=0; i<nevents; i++) {
    //load the i-th event
    Long64_t ientry = LoadTree(i);
    if (ientry < 0) {
      break;
    }
    GetEntry(i);

    //process progress
    if(i!=0 && (i%10000)==0) {
      std::cout << "- processing event " << i << "\r" << std::flush;
    }

    //user code here
    ++eventsRun;

    if(recoMuon_->nMuons < 1) {
      ++skips;
      continue; // Only interested in singlemu events
    }

    for(int recoMu = 0; recoMu < recoMuon_->nMuons; ++recoMu) {
      if(!glmucuts(recoMu)) {
        ++fails;
        continue;
      }
      ++filledTuples;

      Float_t ntupleValues[contDict.size()];

      int gmtMu1, gmtMu2;
      std::pair<bool, bool> diMuMatch = matchDiMuons(recoMu, recoMu, gmtMu1, gmtMu2, GMT, 0, 0.1);

      fillNtuple(recoMu, gmtMu1, gmtMu2, diMuMatch, contDict, ntupleValues);

      ntuple->Fill(ntupleValues);
    }
  }
  out->Write();

  std::cout << "Number of events run over: " << eventsRun << "; filled rows in ntuple: " << filledTuples << std::endl;
  std::cout << "Number of skipped events due to less than one reco muon: " << skips << std::endl;
  std::cout << "Number of skipped muons due to glmuon cuts: " << fails << std::endl;
}

void GMTGhostRateNtupleizer::fillNtuple(int recoMu, int gmtMu1, int gmtMu2, std::pair<bool, bool> diMuMatch, std::vector<std::string> contDict, Float_t ntupleValues[])
{
  for (size_t varIt = 0; varIt < contDict.size(); ++varIt) {
    if (contDict[varIt] == "N_reco") {
      ntupleValues[varIt] = recoMuon_->nMuons;
    }
    if (contDict[varIt] == "pT1_reco") {
      ntupleValues[varIt] = recoMuon_->pt[recoMu];
    }
    if (contDict[varIt] == "Eta1_reco") {
      ntupleValues[varIt] = recoMuon_->eta[recoMu];
    }
    if (contDict[varIt] == "Phi1_reco") {
      ntupleValues[varIt] = recoMuon_->phi[recoMu];
    }

    if (contDict[varIt] == "N_GMT") {
      ntupleValues[varIt] = gmt_->N;
    }

    if (diMuMatch.first == true) {
      if (contDict[varIt] == "pT1_GMT") {
        ntupleValues[varIt] = gmt_->Pt[gmtMu1];
      }
      if (contDict[varIt] == "Eta1_GMT") {
        ntupleValues[varIt] = gmt_->Eta[gmtMu1];
      }
      if (contDict[varIt] == "Phi1_GMT") {
        ntupleValues[varIt] = gmt_->Phi[gmtMu1];
      }
      if (contDict[varIt] == "Qual1_GMT") {
        ntupleValues[varIt] = gmt_->Qual[gmtMu1];
      }
      if (contDict[varIt] == "SubsysID1_GMT") {
        ntupleValues[varIt] =  whichSubsystem(gmtMu1);
      }
    } else {
      if (contDict[varIt] == "pT1_GMT") {
        ntupleValues[varIt] = -1;
      }
      if (contDict[varIt] == "Eta1_GMT") {
        ntupleValues[varIt] = -9999;
      }
      if (contDict[varIt] == "Phi1_GMT") {
        ntupleValues[varIt] = -9999;
      }
      if (contDict[varIt] == "Qual1_GMT") {
        ntupleValues[varIt] = -1;
      }
      if (contDict[varIt] == "SubsysID1_GMT") {
        ntupleValues[varIt] =  -1;
      }
    }
    if (diMuMatch.second == true) {
      if (contDict[varIt] == "pT2_GMT") {
        ntupleValues[varIt] = gmt_->Pt[gmtMu2];
      }
      if (contDict[varIt] == "Eta2_GMT") {
        ntupleValues[varIt] = gmt_->Eta[gmtMu2];
      }
      if (contDict[varIt] == "Phi2_GMT") {
        ntupleValues[varIt] = gmt_->Phi[gmtMu2];
      }
      if (contDict[varIt] == "Qual2_GMT") {
        ntupleValues[varIt] = gmt_->Qual[gmtMu2];
      }
      if (contDict[varIt] == "SubsysID2_GMT") {
        ntupleValues[varIt] =  whichSubsystem(gmtMu2);
      }
    } else {
      if (contDict[varIt] == "pT2_GMT") {
        ntupleValues[varIt] = -1;
      }
      if (contDict[varIt] == "Eta2_GMT") {
        ntupleValues[varIt] = -9999;
      }
      if (contDict[varIt] == "Phi2_GMT") {
        ntupleValues[varIt] = -9999;
      }
      if (contDict[varIt] == "Qual2_GMT") {
        ntupleValues[varIt] = -1;
      }
      if (contDict[varIt] == "SubsysID2_GMT") {
        ntupleValues[varIt] =  -1;
      }
    }

  }
}

muSysEnum GMTGhostRateNtupleizer::whichSubsystem(int mu)
{
  muSysEnum muSys;
  if(gmt_->IdxDTBX[mu] != -1 && gmt_->IdxRPCb[mu] != -1) {
    muSys = eDTRPC;
  } else if (gmt_->IdxCSC[mu] != -1 && gmt_->IdxRPCf[mu] != -1) {
    muSys = eCSCRPC;
  } else if (gmt_->IdxDTBX[mu] != -1) {
    int Ndt= gmt_->Ndt;
    for(int k= 0; k<Ndt; ++k) {
      //Go through Bxdt and search for first muon in triggered beam crossing
      if(gmt_->Bxdt[k]==0) {
        muSys = eDT;
      }
    }
  } else if(gmt_->IdxCSC[mu] != -1) {
    int Ncsc= gmt_->Ncsc;
    for(int k= 0; k<Ncsc; ++k) {
      if(gmt_->Bxcsc[k]==0) {
        muSys = eCSC;
      }
    }
  } else if(gmt_->IdxRPCb[mu] != -1) {
    int Nrpcb= gmt_->Nrpcb;
    for(int k=0; k<Nrpcb; ++k) {
      if(gmt_->Bxrpcb[k]==0) {
        muSys = eRPCb;
      }
    }
  } else if(gmt_->IdxRPCf[mu] != -1) {
    int Nrpcf= gmt_->Nrpcf;
    for(int k=0; k<Nrpcf; ++k) {
      if(gmt_->Bxrpcf[k]==0) {
        muSys = eRPCf;
      }
    }
  } else {
    // This will crash the macro.
    muSys = none;
    std::cout << "This shouldn't happen." << std::endl;
  }

  return muSys;
}
double GMTGhostRateNtupleizer::dphi(int iRecoMu, int iL1Mu, int iL1Sys)
{
  if (recoMuon_->type[iRecoMu] != 0) {
    return -9999;
  } // not a global
  Double_t trigphi = -99999;
  if (iL1Sys==DT) trigphi=gmt_->Phidt[iL1Mu];
  if (iL1Sys==RPCb) trigphi=gmt_->Phirpcb[iL1Mu];
  if (iL1Sys==CSC ) trigphi=gmt_->Phicsc[iL1Mu];
  if (iL1Sys==RPCf) trigphi=gmt_->Phirpcf[iL1Mu];
  if (iL1Sys==GMT)trigphi=gmt_->Phi[iL1Mu];
  Double_t recophi = -99999;
  Double_t recophi2 = -99999;
  if (iL1Sys==DT || iL1Sys==RPCb) {
    recophi=recoMuon_->sa_phi_mb2[iRecoMu]-(pig/144.);
  }
  if ((iL1Sys==CSC|| iL1Sys==RPCf) && recoMuon_->eta[iRecoMu]>=0) recophi=recoMuon_->sa_phi_me2_p[iRecoMu]-(pig/144.);
  if ((iL1Sys==CSC|| iL1Sys==RPCf) && recoMuon_->eta[iRecoMu] <0) recophi=recoMuon_->sa_phi_me2_n[iRecoMu]-(pig/144.);
  if (iL1Sys==GMT) {
    recophi=recoMuon_->sa_phi_mb2[iRecoMu]-(pig/144.);
    if (recoMuon_->eta[iRecoMu]>=0) recophi2=recoMuon_->sa_phi_me2_p[iRecoMu]-(pig/144.);
    else recophi2=recoMuon_->sa_phi_me2_n[iRecoMu]-(pig/144.);
  }
  Double_t newphisep = recophi - trigphi;
  if (newphisep<-pig) newphisep = newphisep + 2*pig;
  if (newphisep> pig) newphisep = newphisep - 2*pig;
  Double_t newphisep2 = recophi2 - trigphi;
  if (newphisep2<-pig) newphisep2 = newphisep2 + 2*pig;
  if (newphisep2> pig) newphisep2 = newphisep2 - 2*pig;
  if (iL1Sys==GMT) {
    if (fabs(newphisep)>fabs(newphisep2)) newphisep = newphisep2;
  }

  if (iL1Sys==RECOPASS) newphisep = 0;
  if (newphisep>1000) return -999;
  return newphisep;
}
double GMTGhostRateNtupleizer::deta(int iRecoMu, int iL1Mu, int iL1Sys)
{
  if (recoMuon_->type[iRecoMu] != 0) {
    return -9999; // not a global
  }
  Double_t trigeta = -99999;
  if (iL1Sys==DT) trigeta=gmt_->Etadt[iL1Mu];
  if (iL1Sys==RPCb) trigeta=gmt_->Etarpcb[iL1Mu];
  if (iL1Sys==CSC ) trigeta=gmt_->Etacsc[iL1Mu];
  if (iL1Sys==RPCf) trigeta=gmt_->Etarpcf[iL1Mu];
  if (iL1Sys==GMT)trigeta=gmt_->Eta[iL1Mu];
  Double_t newetasep = recoMuon_->eta[iRecoMu]-trigeta;
  if (iL1Sys==RECOPASS) newetasep = 0;
  if (newetasep>1000) return -999;
  return newetasep;
}
double GMTGhostRateNtupleizer::bestL1match(int iRecoMu, int &iL1Mu, int iL1Sys, float ptcut, int exclMu)
{
  Double_t bestdeltar = 9999; int cand=9999;
  if (iL1Sys==DT) {
    for (Int_t im=0; im<gmt_->Ndt; im++){
      Double_t deltaphi = dphi(iRecoMu,im,iL1Sys);
      Double_t deltaeta = deta(iRecoMu,im,iL1Sys);
      Double_t deltar = sqrt(deltaphi*deltaphi+deltaeta*deltaeta);
      if (deltar<bestdeltar && gmt_->Bxdt[im]==0 && im != exclMu) {
        bestdeltar=deltar; cand=im;
      }


    }
  }
  if (iL1Sys==RPCb){
    for (Int_t im=0; im<gmt_->Nrpcb; im++){
      Double_t deltaphi = dphi(iRecoMu,im,iL1Sys);
      Double_t deltaeta = deta(iRecoMu,im,iL1Sys);
      Double_t deltar = sqrt(deltaphi*deltaphi+deltaeta*deltaeta);
      if (deltar<bestdeltar && gmt_->Bxrpcb[im]==0 && im != exclMu) {
        bestdeltar=deltar; cand=im;
      }
    }
  }
  if (iL1Sys==RPCf){
    for (Int_t im=0; im<gmt_->Nrpcf; im++){
      Double_t deltaphi = dphi(iRecoMu,im,iL1Sys);
      Double_t deltaeta = deta(iRecoMu,im,iL1Sys);
      Double_t deltar = sqrt(deltaphi*deltaphi+deltaeta*deltaeta);
      if (deltar<bestdeltar && gmt_->Bxrpcf[im]==0 && im != exclMu) {
        bestdeltar=deltar; cand=im;
      }
    }
  }
  if (iL1Sys==CSC){
    for (Int_t im=0; im<gmt_->Ncsc; im++){
      Double_t deltaphi = dphi(iRecoMu,im,iL1Sys);
      Double_t deltaeta = deta(iRecoMu,im,iL1Sys);
      Double_t deltar = sqrt(deltaphi*deltaphi+deltaeta*deltaeta);
      if (deltar<bestdeltar && gmt_->Bxcsc[im]==0 && im != exclMu) {
        bestdeltar=deltar; cand=im;
      }
    }
  }
  if (iL1Sys==GMT){
    for (Int_t im=0; im<gmt_->N; im++){
      Double_t deltaphi = dphi(iRecoMu,im,iL1Sys);
      Double_t deltaeta = deta(iRecoMu,im,iL1Sys);
      Double_t deltar = sqrt(deltaphi*deltaphi+deltaeta*deltaeta);
      if (deltar<bestdeltar && gmt_->CandBx[im]==0 && gmt_->Pt[im]>=ptcut && im != exclMu) {
        bestdeltar=deltar; cand=im;
      }
    }
  }
  if (iL1Sys==RECOPASS){
    bestdeltar=0;
    cand = -1;
  }
  iL1Mu=cand;
  return bestdeltar;
}
std::pair<bool, bool> GMTGhostRateNtupleizer::matchDiMuons(int iRecoMu1, int iRecoMu2, int &L1Mu1, int &L1Mu2, int iL1Sys, float ptcut, float dRmax)
{
  int cand11, cand12, cand21, cand22;
  float dR11, dR12, dR21, dR22, dR1, dR2, dRtot1, dRtot2;
  dR11 = bestL1match(iRecoMu1, cand11, iL1Sys, ptcut, 9999999);
  dR12 = bestL1match(iRecoMu2, cand12, iL1Sys, ptcut, cand11);
  dRtot1 = dR11+dR12;

  dR22 = bestL1match(iRecoMu2, cand22, iL1Sys, ptcut, 9999999);
  dR21 = bestL1match(iRecoMu1, cand21, iL1Sys, ptcut, cand22);
  dRtot2 = dR21+dR22;

  if(dRtot1 < dRtot2) {
    L1Mu1 = cand11;
    L1Mu2 = cand12;
    dR1 = dR11;
    dR2 = dR12;
  } else {
    L1Mu1 = cand21;
    L1Mu2 = cand22;
    dR1 = dR21;
    dR2 = dR22;
  }

  return std::pair<bool, bool>(dR1 < dRmax, dR2 < dRmax);

}

bool GMTGhostRateNtupleizer::glmucuts(int imu){

  return (recoMuon_->type[imu] == 0);

  bool condextrap=false; // true if valid extrapolations exist in relevant eta ranges
  if (recoMuon_->eta[imu] < -10 && condextrap) {
    std::cout << "zeroth" << std::endl;
    std::cout << "eta: " << recoMuon_->eta[imu] << " phi extrapolation: " << recoMuon_->sa_phi_me2_n[imu] << " phi: " << recoMuon_->phi[imu] << " pT: " << recoMuon_->pt[imu] << " pass: " << condextrap << std::endl;
  }
  if (fabs(recoMuon_->eta[imu])<1.0 && recoMuon_->sa_phi_mb2[imu]>-9999) condextrap=true;
  if (recoMuon_->eta[imu] < -10 && condextrap) {
    std::cout << "first" << std::endl;
    std::cout << fabs(recoMuon_->eta[imu]<1.0) << " " << recoMuon_->sa_phi_mb2[imu] << std::endl;
    std::cout << "eta: " << recoMuon_->eta[imu] << " phi extrapolation: " << recoMuon_->sa_phi_me2_n[imu] << " phi: " << recoMuon_->phi[imu] << " pT: " << recoMuon_->pt[imu] << " pass: " << condextrap << std::endl;
  }
  if (recoMuon_->eta[imu]>= 1.0 && recoMuon_->eta[imu]<= 1.2 && recoMuon_->sa_phi_mb2[imu]>-9999 && recoMuon_->sa_phi_me2_p[imu]>-9999) condextrap=true;
  if (recoMuon_->eta[imu] < -10 && condextrap) {
    std::cout << "second" << std::endl;
    std::cout << "eta: " << recoMuon_->eta[imu] << " phi extrapolation: " << recoMuon_->sa_phi_me2_n[imu] << " phi: " << recoMuon_->phi[imu] << " pT: " << recoMuon_->pt[imu] << " pass: " << condextrap << std::endl;
  }
  if (recoMuon_->eta[imu]<=-1.0 && recoMuon_->eta[imu]>=-1.2 && recoMuon_->sa_phi_mb2[imu]>-9999 && recoMuon_->sa_phi_me2_n[imu]>-9999) condextrap=true;
  if (recoMuon_->eta[imu] < -10 && condextrap) {
    std::cout << "third" << std::endl;
    std::cout << "eta: " << recoMuon_->eta[imu] << " phi extrapolation: " << recoMuon_->sa_phi_me2_n[imu] << " phi: " << recoMuon_->phi[imu] << " pT: " << recoMuon_->pt[imu] << " pass: " << condextrap << std::endl;
  }
  if (recoMuon_->eta[imu]>1.2 && recoMuon_->sa_phi_me2_p[imu]>-9999) condextrap=true;
  if (recoMuon_->eta[imu] < -10 && condextrap) {
    std::cout << "fourth" << std::endl;
    std::cout << "eta: " << recoMuon_->eta[imu] << " phi extrapolation: " << recoMuon_->sa_phi_me2_n[imu] << " phi: " << recoMuon_->phi[imu] << " pT: " << recoMuon_->pt[imu] << " pass: " << condextrap << std::endl;
  }
  if (recoMuon_->eta[imu]<-1.2 && recoMuon_->eta[imu]>= -6 && recoMuon_->sa_phi_me2_n[imu]>-9999) condextrap=true;
  if (recoMuon_->eta[imu] < -10 && condextrap) {
    std::cout << "fifth" << std::endl;
    std::cout << "eta: " << recoMuon_->eta[imu] << " phi extrapolation: " << recoMuon_->sa_phi_me2_n[imu] << " phi: " << recoMuon_->phi[imu] << " pT: " << recoMuon_->pt[imu] << " pass: " << condextrap << std::endl;
  }

  bool condpca = false; // true if muon has small extrapolation to IP
  Double_t rr = sqrt(recoMuon_->tr_imp_point_x[imu]*recoMuon_->tr_imp_point_x[imu]+
  recoMuon_->tr_imp_point_y[imu]*recoMuon_->tr_imp_point_y[imu]);
  Double_t zz = recoMuon_->tr_imp_point_z[imu];
  if (rr<.5 && fabs(zz)<10.) condpca=true;

  bool hastrackertrack = ((recoMuon_->howmanytypes[imu]>>4)&0x3); // this is new, means has tight muon (?)

  bool cond = (condextrap && condpca && hastrackertrack);
  // DEBUG
  if (recoMuon_->eta[imu] < -10 && cond) {
    std::cout << "eta: " << recoMuon_->eta[imu] << " phi extrapolation: " << recoMuon_->sa_phi_me2_n[imu] << " phi: " << recoMuon_->phi[imu] << " pT: " << recoMuon_->pt[imu] << " pass: " << cond << std::endl;
  }
  return cond;
}

void GMTGhostRateNtupleizer::toggleBranches()
{
      //Select only needed branches:
    fChain->SetBranchStatus("*",0);

    fChain->SetBranchStatus("Ndt",1);
    fChain->SetBranchStatus("Bxdt",1);

    fChain->SetBranchStatus("Ncsc",1);
    fChain->SetBranchStatus("Bxcsc",1);

    fChain->SetBranchStatus("Nrpcb",1);
    fChain->SetBranchStatus("Bxrpcb",1);

    fChain->SetBranchStatus("Nrpcf",1);
    fChain->SetBranchStatus("Bxrpcf",1);

    fChain->SetBranchStatus("N",1);
    fChain->SetBranchStatus("CandBx",1);
    fChain->SetBranchStatus("Qual",1);
    fChain->SetBranchStatus("Eta",1);
    fChain->SetBranchStatus("Phi",1);
    fChain->SetBranchStatus("Pt",1);
    fChain->SetBranchStatus("Cha",1);
    fChain->SetBranchStatus("IdxDTBX",1);
    fChain->SetBranchStatus("IdxRPCb",1);
    fChain->SetBranchStatus("IdxCSC",1);
    fChain->SetBranchStatus("IdxRPCf",1);

    fChain->SetBranchStatus("nMuons",1);
    fChain->SetBranchStatus("pt",1);
    fChain->SetBranchStatus("eta",1);
    fChain->SetBranchStatus("phi",1);
    fChain->SetBranchStatus("ch",1);
    fChain->SetBranchStatus("howmanytypes",1);
    fChain->SetBranchStatus("type",1);
    fChain->SetBranchStatus("tr_imp_point_x",1);
    fChain->SetBranchStatus("tr_imp_point_y",1);
    fChain->SetBranchStatus("tr_imp_point_z",1);
    fChain->SetBranchStatus("sa_phi_mb2",1);
    fChain->SetBranchStatus("sa_phi_me2_p",1);
    fChain->SetBranchStatus("sa_phi_me2_n",1);

    fChain->SetBranchStatus("hlt",1);
}
