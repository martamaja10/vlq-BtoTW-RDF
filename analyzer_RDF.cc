
// --------------------------------------------------------------------------------------- //
// Implimentation of RDataFrame in C++.					                   //
// Comments on creating a singly produced VLQ search			                   //
// To Run on Command Line:   root -l callRDF.C\(\"Muon(OR)Electron\",\"testNumber\"\,\"root://cmsxrootd.fnal.gov//store/...file.root\")      //
// --------------------------------------------------------------------------------------- //

#define rdf_cxx
#include "analyzer_RDF.h"
#include "lumiMask.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include <TFile.h>
#include <TChain.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH3.h>
#include <algorithm> 
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TVector2.h>
#include <TRandom3.h>
#include <sstream>
#include <chrono> // for high_resolution_clock


using namespace std;
using namespace ROOT::VecOps;


void rdf::analyzer_RDF(TString testNum, TString jesvar)
{
  ROOT::EnableImplicitMT(4);

  cout << "Number of threads: " << ROOT::GetThreadPoolSize() << endl;
  TStopwatch time;
  time.Start();
  string sample = this->sample;
  int samplebin = this->samplebin;
  string year = this->year;
  bool isMC = this->isMC;

  bool debug = false;
  
  cout << "Sample in cc: " << sample << ", bin # " << samplebin << endl;
  cout << "Year in cc: " << year << endl;
  cout << "isMC? " << isMC << ", jesvar = " << jesvar << endl;
  if(!isMC) cout << "Data era = " << era << ", for jec " << jecera << endl;
  else cout << "MC extension tag (blank or ext) = " << era << endl;

  // -------------------------------------------------------
  //               Self-derived corrections
  // -------------------------------------------------------
  
  // polynomials defined in the .h
  TF1 *poly2 = this->poly2;
  TF1 *poly2U = this->poly2U;
  TF1 *poly2D = this->poly2D;
  TF1 *polyHT = this->polyHT;
  TF1 *polyHTU = this->polyHTU;
  TF1 *polyHTD = this->polyHTD;
  auto wjetHTpoly = [poly2,poly2U,poly2D](float &LHE_HT){RVec<double> sf = {poly2->Eval(LHE_HT),poly2U->Eval(LHE_HT),poly2D->Eval(LHE_HT)}; return sf;};
  auto topHTpoly = [polyHT,polyHTU,polyHTD](float &AK4HT){RVec<double> sf = {polyHT->Eval(AK4HT),polyHTU->Eval(AK4HT),polyHTD->Eval(AK4HT)}; return sf;};

  // Lepton scale factors not in correctionLib
  std::vector<std::vector<float>> elecidsfs = this->elIDSF;
  std::vector<std::vector<float>> elecidsfuncs = this->elIDSFUnc;  
  std::vector<std::vector<float>> elecisosfs = this->elISOSF;
  float elecisosfunc = this->elISOSFUnc;  
  std::vector<std::vector<float>> muonisosfs = this->muISOSF;
  std::vector<std::vector<float>> elechltsfs = this->elHLTSF;
  std::vector<std::vector<float>> elechltuncs = this->elHLTSFUnc;

  std::vector<float> elid_pts = {50,100,200,99999}; // from susy double disco UMN paper, UL
  std::vector<float> elid_etas = {-2.5, -2.0, -1.566, -1.442, -0.8, 0.0, 0.8, 1.442, 1.566, 2.0, 2.5};
  std::vector<float> muiso_pts = {50,60,99999};
  std::vector<float> muiso_etas = {0,0.9,1.2,2.1,2.4};
  std::vector<float> elhlt_pts = {50,60,70,80,100,200,300,99999}; // from top charge asym paper, EOY
  std::vector<float> elhlt_etas = {0,0.8,1.442,1.566,2.0,2.5};
  float muonisosfunc = 0.002;

  // DeepJet Loose efficiencies
  std::vector<float> btagpts = this->btagpts;
  std::vector<std::vector<float>> btageffs = this->btageffs;

  // ParticleNet Top and W efficiencies
  std::vector<float> pnetpts = this->pnetpts;
  std::vector<std::vector<float>> pnet_t_t = this->pnet_t_t;
  std::vector<std::vector<float>> pnet_t_W = this->pnet_t_W;
  std::vector<std::vector<float>> pnet_W_t = this->pnet_W_t;
  std::vector<std::vector<float>> pnet_W_W = this->pnet_W_W;
  std::vector<std::vector<float>> pnet_J_t = this->pnet_J_t;
  std::vector<std::vector<float>> pnet_J_W = this->pnet_J_W;
  
  // -------------------------------------------------------
  //               correctionLib corrections
  // -------------------------------------------------------
  
  std::string yrstr, yr, jecyr, jeryr, jecver;
  float deepjetL;
  string mutrig = "TkMu50";
  if(year == "2016APV") {deepjetL = 0.0508; yrstr = "2016preVFP"; yr = "16"; jecyr = "UL16APV"; jeryr = "Summer20UL16APV_JRV3"; jecver = "V7";}
  else if(year == "2016") {deepjetL = 0.0480; yrstr = "2016postVFP"; yr = "16"; jecyr = "UL16"; jeryr = "Summer20UL16_JRV3"; jecver = "V7";}
  else if(year == "2017") {mutrig = "OldMu100_or_TkMu100"; deepjetL = 0.0532; yrstr = "2017"; yr = "17"; jecyr = "UL17"; jeryr = "Summer19UL17_JRV2"; jecver = "V5";}
  else if(year == "2018") {mutrig = "OldMu100_or_TkMu100"; deepjetL = 0.0490; yrstr = "2018"; yr = "18"; jecyr = "UL18"; jeryr = "Summer19UL18_JRV2"; jecver = "V5";}
  else std::cout << "ERROR: Can't parse the year to assign correctionLib json files. Expected 2016, 2016APV, 2017, or 2018. Got: " << year << std::endl;
  
  auto pnetWPs = [year](const RVec<float> &dnnT, const RVec<float> &dnnW){
    float wpT, wpW;
    if(year == "2016APV"){wpT = 0.490; wpW = 0.677;}
    else if(year == "2016"){wpT = 0.495; wpW = 0.668;}
    else if(year == "2017"){wpT = 0.581; wpW = 0.709;}
    else{wpT = 0.580; wpW = 0.700;}
    RVec<int> tag;
    for(int i=0; i<dnnT.size(); i++){
      int tagval = 0; // will stay 0 if not t, W
      if(dnnT[i] > wpT) tagval += 1; // 1 if only top, 3 if top or W
      if(dnnW[i] > wpW) tagval += 2; // 2 if only W, 3 if top or W
      tag.push_back(tagval);
    }
    return tag;
  };

  // -------------------------------------------------------
  //               Open Dataframe + MET Filters
  // -------------------------------------------------------

  auto rdf_input = ROOT::RDataFrame("Events", files); // Initial data
  
  auto METgeneralFilters = rdf_input.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1", "MET Filters")
    .Filter("nJet > 0 && nFatJet > 0", "Event has > 1 AK4 and > 1 AK8");

  auto truth = METgeneralFilters;
  
  // ---------------------------------------------------------
  //                    LEPTON Definitions
  // ---------------------------------------------------------
  
  string elHEMcut = "";
  if(year == "2018") elHEMcut = " && (Electron_eta > -1.479 || (Electron_phi < -1.57 || Electron_phi > -0.87))";
  auto LepDefs = truth.Define("Electron_cutBasedIdNoIso_tight", "Electron_cutBasedIdNoIso_tight(nElectron, Electron_vidNestedWPBitmap)")
    .Define("TPassMu", "abs(Muon_eta)<2.4 && Muon_mediumId==1 && Muon_miniIsoId>=3 && abs(Muon_dz) < 0.5 && Muon_dxy < 0.2")
    .Define("TPassEl", Form("(abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_cutBasedIdNoIso_tight==1 && Electron_miniPFRelIso_all<0.1%s",elHEMcut.c_str()))
    .Define("VetoMu", "TPassMu && (Muon_pt>25)")
    .Define("VetoEl", "TPassEl && (Electron_pt>25)")
    .Define("SignalIsoMu", "TPassMu && (Muon_pt>=55)")
    .Define("SignalIsoEl", "TPassEl && (Electron_pt>=55)")
    .Define("nVetoLep", "(int) (Sum(VetoMu)+Sum(VetoEl))")
    .Define("SMuon_pt", "Muon_pt[SignalIsoMu == true]")
    .Define("SMuon_eta", "Muon_eta[SignalIsoMu == true]")
    .Define("SMuon_phi", "Muon_phi[SignalIsoMu == true]")
    .Define("SMuon_mass", "Muon_mass[SignalIsoMu == true]")
    .Define("SElectron_pt", "Electron_pt[SignalIsoEl == true]")
    .Define("SElectron_eta", "Electron_eta[SignalIsoEl == true]")
    .Define("SElectron_phi", "Electron_phi[SignalIsoEl == true]")
    .Define("SElectron_mass", "Electron_mass[SignalIsoEl == true]")
    .Define("Muon_P4", "fVectorConstructor(Muon_pt,Muon_eta,Muon_phi,Muon_mass)")
    .Define("SMuon_P4", "fVectorConstructor(SMuon_pt,SMuon_eta,SMuon_phi,SMuon_mass)")
    .Define("SElectron_P4", "fVectorConstructor(SElectron_pt,SElectron_eta,SElectron_phi,SElectron_mass)")
    .Define("SMuon_jetIdx", "Muon_jetIdx[SignalIsoMu == true]")
    .Define("SElectron_jetIdx", "Electron_jetIdx[SignalIsoEl]")
    .Define("nSignalIsoMu", "(int) Sum(SignalIsoMu)")
    .Define("nSignalIsoEl", "(int) Sum(SignalIsoEl)")
    .Define("VetoIsoMu", "(VetoMu == true && Muon_pt < 55)")
    .Define("VetoIsoEl", "(VetoEl == true && Electron_pt < 55)")
    .Define("nVetoIsoLep", "(int) (Sum(VetoIsoMu)+Sum(VetoIsoEl))");

  // --------------------------------------------------------
  // 		      LEPTON SELECTION
  // --------------------------------------------------------
  
  string tkmutrig = " || HLT_OldMu100 || HLT_TkMu100";
  string eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon200";
  if(year == "2017" && era == "B"){ tkmutrig = ""; eltrig = "HLT_Ele35_WPTight_Gsf || HLT_Photon200";}
  if(year == "2016" or year == "2016APV"){ tkmutrig = " || HLT_TkMu50"; eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon175";}
  if(year == "2016APV" and (era == "A" || era == "B")){ tkmutrig = ""; eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon175";}
  auto LepSelect = LepDefs.Define("isMu", Form("(nMuon>0) && (HLT_Mu50%s) && (nSignalIsoMu==1) && (nVetoIsoLep==0) && (nElectron == 0 || nSignalIsoEl == 0)",tkmutrig.c_str()))
    .Define("isEl", Form("(nElectron>0) && (%s) && (nSignalIsoEl==1) && (nVetoIsoLep==0) && (nMuon == 0 || nSignalIsoMu == 0)",eltrig.c_str()))
    .Filter("isMu || isEl", "Event is either muon or electron");
  
  auto LepAssign = LepSelect.Define("assignleps", "assign_leps(isMu,isEl,SignalIsoMu,SignalIsoEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_miniPFRelIso_all,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_miniPFRelIso_all)")
    .Define("lepton_pt", "assignleps[0]")
    .Define("lepton_eta", "assignleps[1]")
    .Define("lepton_phi", "assignleps[2]")
    .Define("lepton_mass", "assignleps[3]")
    .Define("lepton_miniIso", "assignleps[4]");
  
  // ---------------------------------------------------------
  //                    MET Selection
  // ---------------------------------------------------------
  
  
  auto METSelect = LepAssign.Filter("MET_pt > 60", "Pass corr MET > 60");
  
  // ---------------------------------------------------------
  // 	  HT Calculation and N Jets cuts
  // ---------------------------------------------------------
  
  auto JetSelect = METSelect.Define("DR_lepJets","DeltaR_VecAndFloat(Jet_eta,Jet_phi,lepton_eta,lepton_phi)")
    .Define("ptrel_lepJets","ptRel(Jet_pt,Jet_eta,Jet_phi,Jet_mass,lepton_pt,lepton_eta,lepton_phi,lepton_mass)")
    .Define("goodJets", "Jet_pt > 30 && abs(Jet_eta) < 2.5 && Jet_jetId > 1 && (DR_lepJets > 0.4 || ptrel_lepJets > 20)")
    .Define("gcJet_HT","Sum(Jet_pt[goodJets == true])")	
    .Define("DR_lepFatJets","DeltaR_VecAndFloat(FatJet_eta,FatJet_phi,lepton_eta,lepton_phi)")
    .Define("goodFatJets", "FatJet_pt > 200 && abs(FatJet_eta) < 2.5 && FatJet_jetId > 1 && (DR_lepFatJets > 0.8)") 
    .Define("NFatJets", "(int) Sum(goodFatJets)")
    .Define("NOS_gcFatJets","(int) Sum(DR_lepFatJets[goodFatJets == true] > TMath::Pi()/2)")
    .Filter("gcJet_HT > 250","Pass HT > 250")				
    .Filter("NFatJets > 0","Pass N good central AK8 > 0")
    .Filter("NOS_gcFatJets > 0","Pass N good central other side AK8 > 0");

  // ---------------------------------------------------------
  // 	  Jet pt ordering, counting, lepton association
  // ---------------------------------------------------------

  auto JetVars = JetSelect.Define("gcHTCorr_top", topHTpoly, {"gcJet_HT"})
    .Define("NJets_central", "(int) Sum(goodJets)")
    .Define("gcJet_pt_unsort", "Jet_pt[goodJets == true]")
    .Define("gcJet_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcJet_pt_unsort))")
    .Define("gcJet_pt","reorder(gcJet_pt_unsort,gcJet_ptargsort)")
    .Define("gcJet_eta", "reorder(Jet_eta[goodJets == true],gcJet_ptargsort)")
    .Define("gcJet_phi", "reorder(Jet_phi[goodJets == true],gcJet_ptargsort)")
    .Define("gcJet_mass", "reorder(Jet_mass[goodJets == true],gcJet_ptargsort)")
    .Define("gcJet_DeepFlav", "reorder(Jet_btagDeepFlavB[goodJets == true],gcJet_ptargsort)")
    .Define("gcJet_DeepFlavL", Form("gcJet_DeepFlav > %f",deepjetL)) 
    .Define("NJets_DeepFlavL", "(int) Sum(gcJet_DeepFlavL)")
    .Define("DR_gcJets_central","reorder(DR_lepJets[goodJets == true],gcJet_ptargsort)")
    .Define("minDR_lepJets","ROOT::VecOps::Min(DR_gcJets_central)")
    .Define("ptrel_atMinDR_lepJets","reorder(ptrel_lepJets[goodJets == true],gcJet_ptargsort)[ROOT::VecOps::ArgMin(DR_gcJets_central)]")
    .Define("OS_gcJets","DR_gcJets_central > TMath::Pi()/2")
    .Define("SS_gcJets","DR_gcJets_central <= TMath::Pi()/2")
    .Define("NOS_gcJets_central","(int) Sum(OS_gcJets)")
    .Define("NSS_gcJets_central","(int) Sum(SS_gcJets)")
    .Define("gcOSJet_pt","gcJet_pt[OS_gcJets == true]")
    .Define("gcOSJet_eta","gcJet_eta[OS_gcJets == true]")
    .Define("gcOSJet_phi","gcJet_phi[OS_gcJets == true]")
    .Define("gcOSJet_mass","gcJet_mass[OS_gcJets == true]")
    .Define("gcOSJet_DeepFlavL","gcJet_DeepFlavL[OS_gcJets == true]")
    .Define("gcSSJet_pt","gcJet_pt[SS_gcJets == true]")
    .Define("gcSSJet_eta","gcJet_eta[SS_gcJets == true]")
    .Define("gcSSJet_phi","gcJet_phi[SS_gcJets == true]")
    .Define("gcSSJet_mass","gcJet_mass[SS_gcJets == true]")
    .Define("gcSSJet_DeepFlavL","gcJet_DeepFlavL[SS_gcJets == true]")
    .Define("NOS_gcJets_DeepFlavL","(int) Sum(gcOSJet_DeepFlavL)")
    .Define("NSS_gcJets_DeepFlavL","(int) Sum(gcSSJet_DeepFlavL)");

  auto ForwardJetVars = JetVars.Define("goodcleanForwardJets", "Jet_pt > 30 && abs(Jet_eta) >= 2.5 && Jet_jetId > 1")
    .Define("NJets_forward", "(int) Sum(goodcleanForwardJets)")
    .Define("gcforwJet_pt_unsort", "Jet_pt[goodcleanForwardJets == true]")
    .Define("gcforwJet_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcforwJet_pt_unsort))")
    .Define("gcforwJet_pt", "reorder(gcforwJet_pt_unsort,gcforwJet_ptargsort)")
    .Define("gcforwJet_eta", "reorder(Jet_eta[goodcleanForwardJets == true],gcforwJet_ptargsort)")
    .Define("gcforwJet_phi", "reorder(Jet_phi[goodcleanForwardJets == true],gcforwJet_ptargsort)")
    .Define("gcforwJet_mass", "reorder(Jet_mass[goodcleanForwardJets == true],gcforwJet_ptargsort)")
    .Define("gcforwJet_DeepFlav", "reorder(Jet_btagDeepFlavB[goodcleanForwardJets == true],gcforwJet_ptargsort)");

  auto FatJetVars = ForwardJetVars.Define("gcFatJet_pt_unsort", "FatJet_pt[goodFatJets == true]")
    .Define("gcFatJet_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcFatJet_pt_unsort))")
    .Define("gcFatJet_pt","reorder(gcFatJet_pt_unsort,gcFatJet_ptargsort)")
    .Define("gcFatJet_eta", "reorder(FatJet_eta[goodFatJets == true],gcFatJet_ptargsort)")
    .Define("gcFatJet_phi", "reorder(FatJet_phi[goodFatJets == true],gcFatJet_ptargsort)")
    .Define("gcFatJet_mass", "reorder(FatJet_mass[goodFatJets == true],gcFatJet_ptargsort)")
    .Define("gcFatJet_sdmass", "reorder(FatJet_msoftdrop[goodFatJets == true],gcFatJet_ptargsort)")
    .Define("DR_gcFatJets", "reorder(DR_lepFatJets[goodFatJets == true],gcFatJet_ptargsort)")
    .Define("minDR_lepFatJets","ROOT::VecOps::Min(DR_gcFatJets)")
    .Define("ptrel_atMinDR_lepFatJets","ptRel(gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass,lepton_pt,lepton_eta,lepton_phi,lepton_mass)[ROOT::VecOps::ArgMin(DR_gcFatJets)]")
    .Define("SS_gcFatJets","DR_gcFatJets <= TMath::Pi()/2")
    .Define("NSS_gcFatJets","(int) Sum(SS_gcFatJets)")
    .Define("OS_gcFatJets","DR_gcFatJets > TMath::Pi()/2")
    .Define("gcOSFatJet_pt","gcFatJet_pt[OS_gcFatJets == true]")
    .Define("gcOSFatJet_eta","gcFatJet_eta[OS_gcFatJets == true]")
    .Define("gcOSFatJet_phi","gcFatJet_phi[OS_gcFatJets == true]")
    .Define("gcOSFatJet_mass","gcFatJet_mass[OS_gcFatJets == true]")
    .Define("gcOSFatJet_sdmass","gcFatJet_sdmass[OS_gcFatJets == true]");
  
  // ---------------------------------------------------------
  // 		Add scale factors and MC jet-based calcs
  // ---------------------------------------------------------
  auto scaleFactors = FatJetVars;
  
  // ---------------------------------------------------------
  // 		JET Tagging variables
  // ---------------------------------------------------------

  auto Taggers = scaleFactors.Define("lepton_lv", "lvConstructor(lepton_pt,lepton_eta,lepton_phi,lepton_mass)")
    .Define("gcJet_ST", "gcJet_HT + lepton_pt + MET_pt")
    .Define("gcFatJet_pNetJ", "reorder(FatJet_particleNet_QCD[goodFatJets == true],gcFatJet_ptargsort)")
    .Define("gcFatJet_pNetTvsQCD", "reorder(FatJet_particleNet_TvsQCD[goodFatJets == true],gcFatJet_ptargsort)")
    .Define("gcFatJet_pNetWvsQCD", "reorder(FatJet_particleNet_WvsQCD[goodFatJets == true],gcFatJet_ptargsort)")
    .Define("gcOSFatJet_pNetJ", "gcFatJet_pNetJ[OS_gcFatJets == true]") 
    .Define("gcOSFatJet_pNetTvsQCD", "gcFatJet_pNetTvsQCD[OS_gcFatJets == true]") 
    .Define("gcOSFatJet_pNetWvsQCD", "gcFatJet_pNetWvsQCD[OS_gcFatJets == true]")
    .Define("gcFatJet_pNetT", "(gcFatJet_pNetTvsQCD * gcFatJet_pNetJ) / (1 - gcFatJet_pNetTvsQCD)")
    .Define("gcFatJet_pNetW", "(gcFatJet_pNetWvsQCD * gcFatJet_pNetJ) / (1 - gcFatJet_pNetWvsQCD)")
    .Define("gcOSFatJet_pNetT", "gcFatJet_pNetT[OS_gcFatJets == true]")
    .Define("gcOSFatJet_pNetW", "gcFatJet_pNetW[OS_gcFatJets == true]")
    .Define("gcFatJet_pNetTag", pnetWPs, {"gcFatJet_pNetTvsQCD", "gcFatJet_pNetWvsQCD"})
    .Define("gcOSFatJet_pNetTag", "gcFatJet_pNetTag[OS_gcFatJets==true]")
    .Define("gcFatJet_nJ", "Sum(gcFatJet_pNetTag == 0)")
    .Define("gcFatJet_nT", "Sum(gcFatJet_pNetTag == 1)")
    .Define("gcFatJet_nW", "Sum(gcFatJet_pNetTag == 2)")
    .Define("gcFatJet_tau21", "reorder((FatJet_tau2 / FatJet_tau1)[goodFatJets == true],gcFatJet_ptargsort)")
    .Define("gcOSFatJet_tau21", "gcFatJet_tau21[OS_gcFatJets == true]")
    .Define("gcFatJet_tau32", "reorder((FatJet_tau3 / FatJet_tau2)[goodFatJets == true],gcFatJet_ptargsort)")
    .Define("gcOSFatJet_tau32", "gcFatJet_tau32[OS_gcFatJets == true]")
    .Define("minDR_leadAK8otherAK8", "minDR_leadJetOtherJet_calc(gcFatJet_eta,gcFatJet_phi)")
    .Define("minDR_leadAK4otherAK4", "minDR_leadJetOtherJet_calc(gcJet_eta,gcJet_phi)")
    .Define("minDR_AK8s_discrete","std::floor(minDR_leadAK8otherAK8/0.5)")
    .Define("minDR_AK4s_discrete","std::floor(minDR_leadAK4otherAK4/0.5)");

  // ---------------------------------------------------------
  // 		W, top, and B reconstruction
  // ---------------------------------------------------------

  auto Reconstruction = Taggers.Define("W_lv", "W_reco(MET_pt,MET_phi,lepton_lv)")
    .Define("W_pt", "W_lv.Pt()")
    .Define("W_eta", "W_lv.Eta()")
    .Define("W_phi", "W_lv.Phi()")
    .Define("W_mass", "W_lv.M()")
    .Define("W_MT", "sqrt(2*lepton_pt*MET_pt*(1-cos(lepton_phi - MET_phi)))")
    .Define("minMlj_output", "minM_lep_jet_calc(gcJet_pt, gcJet_eta, gcJet_phi, gcJet_mass, lepton_lv)")
    .Define("DR_W_lep", "W_lv.DeltaR(lepton_lv)")
    .Define("minM_lep_Jet", "minMlj_output[0]")
    .Define("minM_lep_Jet_jetID", "(int) minMlj_output[1]")
    .Define("minM_lep_Jet_TorW", "isLeptonic_X(minM_lep_Jet)")
    .Define("t_output", "t_reco(minM_lep_Jet_TorW,gcJet_DeepFlavL,gcJet_pt,gcJet_eta,gcJet_phi,gcJet_mass,W_lv,minM_lep_Jet_jetID,NSS_gcJets_DeepFlavL,gcSSJet_DeepFlavL,gcSSJet_pt,gcSSJet_eta,gcSSJet_phi,gcSSJet_mass)")
    .Define("t_pt_SSb", "t_output[5]")
    .Define("t_eta_SSb", "t_output[6]")
    .Define("t_phi_SSb", "t_output[7]")
    .Define("t_mass_SSb", "t_output[8]")
    .Define("DR_W_b_SSb", "t_output[9]")
    .Define("Bprime_output", "BPrime_reco_new(W_lv,NOS_gcJets_DeepFlavL,NSS_gcJets_DeepFlavL,gcSSJet_DeepFlavL,gcOSJet_DeepFlavL,gcOSFatJet_pt,gcOSFatJet_eta,gcOSFatJet_phi,gcOSFatJet_mass,gcOSFatJet_pNetTag,gcOSJet_pt,gcOSJet_eta,gcOSJet_phi,gcOSJet_mass,gcSSJet_pt,gcSSJet_eta,gcSSJet_phi,gcSSJet_mass)")
    .Define("Bprime_mass", "Bprime_output[0]")
    .Define("Bprime_pt", "Bprime_output[1]")
    .Define("Bprime_eta", "Bprime_output[2]")
    .Define("Bprime_phi", "Bprime_output[3]")
    .Define("Bprime_DR", "Bprime_output[4]")
    .Define("Bprime_ptbal", "Bprime_output[5]")
    .Define("Bprime_chi2", "Bprime_output[6]")
    .Define("Bdecay_obs", "Bprime_output[7]")
    .Define("Bprime_chi2_discrete", "Bprime_output[8]");
  
  // -------------------------------------------------
  // 		Save Snapshot to file
  // -------------------------------------------------
  
  cout << "-------------------------------------------------" << endl
       << ">>> Saving " << sample << " Snapshot..." << endl;
  TString finalFile = "RDF_" + sample + era + "_" + year + "_" + testNum.Data() + ".root";
  const char *stdfinalFile = finalFile;
  
  auto ColNames = Reconstruction.GetColumnNames();
  vector<string> snapCol;
  int i = 0;
  for (auto &&ColName : ColNames)
    {
      TString colName = ColName;
      if(colName.Contains("P4") || colName.Contains("Jets") || colName.Contains("FatJets") || colName.Contains("cleanMets") || colName.Contains("Dummy")) continue;
      if(colName.Contains("LHE") && !colName.Contains("Weight") && colName != "LHE_HT" && colName != "LHE_Vpt" && colName != "gcHTCorr_WjetLHE") continue;
      if(colName.BeginsWith("Muon") && !colName.Contains("_tightId") && !colName.Contains("_isPF") && !colName.Contains("tunep") && !colName.Contains("genPartFlav")) continue;
      if(colName.BeginsWith("Electron") && !colName.Contains("genPartFlav")) continue;
      if(colName.BeginsWith("Jet") && !colName.Contains("rawFactor")) continue;
      if(colName.BeginsWith("FatJet") && !colName.Contains("rawFactor")) continue;
      if(colName.BeginsWith("PPS") || colName.BeginsWith("Proton") || colName.BeginsWith("L1_")) continue;
      if(colName.BeginsWith("Gen") || colName.BeginsWith("Soft") || colName.BeginsWith("fixed")) continue;
      if(colName.BeginsWith("Sub") || colName.BeginsWith("RawPuppi") || colName.BeginsWith("Calo") || colName.BeginsWith("Chs")) continue;
      if(colName.BeginsWith("Corr") || colName.BeginsWith("Fsr") || colName.BeginsWith("Iso") || colName.BeginsWith("Tau")) continue;
      if(colName.BeginsWith("SV") || colName.BeginsWith("Puppi") || colName.BeginsWith("Photon") || colName.BeginsWith("Low")) continue;
      if(colName.BeginsWith("HLT") || colName.BeginsWith("HT") || colName.BeginsWith("boosted") || colName.BeginsWith("Deep")) continue;
      if(colName.BeginsWith("Flag") || colName == "Bprime_gen_info" || colName == "t_gen_info" || colName == "W_gen_info" || colName == "metxyoutput") continue;
      if(colName == "assignleps" || colName == "pnetoutput" || colName == "t_output" || colName == "Bprime_output" || colName.BeginsWith("Other")) continue;
      if(colName.BeginsWith("PS") || colName.BeginsWith("PV") || colName.BeginsWith("Tk") || colName.BeginsWith("Trig")) continue;
      if(colName.BeginsWith("nCorr") || colName.BeginsWith("nFsr")) continue;
      if(colName.BeginsWith("nGen") || colName.BeginsWith("nIso") || colName.BeginsWith("nLow")) continue;
      if(colName.BeginsWith("nOther") || colName.BeginsWith("nPS") || colName.BeginsWith("nPhoton")) continue;
      if(colName.BeginsWith("nSV") || colName.BeginsWith("nSub") || colName.BeginsWith("nTau") || colName.BeginsWith("nTrig")) continue;
      if(colName.BeginsWith("nboosted")) continue;

      string name = colName.Data();
      snapCol.push_back(name);
      i++;
    }
  cout << "Number of Columns in Snapshot: " << i << endl;
  
  ROOT::RDF::RSnapshotOptions opts;
  if(jesvar != "Nominal") opts.fMode = "UPDATE";
  Reconstruction.Snapshot("Events_"+jesvar, stdfinalFile, snapCol, opts);
  cout << "Output File: " << finalFile << endl
       << "-------------------------------------------------" << endl;
  
  time.Stop();
  time.Print();
  
  cout << "Cut statistics:" << endl;
  Reconstruction.Report()->Print();

  std::cout << "Got past the report" << std::endl;
  
  if(jesvar == "Nominal"){
    cout << "Adding Counter tree to the file:" << endl;
    auto rdf_runs = ROOT::RDataFrame("Runs", files);
    opts.fMode = "UPDATE";
    rdf_runs.Snapshot("Runs", stdfinalFile, rdf_runs.GetColumnNames(), opts);
  }

  cout << "Done!" << endl;
}