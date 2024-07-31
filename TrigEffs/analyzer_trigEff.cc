
// --------------------------------------------------------------------------------------- //
// Implimentation of RDataFrame in C++.					                   //
// Comments on creating a singly produced VLQ search			                   //
// To Run on Command Line:   root -l callRDF.C\(\"Muon(OR)Electron\",\"testNumber\"\,\"root://cmsxrootd.fnal.gov//store/...file.root\")      //
// --------------------------------------------------------------------------------------- //

#define rdf_cxx
#include "analyzer_trigEff.h"
#include "../lumiMask.h"
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
#include <TH2.h>
#include <algorithm> // std::sort
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TVector2.h>
#include <TRandom3.h>
#include <sstream>
#include <chrono> // for high_resolution_clock
//#include "../correctionlib/include/correction.h"

using namespace std;
using namespace ROOT::VecOps;

void rdf::analyzer_trigEff(TString testNum, TString jesvar)
{
  //  ROOT::EnableImplicitMT();
  TStopwatch time;
  time.Start();
  string sample = this->sample;
  string year = this->year;
  bool isMC = this->isMC;
  bool isSM = this->isSM;
  bool isSE = this->isSE;
  bool debug = false;

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  
  cout << "Sample in cc: " << sample << endl;
  cout << "Year in cc: " << year << endl;
  cout << "isMC? " << isMC << ", jesvar = " << jesvar << endl;
  if(!isMC) cout << "Data era = " << era << ", for jec " << jecera << endl;
  else cout << "MC extension tag (blank or ext) = " << era << endl;

  // -------------------------------------------------------
  //               Golden JSON
  // -------------------------------------------------------
  std::string jsonfile;
  if(year == "2016" or year == "2016APV") jsonfile = "../Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt";
  else if(year == "2017") jsonfile = "../Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt";
  else if(year == "2018") jsonfile = "../Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt";
  else std::cout << "ERROR: Can't parse the year to assign a golden json file. Expected 2016, 2016APV, 2017, or 2018. Got: " << year << std::endl;
  const auto myLumiMask = lumiMask::fromJSON(jsonfile);
  //  std::cout << "Testing the JSON! Known good run/lumi returns: " << myLumiMask.accept(315257, 10) << ", and known bad run returns: " << myLumiMask.accept(315257, 90) << std::endl;
  auto goldenjson = [myLumiMask](unsigned int &run, unsigned int &luminosityBlock){return myLumiMask.accept(run, luminosityBlock);};
  
  // -------------------------------------------------------
  //               correctionLib corrections
  // -------------------------------------------------------
  
  std::string yrstr, yr, jecyr, jeryr, jecver;
  float deepjetL;
  if(year == "2016APV") {deepjetL = 0.0508; yrstr = "2016preVFP"; yr = "16"; jecyr = "UL16APV"; jeryr = "Summer20UL16APV_JRV3"; jecver = "V7";}
  else if(year == "2016") {deepjetL = 0.0480; yrstr = "2016postVFP"; yr = "16"; jecyr = "UL16"; jeryr = "Summer20UL16_JRV3"; jecver = "V7";}
  else if(year == "2017") {deepjetL = 0.0532; yrstr = "2017"; yr = "17"; jecyr = "UL17"; jeryr = "Summer19UL17_JRV2"; jecver = "V5";}
  else if(year == "2018") {deepjetL = 0.0490; yrstr = "2018"; yr = "18"; jecyr = "UL18"; jeryr = "Summer19UL18_JRV2"; jecver = "V5";}
  else std::cout << "ERROR: Can't parse the year to assign correctionLib json files. Expected 2016, 2016APV, 2017, or 2018. Got: " << year << std::endl;

  // -------------------------------------------------------
  //               Open Dataframe + MET Filters
  // -------------------------------------------------------

  auto rdf_input = ROOT::RDataFrame("Events", files); // Initial data
  
  auto METgeneralFilters = rdf_input.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1", "MET Filters");

  auto truth = METgeneralFilters;
  
  // --------------------------------------------------------
  // 	       Golden JSON (Data)
  // --------------------------------------------------------
  
  if(!isMC){ // apply golden json to data
    truth = METgeneralFilters.Define("passesJSON", goldenjson, {"run","luminosityBlock"})
      .Filter("passesJSON == 1", "Data passes Golden JSON");
  }
  
  // ---------------------------------------------------------
  //                    LEPTON Definitions
  // ---------------------------------------------------------
  
  string elHEMcut = "";
  if(year == "2018") elHEMcut = " && (Electron_eta > -1.479 || (Electron_phi < -1.57 || Electron_phi > -0.87))";

  string tageltrig = "HLT_Ele35_WPTight_Gsf == 1 || HLT_Ele38_WPTight_Gsf == 1 || HLT_Photon200";    // tag trig for 2017B,C: only 35 38, no 32
  string tagmutrig = "HLT_Mu50";
  if(year == "2016" || year == "2016APV"){
    tageltrig = "HLT_Ele27_WPTight_Gsf || HLT_Photon175";
    if(year == "2016APV" && era != "A" && era != "B") tagmutrig = "HLT_Mu50 || HLT_TkMu50";
  }
  if((year == "2017" && era != "B") || year == "2018") tagmutrig = "HLT_Mu50 || HLT_OldMu100 || HLT_TkMu100";
   
  string eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon200";
  string tkmutrig = " || HLT_OldMu100 || HLT_TkMu100";
  if(year == "2017" && era == "B"){ tkmutrig = ""; eltrig = "HLT_Ele35_WPTight_Gsf || HLT_Photon200";}
  if(year == "2016" or year == "2016APV"){ tkmutrig = " || HLT_TkMu50"; eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon175";}
  if(year == "2016APV" and (era == "A" || era == "B")){ tkmutrig = ""; eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon175";}

  auto LepDefs = truth.Define("Electron_cutBasedIdNoIso_tight", "Electron_cutBasedIdNoIso_tight(nElectron, Electron_vidNestedWPBitmap)")
    .Define("TagMu","Muon_pt > 55 && abs(Muon_eta)<2.4 && Muon_tightId==1 && Muon_pfIsoId>=3 ") // && abs(Muon_dz) < 0.5 && Muon_dxy < 0.2")
    .Define("ProbeMu","Muon_pt > 30 && abs(Muon_eta)<2.4 && Muon_mediumId==1 && Muon_miniIsoId>=3 && abs(Muon_dz) < 0.5 && Muon_dxy < 0.2")
    .Define("ProbeMunoTag","TagMu == 0 && ProbeMu == 1")
    .Define("TagMunoProbe","ProbeMu == 0 && TagMu == 1")
    .Define("TagEl",Form("Electron_pt > 40 && (abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_cutBased>=4 %s",elHEMcut.c_str()))
    .Define("ProbeEl",Form("(abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_pt > 30 && Electron_cutBasedIdNoIso_tight==1 && Electron_miniPFRelIso_all<0.1%s",elHEMcut.c_str()))
    .Define("ProbeEl2",Form("(abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_pt > 30 && Electron_cutBased>=4 %s",elHEMcut.c_str()))
    .Define("ProbeEl3",Form("(abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_pt > 30 && Electron_mvaFall17V2noIso_WP90 > 0 %s",elHEMcut.c_str()))
    .Define("ProbeEl4",Form("(abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_pt > 30 && Electron_mvaFall17V2noIso_WP90 > 0 && Electron_miniPFRelIso_all < 0.1 %s",elHEMcut.c_str()))    
    .Define("ProbeElnoTag","TagEl == 0 && ProbeEl == 1")
    .Define("ProbeEl2noTag","TagEl == 0 && ProbeEl2 == 1")
    .Define("ProbeEl3noTag","TagEl == 0 && ProbeEl3 == 1")
    .Define("ProbeEl4noTag","TagEl == 0 && ProbeEl4 == 1")
    .Define("TagElnoProbe","ProbeEl == 0 && TagEl == 1")
    .Define("TagElnoProbe2","ProbeEl2 == 0 && TagEl == 1")
    .Define("TagElnoProbe3","ProbeEl3 == 0 && TagEl == 1")
    .Define("TagElnoProbe4","ProbeEl4 == 0 && TagEl == 1")    
    .Define("nTagMu","(int) Sum(TagMu)")
    .Define("nTagEl","(int) Sum(TagEl)")
    .Define("nProbeMu","(int) Sum(ProbeMu)")
    .Define("nProbeEl","(int) Sum(ProbeEl)")
    .Define("nProbeEl2","(int) Sum(ProbeEl2)")
    .Define("nProbeEl3","(int) Sum(ProbeEl3)")
    .Define("nProbeEl4","(int) Sum(ProbeEl4)")
    .Define("nProbeNoTagMu","(int) Sum(ProbeMunoTag)")
    .Define("nProbeNoTagEl","(int) Sum(ProbeElnoTag)")
    .Define("nProbe2NoTagEl","(int) Sum(ProbeEl2noTag)")
    .Define("nProbe3NoTagEl","(int) Sum(ProbeEl3noTag)")
    .Define("nProbe4NoTagEl","(int) Sum(ProbeEl4noTag)")
    .Define("nTagNoProbeMu","(int) Sum(TagMunoProbe)")
    .Define("nTagNoProbeEl","(int) Sum(TagElnoProbe)")
    .Define("nTagNoProbeEl2","(int) Sum(TagElnoProbe2)")
    .Define("nTagNoProbeEl3","(int) Sum(TagElnoProbe3)")
    .Define("nTagNoProbeEl4","(int) Sum(TagElnoProbe4)")
    .Define("passMuTrig",Form("(HLT_Mu50%s)",tkmutrig.c_str()))
    .Define("passElTrig",Form("(%s)",eltrig.c_str()))
    .Define("passTagMuTrig",Form("(%s)",tagmutrig.c_str()))
    .Define("passTagElTrig",Form("(%s)",tageltrig.c_str()))
    .Define("isTagElProbeMu","passTagElTrig && nTagEl == 1 && nProbeMu == 1 && nTagNoProbeMu == 0 && nProbeNoTagEl == 0")
    .Define("isTagMuProbeEl","passTagMuTrig && nTagMu == 1 && nProbeEl == 1 && nTagNoProbeEl == 0 && nProbeNoTagMu == 0")
    .Define("isTagMuProbeEl2","passTagMuTrig && nTagMu == 1 && nProbeEl2 == 1 && nTagNoProbeEl2 == 0 && nProbeNoTagMu == 0")
    .Define("isTagMuProbeEl3","passTagMuTrig && nTagMu == 1 && nProbeEl3 == 1 && nTagNoProbeEl3 == 0 && nProbeNoTagMu == 0")
    .Define("isTagMuProbeEl4","passTagMuTrig && nTagMu == 1 && nProbeEl4 == 1 && nTagNoProbeEl4 == 0 && nProbeNoTagMu == 0");

  // --------------------------------------------------------
  // 		      LEPTON SELECTION
  // --------------------------------------------------------

  auto LepSelect = LepDefs;
  auto LepSelect2 = LepDefs;
  auto LepSelect3 = LepDefs;
  auto LepSelect4 = LepDefs;

  if(isSM){
    
    LepSelect = LepDefs.Filter("nMuon > 0","Event has nMuon > 0")
      .Filter("nTagMu > 0","At least one muon tag")      
      .Filter("nTagMu == 1","Exactly 1 muon tag")
      .Filter("passTagMuTrig == 1","Passes muon tag trig Mu50 OR TkMu variants")
      .Filter("nProbeNoTagMu == 0","No probe muons failing tag criteria")
      .Filter("nElectron > 0","Event has nElectron > 0")
      .Filter("nProbeEl > 0","At least one electron probe")
      .Filter("nProbeEl == 1","Exactly 1 electron probe")
      .Filter("nTagNoProbeEl == 0","No tag electrons failing probe criteria")
      .Filter("isTagMuProbeEl", "Muon tag Electron probe event (eff == 1?)");

    LepSelect2 = LepDefs.Filter("nMuon > 0","Event has nMuon > 0")
      .Filter("nTagMu > 0","At least one muon tag")      
      .Filter("nTagMu == 1","Exactly 1 muon tag")
      .Filter("passTagMuTrig == 1","Passes muon tag trig Mu50 OR TkMu variants")
      .Filter("nProbeNoTagMu == 0","No probe muons failing tag criteria")
      .Filter("nElectron > 0","Event has nElectron > 0")
      .Filter("nProbeEl2 > 0","At least one electron probe")
      .Filter("nProbeEl2 == 1","Exactly 1 electron probe")
      .Filter("nTagNoProbeEl2 == 0","No tag electrons failing probe criteria")
      .Filter("isTagMuProbeEl2", "Muon tag Electron probe event (eff == 1?)");

    LepSelect3 = LepDefs.Filter("nMuon > 0","Event has nMuon > 0")
      .Filter("nTagMu > 0","At least one muon tag")      
      .Filter("nTagMu == 1","Exactly 1 muon tag")
      .Filter("passTagMuTrig == 1","Passes muon tag trig Mu50 OR TkMu variants")
      .Filter("nProbeNoTagMu == 0","No probe muons failing tag criteria")
      .Filter("nElectron > 0","Event has nElectron > 0")
      .Filter("nProbeEl3 > 0","At least one electron probe")
      .Filter("nProbeEl3 == 1","Exactly 1 electron probe")
      .Filter("nTagNoProbeEl3 == 0","No tag electrons failing probe criteria")
      .Filter("isTagMuProbeEl3", "Muon tag Electron probe event (eff == 1?)");

    LepSelect4 = LepDefs.Filter("nMuon > 0","Event has nMuon > 0")
      .Filter("nTagMu > 0","At least one muon tag")      
      .Filter("nTagMu == 1","Exactly 1 muon tag")
      .Filter("passTagMuTrig == 1","Passes muon tag trig Mu50 OR TkMu variants")
      .Filter("nProbeNoTagMu == 0","No probe muons failing tag criteria")
      .Filter("nElectron > 0","Event has nElectron > 0")
      .Filter("nProbeEl4 > 0","At least one electron probe")
      .Filter("nProbeEl4 == 1","Exactly 1 electron probe")
      .Filter("nTagNoProbeEl4 == 0","No tag electrons failing probe criteria")
      .Filter("isTagMuProbeEl4", "Muon tag Electron probe event (eff == 1?)");

  }else if(isSE){

    LepSelect = LepDefs.Filter("nElectron > 0","Event has nElectron > 0")
      .Filter("nTagEl > 0","At least one electron tag")
      .Filter("nTagEl == 1","Exactly 1 electron tag")
      .Filter("passTagElTrig == 1","Passes electron tag trig WPTight OR Photon")
      .Filter("nProbeNoTagEl == 0","No probe electrons failing tag criteria")
      .Filter("nProbeMu > 0","At least one muon probe")
      .Filter("nProbeMu == 1","Exactly 1 muon probe")
      .Filter("nTagNoProbeMu == 0","No tag muons failing probe criteria")
      .Filter("isTagElProbeMu", "Electron tag Muon probe event (eff == 1?)");

  }else{ // MC

    LepSelect = LepDefs.Filter("isTagMuProbeEl || isTagElProbeMu","Event is either muon tag or electron tag");
    LepSelect2 = LepDefs.Filter("isTagMuProbeEl2 || isTagElProbeMu","Event is either muon tag or electron tag");
    LepSelect3 = LepDefs.Filter("isTagMuProbeEl3 || isTagElProbeMu","Event is either muon tag or electron tag");
    LepSelect4 = LepDefs.Filter("isTagMuProbeEl4 || isTagElProbeMu","Event is either muon tag or electron tag");

  }

  auto MetSelect = LepSelect.Filter("MET_pt > 60","MET > 60");
  auto MetSelect2 = LepSelect2.Filter("MET_pt > 60","MET > 60");
  auto MetSelect3 = LepSelect3.Filter("MET_pt > 60","MET > 60");
  auto MetSelect4 = LepSelect4.Filter("MET_pt > 60","MET > 60");

  auto LepAssign = MetSelect.Define("assignleps", "assign_leps(isTagElProbeMu,isTagMuProbeEl,ProbeMu,ProbeEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_jetIdx,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_jetIdx)")
    .Define("probe_pt", "assignleps[0]")
    .Define("probe_eta", "assignleps[1]")
    .Define("probe_phi", "assignleps[2]")
    .Define("probe_mass", "assignleps[3]")
    .Define("probe_jetidx", "assignleps[4]")
    .Define("probe_abseta", "abs(probe_eta)")
    .Define("assigntags", "assign_leps(isTagMuProbeEl,isTagElProbeMu,TagMu,TagEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_jetIdx,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_jetIdx)")
    .Define("tag_pt", "assigntags[0]")
    .Define("tag_eta", "assigntags[1]")
    .Define("tag_phi", "assigntags[2]")
    .Define("tag_mass", "assigntags[3]")
    .Define("tag_jetidx", "assigntags[4]");
  auto LepAssign2 = MetSelect2.Define("assignleps", "assign_leps(isTagElProbeMu,isTagMuProbeEl2,ProbeMu,ProbeEl2,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_jetIdx,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_jetIdx)")
    .Define("probe_pt", "assignleps[0]")
    .Define("probe_eta", "assignleps[1]")
    .Define("probe_phi", "assignleps[2]")
    .Define("probe_mass", "assignleps[3]")
    .Define("probe_jetidx", "assignleps[4]")
    .Define("probe_abseta", "abs(probe_eta)")
    .Define("assigntags", "assign_leps(isTagMuProbeEl2,isTagElProbeMu,TagMu,TagEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_jetIdx,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_jetIdx)")
    .Define("tag_pt", "assigntags[0]")
    .Define("tag_eta", "assigntags[1]")
    .Define("tag_phi", "assigntags[2]")
    .Define("tag_mass", "assigntags[3]")
    .Define("tag_jetidx", "assigntags[4]");
  auto LepAssign3 = MetSelect3.Define("assignleps", "assign_leps(isTagElProbeMu,isTagMuProbeEl3,ProbeMu,ProbeEl3,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_jetIdx,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_jetIdx)")
    .Define("probe_pt", "assignleps[0]")
    .Define("probe_eta", "assignleps[1]")
    .Define("probe_phi", "assignleps[2]")
    .Define("probe_mass", "assignleps[3]")
    .Define("probe_jetidx", "assignleps[4]")
    .Define("probe_abseta", "abs(probe_eta)")
    .Define("assigntags", "assign_leps(isTagMuProbeEl3,isTagElProbeMu,TagMu,TagEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_jetIdx,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_jetIdx)")
    .Define("tag_pt", "assigntags[0]")
    .Define("tag_eta", "assigntags[1]")
    .Define("tag_phi", "assigntags[2]")
    .Define("tag_mass", "assigntags[3]")
    .Define("tag_jetidx", "assigntags[4]");
  auto LepAssign4 = MetSelect4.Define("assignleps", "assign_leps(isTagElProbeMu,isTagMuProbeEl4,ProbeMu,ProbeEl4,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_jetIdx,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_jetIdx)")
    .Define("probe_pt", "assignleps[0]")
    .Define("probe_eta", "assignleps[1]")
    .Define("probe_phi", "assignleps[2]")
    .Define("probe_mass", "assignleps[3]")
    .Define("probe_jetidx", "assignleps[4]")
    .Define("probe_abseta", "abs(probe_eta)")
    .Define("assigntags", "assign_leps(isTagMuProbeEl4,isTagElProbeMu,TagMu,TagEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_jetIdx,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_jetIdx)")
    .Define("tag_pt", "assigntags[0]")
    .Define("tag_eta", "assigntags[1]")
    .Define("tag_phi", "assigntags[2]")
    .Define("tag_mass", "assigntags[3]")
    .Define("tag_jetidx", "assigntags[4]");
  
  // --------------------------------------------------------
  // 		      JET Cleaning and JERC
  // --------------------------------------------------------
  
  auto JetSelect = LepAssign.Filter("nJet > 0", "Event has > 0 AK4")
    //.Filter("nFatJet > 0", "Event has > 0 AK8")
    .Define("DR_tagJets","DeltaR_VecAndFloat(Jet_eta,Jet_phi,tag_eta,tag_phi)")
    .Define("DR_probeJets","DeltaR_VecAndFloat(Jet_eta,Jet_phi,probe_eta,probe_phi)")
    .Define("ptrel_tagJets","ptRel(Jet_pt,Jet_eta,Jet_phi,Jet_mass,tag_pt,tag_eta,tag_phi,tag_mass)")
    .Define("ptrel_probeJets","ptRel(Jet_pt,Jet_eta,Jet_phi,Jet_mass,probe_pt,probe_eta,probe_phi,probe_mass)")
    .Define("goodJets", "Jet_pt > 30 && abs(Jet_eta) < 2.5 && Jet_jetId > 1 && (DR_tagJets > 0.4 || ptrel_tagJets > 20) && (DR_probeJets > 0.4 || ptrel_probeJets > 20)")
    //.Define("nonlepJets","drop_jets(goodJets,probe_jetidx,tag_jetidx,probe_eta,probe_phi,tag_eta,tag_phi)")
    .Define("NgoodJets","(int) Sum(goodJets)")
    .Filter("NgoodJets > 1","Pass 2 AK4s");

  auto JetVars = JetSelect.Define("goodjet_pt","Jet_pt[goodJets == true]")
    .Define("goodjet_HT","Sum(goodjet_pt)")
    .Define("goodJet_DeepFlav", "Jet_btagDeepFlavB[goodJets == true]")
    .Define("goodJet_DeepFlavL", Form("goodJet_DeepFlav > %f",deepjetL))
    .Define("NJets_DeepFlavL", "(int) Sum(goodJet_DeepFlavL)")
    .Define("leadjetpt","goodjet_pt[0]")
    .Define("subleadjetpt","goodjet_pt[1]");

  auto JetSelect2 = LepAssign2.Filter("nJet > 0", "Event has > 0 AK4")
    //.Filter("nFatJet > 0", "Event has > 0 AK8")
    .Define("DR_tagJets","DeltaR_VecAndFloat(Jet_eta,Jet_phi,tag_eta,tag_phi)")
    .Define("DR_probeJets","DeltaR_VecAndFloat(Jet_eta,Jet_phi,probe_eta,probe_phi)")
    .Define("ptrel_tagJets","ptRel(Jet_pt,Jet_eta,Jet_phi,Jet_mass,tag_pt,tag_eta,tag_phi,tag_mass)")
    .Define("ptrel_probeJets","ptRel(Jet_pt,Jet_eta,Jet_phi,Jet_mass,probe_pt,probe_eta,probe_phi,probe_mass)")
    .Define("goodJets", "Jet_pt > 30 && abs(Jet_eta) < 2.5 && Jet_jetId > 1 && (DR_tagJets > 0.4 || ptrel_tagJets > 20) && (DR_probeJets > 0.4 || ptrel_probeJets > 20)")
    //.Define("nonlepJets","drop_jets(goodJets,probe_jetidx,tag_jetidx,probe_eta,probe_phi,tag_eta,tag_phi)")
    .Define("NgoodJets","(int) Sum(goodJets)")
    .Filter("NgoodJets > 1","Pass 2 AK4s");

  auto JetVars2 = JetSelect2.Define("goodjet_pt","Jet_pt[goodJets == true]")
    .Define("goodjet_HT","Sum(goodjet_pt)")
    .Define("goodJet_DeepFlav", "Jet_btagDeepFlavB[goodJets == true]")
    .Define("goodJet_DeepFlavL", Form("goodJet_DeepFlav > %f",deepjetL))
    .Define("NJets_DeepFlavL", "(int) Sum(goodJet_DeepFlavL)")
    .Define("leadjetpt","goodjet_pt[0]")
    .Define("subleadjetpt","goodjet_pt[1]");

  auto JetSelect3 = LepAssign3.Filter("nJet > 0", "Event has > 0 AK4")
    //.Filter("nFatJet > 0", "Event has > 0 AK8")
    .Define("DR_tagJets","DeltaR_VecAndFloat(Jet_eta,Jet_phi,tag_eta,tag_phi)")
    .Define("DR_probeJets","DeltaR_VecAndFloat(Jet_eta,Jet_phi,probe_eta,probe_phi)")
    .Define("ptrel_tagJets","ptRel(Jet_pt,Jet_eta,Jet_phi,Jet_mass,tag_pt,tag_eta,tag_phi,tag_mass)")
    .Define("ptrel_probeJets","ptRel(Jet_pt,Jet_eta,Jet_phi,Jet_mass,probe_pt,probe_eta,probe_phi,probe_mass)")
    .Define("goodJets", "Jet_pt > 30 && abs(Jet_eta) < 2.5 && Jet_jetId > 1 && (DR_tagJets > 0.4 || ptrel_tagJets > 20) && (DR_probeJets > 0.4 || ptrel_probeJets > 20)")
    //.Define("nonlepJets","drop_jets(goodJets,probe_jetidx,tag_jetidx,probe_eta,probe_phi,tag_eta,tag_phi)")
    .Define("NgoodJets","(int) Sum(goodJets)")
    .Filter("NgoodJets > 1","Pass 2 AK4s");

  auto JetVars3 = JetSelect3.Define("goodjet_pt","Jet_pt[goodJets == true]")
    .Define("goodjet_HT","Sum(goodjet_pt)")
    .Define("goodJet_DeepFlav", "Jet_btagDeepFlavB[goodJets == true]")
    .Define("goodJet_DeepFlavL", Form("goodJet_DeepFlav > %f",deepjetL))
    .Define("NJets_DeepFlavL", "(int) Sum(goodJet_DeepFlavL)")
    .Define("leadjetpt","goodjet_pt[0]")
    .Define("subleadjetpt","goodjet_pt[1]");

  auto JetSelect4 = LepAssign4.Filter("nJet > 0", "Event has > 0 AK4")
    //.Filter("nFatJet > 0", "Event has > 0 AK8")
    .Define("DR_tagJets","DeltaR_VecAndFloat(Jet_eta,Jet_phi,tag_eta,tag_phi)")
    .Define("DR_probeJets","DeltaR_VecAndFloat(Jet_eta,Jet_phi,probe_eta,probe_phi)")
    .Define("ptrel_tagJets","ptRel(Jet_pt,Jet_eta,Jet_phi,Jet_mass,tag_pt,tag_eta,tag_phi,tag_mass)")
    .Define("ptrel_probeJets","ptRel(Jet_pt,Jet_eta,Jet_phi,Jet_mass,probe_pt,probe_eta,probe_phi,probe_mass)")
    .Define("goodJets", "Jet_pt > 30 && abs(Jet_eta) < 2.5 && Jet_jetId > 1 && (DR_tagJets > 0.4 || ptrel_tagJets > 20) && (DR_probeJets > 0.4 || ptrel_probeJets > 20)")
    //.Define("nonlepJets","drop_jets(goodJets,probe_jetidx,tag_jetidx,probe_eta,probe_phi,tag_eta,tag_phi)")
    .Define("NgoodJets","(int) Sum(goodJets)")
    .Filter("NgoodJets > 1","Pass 2 AK4s");

  auto JetVars4 = JetSelect4.Define("goodjet_pt","Jet_pt[goodJets == true]")
    .Define("goodjet_HT","Sum(goodjet_pt)")
    .Define("goodJet_DeepFlav", "Jet_btagDeepFlavB[goodJets == true]")
    .Define("goodJet_DeepFlavL", Form("goodJet_DeepFlav > %f",deepjetL))
    .Define("NJets_DeepFlavL", "(int) Sum(goodJet_DeepFlavL)")
    .Define("leadjetpt","goodjet_pt[0]")
    .Define("subleadjetpt","goodjet_pt[1]");

  // -------------------------------------------------
  // 		Store histograms
  // -------------------------------------------------
  
  float ptbinsMu[6] = {55,60,120,200,300,1200}; // more ptbins than POG at high pt, but stats seem ok.
  float ptbinsEl[6] = {55,120,200,1200}; // use same ptbins as POG, except edge at 120 > 115
  float etabinsMu[5] = {0,0.9,1.2,2.1,2.4};
  float etabinsEl[11] = {-2.5, -2.0, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.0, 2.5};

  // auto m01 = JetVars.Filter("isTagElProbeMu == 1")
  //   .Histo2D<float,float>({"TrigEff_Dptbins_Mu",";#mu p_T (GeV);#mu |#eta|",5,ptbinsMu,4,etabinsMu},"probe_pt","probe_abseta");
  // auto m02 = JetVars.Filter("isTagElProbeMu == 1 && passMuTrig == 1")
  //   .Histo2D<float,float>({"TrigEff_Nptbins_Mu",";#mu p_T (GeV);#mu |#eta|",5,ptbinsMu,4,etabinsMu},"probe_pt","probe_abseta");
  // auto m03 = JetVars.Filter("isTagElProbeMu == 1 && goodjet_HT > 250")
  //   .Histo2D<float,float>({"TrigEff_Dptbins_ht250_Mu",";#mu p_T (GeV);#mu |#eta|",5,ptbinsMu,4,etabinsMu},"probe_pt","probe_abseta");
  // auto m04 = JetVars.Filter("isTagElProbeMu == 1 && passMuTrig == 1 && goodjet_HT > 250")
  //   .Histo2D<float,float>({"TrigEff_Nptbins_ht250_Mu",";#mu p_T (GeV);#mu |#eta|",5,ptbinsMu,4,etabinsMu},"probe_pt","probe_abseta");
  // auto m05 = JetVars.Filter("isTagElProbeMu == 1 && goodjet_HT > 250 && NJets_DeepFlavL > 1")
  //   .Histo2D<float,float>({"TrigEff_Dptbins_ht250b2_Mu",";#mu p_T (GeV);#mu |#eta|",5,ptbinsMu,4,etabinsMu},"probe_pt","probe_abseta");
  // auto m06 = JetVars.Filter("isTagElProbeMu == 1 && passMuTrig == 1 && goodjet_HT > 250 && NJets_DeepFlavL > 1")
  //   .Histo2D<float,float>({"TrigEff_Nptbins_ht250b2_Mu",";#mu p_T (GeV);#mu |#eta|",5,ptbinsMu,4,etabinsMu},"probe_pt","probe_abseta");

  auto e01 = JetVars.Filter("isTagMuProbeEl == 1")
    .Histo2D<float,float>({"TrigEff_Dptbins_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e02 = JetVars.Filter("isTagMuProbeEl == 1 && passElTrig == 1")
    .Histo2D<float,float>({"TrigEff_Nptbins_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e03 = JetVars.Filter("isTagMuProbeEl == 1 && leadjetpt > 185")
    .Histo2D<float,float>({"TrigEff_Dptbins_jet185_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e04 = JetVars.Filter("isTagMuProbeEl == 1 && passElTrig == 1 && leadjetpt > 185")
    .Histo2D<float,float>({"TrigEff_Nptbins_jet185_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e05 = JetVars.Filter("isTagMuProbeEl == 1 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Dptbins_ht250_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e06 = JetVars.Filter("isTagMuProbeEl == 1 && passElTrig == 1 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Nptbins_ht250_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e09 = JetVars.Filter("isTagMuProbeEl == 1 && leadjetpt > 185 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Dptbins_jet185ht250_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e10 = JetVars.Filter("isTagMuProbeEl == 1 && passElTrig == 1 && leadjetpt > 185 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Nptbins_jet185ht250_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");

  auto e201 = JetVars2.Filter("isTagMuProbeEl2 == 1")
    .Histo2D<float,float>({"TrigEff_Dptbins_El2",";electron p_T (GeV) ID2;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e202 = JetVars2.Filter("isTagMuProbeEl2 == 1 && passElTrig == 1")
    .Histo2D<float,float>({"TrigEff_Nptbins_El2",";electron p_T (GeV) ID2;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e203 = JetVars2.Filter("isTagMuProbeEl2 == 1 && leadjetpt > 185")
    .Histo2D<float,float>({"TrigEff_Dptbins_jet185_El2",";electron p_T (GeV) ID2;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e204 = JetVars2.Filter("isTagMuProbeEl2 == 1 && passElTrig == 1 && leadjetpt > 185")
    .Histo2D<float,float>({"TrigEff_Nptbins_jet185_El2",";electron p_T (GeV) ID2;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e205 = JetVars2.Filter("isTagMuProbeEl2 == 1 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Dptbins_ht250_El2",";electron p_T (GeV) ID2;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e206 = JetVars2.Filter("isTagMuProbeEl2 == 1 && passElTrig == 1 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Nptbins_ht250_El2",";electron p_T (GeV) ID2;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e209 = JetVars2.Filter("isTagMuProbeEl2 == 1 && leadjetpt > 185 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Dptbins_jet185ht250_El2",";electron p_T (GeV) ID2;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e210 = JetVars2.Filter("isTagMuProbeEl2 == 1 && passElTrig == 1 && leadjetpt > 185 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Nptbins_jet185ht250_El2",";electron p_T (GeV) ID2;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
 
  auto e301 = JetVars3.Filter("isTagMuProbeEl3 == 1")
    .Histo2D<float,float>({"TrigEff_Dptbins_El3",";electron p_T (GeV) ID3;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e302 = JetVars3.Filter("isTagMuProbeEl3 == 1 && passElTrig == 1")
    .Histo2D<float,float>({"TrigEff_Nptbins_El3",";electron p_T (GeV) ID3;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e303 = JetVars3.Filter("isTagMuProbeEl3 == 1 && leadjetpt > 185")
    .Histo2D<float,float>({"TrigEff_Dptbins_jet185_El3",";electron p_T (GeV) ID3;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e304 = JetVars3.Filter("isTagMuProbeEl3 == 1 && passElTrig == 1 && leadjetpt > 185")
    .Histo2D<float,float>({"TrigEff_Nptbins_jet185_El3",";electron p_T (GeV) ID3;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e305 = JetVars3.Filter("isTagMuProbeEl3 == 1 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Dptbins_ht250_El3",";electron p_T (GeV) ID3;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e306 = JetVars3.Filter("isTagMuProbeEl3 == 1 && passElTrig == 1 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Nptbins_ht250_El3",";electron p_T (GeV) ID3;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e309 = JetVars3.Filter("isTagMuProbeEl3 == 1 && leadjetpt > 185 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Dptbins_jet185ht250_El3",";electron p_T (GeV) ID3;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e310 = JetVars3.Filter("isTagMuProbeEl3 == 1 && passElTrig == 1 && leadjetpt > 185 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Nptbins_jet185ht250_El3",";electron p_T (GeV) ID3;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");

  auto e401 = JetVars4.Filter("isTagMuProbeEl4 == 1")
    .Histo2D<float,float>({"TrigEff_Dptbins_El4",";electron p_T (GeV) ID4;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e402 = JetVars4.Filter("isTagMuProbeEl4 == 1 && passElTrig == 1")
    .Histo2D<float,float>({"TrigEff_Nptbins_El4",";electron p_T (GeV) ID4;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e403 = JetVars4.Filter("isTagMuProbeEl4 == 1 && leadjetpt > 185")
    .Histo2D<float,float>({"TrigEff_Dptbins_jet185_El4",";electron p_T (GeV) ID4;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e404 = JetVars4.Filter("isTagMuProbeEl4 == 1 && passElTrig == 1 && leadjetpt > 185")
    .Histo2D<float,float>({"TrigEff_Nptbins_jet185_El4",";electron p_T (GeV) ID4;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e405 = JetVars4.Filter("isTagMuProbeEl4 == 1 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Dptbins_ht250_El4",";electron p_T (GeV) ID4;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e406 = JetVars4.Filter("isTagMuProbeEl4 == 1 && passElTrig == 1 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Nptbins_ht250_El4",";electron p_T (GeV) ID4;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e409 = JetVars4.Filter("isTagMuProbeEl4 == 1 && leadjetpt > 185 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Dptbins_jet185ht250_El4",";electron p_T (GeV) ID4;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  auto e410 = JetVars4.Filter("isTagMuProbeEl4 == 1 && passElTrig == 1 && leadjetpt > 185 && goodjet_HT > 250")
    .Histo2D<float,float>({"TrigEff_Nptbins_jet185ht250_El4",";electron p_T (GeV) ID4;electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");


  // auto e07 = JetVars.Filter("isTagMuProbeEl == 1 && goodjet_HT > 250 && NJets_DeepFlavL > 1")
  //   .Histo2D<float,float>({"TrigEff_Dptbins_ht250b2_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  // auto e08 = JetVars.Filter("isTagMuProbeEl == 1 && passElTrig == 1 && goodjet_HT > 250 && NJets_DeepFlavL > 1")
  //   .Histo2D<float,float>({"TrigEff_Nptbins_ht250b2_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  // auto e11 = JetVars.Filter("isTagMuProbeEl == 1 && leadjetpt > 185 && goodjet_HT > 250 && NJets_DeepFlavL > 1")
  //   .Histo2D<float,float>({"TrigEff_Dptbins_jet185ht250b2_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  // auto e12 = JetVars.Filter("isTagMuProbeEl == 1 && passElTrig == 1 && leadjetpt > 185 && goodjet_HT > 250 && NJets_DeepFlavL > 1")
  //   .Histo2D<float,float>({"TrigEff_Nptbins_jet185ht250b2_El",";electron p_T (GeV);electron #eta",3,ptbinsEl,10,etabinsEl},"probe_pt","probe_eta");
  
  // -------------------------------------------------
  // 		Save histograms to file
  // -------------------------------------------------

  cout << "-------------------------------------------------" << endl
       << ">>> Writing " << sample << " histograms..." << endl;
  TString finalFile = "RDF_" + sample + era + "_" + year + "_" + testNum.Data() + "_TrigEff.root";
  const char *stdfinalFile = finalFile;

  TFile::Open(finalFile,"recreate");
  // m01->Write();
  // m02->Write();
  // m03->Write();
  // m04->Write();
  // m05->Write();
  // m06->Write();

  e01->Write();
  e02->Write();
  e03->Write();
  e04->Write();
  e05->Write();
  e06->Write();
  e09->Write();
  e10->Write();
  // e07->Write();
  // e08->Write();
  // e11->Write();
  // e12->Write();

  e201->Write();
  e202->Write();
  e203->Write();
  e204->Write();
  e205->Write();
  e206->Write();
  e209->Write();
  e210->Write();

  e301->Write();
  e302->Write();
  e303->Write();
  e304->Write();
  e305->Write();
  e306->Write();
  e309->Write();
  e310->Write();

  e401->Write();
  e402->Write();
  e403->Write();
  e404->Write();
  e405->Write();
  e406->Write();
  e409->Write();
  e410->Write();

  time.Stop();
  time.Print();
  
  cout << "Cut statistics:" << endl;
  JetVars.Report()->Print();

  cout << "Done!" << endl;
}
