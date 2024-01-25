// --------------------------------------------------------------------------------------- //
// Implimentation of RDataFrame in C++.					                   //
// Comments on creating a singly produced VLQ search			                   //
// To Run on Command Line:   root -l callRDF.C\(\"Muon(OR)Electron\",\"testNumber\"\,\"root://cmsxrootd.fnal.gov//store/...file.root\")      //
// --------------------------------------------------------------------------------------- //

#define rdf_cxx
#include "analyzer_taggerEff.h"
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
#include <algorithm> // std::sort
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TVector2.h>
#include <TRandom3.h>
#include <sstream>
#include <chrono> // for high_resolution_clock
#include "../correctionlib/include/correction.h"

using namespace std;
using namespace ROOT::VecOps;
using correction::CorrectionSet;

void rdf::analyzer_taggerEff(TString testNum, TString jesvar)
{
  //  ROOT::EnableImplicitMT();
  TStopwatch time;
  time.Start();
  string sample = this->sample;
  string year = this->year;
  bool isMC = this->isMC;

  bool debug = false;
  
  cout << "Sample in cc: " << sample << endl;
  cout << "Year in cc: " << year << endl;
  cout << "isMC? " << isMC << ", jesvar = " << jesvar << endl;
  if(!isMC) cout << "Data era = " << era << ", for jec " << jecera << endl;

  float deepjetL;
  if(year == "2016APV") deepjetL = 0.0508;
  else if(year == "2016") deepjetL = 0.0480;
  else if(year == "2017") deepjetL = 0.0532;
  else if(year == "2018") deepjetL = 0.0490;

  // -------------------------------------------------------
  //               Open Dataframe + GEN info
  // -------------------------------------------------------

  auto rdf_input = ROOT::RDataFrame("Events", files); // Initial data
  
  // auto BprimeGen = rdf_input.Define("Bprime_gen_info", Form("Bprime_gen_info(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status, GenPart_statusFlags)", sample.c_str()))
  //   .Define("Bprime_gen_pt", "Bprime_gen_info[0]")
  //   .Define("Bprime_gen_eta", "(double) Bprime_gen_info[1]")
  //   .Define("Bprime_gen_phi", "(double) Bprime_gen_info[2]")
  //   .Define("Bprime_gen_mass", "Bprime_gen_info[3]")
  //   .Define("Bprime_gen_pdgId", "(int) Bprime_gen_info[4]")
  //   .Define("Bprime_gen_status", "(int) Bprime_gen_info[5]");
  
  // auto truth = BprimeGen.Define("t_gen_info", Form("t_gen_info(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status)", sample.c_str()))
  //   .Define("t_gen_pt", "t_gen_info[0]")
  //   .Define("t_gen_eta", "(double) t_gen_info[1]")
  //   .Define("t_gen_phi", "(double) t_gen_info[2]")
  //   .Define("t_gen_mass", "t_gen_info[3]")
  //   .Define("t_gen_pdgId", "(int) t_gen_info[4]")
  //   .Define("t_gen_status", "(int) t_gen_info[5]")
  //   .Define("daughterb_gen_pt", "t_gen_info[6]")
  //   .Define("daughterb_gen_eta", "(double) t_gen_info[7]")
  //   .Define("daughterb_gen_phi", "(double) t_gen_info[8]")
  //   .Define("daughterb_gen_pdgId", "(int) t_gen_info[9]")
  //   .Define("daughterb_gen_status", "(int) t_gen_info[10]")
  //   .Define("daughterW_gen_pt", "t_gen_info[11]")
  //   .Define("daughterW_gen_eta", "(double) t_gen_info[12]")
  //   .Define("daughterW_gen_phi", "(double) t_gen_info[13]")
  //   .Define("daughterW_gen_mass", "t_gen_info[14]")
  //   .Define("daughterW_gen_pdgId", "(int) t_gen_info[15]")
  //   .Define("daughterW_gen_status", "(int) t_gen_info[16]")
  //   .Define("daughterWb_gen_dr", "DeltaR(daughterW_gen_eta, daughterb_gen_eta, daughterW_gen_phi, daughterb_gen_phi)")
  //   .Define("tDaughter1_gen_pt", "t_gen_info[17]") // e/mu/tau or quark1
  //   .Define("tDaughter1_gen_eta", "(double) t_gen_info[18]")
  //   .Define("tDaughter1_gen_phi", "(double) t_gen_info[19]")
  //   .Define("tDaughter1_gen_mass", "t_gen_info[20]")
  //   .Define("tDaughter1_gen_pdgId", "(int) t_gen_info[21]")
  //   .Define("tDaughter1_gen_status", "(int) t_gen_info[22]")
  //   .Define("tDaughter2_gen_pt", "t_gen_info[23]") // neutrino or quark2
  //   .Define("tDaughter2_gen_eta", "(double) t_gen_info[24]")
  //   .Define("tDaughter2_gen_phi", "(double) t_gen_info[25]")
  //   .Define("tDaughter2_gen_mass", "t_gen_info[26]")
  //   .Define("tDaughter2_gen_pdgId", "(int) t_gen_info[27]")
  //   .Define("tDaughter2_gen_status", "(int) t_gen_info[28]")
  //   .Define("tdaughter12_gen_dr", "DeltaR(tDaughter1_gen_eta, tDaughter2_gen_eta, tDaughter1_gen_phi, tDaughter2_gen_phi)")
  //   .Define("trueLeptonicT", "(int) t_gen_info[29]")
  //   .Define("W_gen_info", Form("W_gen_info(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status, daughterW_gen_pdgId)", sample.c_str()))
  //   .Define("W_gen_pt", "W_gen_info[0]")
  //   .Define("W_gen_eta", "(double) W_gen_info[1]")
  //   .Define("W_gen_phi", "(double) W_gen_info[2]")
  //   .Define("W_gen_mass", "W_gen_info[3]")
  //   .Define("W_gen_pdgId", "(int) W_gen_info[4]")
  //   .Define("W_gen_status", "(int) W_gen_info[5]")
  //   .Define("WDaughter1_gen_pt", "W_gen_info[6]")
  //   .Define("WDaughter1_gen_eta", "(double) W_gen_info[7]")
  //   .Define("WDaughter1_gen_phi", "(double) W_gen_info[8]")
  //   .Define("WDaughter1_gen_mass", "W_gen_info[9]")
  //   .Define("WDaughter1_gen_pdgId", "(int) W_gen_info[10]")
  //   .Define("WDaughter1_gen_status", "(int) W_gen_info[11]")
  //   .Define("WDaughter2_gen_pt", "W_gen_info[12]")
  //   .Define("WDaughter2_gen_eta", "(double) W_gen_info[13]")
  //   .Define("WDaughter2_gen_phi", "(double) W_gen_info[14]")
  //   .Define("WDaughter2_gen_mass", "W_gen_info[15]")
  //   .Define("WDaughter2_gen_pdgId", "(int) W_gen_info[16]")
  //   .Define("WDaughter2_gen_status", "(int) W_gen_info[17]")
  //   .Define("Wdaughter12_gen_dr", "DeltaR(WDaughter1_gen_eta, WDaughter2_gen_eta, WDaughter1_gen_phi, WDaughter2_gen_phi)")
  //   .Define("trueLeptonicW", "(int) W_gen_info[18]")
  //   .Define("trueLeptonicMode", Form("leptonicCheck(\"%s\", trueLeptonicT, trueLeptonicW)", sample.c_str()))
  //   .Define("t_bkg_idx", Form("t_bkg_idx(\"%s\", nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)", sample.c_str()))
  //   .Define("W_bkg_idx", Form("W_bkg_idx(\"%s\", nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags, t_bkg_idx)", sample.c_str()))
  //   .Define("gcFatJet_genmatch", Form("FatJet_matching_bkg(\"%s\", gcFatJet_eta, gcFatJet_phi, NFatJets, gcFatJet_subJetIdx1, nSubJet, SubJet_hadronFlavour, nGenPart, GenPart_pdgId, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, t_bkg_idx, W_bkg_idx)",sample.c_str()));

    
  // ---------------------------------------------------------
  // 	  Jet Selection (must be 1 central)
  // ---------------------------------------------------------
  
  auto JetSelect = rdf_input.Filter("nJet > 0","Pass 1 jet with pt > 15")
    .Define("goodcleanJets", "abs(Jet_eta) < 2.5 && Jet_jetId > 1")
    .Filter("(int) Sum(goodcleanJets) > 0","Pass 1 central jet passing loose ID")
    .Define("gcJet_pt","Jet_pt[goodcleanJets == true]") 
    .Define("gcJet_hflav","Jet_hadronFlavour[goodcleanJets == true]")
    .Define("gcJet_deepjet","Jet_btagDeepFlavB[goodcleanJets == true]")
    .Define("bJet_pt","gcJet_pt[gcJet_hflav == 5]")
    .Define("cJet_pt","gcJet_pt[gcJet_hflav == 4]")
    .Define("lJet_pt","gcJet_pt[gcJet_hflav < 4]")
    .Define("bJet_pt_loose",Form("gcJet_pt[gcJet_hflav == 5 && gcJet_deepjet > %f]",deepjetL))
    .Define("cJet_pt_loose",Form("gcJet_pt[gcJet_hflav == 4 && gcJet_deepjet > %f]",deepjetL))
    .Define("lJet_pt_loose",Form("gcJet_pt[gcJet_hflav < 4 && gcJet_deepjet > %f]",deepjetL))
    .Define("weight","genWeight/abs(genWeight)");

  double ptbins[16] = {15,20,30,50,70,100,150,200,300,400,500,600,800,1000,1200,1500};

  auto h1 = JetSelect.Histo1D({"BEff_Dptbins_b",";Jet p_T (GeV);",15,ptbins},"bJet_pt","weight");
  auto h2 = JetSelect.Histo1D({"BEff_Dptbins_c",";Jet p_T (GeV);",15,ptbins},"cJet_pt","weight");
  auto h3 = JetSelect.Histo1D({"BEff_Dptbins_udsg",";Jet p_T (GeV);",15,ptbins},"lJet_pt","weight");
  auto h4 = JetSelect.Histo1D({"BEffLoose_Nptbins_b",";Jet p_T (GeV);DeepJetL: b efficiency",15,ptbins},"bJet_pt_loose","weight");
  auto h5 = JetSelect.Histo1D({"BEffLoose_Nptbins_c",";Jet p_T (GeV);DeepJetL: c efficiency",15,ptbins},"cJet_pt_loose","weight");
  auto h6 = JetSelect.Histo1D({"BEffLoose_Nptbins_udsg",";Jet p_T (GeV);DeepJetL: udsg efficiency",15,ptbins},"lJet_pt_loose","weight");
  
  // -------------------------------------------------
  // 		Save histos to file
  // -------------------------------------------------
  
  cout << "-------------------------------------------------" << endl
       << ">>> Saving " << sample << " Snapshot..." << endl;
  TString finalFile = "RDF_" + sample + "_" + year + "_DeepJetEff.root";
  const char *stdfinalFile = finalFile;

  TFile::Open(finalFile,"recreate");
  h1->Write();
  h2->Write();
  h3->Write();
  h4->Write();
  h5->Write();
  h6->Write();

  time.Stop();
  time.Print();
  
  cout << "Cut statistics:" << endl;
  JetSelect.Report()->Print();

  cout << "Done!" << endl;
}
