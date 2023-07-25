// --------------------------------------------------------------------------------------- //
// Implimentation of RDataFrame in C++.					                   //
// Comments on creating a singly produced VLQ search			                   //
// To Run on Command Line:   root -l callRDF.C\(\"Muon(OR)Electron\",\"testNumber\"\,\"root://cmsxrootd.fnal.gov//store/...file.root\")      //
// --------------------------------------------------------------------------------------- //

#define rdf_cxx
#include "analyzer_RDF.h"
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
#include <TRandom3.h>
#include <sstream>
#include <chrono> // for high_resolution_clock

using namespace std;
using namespace ROOT::VecOps;

void rdf::analyzer_RDF(TString testNum)
{
     ROOT::EnableImplicitMT();
     TStopwatch time;
     time.Start();
     bool isNominal = this->isNominal;
     string sample = this->sample;
     string year = this->year;

     cout << "Sample in cc: " << sample << endl;
     cout << "Year in cc: " << year << endl;

     // -------------------------------------------------------
     //               Flags and First Filter
     // -------------------------------------------------------
     // Twiki with reccommended ultralegacy values
     auto rdf_input = ROOT::RDataFrame("Events", files); // Initial data

     auto METgeneralFilters = rdf_input.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1", "MET Filters")
                                  .Filter("nJet > 0 && nFatJet > 0", "Event has jets");
     auto truth = METgeneralFilters;



     if (!isSM && !isSE) {
          auto BprimeGen = METgeneralFilters.Define("Bprime_gen_info", Form("Bprime_gen_info(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status, GenPart_statusFlags)", sample.c_str()))
                         .Define("Bprime_gen_pt", "Bprime_gen_info[0]")
                         .Define("Bprime_gen_eta", "(double) Bprime_gen_info[1]")
                         .Define("Bprime_gen_phi", "(double) Bprime_gen_info[2]")
                         .Define("Bprime_gen_mass", "Bprime_gen_info[3]")
                         .Define("Bprime_gen_pdgId", "(int) Bprime_gen_info[4]")
                         .Define("Bprime_gen_status", "(int) Bprime_gen_info[5]");
                         
          truth = BprimeGen.Define("t_gen_info", Form("t_gen_info(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status)", sample.c_str()))
                         .Define("t_gen_pt", "t_gen_info[0]")
                         .Define("t_gen_eta", "(double) t_gen_info[1]")
                         .Define("t_gen_phi", "(double) t_gen_info[2]")
                         .Define("t_gen_mass", "t_gen_info[3]")
                         .Define("t_gen_pdgId", "(int) t_gen_info[4]")
                         .Define("t_gen_status", "(int) t_gen_info[5]")
                         .Define("daughterb_gen_pt", "t_gen_info[6]")
                         .Define("daughterb_gen_eta", "(double) t_gen_info[7]")
                         .Define("daughterb_gen_phi", "(double) t_gen_info[8]")
                         .Define("daughterb_gen_pdgId", "(int) t_gen_info[9]")
                         .Define("daughterb_gen_status", "(int) t_gen_info[10]")
                         .Define("daughterW_gen_pt", "t_gen_info[11]")
                         .Define("daughterW_gen_eta", "(double) t_gen_info[12]")
                         .Define("daughterW_gen_phi", "(double) t_gen_info[13]")
                         .Define("daughterW_gen_mass", "t_gen_info[14]")
                         .Define("daughterW_gen_pdgId", "(int) t_gen_info[15]")
                         .Define("daughterW_gen_status", "(int) t_gen_info[16]")
                         .Define("daughterWb_gen_dr", "DeltaR(daughterW_gen_eta, daughterb_gen_eta, daughterW_gen_phi, daughterb_gen_phi)")
                         .Define("tDaughter1_gen_pt", "t_gen_info[17]") // e/mu/tau or quark1
                         .Define("tDaughter1_gen_eta", "(double) t_gen_info[18]")
                         .Define("tDaughter1_gen_phi", "(double) t_gen_info[19]")
                         .Define("tDaughter1_gen_mass", "t_gen_info[20]")
                         .Define("tDaughter1_gen_pdgId", "(int) t_gen_info[21]")
                         .Define("tDaughter1_gen_status", "(int) t_gen_info[22]")
                         .Define("tDaughter2_gen_pt", "t_gen_info[23]") // neutrino or quark2
                         .Define("tDaughter2_gen_eta", "(double) t_gen_info[24]")
                         .Define("tDaughter2_gen_phi", "(double) t_gen_info[25]")
                         .Define("tDaughter2_gen_mass", "t_gen_info[26]")
                         .Define("tDaughter2_gen_pdgId", "(int) t_gen_info[27]")
                         .Define("tDaughter2_gen_status", "(int) t_gen_info[28]")
                         .Define("tdaughter12_gen_dr", "DeltaR(tDaughter1_gen_eta, tDaughter2_gen_eta, tDaughter1_gen_phi, tDaughter2_gen_phi)")
                         .Define("trueLeptonicT", "(int) t_gen_info[29]")
                         .Define("W_gen_info", Form("W_gen_info(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status, daughterW_gen_pdgId)", sample.c_str()))
                         .Define("W_gen_pt", "W_gen_info[0]")
                         .Define("W_gen_eta", "(double) W_gen_info[1]")
                         .Define("W_gen_phi", "(double) W_gen_info[2]")
                         .Define("W_gen_mass", "W_gen_info[3]")
                         .Define("W_gen_pdgId", "(int) W_gen_info[4]")
                         .Define("W_gen_status", "(int) W_gen_info[5]")
                         .Define("WDaughter1_gen_pt", "W_gen_info[6]")
                         .Define("WDaughter1_gen_eta", "(double) W_gen_info[7]")
                         .Define("WDaughter1_gen_phi", "(double) W_gen_info[8]")
                         .Define("WDaughter1_gen_mass", "W_gen_info[9]")
                         .Define("WDaughter1_gen_pdgId", "(int) W_gen_info[10]")
                         .Define("WDaughter1_gen_status", "(int) W_gen_info[11]")
                         .Define("WDaughter2_gen_pt", "W_gen_info[12]")
                         .Define("WDaughter2_gen_eta", "(double) W_gen_info[13]")
                         .Define("WDaughter2_gen_phi", "(double) W_gen_info[14]")
                         .Define("WDaughter2_gen_mass", "W_gen_info[15]")
                         .Define("WDaughter2_gen_pdgId", "(int) W_gen_info[16]")
                         .Define("WDaughter2_gen_status", "(int) W_gen_info[17]")
                         .Define("Wdaughter12_gen_dr", "DeltaR(WDaughter1_gen_eta, WDaughter2_gen_eta, WDaughter1_gen_phi, WDaughter2_gen_phi)")
                         .Define("trueLeptonicW", "(int) W_gen_info[18]")
                         .Define("trueLeptonicMode", Form("leptonicCheck(\"%s\", trueLeptonicT, trueLeptonicW)", sample.c_str()))
                         .Define("t_bkg_idx", Form("t_bkg_idx(\"%s\", nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)", sample.c_str()))
                         .Define("W_bkg_idx", Form("W_bkg_idx(\"%s\", nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags, t_bkg_idx)", sample.c_str()));
     }

     // ---------------------------------------------------------
     //               Save rdf before any cuts
     // ---------------------------------------------------------
   
     // TString outputFileNC = "RDF_" + sample + "_nocuts_" + testNum + ".root";
     // const char *stdOutputFileNC = outputFileNC;
     // cout << "------------------------------------------------" << endl
     //           << ">>> Saving original Snapshot..." << endl;
     // truth.Snapshot("Events", stdOutputFileNC);
     // cout << "Output File: " << outputFileNC << endl
     // << "-------------------------------------------------" << endl;

     // ---------------------------------------------------------
     //                    MET Filters
     // ---------------------------------------------------------
     

     auto METptFilters = truth.Filter("MET_pt > 60", "Pass MET > 60");

     // ---------------------------------------------------------
     //                    Lepton Filters
     // ---------------------------------------------------------

     auto LepDefs = METptFilters.Define("Electron_cutBasedIdNoIso_tight", Form("Electron_cutBasedIdNoIso_tight(\"%s\", nElectron, Electron_vidNestedWPBitmap)", sample.c_str()))
                        .Define("TPassMu", "abs(Muon_eta)<2.4 && (Muon_highPtId==2) && Muon_miniIsoId>=3")
                        .Define("TPassEl", "(abs(Electron_eta)<2.5) && (Electron_cutBasedIdNoIso_tight==1) && Electron_miniPFRelIso_all<0.1")
                        .Define("VetoMu", "TPassMu && (Muon_pt>25)")
                        .Define("VetoEl", "TPassEl && (Electron_pt>25)")
                        .Define("SignalIsoMu", "TPassMu && (Muon_pt>=55)")
                        .Define("SignalIsoEl", "TPassEl && (Electron_pt>=80)")
                        .Define("nVetoLep", "(int) (Sum(VetoMu)+Sum(VetoEl))")
                        .Define("VMuon_pt", "Muon_pt[VetoMu == true]")
                        .Define("VMuon_eta", "Muon_eta[VetoMu == true]")
                        .Define("VMuon_phi", "Muon_phi[VetoMu == true]")
                        .Define("VMuon_mass", "Muon_mass[VetoMu == true]")
                        .Define("VElectron_pt", "Electron_pt[VetoEl == true]")
                        .Define("VElectron_eta", "Electron_eta[VetoEl == true]")
                        .Define("VElectron_phi", "Electron_phi[VetoEl == true]")
                        .Define("VElectron_mass", "Electron_mass[VetoEl == true]")
                        .Define("VMuon_P4", "fVectorConstructor(VMuon_pt,VMuon_eta,VMuon_phi,VMuon_mass)")
                        .Define("VElectron_P4", "fVectorConstructor(VElectron_pt,VElectron_eta,VElectron_phi,VElectron_mass)")
                        .Define("VMuon_jetIdx", "Muon_jetIdx[VetoMu == true]")
                        .Define("VMuon_miniIsoId", "Muon_miniIsoId[VetoMu]")
                        .Define("VMuon_miniIso", "Muon_miniPFRelIso_all[VetoMu]") // for assign_leps. might not needed
                        .Define("VElectron_jetIdx", "Electron_jetIdx[VetoEl]")
                        .Define("VElectron_miniIso", "Electron_miniPFRelIso_all[VetoEl]")
                        .Define("nSignalIsoMu", "(int) Sum(SignalIsoMu)")
                        .Define("nSignalIsoEl", "(int) Sum(SignalIsoEl)")
                        .Define("VetoIsoMu", "(VMuon_pt<55)")
                        .Define("VetoIsoEl", "(VElectron_pt<80)")
                        .Define("nVetoIsoLep", "(int) (Sum(VetoIsoMu)+Sum(VetoIsoEl))");
     
     auto LepSelect = LepDefs.Define("isMu", "(nMuon>0) && HLT_Mu50 && (nSignalIsoMu==1) && (nVetoIsoLep==0) && (nElectron == 0 || nSignalIsoEl == 0)")
                        .Define("isEl", "(nElectron>0) && (HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165) && (nSignalIsoEl==1) && (nVetoIsoLep==0) && (nMuon == 0 || nSignalIsoMu == 0)")
                        .Filter("isMu || isEl", "Event is either muon or electron");
                        
     auto LepAssign = LepSelect.Define("assignleps", "assign_leps(isMu,isEl,SignalIsoMu,SignalIsoEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_miniPFRelIso_all,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_miniPFRelIso_all)")
                        .Define("lepton_pt", "assignleps[0]")
                        .Define("lepton_eta", "assignleps[1]")
                        .Define("lepton_phi", "assignleps[2]")
                        .Define("lepton_mass", "assignleps[3]")
                        .Define("lepton_miniIso", "assignleps[4]")
                        .Filter("isMu || MET_pt>((130/1.5)*DeltaPhi(lepton_phi, MET_phi)-130)", "Electron Triangle Cut");

     // --------------------------------------------------------
     // 		      JET SELECTION w/ cleaning
     // --------------------------------------------------------

     auto CleanJets = LepAssign.Define("Jet_P4", "fVectorConstructor(Jet_pt,Jet_eta,Jet_phi,Jet_mass)")
                          .Define("cleanJets", "cleanJets(Jet_P4,Jet_rawFactor,VMuon_P4,VMuon_jetIdx,VElectron_P4,VElectron_jetIdx)")
                          .Define("cleanJet_pt", "cleanJets[0]")
                          .Define("cleanJet_eta", "cleanJets[1]")
                          .Define("cleanJet_phi", "cleanJets[2]")
                          .Define("cleanJet_mass", "cleanJets[3]")
                          .Define("cleanJet_rawFactor", "cleanJets[4]")
                          .Define("DR_lepJets","DeltaR_VecAndFloat(cleanJet_eta,cleanJet_phi,lepton_eta,lepton_phi)")
                          .Define("ptrel_lepJets","ptRel(cleanJet_pt,cleanJet_eta,cleanJet_phi,cleanJet_mass,lepton_pt,lepton_eta,lepton_phi,lepton_mass)")
                          .Define("goodcleanJets", "cleanJet_pt > 30 && abs(cleanJet_eta) < 2.4 && Jet_jetId > 1 && (DR_lepJets > 0.4 || ptrel_lepJets > 20)")
                          .Define("NJets_central", "(int) Sum(goodcleanJets)")
                          .Define("gcJet_pt", "cleanJet_pt[goodcleanJets == true]")
                          .Define("gcJet_eta", "cleanJet_eta[goodcleanJets == true]")
                          .Define("gcJet_phi", "cleanJet_phi[goodcleanJets == true]")
                          .Define("gcJet_mass", "cleanJet_mass[goodcleanJets == true]")
                          .Define("gcJet_DeepFlav", "Jet_btagDeepFlavB[goodcleanJets == true]")
                          .Define("gcJet_DeepFlavL", "gcJet_DeepFlav > 0.0490") //0.2783
                          .Define("NJets_DeepFlavL", "(int) Sum(gcJet_DeepFlavL)")
                          .Define("gcJet_DeepFlavL_pt", "gcJet_pt[gcJet_DeepFlavL == true]")
                          .Define("gcJet_DeepFlavL_eta", "gcJet_eta[gcJet_DeepFlavL == true]")
                          .Define("gcJet_DeepFlavL_phi", "gcJet_phi[gcJet_DeepFlavL == true]")
                          .Define("gcJet_DeepFlavL_mass", "gcJet_mass[gcJet_DeepFlavL == true]")
                          .Define("goodcleanForwardJets", "cleanJet_pt > 30 && abs(cleanJet_eta) >= 2.4 && Jet_jetId > 1")
                          .Define("NJets_forward", "(int) Sum(goodcleanForwardJets)")
                          .Define("gcforwJet_pt", "cleanJet_pt[goodcleanForwardJets == true]")
                          .Define("gcforwJet_eta", "cleanJet_eta[goodcleanForwardJets == true]")
                          .Define("gcforwJet_phi", "cleanJet_phi[goodcleanForwardJets == true]")
                          .Define("gcforwJet_mass", "cleanJet_mass[goodcleanForwardJets == true]")
                          .Define("gcforwJet_DeepFlav", "Jet_btagDeepFlavB[goodcleanForwardJets == true]")
                          .Define("DR_lepFatJets","DeltaR_VecAndFloat(FatJet_eta,FatJet_phi,lepton_eta,lepton_phi)")
                          .Define("goodcleanFatJets", "FatJet_pt > 200 && abs(FatJet_eta) < 2.4 && FatJet_jetId > 1 && (DR_lepFatJets > 0.8)")
                          .Define("NFatJets", "(int) Sum(goodcleanFatJets)")
                          .Define("gcFatJet_pt", "FatJet_pt[goodcleanFatJets == true]")
                          .Define("gcFatJet_eta", "FatJet_eta[goodcleanFatJets == true]")
                          .Define("gcFatJet_phi", "FatJet_phi[goodcleanFatJets == true]")
                          .Define("gcFatJet_mass", "FatJet_mass[goodcleanFatJets == true]")
                          .Define("gcFatJet_sdmass", "FatJet_msoftdrop[goodcleanFatJets == true]")
                          .Define("DR_gcJets_central","DeltaR_VecAndFloat(gcJet_eta,gcJet_phi,lepton_eta,lepton_phi)")
                          .Define("DR_gcJets_DeepFlavL","DeltaR_VecAndFloat(gcJet_DeepFlavL_eta,gcJet_DeepFlavL_phi,lepton_eta,lepton_phi)")
                          .Define("DR_gcFatJets","DeltaR_VecAndFloat(gcFatJet_eta,gcFatJet_phi,lepton_eta,lepton_phi)")
                          .Define("OS_gcJets","DR_gcJets_central > TMath::Pi()/2")
                          .Define("SS_gcJets","DR_gcJets_central <= TMath::Pi()/2")
                          .Define("OS_gcJets_DeepFlavL","DR_gcJets_DeepFlavL > TMath::Pi()/2")
                          .Define("SS_gcJets_DeepFlavL","DR_gcJets_DeepFlavL <= TMath::Pi()/2")
                          .Define("OS_gcFatJets","DR_gcFatJets > TMath::Pi()/2")
                          .Define("NOS_gcJets_central","(int) Sum(OS_gcJets)")
                          .Define("NSS_gcJets_central","(int) Sum(SS_gcJets)")
                          .Define("NOS_gcJets_DeepFlavL","(int) Sum(OS_gcJets_DeepFlavL)")
                          .Define("NSS_gcJets_DeepFlavL","(int) Sum(SS_gcJets_DeepFlavL)")
                          .Define("NOS_gcFatJets","(int) Sum(OS_gcFatJets)")
                          .Define("gcOSFatJet_pt","gcFatJet_pt[OS_gcFatJets == true]")
                          .Define("gcOSFatJet_eta","gcFatJet_eta[OS_gcFatJets == true]")
                          .Define("gcOSFatJet_phi","gcFatJet_phi[OS_gcFatJets == true]")
                          .Define("gcOSFatJet_mass","gcFatJet_mass[OS_gcFatJets == true]")
                          .Define("gcOSJet_pt","gcJet_pt[OS_gcJets == true]")
                          .Define("gcOSJet_eta","gcJet_eta[OS_gcJets == true]")
                          .Define("gcOSJet_phi","gcJet_phi[OS_gcJets == true]")
                          .Define("gcOSJet_mass","gcJet_mass[OS_gcJets == true]")
                          .Define("gcSSJet_pt","gcJet_pt[SS_gcJets == true]")
                          .Define("gcSSJet_eta","gcJet_eta[SS_gcJets == true]")
                          .Define("gcSSJet_phi","gcJet_phi[SS_gcJets == true]")
                          .Define("gcSSJet_mass","gcJet_mass[SS_gcJets == true]");

     // ---------------------------------------------------------
     // 	  HT Calculation and Final Preselection Cut
     // ---------------------------------------------------------

     auto HT_calc = CleanJets.Define("Jet_HT","Sum(Jet_pt[goodcleanJets == true])") \
                          .Filter("Jet_HT > 250","Pass HT > 250")						\
                          .Filter("NFatJets > 0","Pass N good central AK8 > 0")
                          .Filter("NOS_gcFatJets > 0","Pass N good central other side AK8 > 0");

     //----------------------------------------------------------
     //       Uncomment from here to the bottom if starting from a preselection file!!
     //----------------------------------------------------------
     // //	auto HT_calc = rdf;

     // ---------------------------------------------------------
     // 		Post Preselection Analysis
     // ---------------------------------------------------------
     auto genttbar = HT_calc;

     if (!isSM && !isSE) {
          genttbar = HT_calc.Define("genttbarMass", Form("genttbarMassCalc(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status)",sample.c_str()));     
     }


     auto postPresel = genttbar.Define("lepton_lv", "lvConstructor(lepton_pt,lepton_eta,lepton_phi,lepton_mass)")
                           .Define("Jet_ST", "Jet_HT + lepton_pt + MET_pt")
                           .Define("FatJet_pt_1", "FatJet_pt[0]")
                           .Define("FatJet_pt_2", "FatJet_pt[1]") 
                           .Define("FatJet_sdMass", "FatJet_msoftdrop[goodcleanFatJets == true]")
                           .Define("FatJet_sdMass_1", "FatJet_sdMass[0]")
                           .Define("FatJet_sdMass_2", "FatJet_sdMass[1]")
                           .Define("dpak8_J", "FatJet_deepTag_QCDothers[goodcleanFatJets == true]")
                           .Define("dpak8_J_1", "dpak8_J[0]")
                           .Define("dpak8_J_2", "dpak8_J[1]")
                           .Define("raw_dpak8_T", "(FatJet_deepTag_TvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_TvsQCD)")
                           .Define("dpak8_T", "raw_dpak8_T[goodcleanFatJets == true]")
                           .Define("dpak8_T_1", "dpak8_T[0]")
                           .Define("dpak8_T_2", "dpak8_T[1]")
                           .Define("raw_dpak8_W", "(FatJet_deepTag_WvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_WvsQCD)")
                           .Define("dpak8_W", "raw_dpak8_W[goodcleanFatJets == true]")
                           .Define("dpak8_W_1", "dpak8_W[0]")
                           .Define("dpak8_W_2", "dpak8_W[1]")
                           .Define("dpak8_tag", "maxFxn(dpak8_J,dpak8_T,dpak8_W)")
                           .Define("dpak8_tag_1", "dpak8_tag[0]")
                           .Define("dpak8_tag_2", "dpak8_tag[1]")
                           .Define("nJ_dpak8", "Sum(dpak8_tag == 0)")
                           .Define("nT_dpak8", "Sum(dpak8_tag == 1)")
                           .Define("nW_dpak8", "Sum(dpak8_tag == 2)")
                           .Define("pNet_J", "FatJet_particleNet_QCD[goodcleanFatJets == true]")
                           .Define("pNet_J_1", "pNet_J[0]")
                           .Define("pNet_J_2", "pNet_J[1]")
                           .Define("raw_pNet_T", "(FatJet_particleNet_TvsQCD * FatJet_particleNet_QCD) / (1 - FatJet_particleNet_TvsQCD)")
                           .Define("pNet_T", "raw_pNet_T[goodcleanFatJets == true]")
                           .Define("pNet_T_alt", "FatJet_particleNet_TvsQCD[goodcleanFatJets == true]")
                           .Define("pNet_T_1", "pNet_T[0]")
                           .Define("pNet_T_2", "pNet_T[1]")
                           .Define("raw_pNet_W", "(FatJet_particleNet_WvsQCD * FatJet_particleNet_QCD) / (1 - FatJet_particleNet_WvsQCD)")
                           .Define("pNet_W", "raw_pNet_W[goodcleanFatJets == true]")
                           .Define("pNet_W_alt", "FatJet_particleNet_WvsQCD[goodcleanFatJets == true]")
                           .Define("pNet_W_1", "pNet_W[0]")
                           .Define("pNet_W_2", "pNet_W[1]")
                           .Define("pNet_tag", "maxFxn(pNet_J,pNet_T,pNet_W)")
                           .Define("OSpNet_tag", "pNet_tag[OS_gcFatJets==true]")
                           .Define("pNet_tag_alt", "JetDiscriminator(pNet_T_alt, pNet_W_alt)")
                           .Define("OSpNet_tag_alt", "pNet_tag_alt[OS_gcFatJets==true]")
                           .Define("pNet_tag_1", "pNet_tag[0]")
                           .Define("pNet_tag_2", "pNet_tag[1]")
                           .Define("nJ_pNet", "Sum(pNet_tag == 0)")
                           .Define("nT_pNet", "Sum(pNet_tag == 1)")
                           .Define("nW_pNet", "Sum(pNet_tag == 2)")
                           .Define("raw_tau21", "(FatJet_tau2 / FatJet_tau1)")
                           .Define("tau21", "raw_tau21[goodcleanFatJets == true]")
                           .Define("tau21_1", "tau21[0]")
                           .Define("tau21_2", "tau21[1]")
                           .Define("minDR_ptRel_lead_lepAK8", "minDR_ptRel_lead_calc(gcFatJet_pt,gcFatJet_eta,gcFatJet_phi, gcFatJet_mass,lepton_lv)")
                           .Define("minDR_lep_FatJet", "minDR_ptRel_lead_lepAK8[0]")
                           .Define("ptRel_lep_FatJet", "minDR_ptRel_lead_lepAK8[1]")
                           .Define("minDR_leadAK8otherAK8", "minDR_ptRel_lead_lepAK8[2]")
                           .Define("minDR_ptRel_lead_lepAK4", "minDR_ptRel_lead_calc(gcJet_pt,gcJet_eta,gcJet_phi, gcJet_mass,lepton_lv)")
                           .Define("minDR_lep_Jet", "minDR_ptRel_lead_lepAK4[0]")
                           .Define("ptRel_lep_Jet", "minDR_ptRel_lead_lepAK4[1]")
                           .Define("W_lv", "W_reco(MET_pt,MET_phi,lepton_lv)")
                           .Define("W_pt", "W_lv.Pt()")
                           .Define("W_eta", "W_lv.Eta()")
                           .Define("W_phi", "W_lv.Phi()")
                           .Define("W_mass", "W_lv.M()")
                           .Define("W_MT", "sqrt(2*lepton_pt*MET_pt*(1-cos(lepton_phi - MET_phi)))")
                           .Define("minMlj_output", Form("minM_lep_jet_calc(\"%s\", gcJet_pt, gcJet_eta, gcJet_phi, gcJet_mass, lepton_lv)",sample.c_str()))
                           .Define("DR_W_lep", "dR_Wt_Calc(W_lv,lepton_lv)")
                           .Define("minM_lep_Jet", "minMlj_output[0]")
                           .Define("minM_lep_Jet_jetID", "(int) minMlj_output[1]")
                           .Define("leptonicParticle", "isLeptonic_X(minM_lep_Jet)")
                           .Define("t_output", "t_reco(leptonicParticle,gcJet_pt,gcJet_eta,gcJet_phi,gcJet_mass, W_lv,minM_lep_Jet,minM_lep_Jet_jetID)")
                           .Define("t_pt", "t_output[0]")
                           .Define("t_eta", "t_output[1]")
                           .Define("t_phi", "t_output[2]")
                           .Define("t_mass", "t_output[3]")
                           .Define("DR_W_b", "t_output[4]")
                           .Define("t_lv", "lvConstructor(t_pt,t_eta,t_phi,t_mass)")
                           .Define("Bprime_output", Form("BPrime_reco_new(t_lv,W_lv,\"%s\",minM_lep_Jet,minM_lep_Jet_jetID,NOS_gcJets_central,NOS_gcJets_DeepFlavL,NSS_gcJets_central,NSS_gcJets_DeepFlavL,SS_gcJets_DeepFlavL,OS_gcJets_DeepFlavL,gcOSFatJet_pt,gcOSFatJet_eta,gcOSFatJet_phi,gcOSFatJet_mass,OSpNet_tag,gcOSJet_pt,gcOSJet_eta,gcOSJet_phi,gcOSJet_mass,gcSSJet_pt,gcSSJet_eta,gcSSJet_phi,gcSSJet_mass)",sample.c_str()))
                           .Define("Bprime_mass", "Bprime_output[0]")
                           .Define("Bprime_pt", "Bprime_output[1]")
                           .Define("Bprime_eta", "Bprime_output[2]")
                           .Define("Bprime_phi", "Bprime_output[3]")
                           .Define("Bprime_DR", "Bprime_output[4]")
                           .Define("Bprime_ptbal", "Bprime_output[5]")
                           .Define("Bprime_chi2", "Bprime_output[6]")
                           //.Define("BPrime_lv", "lvConstructor(Bprime_pt,Bprime_eta,Bprime_phi,Bprime_mass)")
                           .Define("Bdecay_obs", "Bprime_output[7]");
                         
     // -------------------------------------------------
     // 		Save Snapshot to file
     // -------------------------------------------------

     cout << "-------------------------------------------------" << endl
          << ">>> Saving " << sample << " Snapshot..." << endl;
     TString finalFile = "RDF_" + sample + "_finalsel_" + testNum + ".root";
     const char *stdfinalFile = finalFile;

     auto ColNames = postPresel.GetColumnNames();
     vector<string> snapCol;
     int i = 0;
     for (auto &&ColName : ColNames)
     {
          TString colName = ColName;
          if (!colName.Contains("P4") && colName != "cleanJets" && !colName.BeginsWith("PPS") && !colName.BeginsWith("Proton") && !colName.BeginsWith("L1") && !colName.BeginsWith("Gen") && !colName.BeginsWith("Soft") && !colName.BeginsWith("fixed") && !colName.BeginsWith("Sub") && !colName.BeginsWith("LHE") && !colName.BeginsWith("Raw") && !colName.BeginsWith("Calo") && !colName.BeginsWith("Chs") && !colName.BeginsWith("Corr") && !colName.BeginsWith("Fsr") && !colName.BeginsWith("Iso") && !colName.BeginsWith("Tau") && !colName.BeginsWith("SV") && !colName.BeginsWith("Puppi") && !colName.BeginsWith("Jet") && !colName.BeginsWith("FatJet") && !colName.BeginsWith("Photon") && !colName.BeginsWith("Low") && !colName.BeginsWith("HLT") && !colName.BeginsWith("HT") && !colName.BeginsWith("Muon") && !colName.BeginsWith("Electron") && !colName.BeginsWith("boosted") && !colName.BeginsWith("Deep") && !colName.BeginsWith("Flag") && colName != "Bprime_gen_info" && colName != "t_gen_info" && colName != "W_gen_info" && colName != "assignleps" && colName != "t_output" && colName != "Bprime_output" && !colName.BeginsWith("MET") && !colName.BeginsWith("Other") && !colName.BeginsWith("PS") && !colName.BeginsWith("PV") && !colName.BeginsWith("SV") && !colName.BeginsWith("Pile") && !colName.BeginsWith("Tk") && !colName.BeginsWith("Trig") && !colName.BeginsWith("btag") && !colName.BeginsWith("event") && !colName.BeginsWith("gen") && !colName.BeginsWith("run") && !colName.BeginsWith("lum") && !colName.BeginsWith("nCorr") && !colName.BeginsWith("nElectron") && !colName.BeginsWith("nFatJet") && !colName.BeginsWith("nFsr") && !colName.BeginsWith("nGen") && !colName.BeginsWith("nIso") && !colName.BeginsWith("nJet") && !colName.BeginsWith("nLHE") && !colName.BeginsWith("nLow") && !colName.BeginsWith("nMuon") && !colName.BeginsWith("nOther") && !colName.BeginsWith("nPS") && !colName.BeginsWith("nPhoton") && !colName.BeginsWith("nSV") && !colName.BeginsWith("nSub") && !colName.BeginsWith("nTau") && !colName.BeginsWith("nTrig") && !colName.BeginsWith("nboosted"))
          {
               string name = colName.Data();
               snapCol.push_back(name);
               i++;
          }
     }
     cout << "Number of Columns in Snapshot: " << i << endl;

     postPresel.Snapshot("Events", stdfinalFile, snapCol);
     cout << "Output File: " << finalFile << endl
          << "-------------------------------------------------" << endl;

     time.Stop();
     time.Print();

     cout << "Cut statistics:" << endl;
     postPresel.Report()->Print();

     // cout << "Adding Counter tree to the file:" << endl;
     // auto rdf_runs = ROOT::RDataFrame("Runs", files);
     // ROOT::RDF::RSnapshotOptions opts;
     // opts.fMode = "UPDATE";
     // rdf_runs.Snapshot("Runs", stdfinalFile, rdf_runs.GetColumnNames(), opts);

     cout << "Done!" << endl;
}
