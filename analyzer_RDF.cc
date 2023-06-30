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

  // -------------------------------------------------------
  //               Flags and First Filter
  // -------------------------------------------------------
  // Twiki with reccommended ultralegacy values
  auto rdf_input = ROOT::RDataFrame("Events", files); // Initial data

  auto truth = rdf_input.Define("Bprime_gen_info", Form("Bprime_gen_info(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status, GenPart_statusFlags)", sample.c_str()))
                 .Define("Bprime_gen_pt", "Bprime_gen_info[0]")
                 .Define("Bprime_gen_eta", "(double) Bprime_gen_info[1]")
                 .Define("Bprime_gen_phi", "(double) Bprime_gen_info[2]")
                 .Define("Bprime_gen_mass", "Bprime_gen_info[3]")
                 .Define("Bprime_gen_pdgId", "(int) Bprime_gen_info[4]")
                 .Define("Bprime_gen_status", "(int) Bprime_gen_info[5]")
                 .Define("t_gen_info", Form("t_gen_info(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status)", sample.c_str()))
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
                 .Define("trueLeptonicW", "(int) W_gen_info[18]")
                 .Define("trueLeptonicMode", Form("leptonicCheck(\"%s\", trueLeptonicT, trueLeptonicW)", sample.c_str()))
                 .Define("t_bkg_idx", Form("t_bkg_idx(\"%s\", nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)", sample.c_str()))
                 .Define("W_bkg_idx", Form("W_bkg_idx(\"%s\", nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags, t_bkg_idx)", sample.c_str()))
                 .Define("Electron_cutBasedIdNoIso_tight", Form("Electron_cutBasedIdNoIso_tight(\"%s\", nElectron, Electron_vidNestedWPBitmap)", sample.c_str()));
  //  cout << "Number of Events passing Preselection (HT Cut): " << HT_calc.Count().GetValue() << endl;

  // ---------------------------------------------------------
  //               Save rdf before any cuts
  // ---------------------------------------------------------

  // TString outputFileNC = "RDF_" + sample + "_nocuts_" + testNum + ".root";
  // const char *stdOutputFileNC = outputFileNC;
  // cout << "------------------------------------------------" << endl
  //           << ">>> Saving original Snapshot..." << endl;
  // rdf.Snapshot("Events", stdOutputFileNC);
  // cout << "Output File: " << outputFileNC << endl
  //           << "-------------------------------------------------" << endl;

  auto METfilters = truth.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1", "MET Filters")
                        .Filter("MET_pt > 50", "Pass MET > 50")
                        .Filter("nJet > 0 && nFatJet > 0", "Event has jets");

  // ---------------------------------------------------------
  //                    Lepton Filters
  // ---------------------------------------------------------

  auto LepDefs = METfilters.Define("TPassMu", "abs(Muon_eta)<2.4 && (Muon_highPtId==2)")
                     .Define("TPassEl", "(abs(Electron_eta)<2.5) && (abs(Electron_deltaEtaSC+Electron_eta)<2.5) && (Electron_cutBasedIdNoIso_tight==1)")
                     .Define("VetoMu", "TPassMu && (Muon_pt>25)")
                     .Define("VetoEl", "TPassEl && (Electron_pt>25)")
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
                     .Define("VElectron_jetIdx", "Electron_jetIdx[VetoEl]")
                     .Define("VElectron_miniIso", "Electron_miniPFRelIso_all[VetoEl]");

  auto CleanJets = LepDefs.Define("Jet_P4", "fVectorConstructor(Jet_pt,Jet_eta,Jet_phi,Jet_mass)")
                       .Define("cleanJets", "cleanJets(Jet_P4,Jet_rawFactor,VMuon_P4,VMuon_jetIdx,VElectron_P4,VElectron_jetIdx)")
                       .Define("cleanJet_pt", "cleanJets[0]")
                       .Define("cleanJet_eta", "cleanJets[1]")
                       .Define("cleanJet_phi", "cleanJets[2]")
                       .Define("cleanJet_mass", "cleanJets[3]")
                       .Define("cleanJet_rawFactor", "cleanJets[4]")
                       .Define("goodcleanJets", "cleanJet_pt > 30 && abs(cleanJet_eta) < 2.4 && Jet_jetId > 1")
                       .Define("NJets_central", "(int) Sum(goodcleanJets)")
                       .Define("gcJet_pt", "cleanJet_pt[goodcleanJets == true]")
                       .Define("gcJet_eta", "cleanJet_eta[goodcleanJets == true]")
                       .Define("gcJet_phi", "cleanJet_phi[goodcleanJets == true]")
                       .Define("gcJet_mass", "cleanJet_mass[goodcleanJets == true]")
                       .Define("gcJet_DeepFlav", "Jet_btagDeepFlavB[goodcleanJets == true]")
                       .Define("gcJet_DeepFlavM", "gcJet_DeepFlav > 0.2783")
                       .Define("NJets_DeepFlavM", "(int) Sum(gcJet_DeepFlavM)")
                       .Define("goodcleanForwardJets", "cleanJet_pt > 30 && abs(cleanJet_eta) >= 2.4 && Jet_jetId > 1")
                       .Define("NJets_forward", "(int) Sum(goodcleanForwardJets)")
                       .Define("gcforwJet_pt", "cleanJet_pt[goodcleanForwardJets == true]")
                       .Define("gcforwJet_eta", "cleanJet_eta[goodcleanForwardJets == true]")
                       .Define("gcforwJet_phi", "cleanJet_phi[goodcleanForwardJets == true]")
                       .Define("gcforwJet_mass", "cleanJet_mass[goodcleanForwardJets == true]")
                       .Define("gcforwJet_DeepFlav", "Jet_btagDeepFlavB[goodcleanForwardJets == true]")
                       .Define("dR_LIM_AK4", "(float) 0.4")
                       .Define("ptrel25", "25")
                       .Define("VMuon_2Dcut_ptrel25", "cut_ptrel(dR_LIM_AK4, ptrel25, VMuon_P4, NJets_central, gcJet_eta, gcJet_phi, gcJet_pt, gcJet_mass)")
                       .Define("VElectron_2Dcut_ptrel25", "cut_ptrel(dR_LIM_AK4, ptrel25, VElectron_P4, NJets_central, gcJet_eta, gcJet_phi, gcJet_pt, gcJet_mass)")
                       .Define("ptrel40", "40")
                       .Define("VMuon_2Dcut_ptrel40", "cut_ptrel(dR_LIM_AK4, ptrel40, VMuon_P4, NJets_central, gcJet_eta, gcJet_phi, gcJet_pt, gcJet_mass)")
                       .Define("VElectron_2Dcut_ptrel40", "cut_ptrel(dR_LIM_AK4, ptrel40, VElectron_P4, NJets_central, gcJet_eta, gcJet_phi, gcJet_pt, gcJet_mass)");

  // auto LepSelect = CleanJets.Define("SignalMu", "TPassMu && (Muon_pt>55)")			\
  //     .Define("SignalEl", "TPassEl && (Electron_pt>80)")			\
  //     .Define("nSignalMu", "(int) Sum(SignalMu)")			\
  //     .Define("nSignalEl", "(int) Sum(SignalEl)")			\
  //     .Define("isMu","(nMuon>0) && HLT_Mu50 && (nSignalMu==1) && (nVetoLep==0)") \
  //     .Define("isEl","(nElectron>0) && (HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165) && (nSignalEl==1) && (nVetoLep==0)") \
  //   .Filter("isMu || isEl","Event is either muon or electron");

  // auto Lep_df1 = Lep_df0.Define("assignleps","assign_leps(isMu,isEl,TPassMu,TPassEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_miniPFRelIso_all,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_miniPFRelIso_all)") \
  //   .Define("lepton_pt","assignleps[0]")				\
  //   .Define("lepton_eta","assignleps[1]")				\
  //   .Define("lepton_phi","assignleps[2]")				\
  //   .Define("lepton_mass","assignleps[3]")				\
  //   .Define("lepton_miniIso","assignleps[4]");

  // --------------------------------------------------------
  // 		      JET SELECTION w/ cleaning
  // --------------------------------------------------------

  // auto jet_ft0 = Lep_df1.Filter("nJet > 0 && nFatJet > 0","Event has jets");
  // //  cout << "Number of Events with at least one AK4 and AK8 Jet: " << jet_ft0.Count().GetValue() << endl;

  // auto jet_df0 = jet_ft0.Define("goodJets","Jet_pt > 30 && abs(Jet_eta) < 2.4 && Jet_jetId > 1") \
  //   .Define("goodcleanFatJets","cleanJets(FatJet_pt,FatJet_mass,goodFatJets,FatJet_eta,FatJet_phi,\
  // 			      				            lepton_pt,lepton_mass,lepton_eta,lepton_phi,dR_LIM_AK8)")\
  //   .Define("dR_LIM_AK8","(float) 0.8")					\
  //   .Define("goodFatJets","FatJet_jetId > 1 && abs(FatJet_eta) < 2.4 && FatJet_pt > 200") \
  //   .Define("NFatJets","(int) Sum(goodcleanFatJets)")			\
  //   .Define("gcFatJet_pt","FatJet_pt[goodcleanFatJets == true]")	\
  //   .Define("gcFatJet_eta","FatJet_eta[goodcleanFatJets == true]")	\
  //   .Define("gcFatJet_phi","FatJet_phi[goodcleanFatJets == true]")	\
  //   .Define("gcFatJet_mass","FatJet_mass[goodcleanFatJets == true]")	\
  //   .Define("gcFatJet_msoftdrop","FatJet_msoftdrop[goodcleanFatJets == true]");

  // // ---------------------------------------------------------
  // // 	  HT Calculation and Final Preselection Cut
  // // ---------------------------------------------------------
  // auto HT_calc = jet_df0.Define("Jet_HT","Sum(Jet_pt[goodcleanJets == true])") \
  //   .Filter("Jet_HT > 250","Pass HT > 250")						\
  //   .Filter("NFatJets > 0","Pass N good central AK8 > 0");

  // // ---------------------------------------------------------
  // //    Uncomment to save seperate Preselection .root file
  // // ---------------------------------------------------------
  // //TString outputFilePS = "RDF_"+sample+"_presel_"+testNum+".root";
  // //const char* stdOutputFilePS = outputFilePS;
  // //cout << "------------------------------------------------" << endl << ">>> Saving Preselection Snapshot..." << endl;
  // //HT_calc.Snapshot("Events", stdOutputFilePS);
  // //cout << "Output File: " << outputFilePS << endl << "-------------------------------------------------" << endl;
  // // }
  // //----------------------------------------------------------
  // //       Uncomment from here to the bottom if starting from a preselection file!!
  // //----------------------------------------------------------
  // //	auto HT_calc = rdf;

  // // ---------------------------------------------------------
  // // 		Post Preselection Analysis
  // // ---------------------------------------------------------
  // auto postPresel = HT_calc.Define("genttbarMass", Form("genttbarMassCalc(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, \
  // 	GenPart_pt, GenPart_phi, GenPart_eta,			\
  // 	GenPart_genPartIdxMother, GenPart_status)", sample.c_str()))			\
  //   .Define("lepton_lv","lvConstructor(lepton_pt,lepton_eta,lepton_phi,lepton_mass)") \
  //   .Define("Jets_lv","fVectorConstructor(gcJet_pt,gcJet_eta,gcJet_phi,gcJet_mass)") \
  //   .Define("FatJet_lv","fVectorConstructor(gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass)") \
  //   .Define("Jet_ST","Jet_HT + lepton_pt + MET_pt")			\
  //   .Define("FatJet_pt_1","FatJet_pt[0]")				\
  //   .Define("FatJet_pt_2","FatJet_pt[1]")				\
  //   .Define("FatJet_sdMass","FatJet_msoftdrop[goodcleanFatJets == true]") \
  //   .Define("FatJet_sdMass_1","FatJet_sdMass[0]")			\
  //   .Define("FatJet_sdMass_2","FatJet_sdMass[1]")			\
  //   .Define("dpak8_J","FatJet_deepTag_QCDothers[goodcleanFatJets == true]") \
  //   .Define("dpak8_J_1","dpak8_J[0]")					\
  //   .Define("dpak8_J_2","dpak8_J[1]")					\
  //   .Define("raw_dpak8_T","(FatJet_deepTag_TvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_TvsQCD)") \
  //   .Define("dpak8_T","raw_dpak8_T[goodcleanFatJets == true]")		\
  //   .Define("dpak8_T_1","dpak8_T[0]")					\
  //   .Define("dpak8_T_2","dpak8_T[1]")					\
  //   .Define("raw_dpak8_W","(FatJet_deepTag_WvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_WvsQCD)") \
  //   .Define("dpak8_W","raw_dpak8_W[goodcleanFatJets == true]")		\
  //   .Define("dpak8_W_1","dpak8_W[0]")					\
  //   .Define("dpak8_W_2","dpak8_W[1]")					\
  //   .Define("dpak8_tag","maxFxn(dpak8_J,dpak8_T,dpak8_W)")		\
  //   .Define("dpak8_tag_1","dpak8_tag[0]")				\
  //   .Define("dpak8_tag_2","dpak8_tag[1]")				\
  //   .Define("nJ_dpak8","Sum(dpak8_tag == 0)")				\
  //   .Define("nT_dpak8","Sum(dpak8_tag == 1)")				\
  //   .Define("nW_dpak8","Sum(dpak8_tag == 2)")				\
  //   .Define("pNet_J","FatJet_particleNet_QCD[goodcleanFatJets == true]") \
  //   .Define("pNet_J_1","pNet_J[0]")					\
  //   .Define("pNet_J_2","pNet_J[1]")					\
  //   .Define("raw_pNet_T","(FatJet_particleNet_TvsQCD * FatJet_particleNet_QCD) / (1 - FatJet_particleNet_TvsQCD)") \
  //   .Define("pNet_T","raw_pNet_T[goodcleanFatJets == true]")		\
  //   .Define("pNet_T_alt","FatJet_particleNet_TvsQCD[goodcleanFatJets == true]")            \
  //   .Define("pNet_T_1","pNet_T[0]")					\
  //   .Define("pNet_T_2","pNet_T[1]")					\
  //   .Define("raw_pNet_W","(FatJet_particleNet_WvsQCD * FatJet_particleNet_QCD) / (1 - FatJet_particleNet_WvsQCD)") \
  //   .Define("pNet_W","raw_pNet_W[goodcleanFatJets == true]")		\
  //   .Define("pNet_W_alt","FatJet_particleNet_WvsQCD[goodcleanFatJets == true]")
  //   .Define("pNet_W_1","pNet_W[0]")					\
  //   .Define("pNet_W_2","pNet_W[1]")					\
  //   .Define("pNet_tag","maxFxn(pNet_J,pNet_T,pNet_W)")			\
  //   .Define("pNet_tag_alt","JetDiscriminator(pNet_T_alt, pNet_W_alt)")  \
  //   .Define("pNet_tag_1","pNet_tag[0]")					\
  //   .Define("pNet_tag_2","pNet_tag[1]")					\
  //   .Define("nJ_pNet","Sum(pNet_tag == 0)")				\
  //   .Define("nT_pNet","Sum(pNet_tag == 1)")				\
  //   .Define("nW_pNet","Sum(pNet_tag == 2)")				\
  //   .Define("raw_tau21","(FatJet_tau2 / FatJet_tau1)")			\
  //   .Define("tau21","raw_tau21[goodcleanFatJets == true]")		\
  //   .Define("tau21_1","tau21[0]")					\
  //   .Define("tau21_2","tau21[1]")					\
  //   .Define("minDR_ptRel_lead_lepAK8","minDR_ptRel_lead_calc(gcFatJet_pt,gcFatJet_eta,gcFatJet_phi, \
  // 									   gcFatJet_mass,lepton_lv)")\
  //   .Define("minDR_lep_FatJet","minDR_ptRel_lead_lepAK8[0]")		\
  //   .Define("ptRel_lep_FatJet","minDR_ptRel_lead_lepAK8[1]")		\
  //   .Define("minDR_leadAK8otherAK8","minDR_ptRel_lead_lepAK8[2]")	\
  //   .Define("minDR_ptRel_lead_lepAK4","minDR_ptRel_lead_calc(gcJet_pt,gcJet_eta,gcJet_phi, \
  // 					      				   gcJet_mass,lepton_lv)")\
  //   .Define("minDR_lep_Jet","minDR_ptRel_lead_lepAK4[0]")		\
  //   .Define("ptRel_lep_Jet","minDR_ptRel_lead_lepAK4[1]")		\
  //   .Define("DR_lep_FatJets","DR_calc(gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass, \
  // 					      	    lepton_pt,lepton_eta, lepton_phi,lepton_mass)")\
  //   .Define("DR_lep_Jets","DR_calc(gcJet_pt,gcJet_eta,gcJet_phi,gcJet_mass, \
  // 					   	 lepton_pt,lepton_eta,lepton_phi,lepton_mass)")\
  //   .Define("W_lv","W_reco(MET_pt,MET_phi,lepton_lv)")			\
  //   .Define("W_pt", "W_lv.Pt()")
  //   .Define("W_eta", "W_lv.Eta()")
  //   .Define("W_phi", "W_lv.Phi()")
  //   .Define("W_mass", "W_lv.M()")
  //   .Define("W_MT", "sqrt(2*lepton_pt*MET_pt*(1-cos(lepton_phi - MET_phi)))")
  //   .Define("minMlj_output", Form("minM_lep_jet_calc(\"%s\", gcJet_pt, gcJet_eta, gcJet_phi, gcJet_mass, \
  // 	   lepton_lv)", sample.c_str()))							\
  //   .Define("DR_W_lep","dR_Wt_Calc(W_lv,lepton_lv)")			\
  //   .Define("minM_lep_Jet","minMlj_output[0]")				\
  //   .Define("minM_lep_Jet_jetID","(int) minMlj_output[1]")		\
  //   .Define("leptonicParticle","isLeptonic_X(minM_lep_Jet)")		\
  //   .Define("t_output","t_reco(leptonicParticle,gcJet_pt,gcJet_eta,gcJet_phi,gcJet_mass, \
  // 					 W_lv,minM_lep_Jet,minM_lep_Jet_jetID)")\
  //   .Define("t_pt","t_output[0]")					\
  //   .Define("t_eta","t_output[1]")					\
  //   .Define("t_phi","t_output[2]")					\
  //   .Define("t_mass","t_output[3]")					\
  //   .Define("DR_W_b","t_output[4]")					\
  //   .Define("t_lv","lvConstructor(t_pt,t_eta,t_phi,t_mass)")
  //   .Define("Bprime_output","BPrime_reco(t_lv,W_lv,leptonicParticle,\
  // 	  				       gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass,pNet_tag,gcFatJet_msoftdrop)")\
  //   .Define("Bprime_output_alt","BPrime_reco_alt(lepton_lv, t_lv,W_lv,leptonicParticle, gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass,pNet_tag,gcFatJet_msoftdrop)")
  //   .Define("Bprime_mass","Bprime_output[0]")				\
  //   .Define("Bprime_pt","Bprime_output[1]")				\
  //   .Define("Bprime_eta","Bprime_output[2]")				\
  //   .Define("Bprime_phi","Bprime_output[3]")				\
  //   .Define("Bprime_DR","Bprime_output[4]")				\
  //   .Define("Bprime_ptbal","Bprime_output[5]")				\
  //   .Define("Bprime_chi2","Bprime_output[6]")				\
  //   .Define("BPrime_lv","lvConstructor(Bprime_pt,Bprime_eta,Bprime_phi,Bprime_mass)") \
  //   .Define("isValidBDecay","Bprime_output[7]")				\
  //   .Define("taggedWbjetJet","Bprime_output[8]")			\
  //   .Define("taggedTjet","Bprime_output[9]")				\
  //   .Define("taggedWjet","Bprime_output[10]")				\
  //   .Define("dnn_scores",predictMLP,{"pNet_J_1","pNet_T_1","dpak8_T_1","FatJet_pt_1","FatJet_sdMass_1","tau21_1","nT_dpak8","nT_pNet","Jet_HT","Jet_ST","MET_pt","NJets_DeepFlavM","NJets_forward"}) \
  //   .Define("mlp_HT250_WJets","dnn_scores[0]")				\
  //   .Define("mlp_HT250_TTbar","dnn_scores[1]")				\
  //   .Define("mlp_HT250_Bprime","dnn_scores[2]")				\
  //   .Define("mlp_HT500_WJets","dnn_scores[3]")				\
  //   .Define("mlp_HT500_TTbar","dnn_scores[4]")				\
  //   .Define("mlp_HT500_Bprime","dnn_scores[5]")
  //   .Define("genFatJet_matching_sig", FatJet_matching_sig, {"goodcleanFatJets", "gcFatJet_eta", "gcFatJet_phi", "NFatJets", "FatJet_subJetIdx1", "nSubJet", "SubJet_hadronFlavour", "GenPart_pdgId", "daughterb_gen_eta", "daughterb_gen_phi", "tDaughter1_gen_eta", "tDaughter1_gen_phi", "tDaughter1_gen_pdgId", "tDaughter2_gen_eta", "tDaughter2_gen_phi", "tDaughter2_gen_pdgId", "WDaughter1_gen_eta", "WDaughter1_gen_phi", "WDaughter1_gen_pdgId", "WDaughter2_gen_eta", "WDaughter2_gen_phi", "WDaughter2_gen_pdgId"})
  //   .Define("genFatJet_matching_bkg", FatJet_matching_bkg, {"goodcleanFatJets", "gcFatJet_eta", "gcFatJet_phi", "NFatJets", "FatJet_subJetIdx1", "nSubJet", "SubJet_hadronFlavour", "nGenPart", "GenPart_pdgId", "GenPart_phi", "GenPart_eta", "GenPart_genPartIdxMother", "t_bkg_idx", "W_bkg_idx"});

  // // -------------------------------------------------
  // // 		Save Snapshot to file
  // // -------------------------------------------------

  cout << "-------------------------------------------------" << endl
            << ">>> Saving " << sample << " Snapshot..." << endl;
  TString finalFile = "RDF_" + sample + "_finalsel_" + testNum + ".root";
  const char *stdfinalFile = finalFile;

  auto colNames = CleanJets.GetColumnNames();
  vector<string> snapCol;
  int i = 0;
  for (auto &&colName : colNames)
  {
    if (colName != "VMuon_P4" && colName != "VElectron_P4" && colName != "Jet_P4" && colName != "cleanJets")
    {
      snapCol.push_back(colName);
      i++;
    }
  }
  cout << "Number of Columns in Snapshot: " << i << endl;

  CleanJets.Snapshot("Events", stdfinalFile, snapCol);
  cout << "Output File: " << finalFile << endl
            << "-------------------------------------------------" << endl;

  time.Stop();
  time.Print();
  cout << "Cut statistics:" << endl;
  CleanJets.Report()->Print();

  cout << "Adding Counter tree to the file:" << endl;
  auto rdf_runs = ROOT::RDataFrame("Runs", files); 
  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "UPDATE";
  rdf_runs.Snapshot("Runs", stdfinalFile, rdf_runs.GetColumnNames(), opts);

  cout << "Done!" << endl;
}
