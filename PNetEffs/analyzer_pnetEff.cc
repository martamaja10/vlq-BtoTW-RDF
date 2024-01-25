// --------------------------------------------------------------------------------------- //
// Implimentation of RDataFrame in C++.					                   //
// Comments on creating a singly produced VLQ search			                   //
// To Run on Command Line:   root -l callRDF.C\(\"Muon(OR)Electron\",\"testNumber\"\,\"root://cmsxrootd.fnal.gov//store/...file.root\")      //
// --------------------------------------------------------------------------------------- //

#define rdf_cxx
#include "analyzer_pnetEff.h"
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
#include "../../correctionlib/include/correction.h"

using namespace std;
using namespace ROOT::VecOps;
using correction::CorrectionSet;

void rdf::analyzer_pnetEff(TString testNum, TString jesvar)
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
  else cout << "MC extension tag (blank or ext) = " << era << endl;

  // -------------------------------------------------------
  //               correctionLib corrections
  // -------------------------------------------------------
  
  std::string yrstr, yr, jecyr, jeryr, jecver;
  if(year == "2016APV") {yrstr = "2016preVFP"; yr = "16"; jecyr = "UL16APV"; jeryr = "Summer20UL16APV_JRV3"; jecver = "V7";}
  else if(year == "2016") {yrstr = "2016postVFP"; yr = "16"; jecyr = "UL16"; jeryr = "Summer20UL16_JRV3"; jecver = "V7";}
  else if(year == "2017") {yrstr = "2017"; yr = "17"; jecyr = "UL17"; jeryr = "Summer19UL17_JRV2"; jecver = "V5";}
  else if(year == "2018") {yrstr = "2018"; yr = "18"; jecyr = "UL18"; jeryr = "Summer19UL18_JRV2"; jecver = "V5";}
  else std::cout << "ERROR: Can't parse the year to assign correctionLib json files. Expected 2016, 2016APV, 2017, or 2018. Got: " << year << std::endl;
  
  auto ak4corrset = CorrectionSet::from_file("../jsonpog-integration/POG/JME/"+yrstr+"_UL/jet_jerc.json"); 
  auto ak8corrset = CorrectionSet::from_file("../jsonpog-integration/POG/JME/"+yrstr+"_UL/fatJet_jerc.json"); 
  auto ak4ptres = ak4corrset->at(jeryr+"_MC_PtResolution_AK4PFchs"); std::cout << "\t loaded pt res" << std::endl;
  auto ak4jer = ak4corrset->at(jeryr+"_MC_ScaleFactor_AK4PFchs"); std::cout << "\t loaded jer" << std::endl;
  auto ak8corr = ak8corrset->compound().at("Summer19"+jecyr+"_"+jecver+"_MC_L1L2L3Res_AK8PFPuppi"); std::cout << "\t loaded fat jerc MC" << std::endl;
  auto ak8corrUnc = ak8corrset->at("Summer19"+jecyr+"_"+jecver+"_MC_Total_AK8PFPuppi"); std::cout << "\t loaded fat jec unc" << std::endl;

  auto cleanJets = [debug,jesvar,isMC,ak4ptres,ak4jer,ak8corr,ak8corrUnc](const RVec<TLorentzVector> &jt_p4, const RVec<float> &jt_rf, const RVec<float> &jt_murf, const RVec<float> &jt_area, const RVec<float> &jt_em, const RVec<int> &jt_id, const RVec<TLorentzVector> &genjt_p4, const RVec<int> &jt_genidx, const RVec<TLorentzVector> &mu_p4, const RVec<int> mu_jetid, const RVec<TLorentzVector> &el_p4, const RVec<int> &el_jetid, const float &rho, const float &met, const float &phi){
    RVec<float> cleanJetPt(jt_p4.size()), cleanJetEta(jt_p4.size()), cleanJetPhi(jt_p4.size()), cleanJetMass(jt_p4.size()), rawfact(jt_p4.size());
    string jervar = "nom";
    float jesuncmult = 0;
    if(jesvar == "JERup") jervar = "up";
    else if(jesvar == "JERdn") jervar = "down";
    else if(jesvar == "JECup") jesuncmult = 1.0;
    else if(jesvar == "JECdn") jesuncmult = -1.0;
    correction::CompoundCorrection::Ref jescorr;
    correction::Correction::Ref jescorrUnc;
    float drmax = 0.4;
    jescorr = ak8corr; 
    jescorrUnc = ak8corrUnc;
    float metx = met*cos(phi);
    float mety = met*sin(phi);
    if(met > 0 && debug) std::cout<< "Incoming met = " << met << ", phi = " << phi << std::endl;

    for(unsigned int ijet = 0; ijet < jt_p4.size(); ijet++){
      TLorentzVector jet = jt_p4[ijet];
      int jetid = jt_id[ijet];   
      float rf = jt_rf[ijet];

      if(ROOT::VecOps::Mean(jt_area) < 1.0 && met == 0){ // only clean leptons out of AK4 jets (AK8 will require DR > 0.8 from lepton)
	for (unsigned int imu = 0; imu < mu_p4.size(); imu++){
	  if(mu_jetid[imu] != ijet) continue;                      // only consider muons matched to this jet
	  if (jetid < 2 || jet.DeltaR(mu_p4[imu]) > 0.4) continue; // bad jet, or too far from muon
	  jet *= (1 - rf);                                         // first undo the JEC
	  jet -= mu_p4[imu];                                       // subtract muon if it's sensible
	  rf = 0;                                                  // indicate that this is the raw jet
	}
	for (unsigned int iel = 0; iel < el_p4.size(); iel++){     // same for electrons
	  if (el_jetid[iel] != ijet) continue;
	  if (jetid < 2 || jet.DeltaR(el_p4[iel]) > 0.4) continue;	  
	  jet *= (1 - rf); 
	  jet -= el_p4[iel]; 
	  rf = 0; 
	}
      }
      float jes = 1.0; float jesL1 = 1.0; float jer = 1.0; float unc = 1.0;
      jet = jet * (1 - rf);                                                         // rf = 0 if JEC undone above
      if(met > 0 && jt_em[ijet] > 0.9) continue;                                    // not these jets for MET	
      if(met > 0) jet *= (1 - jt_murf[ijet]);                                       // further correct raw to muon-substracted raw for T1.
      float rawpt = jet.Pt();
      jes = jescorr->evaluate({jt_area[ijet],jet.Eta(),rawpt,rho});                 // Data & MC get jes
      if(isMC){
	float res = ak4ptres->evaluate({jet.Eta(),rawpt*jes,rho});
	float sf = ak4jer->evaluate({jet.Eta(),jervar});
	bool smeared = false;                                                       // MC only gets a JER smear, one of 2 methods below:
	if(jt_genidx[ijet] > -1 && genjt_p4[jt_genidx[ijet]].Pt() > 0){	  
	  double dPt = fabs(genjt_p4[jt_genidx[ijet]].Pt() - rawpt*jes);
	  double dR = genjt_p4[jt_genidx[ijet]].DeltaR(jet);
	  if(dR < drmax && dPt < 3*rawpt*jes*res){
	    jer = max(0.0, 1.0 + (sf - 1.0)*(rawpt*jes - genjt_p4[jt_genidx[ijet]].Pt())/(rawpt*jes));
	    smeared = true;
	  }
	}
	if(!smeared){
	  TRandom3 rand(abs(static_cast<int>(jet.Phi()*1e4)));
	  jer = max(0.0, 1.0 + rand.Gaus(0, res)*sqrt(max(0.0, sf*sf - 1.0)));
	}
	unc = 1.0 + jesuncmult*(jescorrUnc->evaluate({jet.Eta(),rawpt*jes*jer}));   // MC gets JEC unc
      }	
      TLorentzVector jetL1 = jet*jesL1*jer*unc;
      jet = jet*jes*jer*unc;                                                        // evals to jes*1 for data.
      rf = 1.0 - 1.0/(jes*jer*unc);	
      if(jet.Pt() > 15){
	metx += (jetL1 - jet).Px();
	mety += (jetL1 - jet).Py();
      }
      cleanJetPt[ijet] = jet.Pt();
      cleanJetEta[ijet] = jet.Eta();
      cleanJetPhi[ijet] = jet.Phi();
      cleanJetMass[ijet] = jet.M();
      rawfact[ijet] = rf;
    }
    TVector2 corrmet(metx,mety);
    RVec<float> corrmets = {float(corrmet.Mod()),float(TVector2::Phi_mpi_pi(corrmet.Phi()))};
    
    RVec<RVec<float>> output = {cleanJetPt, cleanJetEta, cleanJetPhi, cleanJetMass, rawfact, corrmets};
    return output;
  };

  auto pnetWPs = [year](const RVec<float> &dnnT, const RVec<float> &dnnW){
    float wpT = 0.490; float wpW = 0.677;
    if(year == "2016"){wpT = 0.495; wpW = 0.668;}
    else if(year == "2017"){wpT = 0.581; wpW = 0.709;}
    else if(year == "2018"){wpT = 0.580; wpW = 0.700;}
    RVec<int> tag;
    for(int i=0; i<dnnT.size(); i++){
      if(dnnT[i] > wpT) tag.push_back(1);
      else if(dnnW[i] > wpW) tag.push_back(2);
      else tag.push_back(0);
    }
    return tag;
  };

  // -------------------------------------------------------
  //               Open Dataframe + MET Filters
  // -------------------------------------------------------

  auto rdf_input = ROOT::RDataFrame("Events", files); // Initial data
  
  auto METgeneralFilters = rdf_input.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1", "MET Filters")
    .Filter("nFatJet > 0", "Event > 1 AK8");

  // --------------------------------------------------------
  // 	       Golden JSON (Data) || GEN Info (MC)
  // --------------------------------------------------------
   
  auto truth = METgeneralFilters.Define("t_gen_info", Form("t_gen_info(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status)", sample.c_str()))
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
  
  auto LepSelect = LepDefs.Define("isMu", "(nMuon>0) && (nSignalIsoMu==1) && (nVetoIsoLep==0) && (nElectron == 0 || nSignalIsoEl == 0)")
    .Define("isEl", "(nElectron>0) && (nSignalIsoEl==1) && (nVetoIsoLep==0) && (nMuon == 0 || nSignalIsoMu == 0)")
    .Filter("isMu || isEl", "Event is either muon or electron");
  
  auto LepAssign = LepSelect.Define("assignleps", "assign_leps(isMu,isEl,SignalIsoMu,SignalIsoEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_miniPFRelIso_all,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_miniPFRelIso_all)")
    .Define("lepton_eta", "assignleps[1]")
    .Define("lepton_phi", "assignleps[2]");
  
  // --------------------------------------------------------
  // 		      JET Cleaning and JERC
  // --------------------------------------------------------
  
  auto Jet4vecs = LepAssign.Define("FatJet_P4", "fVectorConstructor(FatJet_pt,FatJet_eta,FatJet_phi,FatJet_mass)")
    .Define("DummyZero","float(0.0)");

  auto CleanJets = Jet4vecs.Define("GenJetAK8_P4", "fVectorConstructor(GenJetAK8_pt,GenJetAK8_eta,GenJetAK8_phi,GenJetAK8_mass)")
    .Define("cleanFatJets", cleanJets, {"FatJet_P4","FatJet_rawFactor","FatJet_rawFactor","FatJet_area","FatJet_area","FatJet_jetId","GenJetAK8_P4","FatJet_genJetAK8Idx","SMuon_P4","SMuon_jetIdx","SElectron_P4","SElectron_jetIdx","fixedGridRhoFastjetAll","DummyZero","DummyZero"}); // args 2 and 4 are dummies
  
  auto JetAssign = CleanJets.Define("cleanFatJet_pt", "cleanFatJets[0]")
    .Define("cleanFatJet_eta", "cleanFatJets[1]")
    .Define("cleanFatJet_phi", "cleanFatJets[2]")
    .Define("cleanFatJet_mass", "cleanFatJets[3]")
    .Define("cleanFatJet_rawFactor", "cleanFatJets[4]");

  // ---------------------------------------------------------
  // 	  HT Calculation and N Jets cuts
  // ---------------------------------------------------------
  
  auto JetSelect = JetAssign.Define("DR_lepFatJets","DeltaR_VecAndFloat(cleanFatJet_eta,cleanFatJet_phi,lepton_eta,lepton_phi)")
    .Define("goodcleanFatJets", "cleanFatJet_pt > 200 && abs(cleanFatJet_eta) < 2.5 && FatJet_jetId > 1 && (DR_lepFatJets > 0.8)") 
    .Define("NFatJets", "(int) Sum(goodcleanFatJets)")
    .Define("NOS_gcFatJets","(int) Sum(DR_lepFatJets[goodcleanFatJets == true] > TMath::Pi()/2)")
    .Filter("NFatJets > 0","Pass N good central AK8 > 0")
    .Filter("NOS_gcFatJets > 0","Pass N good central other side AK8 > 0");

  auto FatJetVars = JetSelect.Define("gcFatJet_pt_unsort", "FatJet_pt[goodcleanFatJets == true]")
    .Define("gcFatJet_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcFatJet_pt_unsort))")
    .Define("gcFatJet_pt","reorder(gcFatJet_pt_unsort,gcFatJet_ptargsort)")
    .Define("gcFatJet_eta", "reorder(FatJet_eta[goodcleanFatJets == true],gcFatJet_ptargsort)")
    .Define("gcFatJet_phi", "reorder(FatJet_phi[goodcleanFatJets == true],gcFatJet_ptargsort)")
    .Define("DR_gcFatJets", "reorder(DR_lepFatJets[goodcleanFatJets == true],gcFatJet_ptargsort)")
    .Define("OS_gcFatJets","DR_gcFatJets > TMath::Pi()/2")
    .Define("gcOSFatJet_pt","gcFatJet_pt[OS_gcFatJets == true]")
    .Define("gcOSFatJet_eta","gcFatJet_eta[OS_gcFatJets == true]")
    .Define("gcOSFatJet_phi","gcFatJet_phi[OS_gcFatJets == true]");
  
  // ---------------------------------------------------------
  // 		Add scale factors and MC jet-based calcs
  // ---------------------------------------------------------
  auto scaleFactors = FatJetVars.Define("gcFatJet_subJetIdx1","reorder(FatJet_subJetIdx1[goodcleanFatJets == true],gcFatJet_ptargsort)")
    .Define("gcFatJet_genmatch", Form("FatJet_matching_bkg(\"%s\", gcFatJet_eta, gcFatJet_phi, NFatJets, gcFatJet_subJetIdx1, nSubJet, SubJet_hadronFlavour, nGenPart, GenPart_pdgId, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, t_bkg_idx, W_bkg_idx)",sample.c_str())) // pt-ordered
    .Define("gcOSFatJet_genmatch", "gcFatJet_genmatch[OS_gcFatJets == true]");

  // ---------------------------------------------------------
  // 		JET Tagging variables
  // ---------------------------------------------------------

  auto Taggers = scaleFactors.Define("gcFatJet_pNetJ", "reorder(FatJet_particleNet_QCD[goodcleanFatJets == true],gcFatJet_ptargsort)")
    .Define("gcFatJet_pNetTvsQCD", "reorder(FatJet_particleNet_TvsQCD[goodcleanFatJets == true],gcFatJet_ptargsort)")
    .Define("gcFatJet_pNetWvsQCD", "reorder(FatJet_particleNet_WvsQCD[goodcleanFatJets == true],gcFatJet_ptargsort)")
    .Define("gcFatJet_pNetT", "(gcFatJet_pNetTvsQCD * gcFatJet_pNetJ) / (1 - gcFatJet_pNetTvsQCD)")
    .Define("gcFatJet_pNetW", "(gcFatJet_pNetWvsQCD * gcFatJet_pNetJ) / (1 - gcFatJet_pNetWvsQCD)")
    .Define("gcFatJet_pNetTag", "maxFxn(gcFatJet_pNetJ,gcFatJet_pNetT,gcFatJet_pNetW)")
    .Define("gcFatJet_pNetTag_alt", pnetWPs, {"gcFatJet_pNetTvsQCD", "gcFatJet_pNetWvsQCD"})
    .Define("gcOSFatJet_pNetTag", "gcFatJet_pNetTag[OS_gcFatJets==true]")
    .Define("gcOSFatJet_pNetTag_alt", "gcFatJet_pNetTag_alt[OS_gcFatJets==true]");

    
  // ---------------------------------------------------------
  // 	  Efficiency analysis
  //      Reminder: from genmatch, -6 = leptonic top and -24 = leptonic W!
  // ---------------------------------------------------------
  
  auto Efficiency = Taggers.Define("t_match","gcOSFatJet_pt[gcOSFatJet_genmatch == 6]")
    .Define("W_match","gcOSFatJet_pt[gcOSFatJet_genmatch == 24]")
    .Define("J_match","gcOSFatJet_pt[abs(gcOSFatJet_genmatch) < 6]")
    .Define("t_match_t_tag","gcOSFatJet_pt[gcOSFatJet_genmatch == 6 && gcOSFatJet_pNetTag_alt == 1]")
    .Define("t_match_W_tag","gcOSFatJet_pt[gcOSFatJet_genmatch == 6 && gcOSFatJet_pNetTag_alt == 2]")
    .Define("t_match_J_tag","gcOSFatJet_pt[gcOSFatJet_genmatch == 6 && gcOSFatJet_pNetTag_alt == 0]")
    .Define("W_match_t_tag","gcOSFatJet_pt[gcOSFatJet_genmatch == 24 && gcOSFatJet_pNetTag_alt == 1]")
    .Define("W_match_W_tag","gcOSFatJet_pt[gcOSFatJet_genmatch == 24 && gcOSFatJet_pNetTag_alt == 2]")
    .Define("W_match_J_tag","gcOSFatJet_pt[gcOSFatJet_genmatch == 24 && gcOSFatJet_pNetTag_alt == 0]")
    .Define("J_match_t_tag","gcOSFatJet_pt[abs(gcOSFatJet_genmatch) < 6 && gcOSFatJet_pNetTag_alt == 1]")
    .Define("J_match_J_tag","gcOSFatJet_pt[abs(gcOSFatJet_genmatch) < 6 && gcOSFatJet_pNetTag_alt == 0]")
    .Define("weight","genWeight/abs(genWeight)");

  double ptbins[10] = {200,300,400,500,600,800,1000,1200,1500,2000};

  auto h01 = Efficiency.Histo1D({"PNEff_Dptbins_t",";Jet p_T (GeV);",9,ptbins},"t_match","weight");
  auto h02 = Efficiency.Histo1D({"PNEff_Dptbins_W",";Jet p_T (GeV);",9,ptbins},"W_match","weight");
  auto h03 = Efficiency.Histo1D({"PNEff_Dptbins_J",";Jet p_T (GeV);",9,ptbins},"J_match","weight");
  auto h04 = Efficiency.Histo1D({"PNEff_Nptbins_t_t",";Jet p_T (GeV);t efficiency",9,ptbins},"t_match_t_tag","weight");
  auto h05 = Efficiency.Histo1D({"PNEff_Nptbins_t_W",";Jet p_T (GeV);t-as-W rate",9,ptbins},"t_match_W_tag","weight");
  auto h06 = Efficiency.Histo1D({"PNEff_Nptbins_t_J",";Jet p_T (GeV);t-as-QCD rate",9,ptbins},"t_match_J_tag","weight");
  auto h07 = Efficiency.Histo1D({"PNEff_Nptbins_W_t",";Jet p_T (GeV);W-as-t rate",9,ptbins},"W_match_t_tag","weight");
  auto h08 = Efficiency.Histo1D({"PNEff_Nptbins_W_W",";Jet p_T (GeV);W efficiency",9,ptbins},"W_match_W_tag","weight");
  auto h09 = Efficiency.Histo1D({"PNEff_Nptbins_W_J",";Jet p_T (GeV);W-as-QCD rate",9,ptbins},"W_match_J_tag","weight");
  auto h10 = Efficiency.Histo1D({"PNEff_Nptbins_J_t",";Jet p_T (GeV);QCD-as-t rate",9,ptbins},"J_match_t_tag","weight");
  auto h11 = Efficiency.Histo1D({"PNEff_Nptbins_J_J",";Jet p_T (GeV);QCD efficiency",9,ptbins},"J_match_J_tag","weight");
  
  // -------------------------------------------------
  // 		Save histos to file
  // -------------------------------------------------
  
  cout << "-------------------------------------------------" << endl
       << ">>> Saving " << sample << " Snapshot..." << endl;
  TString finalFile = "RDF_" + sample + "_" + year + "_PNetEff.root";
  const char *stdfinalFile = finalFile;

  TFile::Open(finalFile,"recreate");
  h01->Write();
  h02->Write();
  h03->Write();
  h04->Write();
  h05->Write();
  h06->Write();
  h07->Write();
  h08->Write();
  h09->Write();
  h10->Write();
  h11->Write();

  time.Stop();
  time.Print();
  
  cout << "Cut statistics:" << endl;
  Efficiency.Report()->Print();

  cout << "Done!" << endl;
}
