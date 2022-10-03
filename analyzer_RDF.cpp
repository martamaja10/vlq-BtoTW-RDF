// --------------------------------------------------------------------------------------- //
// Implimentation of RDataFrame in C++.					                   //
// Comments on creating a singly produced VLQ search			                   //
// To Run on Command Line:   root -l callRDF.C\(\"Muon(OR)Electron\",\"testNumber\"\)      //
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
//#include "lwtnn/lwtnn/interface/parse_json.hh"
//#include "lwtnn/lwtnn/interface/LightweightNeuralNetwork.hh"

using namespace ROOT::VecOps;
void rdf::analyzer_RDF(std::string sample, TString chan, TString testNum, int year)
{
	ROOT::EnableImplicitMT();
//	const std::string samplesBasePath = "root://cmsxrootd.fnal.gov/";
	TStopwatch time;
	time.Start();
// Leptonic decay (1Lep)
	TString nChan = "n"+chan;
	const char* stdNChan = nChan;
	std::cout << "Channel Type: " << chan << std::endl;
	bool isNominal = isNominal;
// Regions will be signal and TTTosemileptonic which is background
	TString region = "Signal";
	if(isSig == true){region = "Signal";} // BPrimeBPrime
	else if(isTT == true){region = "TTToSemiLeptonic";} // TTToSemiLeptonic
	std::cout<< "Region: " << region << std::endl;
	// -------------------------------------------------------
	// 			DNN Stuff
	// -------------------------------------------------------
// Need new json files
//	std::string dnnFileBB = "vlq_mlp_Sept09_BB_3arch.json";
	// std::string dnnFileBB = "vlq_mlp_June_08_20_BB_3arch.json";
	// std::ifstream input_cfgBB( dnnFileBB );
	// lwt::JSONConfig cfgBB = lwt::parse_json(input_cfgBB);
	// lwt::LightweightNeuralNetwork* lwtnnBB = new lwt::LightweightNeuralNetwork(cfgBB.inputs, cfgBB.layers, cfgBB.outputs);

	// --------------------------------------------------------------------------------------------------------------------
	// 							LAMBDA FXNS
	// --------------------------------------------------------------------------------------------------------------------
	
	// ----------------------------------------------------
	//   		DEACY/GENTTBAR CALCULATOR:
	// ----------------------------------------------------
// Single decay B-->tW
	auto decayModeSelection_genTTbarMassCalc = [region](unsigned int nGenPart, ROOT::VecOps::RVec<int>& GenPart_pdgId, ROOT::VecOps::RVec<float>& GenPart_mass, ROOT::VecOps::RVec<float>& GenPart_pt, ROOT::VecOps::RVec<float>& GenPart_phi, ROOT::VecOps::RVec<float>& GenPart_eta, ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother, ROOT::VecOps::RVec<int>& GenPart_status)
	{
		int returnVar = 0;
		if(region == "TTToSemiLeptonic")
		{
			int genTTbarMass = -999;
			double topPtWeight = 1.0;
			TLorentzVector top, antitop;
			bool gottop = false;
			bool gotantitop = false;
			bool gottoppt = false;
			bool gotantitoppt = false;
			float toppt, antitoppt;
			for(unsigned int p = 0; p < nGenPart; p++)
			{
				int id = GenPart_pdgId[p];
				if (abs(id) != 6){continue;}
				if (GenPart_mass[p] < 10){continue;}
				int motherid = GenPart_pdgId[GenPart_genPartIdxMother[p]];
				if(abs(motherid) != 6)
				{
					if (!gottop && id == 6)
					{
						top.SetPtEtaPhiM(GenPart_pt[p], GenPart_eta[p], GenPart_phi[p], GenPart_mass[p]);
						gottop = true;
					}
					if (!gotantitop && id == -6)
					{
						antitop.SetPtEtaPhiM(GenPart_pt[p], GenPart_eta[p], GenPart_phi[p], GenPart_mass[p]);
						gotantitop = true;
					}
				}
				if(GenPart_status[p] == 62)
				{
					if (!gottoppt && id == 6)
					{
						toppt = GenPart_pt[p];
						gottoppt = true;
					}
					if (!gotantitoppt && id == -6)
					{
						antitoppt = GenPart_pt[p];
						gotantitoppt = true;
					}
				}
			}
			if(gottop && gotantitop){genTTbarMass = (top+antitop).M();}
			if(gottoppt && gotantitoppt)
			{
				float SFtop = TMath::Exp(0.0615-0.0005*toppt);
				float SFantitop = TMath::Exp(0.0615-0.0005*antitoppt);
				topPtWeight = TMath::Sqrt(SFtop*SFantitop);
			}
			returnVar = genTTbarMass;
		}
		return returnVar;
	};

	// ----------------------------------------------------
	//     minMleppJet VECTOR RETURN + NJETSDEEPFLAV
	// ----------------------------------------------------
// Minimum mass of lepton plus jet
	auto minMleppJet_calc = [isNominal](ROOT::VecOps::RVec<float>& jet_pt, ROOT::VecOps::RVec<float>& jet_eta, ROOT::VecOps::RVec<float>& jet_phi, ROOT::VecOps::RVec<float>& jet_mass, TLorentzVector lepton_lv, ROOT::VecOps::RVec<float> Jet_btagDeepFlavB_GCJ)
	{
		float ind_MinMlj = -1; // This gets changed into int in .Define()
		float minMleppJet = 1e8;
		ROOT::VecOps::RVec<int> theJetBTag_JetSubCalc_PtOrdered (jet_pt.size(),0);
		float NJetsDeepFlavwithSF_JetSubCalc = 0;
		TLorentzVector jet_lv;
		
		for(unsigned int ijet=0; ijet < jet_pt.size(); ijet++)
		{
			jet_lv.SetPtEtaPhiM(jet_pt.at(ijet),jet_eta.at(ijet),jet_phi.at(ijet),jet_mass.at(ijet));
			if(Jet_btagDeepFlavB_GCJ[ijet] > 0.3033){theJetBTag_JetSubCalc_PtOrdered.at(ijet) = 1;} // BTagged or not
			else if(Jet_btagDeepFlavB_GCJ[ijet] < 0.3033){theJetBTag_JetSubCalc_PtOrdered.at(ijet) = 0;}
			
			if((lepton_lv + jet_lv).M() < minMleppJet)
			{
				minMleppJet = fabs((lepton_lv + jet_lv).M());
				ind_MinMlj = ijet;
			}
			if(isNominal && theJetBTag_JetSubCalc_PtOrdered.at(ijet) == 1){NJetsDeepFlavwithSF_JetSubCalc += 1;}
		}
		ROOT::VecOps::RVec<float> minMlj = {minMleppJet,ind_MinMlj,NJetsDeepFlavwithSF_JetSubCalc};
		return minMlj;
	};

	// -----------------------------------------------
	//   LWTNN IMPLIMENTATION AND MYMAP CALCULATION
	// -----------------------------------------------
// jec_ TBD
// Slim down for only B version 

 	auto predictMLP = [lwtnnBB](float corr_met_MultiLepCalc, float AK4HT, int NJets_JetSubCalc, int NJetsAK8_JetSubCalc, float AK4HTpMETpLepPt, float jetPt_1, float jetPt_2, float jetPt_3, float sdMass_1, float sdMass_2, float sdMass_3, float dnnJ_1, float dnnJ_2, float dnnJ_3, float dnnT_1, float dnnT_2, float dnnT_3, float dnnH_1, float dnnH_2, float dnnH_3, float dnnZ_1, float dnnZ_2, float dnnZ_3, float dnnW_1, float dnnW_2, float dnnW_3, float dnnB_1, float dnnB_2, float dnnB_3, int dnnLargest_1, int dnnLargest_2, int dnnLargest_3, int nJ_DeepAK8, int nT_DeepAK8, int nH_DeepAK8, int nZ_DeepAK8, int nW_DeepAK8, int nB_DeepAK8, float tau21_1, float tau21_2, float tau21_3, float minDR_leadAK8otherAK8)
	{
	  ROOT::VecOps::RVec<float> dnn_SigWjetTtbar (6,0);
	  std::map<std::string,double> varMapBB;
	  std::map<std::string,double> myMapBB;
		
	  myMapBB = {
	    {"WjetsBB",-999},
	    {"ttbarBB",-999},
	    {"Bprime",-999},
	  };
	  
	  varMapBB = {
	    {"corr_met_MultiLepCalc", corr_met_MultiLepCalc},
	    {"AK4HTpMETpLepPt", AK4HTpMETpLepPt},
	    {"AK4HT", AK4HT},
	    {"NJets_JetSubCalc", NJets_JetSubCalc},
	    {"NJetsAK8_JetSubCalc", NJetsAK8_JetSubCalc},
	    {"minDR_leadAK8otherAK8", minDR_leadAK8otherAK8},
	    {"nH_DeepAK8", nH_DeepAK8},
	    {"nT_DeepAK8", nT_DeepAK8},
	    {"jetPt_1", jetPt_1},
	    {"jetPt_2", jetPt_2},
	    {"jetPt_3", jetPt_3},
	    {"sdMass_1", sdMass_1},
	    {"sdMass_2", sdMass_2},
	    {"sdMass_3", sdMass_3},
	    {"dnnLargest_1", dnnLargest_1},
	    {"dnnLargest_2", dnnLargest_2},
	    {"dnnLargest_3", dnnLargest_3},
	    {"dnnJ_1", dnnJ_1},
	    {"dnnJ_2", dnnJ_2},
	    {"dnnJ_3", dnnJ_3},
	    {"dnnH_2", dnnH_2},
	    {"dnnH_3", dnnH_3},
	    {"dnnT_1", dnnT_1},
	  };
	  
	  myMapBB = lwtnnBB->compute(varMapBB);
	  dnn_SigWjetTtbar[0] = myMapBB["Bprime"];
	  dnn_SigWjetTtbar[1] = myMapBB["WjetsBB"];
	  dnn_SigWjetTtbar[2] = myMapBB["ttbarBB"];
	  
	  return dnn_SigWjetTtbar;
	};

	// -------------------------------------------------------
	//               Flags and First Filter 
	// -------------------------------------------------------
// Twiki with reccommended ultralegacy values
	auto rdf = ROOT::RDataFrame("Events",sample); // Initial data
	std::cout << "Number of Events: " << rdf.Count().GetValue() << std::endl;

	//	auto largeFlags = rdf.Filter("Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalTightHalo2016Filter == 1")

	auto largeFlags = rdf.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1")
			     .Define("corr_met_MultiLepCalc","MET_pt")\
			     .Define("corr_met_phi_MultiLepCalc","MET_phi")\
			     .Filter("corr_met_MultiLepCalc > 50");
	std::cout << "Number of Events post Flags: " << largeFlags.Count().GetValue() << std::endl;
	auto Lep_df0 = largeFlags.Define("nLeptons", stdNChan)\
				 .Filter("nLeptons > 0");  // REQUIRED to come before any .Define operators for channels
	std::cout << "Number of Events with Leptons: " << Lep_df0.Count().GetValue() << std::endl;

	// --------------------------------------------------------
	//      Initialize Stuff - careful where you put these
	// --------------------------------------------------------
	auto xLep = Lep_df0;

	// ---------------------------------------------------------
	//                    Lepton Filters
	// ---------------------------------------------------------
	//Compare triggers to gridpacks and cut out unneccessary or bad ones
	if(chan == "Muon")
	{
		auto HLTTriggers = Lep_df0.Filter("HLT_Mu50 == 1 || HLT_Mu15_IsoVVVL_PFHT450 == 1 || \
						   HLT_Mu15_IsoVVVL_PFHT600 == 1 || HLT_Mu50_IsoVVVL_PFHT450 == 1");
		auto Lep_df1 = HLTTriggers.Define("LPassMu","Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_miniIsoId >= 1 && Muon_looseId == true")\
					  .Define("nLPassMu","(int) Sum(LPassMu)")\
					  .Filter("nLPassMu == 1");
		std::cout << "Number of Events with Loose Muons: " << Lep_df1.Count().GetValue() << std::endl;
		auto Lep_df2 = Lep_df1.Define("TPassMu","Muon_pt > 50 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_miniIsoId >= 3")\
				      .Define("nTPassMu","(int) Sum(TPassMu)")\
				      .Filter("nTPassMu == 1");
		std::cout << "Number of Events with Tight Muons: " << Lep_df2.Count().GetValue() << std::endl;
		auto Lep_dfEl = Lep_df2.Define("passEl","Electron_pt > 10 && Electron_miniPFRelIso_all < 0.4 && \
							 Electron_mvaFall17V2noIso_WPL == true && abs(Electron_eta) < 2.5")\
				       .Define("nPassEl","(int) Sum(passEl)")\
				       .Filter("nElectron == 0 || (nElectron > 0 && nPassEl == 0)","Good Muon events if no electrons existed or if ");
		std::cout << "Good Muon Events with or without the existance of Electrons: " << Lep_dfEl.Count().GetValue() << std::endl;
		xLep = Lep_dfEl.Define("leptonPt_MultiLepCalc","Muon_pt[TPassMu == true]")\
			       .Define("leptonEta_MultiLepCalc","Muon_eta[TPassMu == true]")\
			       .Define("leptonPhi_MultiLepCalc","Muon_phi[TPassMu == true]")\
			       .Define("leptonMass_MultiLepCalc","Muon_mass[TPassMu == true]")\
			       .Define("leptonMiniIso_MultiLepCalc","Muon_miniPFRelIso_all[TPassMu == true]");
	}
	else if(chan == "Electron")
	{
		auto HLTTriggers = Lep_df0.Filter("HLT_Ele38_WPTight_Gsf == 1 || HLT_Ele35_WPTight_Gsf == 1 || \
						   HLT_Ele15_IsoVVVL_PFHT450 == 1 || HLT_Ele50_IsoVVVL_PFHT450 == 1");
		auto Lep_df1 = HLTTriggers.Define("LPassEl","Electron_pt > 10 && Electron_miniPFRelIso_all < 0.4 && \
							     Electron_mvaFall17V2noIso_WPL == true && abs(Electron_eta) < 2.5")\
					  .Define("nLPassEl","(int) Sum(LPassEl)")\
					  .Filter("nLPassEl == 1");
		std::cout << "Number of Events with Loose Electrons: " << Lep_df1.Count().GetValue() << std::endl;
		auto Lep_df2 = Lep_df1.Define("TPassEl","Electron_pt > 30 && Electron_mvaFall17V2noIso_WP90 == true && \
							 Electron_miniPFRelIso_all < 0.1 && abs(Electron_eta) < 2.5")\
				      .Define("nTPassEl","(int) Sum(TPassEl)")\
				      .Filter("nTPassEl == 1");
		std::cout << "Number of Events with Tight Electrons: " << Lep_df2.Count().GetValue() << std::endl;
		auto Lep_dfMu = Lep_df2.Define("passMu","Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_miniIsoId >= 1 && Muon_looseId == true")\
				       .Define("nPassMu","(int) Sum(passMu)")\
				       .Filter("nMuon == 0 || (nMuon > 0 && nPassMu == 0)","Good Electron events if no Muons existed");
		std::cout << "Good Electron Events with or without the existance of Muons: " << Lep_dfMu.Count().GetValue() << std::endl;
		xLep = Lep_dfMu.Define("leptonPt_MultiLepCalc","Electron_pt[TPassEl == true]")\
			       .Define("leptonEta_MultiLepCalc","Electron_eta[TPassEl == true]")\
			       .Define("leptonPhi_MultiLepCalc","Electron_phi[TPassEl == true]")\
			       .Define("leptonMass_MultiLepCalc","Electron_mass[TPassEl == true]")\
			       .Define("leptonMiniIso_MultiLepCalc","Electron_miniPFRelIso_all[TPassEl == true]");
	}
	// --------------------------------------------------------
	// 		     AK4 JETS w/ cleaning
	// --------------------------------------------------------
	// Supposedly, this is where the particleNet variables will go
	auto jet_ft0 = xLep.Filter("nJet > 0 && nFatJet > 0");
	std::cout << "Number of Events with at least one AK4 and AK8 Jet: " << jet_ft0.Count().GetValue() << std::endl;

	auto jet_df0 = jet_ft0.Define("goodAK4Jets","Jet_pt > 30 && abs(Jet_eta) < 2.4 && Jet_jetId > 1")\
			      .Define("dR_LIM_AK4","(float) 0.4")\
			      .Define("GCAK4Jets","cleanJets(Jet_pt,Jet_mass,goodAK4Jets,Jet_eta,Jet_phi,\
			      					      leptonPt_MultiLepCalc,leptonMass_MultiLepCalc,\
			      					      leptonEta_MultiLepCalc,leptonPhi_MultiLepCalc,dR_LIM_AK4)")\
                             .Define("NJets_JetSubCalc","(int) Sum(GCAK4Jets)")\
			      .Define("theJetPt_JetSubCalc_PtOrdered","Jet_pt[GCAK4Jets == true]")\
			      .Define("theJetEta_JetSubCalc_PtOrdered","Jet_eta[GCAK4Jets == true]")\
			      .Define("theJetPhi_JetSubCalc_PtOrdered","Jet_phi[GCAK4Jets == true]")\
			      .Define("theJetMass_JetSubCalc_PtOrdered","Jet_mass[GCAK4Jets == true]")\
			      .Define("Jet_btagDeepFlavB_GCJ","Jet_btagDeepFlavB[GCAK4Jets == true]")\
			      .Define("goodAK8Jets","FatJet_jetId > 1 && abs(FatJet_eta) < 2.4 && FatJet_pt > 200")\
			      .Define("dR_LIM_AK8","(float) 0.8")\
			      .Define("GCAK8Jets","cleanJets(FatJet_pt,FatJet_mass,goodAK8Jets,FatJet_eta,FatJet_phi,\
			      				     leptonPt_MultiLepCalc,leptonMass_MultiLepCalc,leptonEta_MultiLepCalc,\
			      				     leptonPhi_MultiLepCalc,dR_LIM_AK8)")\
			      .Define("NJetsAK8_JetSubCalc","(int) Sum(GCAK8Jets)")\
			      .Define("theJetAK8Pt_JetSubCalc_PtOrdered","FatJet_pt[GCAK8Jets == true]")\
			      .Define("theJetAK8Eta_JetSubCalc_PtOrdered","FatJet_eta[GCAK8Jets == true]")\
			      .Define("theJetAK8Phi_JetSubCalc_PtOrdered","FatJet_phi[GCAK8Jets == true]")\
			      .Define("theJetAK8Mass_JetSubCalc_PtOrdered","FatJet_mass[GCAK8Jets == true]")\
			      .Define("theJetAK8SoftDropCorr_PtOrdered","FatJet_msoftdrop[GCAK8Jets == true]");

	// ---------------------------------------------------------
	// 	  HT Calculation and Final Preselection Cut
	// ---------------------------------------------------------
	auto HT_calc = jet_df0.Define("AK4HT","Sum(Jet_pt[GCAK4Jets == true])")\
			      .Filter("AK4HT > 510");
	std::cout << "Number of Events passing Preselection (HT Cut): " << HT_calc.Count().GetValue() << std::endl;

	// ---------------------------------------------------------
	//    Uncomment to save seperate Preselection .root file
	// ---------------------------------------------------------
// 	TString outputFilePS = "RDF_"+region+"_FullPreselectionHFIles_"+chan+"_PS.root";
// 	const char* stdOutputFilePS = outputFilePS;
// 	std::cout << "------------------------------------------------" << std::endl << ">>> Saving Preselection Snapshot..." << std::endl;
// 	HT_calc.Snapshot("Events", stdOutputFilePS);
// 	std::cout << "Output File: " << outputFilePS << std::endl << "-------------------------------------------------" << std::endl;
// }
	//----------------------------------------------------------
	//       Uncomment if using a preselection file
	//----------------------------------------------------------
	//	auto HT_calc = rdf;	


	// ---------------------------------------------------------
	// 		Post Preselection Analysis
	// ---------------------------------------------------------
//      Look for particle net instead of deepAK8 (look at the large chart with variable names)
	auto postPresel = HT_calc.Define("decayMode_or_genTTbarMass",decayModeSelection_genTTbarMassCalc,{"nGenPart","GenPart_pdgId","GenPart_mass", \
	      "GenPart_pt","GenPart_phi","GenPart_eta",			\
	      "GenPart_genPartIdxMother","GenPart_status"})		\
	  .Filter("NJetsAK8_JetSubCalc > 0")				\
	  .Define("lepton_lv","fVectorConstructor(leptonPt_MultiLepCalc,leptonEta_MultiLepCalc,\
				      					 leptonPhi_MultiLepCalc,leptonMass_MultiLepCalc)")\
	  .Define("AK4jet_lv","fVectorConstructor(theJetPt_JetSubCalc_PtOrdered,theJetEta_JetSubCalc_PtOrdered,\
				      					 theJetPhi_JetSubCalc_PtOrdered,theJetMass_JetSubCalc_PtOrdered)")\
	  .Define("AK8jet_lv","fVectorConstructor(theJetAK8Pt_JetSubCalc_PtOrdered,theJetAK8Eta_JetSubCalc_PtOrdered,\
				      					 theJetAK8Phi_JetSubCalc_PtOrdered,theJetAK8Mass_JetSubCalc_PtOrdered)")\
	  .Define("AK4HTpMETpLepPt","AK4HT + leptonPt_MultiLepCalc[0] + corr_met_MultiLepCalc")	\
	  .Define("jetPt_1","theJetAK8Pt_JetSubCalc_PtOrdered[0]")	\
	  .Define("jetPt_2","theJetAK8Pt_JetSubCalc_PtOrdered[1]")	\
	  .Define("jetPt_3","theJetAK8Pt_JetSubCalc_PtOrdered[2]")	\
	  .Define("sdMass","FatJet_msoftdrop[GCAK8Jets == true]")	\
	  .Define("sdMass_1","sdMass[0]")				\
	  .Define("sdMass_2","sdMass[1]")				\
	  .Define("sdMass_3","sdMass[2]")				\
	  .Define("dnnJ","FatJet_deepTag_QCDothers[GCAK8Jets == true]")	\
	  .Define("dnnJ_1","dnnJ[0]")					\
	  .Define("dnnJ_2","dnnJ[1]")					\
	  .Define("dnnJ_3","dnnJ[2]")					\
	  .Define("int_dnnT","(FatJet_deepTag_TvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_TvsQCD)") \
	  .Define("dnnT","int_dnnT[GCAK8Jets == true]")			\
	  .Define("dnnT_1","dnnT[0]")					\
	  .Define("dnnT_2","dnnT[1]")					\
	  .Define("dnnT_3","dnnT[2]")					\
	  .Define("dnnH","FatJet_deepTag_H[GCAK8Jets == true]")		\
	  .Define("dnnH_1","dnnH[0]")					\
	  .Define("dnnH_2","dnnH[1]")					\
	  .Define("dnnH_3","dnnH[2]")					\
	  .Define("int_dnnZ","(FatJet_deepTag_ZvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_ZvsQCD)") \
	  .Define("dnnZ","int_dnnZ[GCAK8Jets == true]")			\
	  .Define("dnnZ_1","dnnZ[0]")					\
	  .Define("dnnZ_2","dnnZ[1]")					\
	  .Define("dnnZ_3","dnnZ[2]")					\
	  .Define("int_dnnW","(FatJet_deepTag_WvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_WvsQCD)") \
	  .Define("dnnW","int_dnnW[GCAK8Jets == true]")			\
	  .Define("dnnW_1","dnnW[0]")					\
	  .Define("dnnW_2","dnnW[1]")					\
	  .Define("dnnW_3","dnnW[2]")					\
	  .Define("int_dnnB","(FatJet_deepTag_QCD - FatJet_deepTag_QCDothers)")	\
	  .Define("dnnB","int_dnnB[GCAK8Jets == true]")			\
	  .Define("dnnB_1","dnnB[0]")					\
	  .Define("dnnB_2","dnnB[1]")					\
	  .Define("dnnB_3","dnnB[2]")					\
	  .Define("dnnLargest","maxFxn(dnnJ,dnnT,dnnH,dnnZ,dnnW,dnnB)")	\
	  .Define("dnnLargest_1","dnnLargest[0]")			\
	  .Define("dnnLargest_2","dnnLargest[1]")			\
	  .Define("dnnLargest_3","dnnLargest[2]")			\
	  .Define("nJ_DeepAK8","Sum(dnnLargest == 0)")			\
	  .Define("nT_DeepAK8","Sum(dnnLargest == 1)")			\
	  .Define("nH_DeepAK8","Sum(dnnLargest == 2)")			\
	  .Define("nZ_DeepAK8","Sum(dnnLargest == 3)")			\
	  .Define("nW_DeepAK8","Sum(dnnLargest == 4)")			\
	  .Define("nB_DeepAK8","Sum(dnnLargest == 5)")			\
	  .Define("int_tau21","(FatJet_tau2 / FatJet_tau1)")		\
	  .Define("tau21","int_tau21[GCAK8Jets == true]")		\
	  .Define("tau21_1","tau21[0]")					\
	  .Define("tau21_2","tau21[1]")					\
	  .Define("tau21_3","tau21[2]")					\
	  .Define("minDR_ptRel_lead_lepAK8","minDR_ptRel_lead_calc(theJetAK8Pt_JetSubCalc_PtOrdered,\
											  theJetAK8Eta_JetSubCalc_PtOrdered, \
				      							  theJetAK8Phi_JetSubCalc_PtOrdered, \
											  theJetAK8Mass_JetSubCalc_PtOrdered, \
				      							  lepton_lv)")\
	  .Define("minDR_lepAK8","minDR_ptRel_lead_lepAK8[0]")		\
	  .Define("ptRel_lepAK8","minDR_ptRel_lead_lepAK8[1]")		\
	  .Define("minDR_leadAK8otherAK8","minDR_ptRel_lead_lepAK8[2]")	\
	  .Define("minDR_ptRel_lead_lepJets","minDR_ptRel_lead_calc(theJetPt_JetSubCalc_PtOrdered,\
				      							   theJetEta_JetSubCalc_PtOrdered, \
				      							   theJetPhi_JetSubCalc_PtOrdered, \
				      							   theJetMass_JetSubCalc_PtOrdered, \
				      							   lepton_lv)")\
	  .Define("minDR_lepJet","minDR_ptRel_lead_lepJets[0]")		\
	  .Define("ptRel_lepJet","minDR_ptRel_lead_lepJets[1]")		\
	  .Define("DR_lepAK8s","DR_calc(theJetAK8Pt_JetSubCalc_PtOrdered,theJetAK8Eta_JetSubCalc_PtOrdered,\
				      			       theJetAK8Phi_JetSubCalc_PtOrdered,theJetAK8Mass_JetSubCalc_PtOrdered, \
				      			       leptonPt_MultiLepCalc,leptonEta_MultiLepCalc, \
				      			       leptonPhi_MultiLepCalc,leptonMass_MultiLepCalc)")\
	  .Define("DR_lepJets","DR_calc(theJetPt_JetSubCalc_PtOrdered,theJetEta_JetSubCalc_PtOrdered,\
				      			       theJetPhi_JetSubCalc_PtOrdered,theJetMass_JetSubCalc_PtOrdered, \
				      			       leptonPt_MultiLepCalc,leptonEta_MultiLepCalc, \
				      			       leptonPhi_MultiLepCalc,leptonMass_MultiLepCalc)")\
	  .Define("W_lv","lpNu_WCalc(corr_met_MultiLepCalc,corr_met_phi_MultiLepCalc,lepton_lv)") \
	  .Define("minMlj",minMleppJet_calc,{"theJetPt_JetSubCalc_PtOrdered","theJetEta_JetSubCalc_PtOrdered", \
		"theJetPhi_JetSubCalc_PtOrdered","theJetMass_JetSubCalc_PtOrdered", \
		"lepton_lv","Jet_btagDeepFlavB_GCJ"})			\
	  .Define("W_dRLep","dR_Wt_Calc(W_lv,lepton_lv)")		\
	  .Define("minMleppJet","minMlj[0]")				\
	  .Define("ind_MinMlj","(int) minMlj[1]")			\
	  .Define("NJetsDeepFlavwithSF_JetSubCalc","(int) minMlj[2]")	\
	  .Define("isLeptonic","isLeptonic_X(minMleppJet)")		\
	  .Define("t_lv","lpNu_t_Calc(isLeptonic,theJetPt_JetSubCalc_PtOrdered,theJetEta_JetSubCalc_PtOrdered,\
				   			     theJetPhi_JetSubCalc_PtOrdered,theJetMass_JetSubCalc_PtOrdered, \
				   			     W_lv,minMleppJet,ind_MinMlj)")\
	  .Define("t_pt","t_lv[0]")					\
	  .Define("t_eta","t_lv[1]")					\
	  .Define("t_phi","t_lv[2]")					\
	  .Define("t_mass","t_lv[3]")					\
	  .Define("t_dRWb","t_lv[4]")					\
	  .Define("top_lv","top_lvConstructor(t_pt,t_eta,t_phi,t_mass)")
	  .Define("tj_vec","three_jet(top_lv,W_lv,isLeptonic,theJetAK8Pt_JetSubCalc_PtOrdered,theJetAK8Eta_JetSubCalc_PtOrdered,\
	  			      			     theJetAK8Phi_JetSubCalc_PtOrdered,theJetAK8Mass_JetSubCalc_PtOrdered,dnnT,dnnH,dnnZ, \
	  			      			     dnnW,dnnB,pnLargest,theJetAK8SoftDropCorr_PtOrdered)")\
	  .Define("Bprime1_DeepAK8_Mass","tj_vec[0]")			\
	  .Define("Bprime1_DeepAK8_Pt","tj_vec[1]")			\
	  .Define("Bprime1_DeepAK8_Eta","tj_vec[2]")			\
	  .Define("Bprime1_DeepAK8_Phi","tj_vec[3]")			\
	  .Define("Bprime1_DeepAK8_deltaR","tj_vec[4]")			\
	  .Define("BPrime1_lv","top_lvConstructor(Bprime1_DeepAK8_Pt,Bprime1_DeepAK8_Eta,Bprime1_DeepAK8_Phi,Bprime1_DeepAK8_Mass)")\
	  .Define("isValidBDecay","tj_vec[5]")\
	  .Define("taggedWbjetJet","tj_vec[6]")\
	  .Define("taggedTjet","tj_vec[7]")\
	  .Define("taggedWjet","tj_vec[8]")\
	  .Define("int_pnT","(FatJet_particleNet_TvsQCD * FatJet_particleNet_QCD) / (1 - FatJet_particleNet_TvsQCD)") \
	  .Define("pnT","int_pnT[GCAK8Jets == true]")			\
	  .Define("pnT_1","pnT[0]")					\
	  .Define("pnT_2","pnT[1]")					\
	  .Define("pnT_3","pnT[2]")					\
	  .Define("pnH","FatJet_particleNet_H4qvsQCD[GCAK8Jets == true]")		\ //I'm not sure about this function
	  .Define("pnH_1","pnH[0]")					\
	  .Define("pnH_2","pnH[1]")					\
	  .Define("pnH_3","pnH[2]")					\
	  .Define("int_pnZ","(FatJet_particleNet_ZvsQCD * FatJet_particleNet_QCD) / (1 - FatJet_particleNet_ZvsQCD)") \
	  .Define("pnZ","int_pnZ[GCAK8Jets == true]")			\
	  .Define("pnZ_1","pnZ[0]")					\
	  .Define("pnZ_2","pnZ[1]")					\
	  .Define("pnZ_3","pnZ[2]")					\
	  .Define("int_pnW","(FatJet_particleNet_WvsQCD * FatJet_particleNet_QCD) / (1 - FatJet_particleNet_WvsQCD)") \
	  .Define("pnW","int_pnW[GCAK8Jets == true]")			\
	  .Define("pnW_1","pnW[0]")					\
	  .Define("pnW_2","pnW[1]")					\
	  .Define("pnW_3","pnW[2]")					\
	  .Define("int_dnnB","(FatJet_deepTag_QCD - FatJet_deepTag_QCDothers)") //Can't find corrolary for QCDothers
	//   .Define("dnnB","int_dnnB[GCAK8Jets == true]")			\
	//   .Define("dnnB_1","dnnB[0]")					\
	//   .Define("dnnB_2","dnnB[1]")					\
	//   .Define("dnnB_3","dnnB[2]")					\
	//   .Define("dnnLargest","maxFxn(dnnJ,dnnT,dnnH,dnnZ,dnnW,dnnB)")	\
	//   .Define("dnnLargest_1","dnnLargest[0]")			\
	//   .Define("dnnLargest_2","dnnLargest[1]")			\
	//   .Define("dnnLargest_3","dnnLargest[2]")
	  ;
	// -------------------------------------------------
	// 		Save Snapshot to file
	// -------------------------------------------------
	// Change final snapshot file under new name so we don't get confused
	std::cout << "Number of Events passing post preselection: " << postPresel.Count().GetValue() << std::endl;
	auto dfFinal = postPresel;
	std::cout << "Number of Events passing post preselection(dfFinal Check): " << dfFinal.Count().GetValue() << std::endl;

	std::cout << "-------------------------------------------------" << std::endl << ">>> Saving " << region << " Snapshot..." << std::endl;
	TString outputFileFNAL = "RDF_"+region+"_seniorResearch_"+chan+"_HnanoAOD"+testNum+".root";
	const char* stdOutputFileFNAL = outputFileFNAL;
	dfFinal.Snapshot("Events", stdOutputFileFNAL);
	std::cout << "Output File: " << outputFileFNAL << std::endl << "-------------------------------------------------" << std::endl;
	time.Stop();
	time.Print();
}
