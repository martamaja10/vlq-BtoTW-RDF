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
#include "../lwtnn/include/lwtnn/parse_json.hh"
#include "../lwtnn/include/lwtnn/LightweightNeuralNetwork.hh"

using namespace ROOT::VecOps;
void rdf::analyzer_RDF(std::string filename, TString testNum, int year)
{
  ROOT::EnableImplicitMT();
  TStopwatch time;
  time.Start();
  // TString nChan = "n"+chan;
  // const char* stdNChan = nChan;
  // std::cout << "Channel Type: " << chan << std::endl;
  bool isNominal = isNominal;
  
  // Samples will be signal, ttbar (a background), and possibly others...
  TString sample = "singletop";
  if(isTT == true){sample = "ttbar";} 
  else if(isSig == true){sample = "Bprime";}
  else if(isMadgraphBkg == true){sample = "wjets";}
  std::cout<< "Sample: " << sample << std::endl;
  
  // -------------------------------------------------------
  // 			DNN Stuff
  // -------------------------------------------------------
  // Need new json files for B -> tW!
  std::ifstream input_mlpN("mlp_new.json");
  std::ifstream input_mlpO("mlp_old.json");
  lwt::JSONConfig mlpN = lwt::parse_json(input_mlpN);
  lwt::JSONConfig mlpO = lwt::parse_json(input_mlpO);
  lwt::LightweightNeuralNetwork* lwtnn_mlpN = new lwt::LightweightNeuralNetwork(mlpN.inputs, mlpN.layers, mlpN.outputs);
  lwt::LightweightNeuralNetwork* lwtnn_mlpO = new lwt::LightweightNeuralNetwork(mlpO.inputs, mlpO.layers, mlpO.outputs);
  
  // --------------------------------------------------------------------------------------------------------------------
  // 							LAMBDA FXNS
  // --------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------           
  //           Bprime truth extraction:                      
  // ---------------------------------------------------- 
  auto Bprime_gen_info=[sample](unsigned int nGenPart, ROOT::VecOps::RVec<int>& GenPart_pdgId, ROOT::VecOps::RVec<float>& GenPart_mass, ROOT::VecOps::RVec<float>& GenPart_pt, ROOT::VecOps::RVec<float>& GenPart_phi, ROOT::VecOps::RVec<float>& GenPart_eta, ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother, ROOT::VecOps::RVec<int>& GenPart_status){
    ROOT::VecOps::RVec<int> BPrimeInfo (6,-999);

    if(sample!="Bprime"){return BPrimeInfo;}

    for(unsigned int p=0; p<nGenPart; p++){
      int id=GenPart_pdgId[p];
      if(abs(id) != 6000007){continue;}
      int motherid = GenPart_pdgId[GenPart_genPartIdxMother[p]];
      if(motherid != id){ // takes the second B'                         
        BPrimeInfo[0] = GenPart_pt[p];
        BPrimeInfo[1] = GenPart_eta[p];
        BPrimeInfo[2] = GenPart_phi[p];
        BPrimeInfo[3] = GenPart_mass[p];
        BPrimeInfo[4] = GenPart_pdgId[p];
        BPrimeInfo[5] = GenPart_status[p];
      }
    }
    return BPrimeInfo;
  };

  // ----------------------------------------------------                        
  //           t truth extraction:    
  // ---------------------------------------------------- 
  auto t_daughter_gen=[sample](unsigned int nGenPart, ROOT::VecOps::RVec<int>& GenPart_pdgId, ROOT::VecOps::RVec<float>& GenPart_mass, ROOT::VecOps::RVec<float>& GenPart_pt, ROOT::VecOps::RVec<float>& GenPart_phi, ROOT::VecOps::RVec<float>& GenPart_eta, ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother, ROOT::VecOps::RVec<int>& GenPart_status){
    ROOT::VecOps::RVec<int> t_daughter_gen_info(21,-999);

    if(sample!="Bprime"){return t_daughter_gen_info;}

    int t_motherIdx = -1, W_motherIdx = -1; // if not t->t' or W->W', no need to store the mother idx info
    int trueLeptonicT = 0; // 0-false 1-true                                                                                     

    for(unsigned int i=0; i<nGenPart; i++){
      int id = GenPart_pdgId[i];
      int motherid = GenPart_pdgId[GenPart_genPartIdxMother[i]];

      if(abs(motherid)!=6000007){continue;}
      if(abs(id)!=6){continue;} // B->tW, pick out t              

      int igen = i;

      for(unsigned int j=i; j<nGenPart; j++){
        if((GenPart_genPartIdxMother[j])!=i){continue;} // pick out the daughter of t      
        if((GenPart_pdgId[j]) != id){continue;}
        t_motherIdx = j; // store idx of t': t->t'    
        igen = j;
      }

      t_daughter_gen_info[0] = GenPart_pt[igen];
      t_daughter_gen_info[1] = GenPart_eta[igen];
      t_daughter_gen_info[2] = GenPart_phi[igen];
      t_daughter_gen_info[3] = GenPart_mass[igen];
      t_daughter_gen_info[4] = GenPart_pdgId[igen];
      t_daughter_gen_info[5] = GenPart_status[igen];
      t_daughter_gen_info[6] = t_motherIdx;

      for(unsigned int j=igen; j<nGenPart; j++){
        int j_id = GenPart_pdgId[j]; // look for W in t->Wb                                                                       
        int j_mother = GenPart_genPartIdxMother[j];
        if(j_mother!=igen){continue;}
        if(abs(j_id)!=24){continue;}

        int jgen = j;
        for(unsigned int k = igen; k<nGenPart; k++){
          int k_id = GenPart_pdgId[k];
          int k_mother = GenPart_genPartIdxMother[k];
          if(k_mother!=j){continue;}
          if(k_id == j_id){
            W_motherIdx = k_id; // store idx of W': W->W'                          
            jgen = k;
          }
        }

	t_daughter_gen_info[7] = GenPart_pt[jgen];
        t_daughter_gen_info[8] = GenPart_eta[jgen];
        t_daughter_gen_info[9] = GenPart_phi[jgen];
        t_daughter_gen_info[10] = GenPart_mass[jgen];
        t_daughter_gen_info[11] = GenPart_pdgId[jgen];
        t_daughter_gen_info[12] = GenPart_status[jgen];
        t_daughter_gen_info[13] = W_motherIdx;

	for(unsigned int p=jgen; p<nGenPart; p++){
          int p_id = GenPart_pdgId[p];
          int p_mother = GenPart_genPartIdxMother[p];
          if(p_mother!=jgen){continue;}
          if(abs(p_id)==11 | abs(p_id)==13){ // why are we not taking tau?                     
            trueLeptonicT = 1;
            t_daughter_gen_info[14] = GenPart_pt[p];
            t_daughter_gen_info[15] = GenPart_eta[p];
            t_daughter_gen_info[16] = GenPart_phi[p];
            t_daughter_gen_info[17] = GenPart_mass[p];
            t_daughter_gen_info[18] = GenPart_pdgId[p];
            t_daughter_gen_info[19] = GenPart_status[p];
            t_daughter_gen_info[20] = trueLeptonicT;
          }
        }
      }
    }
    return t_daughter_gen_info;
  };

  // ----------------------------------------------------           
  //           W truth extraction: 
  // ---------------------------------------------------- 
  auto W_daughter_gen=[sample](unsigned int nGenPart, ROOT::VecOps::RVec<int>& GenPart_pdgId, ROOT::VecOps::RVec<float>& GenPart_mass, ROOT::VecOps::RVec<float>& GenPart_pt, ROOT::VecOps::RVec<float>& GenPart_phi, ROOT::VecOps::RVec<float>& GenPart_eta, ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother, ROOT::VecOps::RVec<int>& GenPart_status){
    ROOT::VecOps::RVec<int> W_daughter_gen_info(14,-999);

    if(sample!="Bprime"){return W_daughter_gen_info;}

    int W_motherIdx = -1; // if not W->W', no need to store the mother idx info                             
    int trueW_decayMode = 0; // 0-did not decay, 1-leptonic, 2-hadronic/tau                                 

    for(unsigned int i=0; i<nGenPart; i++){
      int id = GenPart_pdgId[i];
      int motherid = GenPart_pdgId[GenPart_genPartIdxMother[i]];

      if(abs(motherid)!=6000007){continue;}
      if(abs(id)!=24){continue;} // B->tW, pick out W                                                       

      int igen = i;

      for(unsigned int j=i; j<nGenPart; j++){
        if((GenPart_genPartIdxMother[j])!=i){continue;} // pick out the daughter of W                                           
        if((GenPart_pdgId[j]) != id){continue;}
        W_motherIdx = j; // store idx of W': W->W'                                                                              
        igen = j;
      }

      W_daughter_gen_info[0] = GenPart_pt[igen];
      W_daughter_gen_info[1] = GenPart_eta[igen];
      W_daughter_gen_info[2] = GenPart_phi[igen];
      W_daughter_gen_info[3] = GenPart_mass[igen];
      W_daughter_gen_info[4] = GenPart_pdgId[igen];
      W_daughter_gen_info[5] = GenPart_status[igen];
      W_daughter_gen_info[6] = W_motherIdx;

      for(unsigned int p=igen; p<nGenPart; p++){
        int p_id = GenPart_pdgId[p];
        int p_mother = GenPart_genPartIdxMother[p];
        if(p_mother!=igen){continue;}
        if(abs(p_id)==11 | abs(p_id)==13){ // why are we not taking tau?                                                        
          trueW_decayMode = 1;
          W_daughter_gen_info[7] = GenPart_pt[p];
          W_daughter_gen_info[8] = GenPart_eta[p];
          W_daughter_gen_info[9] = GenPart_phi[p];
          W_daughter_gen_info[10] = GenPart_mass[p];
          W_daughter_gen_info[11] = GenPart_pdgId[p];
          W_daughter_gen_info[12] = GenPart_status[p];
          W_daughter_gen_info[13] = trueW_decayMode;
        }
        else{trueW_decayMode = 2;}
      }
    }
    return W_daughter_gen_info;
  };

  // ----------------------------------------------------
  //   		ttbar background mass CALCULATOR:
  // ----------------------------------------------------
  
  auto genttbarMassCalc = [sample](unsigned int nGenPart, ROOT::VecOps::RVec<int>& GenPart_pdgId, ROOT::VecOps::RVec<float>& GenPart_mass, ROOT::VecOps::RVec<float>& GenPart_pt, ROOT::VecOps::RVec<float>& GenPart_phi, ROOT::VecOps::RVec<float>& GenPart_eta, ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother, ROOT::VecOps::RVec<int>& GenPart_status)
    {
      int returnVar = 0;
      if(sample == "ttbar")
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
  //     minM_lep_jet VECTOR RETURN + NJETSDEEPFLAV
  // ----------------------------------------------------
  
  auto minM_lep_jet_calc = [isNominal](ROOT::VecOps::RVec<float>& jet_pt, ROOT::VecOps::RVec<float>& jet_eta, ROOT::VecOps::RVec<float>& jet_phi, ROOT::VecOps::RVec<float>& jet_mass, TLorentzVector lepton_lv)
    {
      float ind_MinMlj = -1; // This gets changed into int in .Define()
      float minMleppJet = 1e8;
      TLorentzVector jet_lv;
      
      for(unsigned int ijet=0; ijet < jet_pt.size(); ijet++)
	{
	  jet_lv.SetPtEtaPhiM(jet_pt.at(ijet),jet_eta.at(ijet),jet_phi.at(ijet),jet_mass.at(ijet));			
	  if((lepton_lv + jet_lv).M() < minMleppJet)
	    {
	      minMleppJet = fabs((lepton_lv + jet_lv).M());
	      ind_MinMlj = ijet;
	    }
	}
      ROOT::VecOps::RVec<float> minMlj = {minMleppJet,ind_MinMlj};
      return minMlj;
    };
  
  // -----------------------------------------------
  //   LWTNN IMPLIMENTATION AND MYMAP CALCULATION
  // -----------------------------------------------	
  
  auto predictMLP = [lwtnn_mlpN,lwtnn_mlpO](float pNet_J_1, float pNet_T_1, float dpak8_T_1, float FatJet_pt_1, float FatJet_sdMass_1, float tau21_1, int nT_dpak8, int nT_pNet, float Jet_HT, float Jet_ST, float MET_pt, int NJets_DeepFlavM, int NJets_forward)
  {
    ROOT::VecOps::RVec<float> mlp_scores (6,-1);
    std::map<std::string,double> varMapB;
    std::map<std::string,double> myMapB;
  
    varMapB = {
      {"pNet_J_1",	 pNet_J_1},	 
      {"pNet_T_1",	 pNet_T_1},	 
      {"dpak8_T_1",	 dpak8_T_1},	 
      {"FatJet_pt_1",	 FatJet_pt_1},
      {"FatJet_sdMass_1",FatJet_sdMass_1},
      {"tau21_1",	 tau21_1},
      {"nT_dpak8",	 nT_dpak8},
      {"nT_pNet",	 nT_pNet},
      {"Jet_HT",	 Jet_HT},
      {"Jet_ST",	 Jet_ST},
      {"MET_pt",	 MET_pt},
      {"NJets_DeepFlavM",NJets_DeepFlavM},
      {"NJets_forward",  NJets_forward},
    };

    myMapB = {{"WJets",-999},{"TTbar",-999},{"Bprime",-999}};    
    myMapB = lwtnn_mlpN->compute(varMapB);
    mlp_scores[0] = myMapB["WJets"];
    mlp_scores[1] = myMapB["TTbar"];
    mlp_scores[2] = myMapB["Bprime"];

    myMapB = {{"WJets",-999},{"TTbar",-999},{"Bprime",-999}};    
    myMapB = lwtnn_mlpO->compute(varMapB);
    mlp_scores[3] = myMapB["WJets"];
    mlp_scores[4] = myMapB["TTbar"];
    mlp_scores[5] = myMapB["Bprime"];

    return mlp_scores;
  };
  
  // -------------------------------------------------------
  //               Flags and First Filter 
  // -------------------------------------------------------
  // Twiki with reccommended ultralegacy values
  auto rdf = ROOT::RDataFrame("Events",filename); // Initial data
  //  std::cout << "Number of Events: " << rdf.Count().GetValue() << std::endl;
  
  auto METfilters = rdf.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1","MET Filters")
    .Filter("MET_pt > 50","Pass MET > 50");
  //  std::cout << "Number of Events post MET filters: " << METfilters.Count().GetValue() << std::endl;

  // ---------------------------------------------------------
  //                    Lepton Filters
  // ---------------------------------------------------------

  auto Lep_df0 = METfilters.Define("TPassMu","Muon_pt > 50 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_miniIsoId >= 3") \
    .Define("nTPassMu","(int) Sum(TPassMu)")				\
    .Define("TPassEl","Electron_pt > 50 && Electron_mvaFall17V2noIso_WP90 == true && \
							 Electron_miniPFRelIso_all < 0.1 && abs(Electron_eta) < 2.5")\
    .Define("nTPassEl","(int) Sum(TPassEl)")				\
    .Define("isMu","(nMuon > 0 && nTPassMu == 1 && HLT_Mu50 == 1 && (nElectron == 0 || (nElectron > 0 && nTPassEl == 0)))") \
    .Define("isEl","(nElectron > 0 && nTPassEl == 1 && (HLT_Ele38_WPTight_Gsf == 1 || HLT_Ele35_WPTight_Gsf == 1) && (nMuon == 0 || (nMuon > 0 && nTPassMu == 0)))") \
    .Filter("isMu || isEl","Event is either muon or electron");
  
  auto Lep_df1 = Lep_df0.Define("assignleps","assign_leps(isMu,isEl,TPassMu,TPassEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_miniPFRelIso_all,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_miniPFRelIso_all)") \
    .Define("lepton_pt","assignleps[0]")				\
    .Define("lepton_eta","assignleps[1]")				\
    .Define("lepton_phi","assignleps[2]")				\
    .Define("lepton_mass","assignleps[3]")				\
    .Define("lepton_miniIso","assignleps[4]");
  
    
  // --------------------------------------------------------
  // 		      JET SELECTION w/ cleaning
  // --------------------------------------------------------
  
  auto jet_ft0 = Lep_df1.Filter("nJet > 0 && nFatJet > 0","Event has jets");
  //  std::cout << "Number of Events with at least one AK4 and AK8 Jet: " << jet_ft0.Count().GetValue() << std::endl;
  
  auto jet_df0 = jet_ft0.Define("goodJets","Jet_pt > 30 && abs(Jet_eta) < 2.4 && Jet_jetId > 1") \
    .Define("dR_LIM_AK4","(float) 0.4")					\
    .Define("goodcleanJets","cleanJets(Jet_pt,Jet_mass,goodJets,Jet_eta,Jet_phi,\
			      					 lepton_pt,lepton_mass,lepton_eta,lepton_phi,dR_LIM_AK4)")\
    .Define("NJets_central","(int) Sum(goodcleanJets)")			\
    .Define("gcJet_pt","Jet_pt[goodcleanJets == true]")			\
    .Define("gcJet_eta","Jet_eta[goodcleanJets == true]")		\
    .Define("gcJet_phi","Jet_phi[goodcleanJets == true]")		\
    .Define("gcJet_mass","Jet_mass[goodcleanJets == true]")		\
    .Define("gcJet_DeepFlav","Jet_btagDeepFlavB[goodcleanJets == true]") \
    .Define("gcJet_DeepFlavM","gcJet_DeepFlav > 0.2783")		\
    .Define("NJets_DeepFlavM","(int) Sum(gcJet_DeepFlavM)")		\
    .Define("forwardJets","Jet_pt > 30 && abs(Jet_eta) >= 2.4 && Jet_jetId > 1") \
    .Define("goodcleanForwardJets","cleanJets(Jet_pt,Jet_mass,forwardJets,Jet_eta,Jet_phi,\
			      					 lepton_pt,lepton_mass,lepton_eta,lepton_phi,dR_LIM_AK4)")\
    .Define("NJets_forward","(int) Sum(goodcleanForwardJets)")		\
    .Define("gcforwJet_pt","Jet_pt[goodcleanForwardJets == true]")	\
    .Define("gcforwJet_eta","Jet_eta[goodcleanForwardJets == true]")	\
    .Define("gcforwJet_phi","Jet_phi[goodcleanForwardJets == true]")	\
    .Define("gcforwJet_mass","Jet_mass[goodcleanForwardJets == true]")	\
    .Define("gcforwJet_DeepFlav","Jet_btagDeepFlavB[goodcleanForwardJets == true]") \
    .Define("goodFatJets","FatJet_jetId > 1 && abs(FatJet_eta) < 2.4 && FatJet_pt > 200") \
    .Define("dR_LIM_AK8","(float) 0.8")					\
    .Define("goodcleanFatJets","cleanJets(FatJet_pt,FatJet_mass,goodFatJets,FatJet_eta,FatJet_phi,\
			      				            lepton_pt,lepton_mass,lepton_eta,lepton_phi,dR_LIM_AK8)")\
    .Define("NFatJets","(int) Sum(goodcleanFatJets)")			\
    .Define("gcFatJet_pt","FatJet_pt[goodcleanFatJets == true]")	\
    .Define("gcFatJet_eta","FatJet_eta[goodcleanFatJets == true]")	\
    .Define("gcFatJet_phi","FatJet_phi[goodcleanFatJets == true]")	\
    .Define("gcFatJet_mass","FatJet_mass[goodcleanFatJets == true]")	\
    .Define("gcFatJet_msoftdrop","FatJet_msoftdrop[goodcleanFatJets == true]");
  
  // ---------------------------------------------------------
  // 	  HT Calculation and Final Preselection Cut
  // ---------------------------------------------------------
  auto HT_calc = jet_df0.Define("Jet_HT","Sum(Jet_pt[goodcleanJets == true])") \
    .Filter("Jet_HT > 250","Pass HT > 250")						\
    .Filter("NFatJets > 0","Pass N good central AK8 > 0");				        
  //  std::cout << "Number of Events passing Preselection (HT Cut): " << HT_calc.Count().GetValue() << std::endl;
  
  // ---------------------------------------------------------
  //    Uncomment to save seperate Preselection .root file
  // ---------------------------------------------------------
  // TString outputFilePS = "RDF_"+sample+"_presel_"+chan+"_"+testnum+".root";
  // const char* stdOutputFilePS = outputFilePS;
  // std::cout << "------------------------------------------------" << std::endl << ">>> Saving Preselection Snapshot..." << std::endl;
  // HT_calc.Snapshot("Events", stdOutputFilePS);
  // std::cout << "Output File: " << outputFilePS << std::endl << "-------------------------------------------------" << std::endl;
  // }
  //----------------------------------------------------------
  //       Uncomment from here to the bottom if starting from a preselection file!!
  //----------------------------------------------------------
  //	auto HT_calc = rdf;	
  
  // ---------------------------------------------------------
  // 		Post Preselection Analysis
  // ---------------------------------------------------------
  auto postPresel = HT_calc.Define("genttbarMass",genttbarMassCalc,{"nGenPart","GenPart_pdgId","GenPart_mass", \
	"GenPart_pt","GenPart_phi","GenPart_eta",			\
	"GenPart_genPartIdxMother","GenPart_status"})			\
    .Define("lepton_lv","lvConstructor(lepton_pt,lepton_eta,lepton_phi,lepton_mass)") \
    .Define("Jets_lv","fVectorConstructor(gcJet_pt,gcJet_eta,gcJet_phi,gcJet_mass)") \
    .Define("FatJet_lv","fVectorConstructor(gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass)") \
    .Define("Jet_ST","Jet_HT + lepton_pt + MET_pt")			\
    .Define("FatJet_pt_1","FatJet_pt[0]")				\
    .Define("FatJet_pt_2","FatJet_pt[1]")				\
    .Define("FatJet_sdMass","FatJet_msoftdrop[goodcleanFatJets == true]") \
    .Define("FatJet_sdMass_1","FatJet_sdMass[0]")			\
    .Define("FatJet_sdMass_2","FatJet_sdMass[1]")			\
    .Define("dpak8_J","FatJet_deepTag_QCDothers[goodcleanFatJets == true]") \
    .Define("dpak8_J_1","dpak8_J[0]")					\
    .Define("dpak8_J_2","dpak8_J[1]")					\
    .Define("raw_dpak8_T","(FatJet_deepTag_TvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_TvsQCD)") \
    .Define("dpak8_T","raw_dpak8_T[goodcleanFatJets == true]")		\
    .Define("dpak8_T_1","dpak8_T[0]")					\
    .Define("dpak8_T_2","dpak8_T[1]")					\
    .Define("raw_dpak8_W","(FatJet_deepTag_WvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_WvsQCD)") \
    .Define("dpak8_W","raw_dpak8_W[goodcleanFatJets == true]")		\
    .Define("dpak8_W_1","dpak8_W[0]")					\
    .Define("dpak8_W_2","dpak8_W[1]")					\
    .Define("dpak8_tag","maxFxn(dpak8_J,dpak8_T,dpak8_W)")		\
    .Define("dpak8_tag_1","dpak8_tag[0]")				\
    .Define("dpak8_tag_2","dpak8_tag[1]")				\
    .Define("nJ_dpak8","Sum(dpak8_tag == 0)")				\
    .Define("nT_dpak8","Sum(dpak8_tag == 1)")				\
    .Define("nW_dpak8","Sum(dpak8_tag == 2)")				\
    .Define("pNet_J","FatJet_particleNet_QCD[goodcleanFatJets == true]") \
    .Define("pNet_J_1","pNet_J[0]")					\
    .Define("pNet_J_2","pNet_J[1]")					\
    .Define("raw_pNet_T","(FatJet_particleNet_TvsQCD * FatJet_particleNet_QCD) / (1 - FatJet_particleNet_TvsQCD)") \
    .Define("pNet_T","raw_pNet_T[goodcleanFatJets == true]")		\
    .Define("pNet_T_1","pNet_T[0]")					\
    .Define("pNet_T_2","pNet_T[1]")					\
    .Define("raw_pNet_W","(FatJet_particleNet_WvsQCD * FatJet_particleNet_QCD) / (1 - FatJet_particleNet_WvsQCD)") \
    .Define("pNet_W","raw_pNet_W[goodcleanFatJets == true]")		\
    .Define("pNet_W_1","pNet_W[0]")					\
    .Define("pNet_W_2","pNet_W[1]")					\
    .Define("pNet_tag","maxFxn(pNet_J,pNet_T,pNet_W)")			\
    .Define("pNet_tag_1","pNet_tag[0]")					\
    .Define("pNet_tag_2","pNet_tag[1]")					\
    .Define("nJ_pNet","Sum(pNet_tag == 0)")				\
    .Define("nT_pNet","Sum(pNet_tag == 1)")				\
    .Define("nW_pNet","Sum(pNet_tag == 2)")				\
    .Define("raw_tau21","(FatJet_tau2 / FatJet_tau1)")			\
    .Define("tau21","raw_tau21[goodcleanFatJets == true]")		\
    .Define("tau21_1","tau21[0]")					\
    .Define("tau21_2","tau21[1]")					\
    .Define("minDR_ptRel_lead_lepAK8","minDR_ptRel_lead_calc(gcFatJet_pt,gcFatJet_eta,gcFatJet_phi, \
									   gcFatJet_mass,lepton_lv)")\
    .Define("minDR_lep_FatJet","minDR_ptRel_lead_lepAK8[0]")		\
    .Define("ptRel_lep_FatJet","minDR_ptRel_lead_lepAK8[1]")		\
    .Define("minDR_leadAK8otherAK8","minDR_ptRel_lead_lepAK8[2]")	\
    .Define("minDR_ptRel_lead_lepAK4","minDR_ptRel_lead_calc(gcJet_pt,gcJet_eta,gcJet_phi, \
					      				   gcJet_mass,lepton_lv)")\
    .Define("minDR_lep_Jet","minDR_ptRel_lead_lepAK4[0]")		\
    .Define("ptRel_lep_Jet","minDR_ptRel_lead_lepAK4[1]")		\
    .Define("DR_lep_FatJets","DR_calc(gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass, \
					      	    lepton_pt,lepton_eta, lepton_phi,lepton_mass)")\
    .Define("DR_lep_Jets","DR_calc(gcJet_pt,gcJet_eta,gcJet_phi,gcJet_mass, \
					   	 lepton_pt,lepton_eta,lepton_phi,lepton_mass)")\
    .Define("W_lv","W_reco(MET_pt,MET_phi,lepton_lv)")			\
    .Define("minMlj_output",minM_lep_jet_calc,{"gcJet_pt","gcJet_eta", "gcJet_phi","gcJet_mass", \
	  "lepton_lv"})							\
    .Define("DR_W_lep","dR_Wt_Calc(W_lv,lepton_lv)")			\
    .Define("minM_lep_Jet","minMlj_output[0]")				\
    .Define("minM_lep_Jet_jetID","(int) minMlj_output[1]")		\
    .Define("leptonicParticle","isLeptonic_X(minM_lep_Jet)")		\
    .Define("t_output","t_reco(leptonicParticle,gcJet_pt,gcJet_eta,gcJet_phi,gcJet_mass, \
					 W_lv,minM_lep_Jet,minM_lep_Jet_jetID)")\
    .Define("t_pt","t_output[0]")					\
    .Define("t_eta","t_output[1]")					\
    .Define("t_phi","t_output[2]")					\
    .Define("t_mass","t_output[3]")					\
    .Define("DR_W_b","t_output[4]")					\
    .Define("t_lv","lvConstructor(t_pt,t_eta,t_phi,t_mass)")
    .Define("Bprime_output","BPrime_reco(t_lv,W_lv,leptonicParticle,\
	  				       gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass,pNet_tag,gcFatJet_msoftdrop)")\
    .Define("Bprime_mass","Bprime_output[0]")				\
    .Define("Bprime_pt","Bprime_output[1]")				\
    .Define("Bprime_eta","Bprime_output[2]")				\
    .Define("Bprime_phi","Bprime_output[3]")				\
    .Define("Bprime_DR","Bprime_output[4]")				\
    .Define("Bprime_ptbal","Bprime_output[5]")				\
    .Define("Bprime_chi2","Bprime_output[6]")				\
    .Define("BPrime_lv","lvConstructor(Bprime_pt,Bprime_eta,Bprime_phi,Bprime_mass)") \
    .Define("isValidBDecay","Bprime_output[7]")				\
    .Define("taggedWbjetJet","Bprime_output[8]")			\
    .Define("taggedTjet","Bprime_output[9]")				\
    .Define("taggedWjet","Bprime_output[10]")				\
    .Define("dnn_scores",predictMLP,{"pNet_J_1","pNet_T_1","dpak8_T_1","FatJet_pt_1","FatJet_sdMass_1","tau21_1","nT_dpak8","nT_pNet","Jet_HT","Jet_ST","MET_pt","NJets_DeepFlavM","NJets_forward"}) \
    .Define("mlp_HT250_WJets","dnn_scores[0]")				\
    .Define("mlp_HT250_TTbar","dnn_scores[1]")				\
    .Define("mlp_HT250_Bprime","dnn_scores[2]")				\
    .Define("mlp_HT500_WJets","dnn_scores[3]")				\
    .Define("mlp_HT500_TTbar","dnn_scores[4]")				\
    .Define("mlp_HT500_Bprime","dnn_scores[5]")
    .Define("Bprime_gen_info", Bprime_gen_info, {"nGenPart", "GenPart_pdgId\
", "GenPart_mass", "GenPart_pt", "GenPart_phi", "GenPart_eta", "GenPart_genPartIdxMother", "GenPart_status"})
    .Define("Bprime_gen_pt", "Bprime_gen_info[0]")
    .Define("Bprime_gen_eta", "Bprime_gen_info[1]")
    .Define("Bprime_gen_phi", "Bprime_gen_info[2]")
    .Define("Bprime_gen_mass", "Bprime_gen_info[3]")
    .Define("Bprime_gen_pdgId", "Bprime_gen_info[4]")
    .Define("Bprime_gen_status", "Bprime_gen_info[5]")
    .Define("t_daughter_gen_info", t_daughter_gen, {"nGenPart", "GenPart_pdgId", "GenPart_mass", "GenPart_pt", "GenP\
art_phi", "GenPart_eta", "GenPart_genPartIdxMother", "GenPart_status"})
    .Define("t_gen_pt", "t_daughter_gen_info[0]")
    .Define("t_gen_eta", "t_daughter_gen_info[1]")
    .Define("t_gen_phi", "t_daughter_gen_info[2]")
    .Define("t_gen_mass", "t_daughter_gen_info[3]")
    .Define("t_gen_pdgId", "t_daughter_gen_info[4]")
    .Define("t_gen_status", "t_daughter_gen_info[5]")
    .Define("t_motherIdx", "t_daughter_gen_info[6]")
    .Define("daughterW_gen_pt", "t_daughter_gen_info[7]")
    .Define("daughterW_gen_eta", "t_daughter_gen_info[8]")
    .Define("daughterW_gen_phi", "t_daughter_gen_info[9]")
    .Define("daughterW_gen_mass", "t_daughter_gen_info[10]")
    .Define("daughterW_gen_pdgId", "t_daughter_gen_info[11]")
    .Define("daughterW_gen_status", "t_daughter_gen_info[12]")
    .Define("daughterW_motherIdx", "t_daughter_gen_info[13]")
    .Define("Tlepton_gen_pt", "t_daughter_gen_info[14]")
    .Define("Tlepton_gen_eta", "t_daughter_gen_info[15]")
    .Define("Tlepton_gen_phi", "t_daughter_gen_info[16]")
    .Define("Tlepton_gen_mass", "t_daughter_gen_info[17]")
    .Define("Tlepton_gen_pdgId", "t_daughter_gen_info[18]")
    .Define("Tlepton_gen_status", "t_daughter_gen_info[19]")
    .Define("trueLeptonicT", "t_daughter_gen_info[20]")
    .Define("W_daughter_gen_info", W_daughter_gen, {"nGenPart", "GenPart_pdgId", "GenPart_mass", "GenPart_pt", "GenPart_phi", "GenPart_eta", "GenPart_genPartIdxMother", "GenPart_status"})
    .Define("W_gen_pt", "W_daughter_gen_info[0]")
    .Define("W_gen_eta", "W_daughter_gen_info[1]")
    .Define("W_gen_phi", "W_daughter_gen_info[2]")
    .Define("W_gen_mass", "W_daughter_gen_info[3]")
    .Define("W_gen_pdgId", "W_daughter_gen_info[4]")
    .Define("W_gen_status", "W_daughter_gen_info[5]")
    .Define("W_motherIdx", "W_daughter_gen_info[6]")
    .Define("Wlepton_gen_pt", "W_daughter_gen_info[7]")
    .Define("wlepton_gen_eta", "W_daughter_gen_info[8]")
    .Define("Wlepton_gen_phi", "W_daughter_gen_info[9]")
    .Define("Wlepton_gen_mass", "W_daughter_gen_info[10]")
    .Define("Wlepton_gen_pdgId", "W_daughter_gen_info[11]")
    .Define("Wlepton_gen_status", "W_daughter_gen_info[12]")
    .Define("trueW_decayMode", "W_daughter_gen_info[13]");
  // -------------------------------------------------
  // 		Save Snapshot to file
  // -------------------------------------------------
  
  std::cout << "-------------------------------------------------" << std::endl << ">>> Saving " << sample << " Snapshot..." << std::endl;
  TString finalFile = "RDF_"+sample+"_finalsel_"+testNum+".root";
  const char* stdfinalFile = finalFile;
  postPresel.Snapshot("Events", stdfinalFile);
  std::cout << "Output File: " << finalFile << std::endl << "-------------------------------------------------" << std::endl;
  time.Stop();
  time.Print();
  std::cout << "Cut statistics:" << std::endl;
  postPresel.Report()->Print();
}
