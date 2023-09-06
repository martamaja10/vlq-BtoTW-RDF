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
#include <algorithm> // std::sort
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <sstream>
#include <chrono> // for high_resolution_clock
#include "../correctionlib/include/correction.h"

using namespace std;
using namespace ROOT::VecOps;
using correction::CorrectionSet;

void rdf::analyzer_RDF(TString testNum)
{
  ROOT::EnableImplicitMT();
  TStopwatch time;
  time.Start();
  string sample = this->sample;
  string year = this->year;
  
  cout << "Sample in cc: " << sample << endl;
  cout << "Year in cc: " << year << endl;
  cout << "isMC? " << isMC << endl;
  if(!isMC) cout << "Data era = " << era << endl;
  
  // -------------------------------------------------------
  //               Golden JSON
  // -------------------------------------------------------
  std::string jsonfile;
  if(year == "2016" or year == "2016APV") jsonfile = "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt";
  else if(year == "2017") jsonfile = "Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt";
  else if(year == "2018") jsonfile = "Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt";
  else std::cout << "ERROR: Can't parse the year to assign a golden json file. Expected 2016, 2016APV, 2017, or 2018. Got: " << year << std::endl;
  const auto myLumiMask = lumiMask::fromJSON(jsonfile);
  //  std::cout << "Testing the JSON! Known good run/lumi returns: " << myLumiMask.accept(315257, 10) << ", and known bad run returns: " << myLumiMask.accept(315257, 90) << std::endl;
  auto goldenjson = [myLumiMask](unsigned int &run, unsigned int &luminosityBlock){return myLumiMask.accept(run, luminosityBlock);};
  
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

  // Electron SF for TightID without iso
  std::vector<float> elid_pts;
  std::vector<float> elid_etas = {-2.5, -2.0, -1.566, -1.442, -0.8, 0.0, 0.8, 1.442, 1.566, 2.0, 2.5};
  std::vector<std::vector<float>> elecidsfs;
  std::vector<std::vector<float>> elecidsfuncs;  

  // Electron SF for Triggers
  std::vector<float> elhlt_pts = {50,60,70,80,100,200,300,99999};
  std::vector<float> elhlt_etas = {0,0.8,1.442,1.566,2.0,2.5};
  std::vector<std::vector<float>> elechltsfs;
  std::vector<std::vector<float>> elechltuncs;
  if(year == "2016" or year == "2016APV"){
    elid_pts = {50,90,99999};
    elecidsfs = {
      {0.983, 0.989, 0.978, 0.974, 0.966, 0.985, 0.985, 0.951, 0.987, 0.973},
      {1.002, 0.989, 1.094, 0.983, 0.976, 0.998, 1.018, 1.062, 0.988, 0.995}
    };
    elecidsfuncs = {
      {0.009, 0.007, 0.120, 0.004, 0.004, 0.004, 0.004, 0.120, 0.007, 0.009},
      {0.019, 0.017, 0.061, 0.016, 0.007, 0.007, 0.016, 0.061, 0.017, 0.018},
    };
    elechltsfs = {
      {0.969994,0.999272,1.000000,0.955299,0.972568},
      {0.979406,0.976827,0.948699,0.925570,0.904762},
      {0.998922,0.948318,1.086539,0.935962,0.990274},
      {0.981869,0.979395,1.035172,0.991033,1.042024},
      {0.981854,0.999758,1.023505,1.006556,0.980000},
      {0.966784,0.990014,1.000000,1.029556,1.000000},
      {0.986579,1.000000,1.000000,1.000000,1.000000}
    };
    elechltuncs = {
      {0.113757,0.104129,0.000000,0.083539,0.074347},
      {0.087648,0.097016,0.157331,0.128106,0.095238},
      {0.078416,0.075684,0.094028,0.118257,0.080618},
      {0.078566,0.042963,0.036409,0.034765,0.043790},
      {0.060464,0.025592,0.024057,0.006599,0.020000},
      {0.071746,0.016706,0.000000,0.030430,0.000000},
      {0.043074,0.000000,1.000000,0.000000,1.000000}
    };
  }else if(year == "2017"){
    elid_pts = {50,100,200,99999};
    elecidsfs = {
      {0.936, 0.964, 0.949, 0.968, 0.973, 0.971, 0.968, 0.916, 0.963, 0.934},
      {0.979, 0.981, 1.046, 0.987, 0.984, 0.987, 0.989, 0.983, 0.975, 0.962},
      {0.937, 0.969, 0.930, 0.966, 0.971, 0.967, 1.005, 0.912, 0.952, 1.073}
    };
    elecidsfuncs = {
      {0.021, 0.007, 0.025, 0.003, 0.005, 0.005, 0.003, 0.025, 0.007, 0.021},
      {0.020, 0.020, 0.050, 0.010, 0.006, 0.006, 0.010, 0.050, 0.020, 0.020},
      {0.077, 0.037, 0.179, 0.024, 0.020, 0.020, 0.025, 0.186, 0.038, 0.081}
    };
    elechltsfs = {
      {0.958734,0.879279,0.800000,0.930872,0.819886},
      {0.972084,0.890878,1.091776,0.814388,0.934196},
      {0.936706,0.887334,1.062766,0.932065,0.935490},
      {0.939349,0.976644,0.868087,0.989084,1.057438},
      {0.993655,0.979765,0.940140,0.972614,0.955556},
      {0.985202,0.995961,1.000000,0.909091,1.000000},
      {0.977273,0.942952,1.000000,1.000000,1.000000}
    };
    elechltuncs = {
      {0.131348,0.166195,0.687125,0.291188,0.370062},
      {0.098831,0.124641,0.772932,0.203191,0.347070},
      {0.101623,0.145139,0.693954,0.236831,0.382964},
      {0.083849,0.121101,0.499959,0.204401,0.338514},
      {0.060334,0.083162,0.354400,0.153739,0.239221},
      {0.126580,0.187983,0.764234,0.345859,0.712837},
      {0.251020,0.428866,1.000000,0.851152,1.000000}
    };
  }else{
    elid_pts = {50,100,200,99999};
    elecidsfs = {
      {1.036, 0.982, 0.989, 0.970, 0.973, 0.974, 0.974, 0.980, 0.991, 1.013},
      {1.024, 0.993, 0.944, 0.980, 0.975, 0.985, 0.987, 1.079, 1.005, 1.004},
      {0.982, 1.026, 0.987, 0.956, 1.021, 0.972, 0.986, 1.085, 0.967, 1.064}
    };
    elecidsfuncs = {
      {0.028, 0.004, 0.040, 0.005, 0.005, 0.005, 0.005, 0.040, 0.004, 0.028},
      {0.026, 0.024, 0.074, 0.016, 0.009, 0.010, 0.016, 0.075, 0.024, 0.026},
      {0.083, 0.047, 0.159, 0.038, 0.027, 0.024, 0.038, 0.163, 0.042, 0.084}
    };
    elechltsfs = {
      {0.948653,0.946580,1.030069,0.789884,0.914935},
      {0.959176,0.975594,0.982969,0.878223,0.978510},
      {0.914667,0.921914,0.859925,0.907098,0.856968},
      {0.936572,0.964878,0.790190,0.930776,1.000734},
      {0.977059,0.974735,0.939183,0.934640,0.870107},
      {0.988220,0.981175,1.000000,0.956091,1.011764},
      {1.002034,0.974624,0.000000,0.857143,1.000000}
    };
    elechltuncs = {
      {0.089153,0.129760,0.493564,0.183087,0.280660},
      {0.068586,0.097296,0.411438,0.140571,0.206748},
      {0.071284,0.102117,0.422389,0.154512,0.228542},
      {0.059381,0.082177,0.313267,0.126428,0.191029},
      {0.042942,0.058278,0.269807,0.094799,0.145511},
      {0.090644,0.142771,1.000000,0.234939,0.530056},
      {0.187096,0.275891,0.000000,0.505971,1.000000}
    };
  }    
  
  // -------------------------------------------------------
  //               correctionLib corrections
  // -------------------------------------------------------
  
  std::string yrstr, yr, jecyr, jeryr, jecver;
  float deepjetL;
  if(year == "2016") {deepjetL = 0.0508; yrstr = "2016preVFP_UL"; yr = "16"; jecyr = "UL16"; jeryr = "Summer20UL16_JRV3"; jecver = "V7";}
  else if(year == "2016APV") {deepjetL = 0.0480; yrstr = "2016postVFP_UL"; yr = "16"; jecyr = "UL16APV"; jeryr = "Summer20UL16APV_JRV3"; jecver = "V7";}
  else if(year == "2017") {deepjetL = 0.0532; yrstr = "2017_UL"; yr = "17"; jecyr = "UL17"; jeryr = "Summer19UL17_JRV2"; jecver = "V5";}
  else if(year == "2018") {deepjetL = 0.0490; yrstr = "2018_UL"; yr = "18"; jecyr = "UL18"; jeryr = "Summer19UL18_JRV2"; jecver = "V5";}
  else std::cout << "ERROR: Can't parse the year to assign correctionLib json files. Expected 2016, 2016APV, 2017, or 2018. Got: " << year << std::endl;
  
  auto pileupcorrset = CorrectionSet::from_file("jsonpog-integration/POG/LUM/"+yrstr+"/puWeights.json");
  auto electroncorrset = CorrectionSet::from_file("jsonpog-integration/POG/EGM/"+yrstr+"/electron.json");
  auto muoncorrset = CorrectionSet::from_file("jsonpog-integration/POG/MUO/"+yrstr+"/muon_Z.json");
  auto btagcorrset = CorrectionSet::from_file("jsonpog-integration/POG/BTV/"+yrstr+"/btagging.json");
  auto ak4corrset = CorrectionSet::from_file("jsonpog-integration/POG/JME/"+yrstr+"/jet_jerc.json");
  //  auto ak8corrset = CorrectionSet::from_file("jsonpog-integration/POG/JME/"+yrstr+"/fatjet_jerc.json");
  auto metcorrset = CorrectionSet::from_file("jsonpog-integration/POG/JME/"+yrstr+"/met.json");
  auto jetvetocorrset = CorrectionSet::from_file("jsonpog-integration/POG/JME/"+yrstr+"/jetvetomaps.json");

  auto pileupcorr = pileupcorrset->at("Collisions"+yr+"_UltraLegacy_goldenJSON");
  auto electroncorr = electroncorrset->at("UL-Electron-ID-SF");
  auto muonidcorr = muoncorrset->at("NUM_HighPtID_DEN_genTracks");
  auto muonhltcorr = muoncorrset->at("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose");
  auto btagcorr = btagcorrset->at("deepJet_shape");
  if(isMC) std::cout << "loading jec Summer" << jecyr << "_" << jecver << "_MC_L1L2L3Res_AK4PFchs "<< std::endl;
  else std::cout << "loading jec Summer19" << jecyr << "_Run" << era << "_" << jecver << "_DATA_L1L2L3Res_AK4PFchs" << std::endl;
  auto ak4corr = ak4corrset->compound().at("Summer19"+jecyr+"_"+jecver+"_MC_L1L2L3Res_AK4PFchs");
  if(!isMC) ak4corr = ak4corrset->compound().at("Summer19"+jecyr+"_Run"+era+"_"+jecver+"_DATA_L1L2L3Res_AK4PFchs"); // call evaluate with area, eta, py, rho
  auto ak4corrUnc = ak4corrset->at("Summer19"+jecyr+"_"+jecver+"_MC_Total_AK4PFchs"); // call evaluate with eta. pt
  std:: cout << "loading jer" << std::endl;
  auto ak4ptres = ak4corrset->at(jeryr+"_MC_PtResolution_AK4PFchs"); // call evaluate with eta, pt, rho
  auto ak4jer = ak4corrset->at(jeryr+"_MC_ScaleFactor_AK4PFchs"); // call evaluate with eta, shift --> figure out the smear as before
  //  auto ak8corr = ak8corrset->at("");
  auto metcorr_ptdata = metcorrset->at("pt_metphicorr_pfmet_data");
  auto metcorr_phidata = metcorrset->at("phi_metphicorr_pfmet_data");
  auto metcorr_ptmc = metcorrset->at("pt_metphicorr_pfmet_mc");
  auto metcorr_phimc = metcorrset->at("phi_metphicorr_pfmet_mc");
  auto jetvetocorr = jetvetocorrset->at("Summer19UL"+yr+"_V1");
  
  auto pufunc = [pileupcorr](const float &numTrueInt){
    RVec<double> pu = {pileupcorr->evaluate({numTrueInt, "nominal"}), pileupcorr->evaluate({numTrueInt, "up"}), pileupcorr->evaluate({numTrueInt, "down"})};
    return pu;
  };
  auto recofunc = [electroncorr, year](const float &pt, const float &eta, const bool &isEl){
    if(isEl == 0) { RVec<double> mu = {1.0, 1.0, 1.0}; return mu; }
    else{
      RVec<double> el = {electroncorr->evaluate({year,"sf","RecoAbove20",eta,pt}), 
			 electroncorr->evaluate({year,"sfup","RecoAbove20",eta,pt}), 
      			 electroncorr->evaluate({year,"sfdown","RecoAbove20",eta,pt})};
      return el;
    }
  }; 
  auto idfunc = [muonidcorr,elid_pts,elid_etas,elecidsfs,elecidsfuncs,yrstr](const float &pt, const float &eta, const bool &isEl){
    RVec<double> id;
    if(isEl > 0){
      int ptbin = (std::upper_bound(elid_pts.begin(), elid_pts.end(), pt) - elid_pts.begin())-1;
      int etabin = (std::upper_bound(elid_etas.begin(), elid_etas.end(), eta) - elid_etas.begin())-1;
      id = {elecidsfs[ptbin][etabin], elecidsfuncs[ptbin][etabin]};      
    }else{
      id = {muonidcorr->evaluate({yrstr,abs(eta),pt,"sf"}), 
	    muonidcorr->evaluate({yrstr,abs(eta),pt,"systup"}), 
	    muonidcorr->evaluate({yrstr,abs(eta),pt,"systdown"})};
    }
    return id;
  }; 
  // // FIXME: calculate MiniIso SFs for UL. Put in an "isofunc". They are usually == 1.
  auto hltfunc = [muonhltcorr,elhlt_pts,elhlt_etas,elechltsfs,elechltuncs,year,yrstr](const float &pt, const float &eta, const bool &isEl){
    RVec<double> hlt;
    if(isEl > 0){
      int ptbin = (std::upper_bound(elhlt_pts.begin(), elhlt_pts.end(), pt) - elhlt_pts.begin())-1;
      int etabin = (std::upper_bound(elhlt_etas.begin(), elhlt_etas.end(), eta) - elhlt_etas.begin())-1;
      hlt = {elechltsfs[ptbin][etabin], elechltuncs[ptbin][etabin]};      
    }else{
      hlt = {muonhltcorr->evaluate({yrstr,abs(eta),pt,"sf"}), 
	     muonhltcorr->evaluate({yrstr,abs(eta),pt,"systup"}), 
	     muonhltcorr->evaluate({yrstr,abs(eta),pt,"systdown"})};
    }
    return hlt;
  }; 
  auto metfunc = [metcorr_ptdata,metcorr_ptmc](const float &met, const float &metphi, const float &npv, const float &run){
    if(run > 100000) return metcorr_ptdata->evaluate({met, metphi, npv, run});
    else return metcorr_ptmc->evaluate({met, metphi, npv, run});
  };
  auto metphifunc = [metcorr_phidata,metcorr_phimc](const float &met, const float &metphi, const float &npv, const float &run){
    if(run > 100000) return metcorr_phidata->evaluate({met, metphi, npv, run});
    else return metcorr_phimc->evaluate({met, metphi, npv, run});
  };
  auto jetvetofunc = [jetvetocorr](const RVec<float> &eta, const RVec<float> &phi){
    RVec<double> map;
    for(unsigned int ijet = 0; ijet < eta.size(); ijet++){
      float phitemp = phi.at(ijet);
      if(phitemp < -3.14159) phitemp = -3.14159;
      else if(phitemp > 3.14159) phitemp = 3.14159;
      map.push_back(jetvetocorr->evaluate({"jetvetomap",eta.at(ijet),phitemp}));
    }
    return map;
  };
  
  std::string nominal = "central";
  if(jesvar == "JECup") nominal = "up_jes";
  else if (jesvar == "JECdn") nominal = "down_jes";
  auto btagshapefunc = [btagcorr,nominal](const RVec<float> &pt, const RVec<float> &eta, const RVec<float> &disc, const RVec<int> &flav){
    RVec<float> weights(17, 1.0); // collect product of SFs over jets
    for(unsigned int ijet = 0; ijet < eta.size(); ijet++){
      weights[0] *= btagcorr->evaluate({nominal, flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
      if(flav.at(ijet) != 4){
	weights[1] *= btagcorr->evaluate({"up_hf", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[2] *= btagcorr->evaluate({"down_hf", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[3] *= btagcorr->evaluate({"up_lf", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[4] *= btagcorr->evaluate({"down_lf", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[5] *= btagcorr->evaluate({"up_hfstats1", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[6] *= btagcorr->evaluate({"down_hfstats1", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[7] *= btagcorr->evaluate({"up_hfstats2", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[8] *= btagcorr->evaluate({"down_hfstats2", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[9] *= btagcorr->evaluate({"up_lfstats1", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[10] *= btagcorr->evaluate({"down_lfstats1", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[11] *= btagcorr->evaluate({"up_lfstats2", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[12] *= btagcorr->evaluate({"down_lfstats2", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
      }
      else{
	weights[13] *= btagcorr->evaluate({"up_cferr1", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[14] *= btagcorr->evaluate({"down_cferr1", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[15] *= btagcorr->evaluate({"up_cferr2", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
	weights[16] *= btagcorr->evaluate({"down_cferr2", flav.at(ijet), abs(eta.at(ijet)), pt.at(ijet), disc.at(ijet)});
      }
    }
    return weights;      
  };
  
  // FIXME: figure out the JERC...

  // -------------------------------------------------------
  //               Open Dataframe + MET Filters
  // -------------------------------------------------------

  auto rdf_input = ROOT::RDataFrame("Events", files); // Initial data
  
  auto METgeneralFilters = rdf_input.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1", "MET Filters")
    .Filter("nJet > 0 && nFatJet > 0", "Event has > 1 AK4 and > 1 AK8");
  auto truth = METgeneralFilters;
  
  // --------------------------------------------------------
  // 	       Golden JSON (Data) || GEN Info (MC)
  // --------------------------------------------------------
  
  if(!isMC){ // apply golden json to data
    truth = METgeneralFilters.Define("passesJSON", goldenjson, {"run","luminosityBlock"})
      .Filter("passesJSON == 1", "Data passes Golden JSON");
  }
  
  else{ // calculate some gen-level stuff
    auto BprimeGen = METgeneralFilters.Define("Bprime_gen_info", Form("Bprime_gen_info(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status, GenPart_statusFlags)", sample.c_str()))
      .Define("Bprime_gen_pt", "Bprime_gen_info[0]")
      .Define("Bprime_gen_eta", "(double) Bprime_gen_info[1]")
      .Define("Bprime_gen_phi", "(double) Bprime_gen_info[2]")
      .Define("Bprime_gen_mass", "Bprime_gen_info[3]")
      .Define("Bprime_gen_pdgId", "(int) Bprime_gen_info[4]")
      .Define("Bprime_gen_status", "(int) Bprime_gen_info[5]");
    
    auto wjetHT = BprimeGen;
    if(!isVV) wjetHT = BprimeGen.Define("HTCorr_WjetLHE", wjetHTpoly, {"LHE_HT"}); // LHE_HT missing for Pythia-only samples
    
    truth = wjetHT.Define("t_gen_info", Form("t_gen_info(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status)", sample.c_str()))
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
      .Define("W_bkg_idx", Form("W_bkg_idx(\"%s\", nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags, t_bkg_idx)", sample.c_str()))
      .Define("PileupWeights", pufunc, {"Pileup_nTrueInt"});

  }
  
  // ---------------------------------------------------------
  //                    MET Filters
  // ---------------------------------------------------------
  
  
  auto METptFilters = truth.Filter("MET_pt > 60", "Pass MET > 60");
  
  // ---------------------------------------------------------
  //                    LEPTON Definitions
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
    .Define("VElectron_jetIdx", "Electron_jetIdx[VetoEl]")
    .Define("nSignalIsoMu", "(int) Sum(SignalIsoMu)")
    .Define("nSignalIsoEl", "(int) Sum(SignalIsoEl)")
    .Define("VetoIsoMu", "(VMuon_pt<55)")
    .Define("VetoIsoEl", "(VElectron_pt<80)")
    .Define("nVetoIsoLep", "(int) (Sum(VetoIsoMu)+Sum(VetoIsoEl))");

  // --------------------------------------------------------
  // 		      LEPTON SELECTION
  // --------------------------------------------------------
  
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
    .Define("goodcleanJets", "cleanJet_pt > 30 && abs(cleanJet_eta) < 2.5 && Jet_jetId > 1 && (DR_lepJets > 0.4 || ptrel_lepJets > 20)")
    .Define("NJets_central", "(int) Sum(goodcleanJets)")
    .Define("gcJet_pt", "cleanJet_pt[goodcleanJets == true]")
    .Define("gcJet_eta", "cleanJet_eta[goodcleanJets == true]")
    .Define("gcJet_phi", "cleanJet_phi[goodcleanJets == true]")
    .Define("gcJet_mass", "cleanJet_mass[goodcleanJets == true]")
    .Define("gcJet_vetomap", jetvetofunc, {"gcJet_eta","gcJet_phi"})
    .Define("gcJet_DeepFlav", "Jet_btagDeepFlavB[goodcleanJets == true]")
    .Define("gcJet_DeepFlavL", Form("gcJet_DeepFlav > %f",deepjetL)) 
    .Define("NJets_DeepFlavL", "(int) Sum(gcJet_DeepFlavL)")
    .Define("gcJet_DeepFlavL_pt", "gcJet_pt[gcJet_DeepFlavL == true]")
    .Define("gcJet_DeepFlavL_eta", "gcJet_eta[gcJet_DeepFlavL == true]")
    .Define("gcJet_DeepFlavL_phi", "gcJet_phi[gcJet_DeepFlavL == true]")
    .Define("gcJet_DeepFlavL_mass", "gcJet_mass[gcJet_DeepFlavL == true]")
    .Define("goodcleanForwardJets", "cleanJet_pt > 30 && abs(cleanJet_eta) >= 2.5 && Jet_jetId > 1")
    .Define("NJets_forward", "(int) Sum(goodcleanForwardJets)")
    .Define("gcforwJet_pt", "cleanJet_pt[goodcleanForwardJets == true]")
    .Define("gcforwJet_eta", "cleanJet_eta[goodcleanForwardJets == true]")
    .Define("gcforwJet_phi", "cleanJet_phi[goodcleanForwardJets == true]")
    .Define("gcforwJet_mass", "cleanJet_mass[goodcleanForwardJets == true]")
    .Define("gcforwJet_DeepFlav", "Jet_btagDeepFlavB[goodcleanForwardJets == true]")
    .Define("DR_lepFatJets","DeltaR_VecAndFloat(FatJet_eta,FatJet_phi,lepton_eta,lepton_phi)")
    .Define("goodcleanFatJets", "FatJet_pt > 200 && abs(FatJet_eta) < 2.5 && FatJet_jetId > 1 && (DR_lepFatJets > 0.8)")
    .Define("NFatJets", "(int) Sum(goodcleanFatJets)")
    .Define("gcFatJet_pt", "FatJet_pt[goodcleanFatJets == true]")
    .Define("gcFatJet_eta", "FatJet_eta[goodcleanFatJets == true]")
    .Define("gcFatJet_phi", "FatJet_phi[goodcleanFatJets == true]")
    .Define("gcFatJet_mass", "FatJet_mass[goodcleanFatJets == true]")
    .Define("gcFatJet_sdmass", "FatJet_msoftdrop[goodcleanFatJets == true]")
    .Define("gcFatJet_vetomap", jetvetofunc, {"gcFatJet_eta","gcFatJet_phi"})
    .Define("DR_gcJets_central","DR_lepJets[goodcleanJets == true]")
    .Define("minDR_lepJets","ROOT::VecOps::Min(DR_gcJets_central)")
    .Define("ptrel_atMinDR_lepJets","ptrel_lepJets[goodcleanJets == true][ROOT::VecOps::ArgMin(DR_gcJets_central)]")
    .Define("DR_gcJets_DeepFlavL","DR_gcJets_central[gcJet_DeepFlavL == true]")
    .Define("DR_gcFatJets","DR_lepFatJets[goodcleanFatJets == true]")
    .Define("minDR_lepFatJets","ROOT::VecOps::Min(DR_gcFatJets)")
    .Define("ptrel_atMinDR_lepFarJets","ptRel(gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass,lepton_pt,lepton_eta,lepton_phi,lepton_mass)[ROOT::VecOps::ArgMin(DR_gcFatJets)]")
    .Define("OS_gcJets","DR_gcJets_central > TMath::Pi()/2")
    .Define("SS_gcJets","DR_gcJets_central <= TMath::Pi()/2")
    .Define("OS_gcJets_DeepFlavL","DR_gcJets_DeepFlavL > TMath::Pi()/2")
    .Define("SS_gcJets_DeepFlavL","DR_gcJets_DeepFlavL <= TMath::Pi()/2")
    .Define("OS_gcFatJets","DR_gcFatJets > TMath::Pi()/2")
    .Define("SS_gcFatJets","DR_gcFatJets <= TMath::Pi()/2")
    .Define("NOS_gcJets_central","(int) Sum(OS_gcJets)")
    .Define("NSS_gcJets_central","(int) Sum(SS_gcJets)")
    .Define("NOS_gcJets_DeepFlavL","(int) Sum(OS_gcJets_DeepFlavL)")
    .Define("NSS_gcJets_DeepFlavL","(int) Sum(SS_gcJets_DeepFlavL)")
    .Define("NOS_gcFatJets","(int) Sum(OS_gcFatJets)")
    .Define("NSS_gcFatJets","(int) Sum(SS_gcFatJets)")
    .Define("gcOSFatJet_pt","gcFatJet_pt[OS_gcFatJets == true]")
    .Define("gcOSFatJet_eta","gcFatJet_eta[OS_gcFatJets == true]")
    .Define("gcOSFatJet_phi","gcFatJet_phi[OS_gcFatJets == true]")
    .Define("gcOSFatJet_mass","gcFatJet_mass[OS_gcFatJets == true]")
    .Define("gcOSFatJet_sdmass","gcFatJet_sdmass[OS_gcFatJets == true]")
    .Define("gcOSJet_pt","gcJet_pt[OS_gcJets == true]")
    .Define("gcOSJet_eta","gcJet_eta[OS_gcJets == true]")
    .Define("gcOSJet_phi","gcJet_phi[OS_gcJets == true]")
    .Define("gcOSJet_mass","gcJet_mass[OS_gcJets == true]")
    .Define("gcSSJet_pt","gcJet_pt[SS_gcJets == true]")
    .Define("gcSSJet_eta","gcJet_eta[SS_gcJets == true]")
    .Define("gcSSJet_phi","gcJet_phi[SS_gcJets == true]")
    .Define("gcSSJet_mass","gcJet_mass[SS_gcJets == true]");
  
  // ---------------------------------------------------------
  // 	  HT Calculation and initial Cuts
  // ---------------------------------------------------------
  
  auto HT_calc = CleanJets.Define("gcJet_HT","Sum(Jet_pt[goodcleanJets == true])") \
    .Filter("gcJet_HT > 250","Pass HT > 250")				\
    .Define("HTCorr_top", topHTpoly, {"gcJet_HT"})
    .Filter("NFatJets > 0","Pass N good central AK8 > 0")
    .Filter("NOS_gcFatJets > 0","Pass N good central other side AK8 > 0");
  
  
  // ---------------------------------------------------------
  // 		JET Tagging and RECONSTRUCTION
  // ---------------------------------------------------------
  auto genttbarJets = HT_calc;
  
  if (isMC) {
    genttbarJets = HT_calc.Define("genttbarMass", Form("genttbarMassCalc(\"%s\", nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status)",sample.c_str()))
      .Define("gcFatJet_genmatch", Form("FatJet_matching_bkg(\"%s\", goodcleanFatJets, gcFatJet_eta, gcFatJet_phi, NFatJets, FatJet_subJetIdx1, nSubJet, SubJet_hadronFlavour, nGenPart, GenPart_pdgId, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, t_bkg_idx, W_bkg_idx)",sample.c_str()))
      .Define("elRecoSF", recofunc, {"lepton_pt","lepton_eta","isEl"})
      .Define("leptonIDSF", idfunc, {"lepton_pt","lepton_eta","isEl"})
      .Define("leptonHLTSF", hltfunc, {"lepton_pt","lepton_eta","isEl"})
      .Define("gcJet_hflav","Jet_hadronFlavour[goodcleanJets == true]")
      .Define("btagWeights",btagshapefunc, {"gcJet_pt", "gcJet_eta", "gcJet_DeepFlav", "gcJet_hflav"});
    // FIXME add iso 
  }
  auto postPresel = genttbarJets.Define("lepton_lv", "lvConstructor(lepton_pt,lepton_eta,lepton_phi,lepton_mass)")
    .Define("gcJet_ST", "gcJet_HT + lepton_pt + MET_pt")
    .Define("gcFatJet_pNetJ", "FatJet_particleNet_QCD[goodcleanFatJets == true]")
    .Define("gcOSFatJet_pNetJ", "gcFatJet_pNetJ[OS_gcFatJets == true]") 
    .Define("gcFatJet_pNetT", "((FatJet_particleNet_TvsQCD * FatJet_particleNet_QCD) / (1 - FatJet_particleNet_TvsQCD))[goodcleanFatJets == true]")
    .Define("gcFatJet_pNetTvsQCD", "FatJet_particleNet_TvsQCD[goodcleanFatJets == true]")
    .Define("gcOSFatJet_pNetT", "gcFatJet_pNetT[OS_gcFatJets == true]")
    .Define("gcOSFatJet_pNetTvsQCD", "gcFatJet_pNetTvsQCD[OS_gcFatJets == true]") 
    .Define("gcFatJet_pNetW", "((FatJet_particleNet_WvsQCD * FatJet_particleNet_QCD) / (1 - FatJet_particleNet_WvsQCD))[goodcleanFatJets == true]")
    .Define("gcFatJet_pNetWvsQCD", "FatJet_particleNet_WvsQCD[goodcleanFatJets == true]")
    .Define("gcOSFatJet_pNetW", "gcFatJet_pNetW[OS_gcFatJets == true]")
    .Define("gcOSFatJet_pNetWvsQCD", "gcFatJet_pNetWvsQCD[OS_gcFatJets == true]")
    .Define("gcFatJet_pNetTag", "maxFxn(gcFatJet_pNetJ,gcFatJet_pNetT,gcFatJet_pNetW)")
    .Define("gcFatJet_pNetTag_alt", "JetDiscriminator(gcFatJet_pNetTvsQCD, gcFatJet_pNetWvsQCD)")
    .Define("gcOSFatJet_pNetTag", "gcFatJet_pNetTag[OS_gcFatJets==true]")
    .Define("gcOSFatJet_pNetTag_alt", "gcFatJet_pNetTag_alt[OS_gcFatJets==true]")
    .Define("gcFatJet_nJ", "Sum(gcFatJet_pNetTag == 0)")
    .Define("gcFatJet_nT", "Sum(gcFatJet_pNetTag == 1)")
    .Define("gcFatJet_nW", "Sum(gcFatJet_pNetTag == 2)")
    .Define("gcFatJet_tau21", "(FatJet_tau2 / FatJet_tau1)[goodcleanFatJets == true]")
    .Define("gcOSFatJet_tau21", "gcFatJet_tau21[OS_gcFatJets == true]")
    .Define("gcFatJet_tau32", "(FatJet_tau3 / FatJet_tau2)[goodcleanFatJets == true]")
    .Define("gcOSFatJet_tau32", "gcFatJet_tau32[OS_gcFatJets == true]")
    .Define("minDR_leadAK8otherAK8", "minDR_leadJetOtherJet_calc(gcFatJet_eta,gcFatJet_phi)")
    .Define("minDR_leadAK4otherAK4", "minDR_leadJetOtherJet_calc(gcJet_eta,gcJet_phi)")
    .Define("minDR_AK8s_discrete","std::floor(minDR_leadAK8otherAK8/0.5)")
    .Define("minDR_AK4s_discrete","std::floor(minDR_leadAK4otherAK4/0.5)")
    .Define("W_lv", "W_reco(MET_pt,MET_phi,lepton_lv)")
    .Define("W_pt", "W_lv.Pt()")
    .Define("W_eta", "W_lv.Eta()")
    .Define("W_phi", "W_lv.Phi()")
    .Define("W_mass", "W_lv.M()")
    .Define("W_MT", "sqrt(2*lepton_pt*MET_pt*(1-cos(lepton_phi - MET_phi)))")
    .Define("minMlj_output", "minM_lep_jet_calc(gcJet_pt, gcJet_eta, gcJet_phi, gcJet_mass, lepton_lv)")
    .Define("DR_W_lep", "dR_Wt_Calc(W_lv,lepton_lv)")
    .Define("minM_lep_Jet", "minMlj_output[0]")
    .Define("minM_lep_Jet_jetID", "(int) minMlj_output[1]")
    .Define("minM_lep_Jet_TorW", "isLeptonic_X(minM_lep_Jet)")
    .Define("t_output", "t_reco(minM_lep_Jet_TorW,gcJet_pt,gcJet_eta,gcJet_phi,gcJet_mass, W_lv,minM_lep_Jet,minM_lep_Jet_jetID)")
    .Define("t_pt", "t_output[0]")
    .Define("t_eta", "t_output[1]")
    .Define("t_phi", "t_output[2]")
    .Define("t_mass", "t_output[3]")
    .Define("DR_W_b", "t_output[4]")
    .Define("t_lv", "lvConstructor(t_pt,t_eta,t_phi,t_mass)")
    .Define("Bprime_output", Form("BPrime_reco_new(t_lv,W_lv,\"%s\",minM_lep_Jet,minM_lep_Jet_jetID,NOS_gcJets_central,NOS_gcJets_DeepFlavL,NSS_gcJets_central,NSS_gcJets_DeepFlavL,SS_gcJets_DeepFlavL,OS_gcJets_DeepFlavL,gcOSFatJet_pt,gcOSFatJet_eta,gcOSFatJet_phi,gcOSFatJet_mass,gcOSFatJet_pNetTag,gcOSJet_pt,gcOSJet_eta,gcOSJet_phi,gcOSJet_mass,gcSSJet_pt,gcSSJet_eta,gcSSJet_phi,gcSSJet_mass)",sample.c_str()))
    .Define("Bprime_mass", "Bprime_output[0]")
    .Define("Bprime_pt", "Bprime_output[1]")
    .Define("Bprime_eta", "Bprime_output[2]")
    .Define("Bprime_phi", "Bprime_output[3]")
    .Define("Bprime_DR", "Bprime_output[4]")
    .Define("Bprime_ptbal", "Bprime_output[5]")
    .Define("Bprime_chi2", "Bprime_output[6]")
    .Define("Bdecay_obs", "Bprime_output[7]");
  
  // -------------------------------------------------
  // 		Save Snapshot to file
  // -------------------------------------------------
  
  cout << "-------------------------------------------------" << endl
       << ">>> Saving " << sample << " Snapshot..." << endl;
  TString finalFile = "RDF_" + sample + "_finalsel_" + testNum.Data() + ".root";
  const char *stdfinalFile = finalFile;
  
  auto ColNames = postPresel.GetColumnNames();
  vector<string> snapCol;
  int i = 0;
  for (auto &&ColName : ColNames)
    {
      TString colName = ColName;
      if(colName.Contains("P4") || colName.Contains("cleanJets")) continue;
      if(colName.Contains("LHE") && !colName.Contains("Weight") && colName != "LHE_HT" && colName != "LHE_Vpt") continue;
      if(colName.BeginsWith("Muon") && !colName.Contains("_tightId") && !colName.Contains("_isPF") && !colName.Contains("tunep") && !colName.Contains("genPartFlav")) continue;
      if(colName.BeginsWith("Electron") && !colName.Contains("genPartFlav")) continue;
      if(colName.BeginsWith("Jet") && !colName.Contains("rawFactor")) continue;
      if(colName.BeginsWith("FatJet") && !colName.Contains("rawFactor")) continue;
      if(colName.BeginsWith("PPS") || colName.BeginsWith("Proton") || colName.BeginsWith("L1_")) continue;
      if(colName.BeginsWith("Gen") || colName.BeginsWith("Soft") || colName.BeginsWith("fixed")) continue;
      if(colName.BeginsWith("Sub") || colName.BeginsWith("Raw") || colName.BeginsWith("Calo") || colName.BeginsWith("Chs")) continue;
      if(colName.BeginsWith("Corr") || colName.BeginsWith("Fsr") || colName.BeginsWith("Iso") || colName.BeginsWith("Tau")) continue;
      if(colName.BeginsWith("SV") || colName.BeginsWith("Puppi") || colName.BeginsWith("Photon") || colName.BeginsWith("Low")) continue;
      if(colName.BeginsWith("HLT") || colName.BeginsWith("HT") || colName.BeginsWith("boosted") || colName.BeginsWith("Deep")) continue;
      if(colName.BeginsWith("Flag") || colName == "Bprime_gen_info" || colName == "t_gen_info" || colName == "W_gen_info") continue;
      if(colName == "assignleps" || colName == "t_output" || colName == "Bprime_output" || colName.BeginsWith("Other")) continue;
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
  
  postPresel.Snapshot("Events", stdfinalFile, snapCol);
  cout << "Output File: " << finalFile << endl
       << "-------------------------------------------------" << endl;
  
  time.Stop();
  time.Print();
  
  cout << "Cut statistics:" << endl;
  postPresel.Report()->Print();
  
  cout << "Adding Counter tree to the file:" << endl;
  auto rdf_runs = ROOT::RDataFrame("Runs", files);
  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "UPDATE";
  rdf_runs.Snapshot("Runs", stdfinalFile, rdf_runs.GetColumnNames(), opts);

  delete poly2; delete poly2U; delete poly2D;
  delete polyHT; delete polyHTU; delete polyHTD;
  
  cout << "Done!" << endl;
}
