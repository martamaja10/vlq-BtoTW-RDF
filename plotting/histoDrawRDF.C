//To run on command line root -l histoDrawRDF.C

#define rdf_cxx
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

using namespace ROOT::VecOps;
void histoDrawRDF()
{
//	ROOT::EnableImplicitMT();
	TFile *outputFile=new TFile("seniorResearchHistograms_withDataObs.root","RECREATE");

	// ------------------------------------------------------------------------------------------------
	// Signal Sample Mass 800 GeV
	// ------------------------------------------------------------------------------------------------
	int eventsNumber_sig800 = 99800;
	ROOT::RDataFrame a("Events", "RDF_Signal_seniorResearch_Muon__m800.root");
	auto b = a.Filter("Bprime1_DeepAK8_Mass && isValidBDecay == true").Define("weights_sig800","Generator_weight*137000*1/(abs(Generator_weight)*99800)");
	auto c = b.Filter("taggedWbjetJet == true");
	auto d = b.Filter("taggedTjet == true");
	auto e = b.Filter("taggedWjet == true");

	auto sig800_Wb = c.Histo1D({"sig800_Wb",";Mass (GeV);Events / bin", 25, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_sig800");
	sig800_Wb->Write();
	auto sig800_W = d.Histo1D({"sig800_W",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_sig800");
	sig800_W->Write();
	auto sig800_T = e.Histo1D({"sig800_T",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_sig800");
	sig800_T->Write();
	// ------------------------------------------------------------------------------------------------

	// ------------------------------------------------------------------------------------------------
	// Signal Sample Mass 1400 GeV
	// ------------------------------------------------------------------------------------------------
	int eventsNumber_sig1400 = 99600;
	ROOT::RDataFrame i("Events", "RDF_Signal_seniorResearch_Muon__m1400.root");
	auto j = i.Filter("Bprime1_DeepAK8_Mass && isValidBDecay == true").Define("weights_sig1400","Generator_weight*137000*1/(abs(Generator_weight)*99600)");
	auto k = j.Filter("taggedWbjetJet == true");
	auto l = j.Filter("taggedTjet == true");
	auto m = j.Filter("taggedWjet == true");
	
	auto sig1400_Wb = k.Histo1D({"sig1400_Wb",";Mass (GeV);Events / bin", 25, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_sig1400");
	sig1400_Wb->Write();
	auto sig1400_W = l.Histo1D({"sig1400_W",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_sig1400");
	sig1400_W->Write();
	auto sig1400_T = m.Histo1D({"sig1400_T",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_sig1400");
	sig1400_T->Write();
	// ------------------------------------------------------------------------------------------------

	// ------------------------------------------------------------------------------------------------
	// Signal Sample Mass 2000 GeV
	// ------------------------------------------------------------------------------------------------
	int eventsNumber_sig2000 = 99600;
	ROOT::RDataFrame q("Events", "RDF_Signal_seniorResearch_Muon__m2000.root");
	auto r = q.Filter("Bprime1_DeepAK8_Mass && isValidBDecay == true").Define("weights_sig2000","Generator_weight*137000*1/(abs(Generator_weight)*99600)");
	auto s = r.Filter("taggedWbjetJet == true");
	auto t = r.Filter("taggedTjet == true");
	auto u = r.Filter("taggedWjet == true");
	
	auto sig2000_Wb = s.Histo1D({"sig2000_Wb",";Mass (GeV);Events / bin", 25, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_sig2000");
	sig2000_Wb->Write();
	auto sig2000_W = t.Histo1D({"sig2000_W",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_sig2000");
	sig2000_W->Write();
	auto sig2000_T = u.Histo1D({"sig2000_T",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_sig2000");
	sig2000_T->Write();
	// ------------------------------------------------------------------------------------------------

	// ------------------------------------------------------------------------------------------------
	// TTToSemiLeptonic Background File
	// ------------------------------------------------------------------------------------------------
	//int eventsNumber_TTToSemiLeptonic = 6544000;
	int eventsNumber_TTToSemiLeptonic = 302990339;
	//ROOT::RDataFrame y("Events", "RDF_TTToSemiLeptonic_seniorResearch_Muon__TTToSemiLeptonic05.root");
	ROOT::RDataFrame y("Events", "root://cmseos.fnal.gov//store/user/jmanagan/SingleProd_RDF_TTToSemiLeptonic_hadd.root");
	auto z = y.Filter("Bprime1_DeepAK8_Mass && isValidBDecay == true").Define("weights_TTToSemiLeptonic","Generator_weight*137000*831.76/(abs(Generator_weight)*302990339)");
	auto aa = z.Filter("taggedWbjetJet == true");
	auto bb = z.Filter("taggedTjet == true");
	auto cc = z.Filter("taggedWjet == true");
	
	auto TTToSemiLeptonic_Wb = aa.Histo1D({"TTToSemiLeptonic_Wb",";Mass (GeV);Events / bin", 25, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_TTToSemiLeptonic");
	TTToSemiLeptonic_Wb->Write();
	auto TTToSemiLeptonic_W = bb.Histo1D({"TTToSemiLeptonic_W",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_TTToSemiLeptonic");
	TTToSemiLeptonic_W->Write();
	auto TTToSemiLeptonic_T = cc.Histo1D({"TTToSemiLeptonic_T",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_TTToSemiLeptonic");
	TTToSemiLeptonic_T->Write();
	// ------------------------------------------------------------------------------------------------

	// ------------------------------------------------------------------------------------------------	
	// WToLNu MadGraoh Background	
	// ------------------------------------------------------------------------------------------------
	//int eventsNumber_WJetsToLNu = 4157920;
	int eventsNumber_WJetsToLNu = 47876094;
	//ROOT::RDataFrame gg("Events", "RDF_MadgraphBkg_seniorResearch_Muon_WJetsToLNu04.root");
	ROOT::RDataFrame gg("Events", "root://cmseos.fnal.gov//store/user/jmanagan/SingleProd_RDF_WJetsToLNu_hadd.root");
	auto hh = gg.Filter("Bprime1_DeepAK8_Mass && isValidBDecay == true").Define("weights_WJetsToLNu","Generator_weight*137000*61526.7/(abs(Generator_weight)*47876094)");
	auto ii = hh.Filter("taggedWbjetJet == true");
	auto jj = hh.Filter("taggedTjet == true");
	auto kk = hh.Filter("taggedWjet == true");
	
	auto WJetsToLNu_Wb = ii.Histo1D({"WJetsToLNu_Wb",";Mass (GeV);Events / bin", 25, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_WJetsToLNu");
	WJetsToLNu_Wb->Write();
	auto WJetsToLNu_W = jj.Histo1D({"WJetsToLNu_W",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_WJetsToLNu");
	WJetsToLNu_W->Write();
	auto WJetsToLNu_T = kk.Histo1D({"WJetsToLNu_T",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_WJetsToLNu");
	WJetsToLNu_T->Write();

	auto data_obs_Wb = ii.Histo1D({"data_obs_Wb",";Mass (GeV);Events / bin", 25, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_WJetsToLNu");
	data_obs_Wb->Write();
	auto data_obs_W = jj.Histo1D({"data_obs_W",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_WJetsToLNu");
	data_obs_W->Write();
	auto data_obs_T = kk.Histo1D({"data_obs_T",";Mass (GeV);Events / bin", 50, 0.25, 3000}, "Bprime1_DeepAK8_Mass","weights_WJetsToLNu");
	data_obs_T->Write();

	outputFile->Close();
	// ------------------------------------------------------------------------------------------------
};
