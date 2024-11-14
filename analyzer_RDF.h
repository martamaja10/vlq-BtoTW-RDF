//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 16 20:39:17 2015 by ROOT version 6.02/05
// from TTree ljmet/ljmet
// found on file: /eos/uscms/store/user/lpcljm/LJMet_1lep_110915/X53X53_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns/X53X53_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns_1.root
//////////////////////////////////////////////////////////

#ifndef rdf_h
#define rdf_h
#pragma once

#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1.h"
#include "TF1.h"

// Header file for the classes stored in the TTree if any.
#include "vector"
#include <string>
#include "TLorentzVector.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "/Users/martaczurylo/miniconda3/lib/python3.10/site-packages/correctionlib/include/correction.h"
R__LOAD_LIBRARY(/Users/martaczurylo/miniconda3/lib/python3.10/site-packages/correctionlib/lib/libcorrectionlib.dylib);


enum shift : char;

using namespace std;
using namespace ROOT::VecOps;
using correction::CorrectionSet;

class rdf
{
private:
  TTree *inputTree; //! pointer to the analyzed TTree or TChain
  TFile *inputFile;
  TString psOutName, fsOutName;
  Int_t fCurrent; //! current Tree number in a TChain
  vector<string> files;

  // Polynominals for LHE-based WJets HT scaling
  TF1 *poly2 = new TF1("poly2","max(0.402806,0.998174 + (8.40861e-05)*x + (-6.63274e-07)*x*x + (4.09272e-10)*x*x*x + (-9.50233e-14)*x*x*x*x + (7.59648e-18)*x*x*x*x*x)",100,5000);  
  TF1 *poly2U = new TF1("poly2U","max([6],[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)",100,5000);   
  TF1 *poly2D = new TF1("poly2D","max([6],[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)",100,5000);
  // Polynomials for data-derived TOP HT scaling
  TF1 *polyHT = new TF1("polyHT","min(1.0,max([0] + [1]*x,[2]))",700,5000);
  TF1 *polyHTU = new TF1("polyHTU","min(1.0,max([0] + [1]*x + sqrt([3] + 2*x*[4] + x*x*[5]),[2] + [6]))",700,5000);
  TF1 *polyHTD = new TF1("polyHTD","min(1.0,max([0] + [1]*x - sqrt([3] + 2*x*[4] + x*x*[5]),[2] - [6]))",700,5000);

  std::vector<std::vector<float>> elIDSF;
  std::vector<std::vector<float>> elIDSFUnc;  
  std::vector<std::vector<float>> elISOSF;
  float elISOSFUnc;  
  std::vector<std::vector<float>> muISOSF;
  std::vector<std::vector<float>> elHLTSF;
  std::vector<std::vector<float>> elHLTSFUnc;
  std::vector<std::vector<float>> muHLTSF;
  std::vector<std::vector<float>> muHLTSFUnc;

  std::vector<float> btagpts;
  std::vector<std::vector<float>> btageffs;

  std::vector<float> pnetpts;
  std::vector<std::vector<float>> pnet_t_t;
  std::vector<std::vector<float>> pnet_t_W;
  std::vector<std::vector<float>> pnet_W_t;
  std::vector<std::vector<float>> pnet_W_W;
  std::vector<std::vector<float>> pnet_J_t;
  std::vector<std::vector<float>> pnet_J_W;

  Bool_t isSig;
  Bool_t isTOP;
  Bool_t isMadgraphBkg; // W, Z, QCD
  Bool_t isMC;
  Bool_t isSM;
  Bool_t isSE;
  Bool_t isTT;
  Bool_t isVV;
  Bool_t isTTincMtt0to1000;
  Bool_t isTTincMtt0to700;
  Bool_t isTTincMtt700to1000;
  Bool_t isTTincMtt1000toInf;
  int    samplebin;
  string sample;
  string year;
  string era;
  string jecera;

  // Fixed size dimensions of array or collections stored in the TTree if any.
public:
  // Main Methods
  rdf(string inputFileName, string testNum1, string testNum2, string yearIn);
  RVec<RVec<float>> cleanJets(RVec<TLorentzVector> &jt_p4, RVec<float> &jt_rf, RVec<TLorentzVector> &mu_p4, RVec<int> mu_jetid, RVec<TLorentzVector> &el_p4, RVec<int> &el_jetid);
  virtual ~rdf();
  virtual void analyzer_RDF(TString testNum, TString jesvar);
};

#endif

#ifdef rdf_cxx

rdf::rdf(string inputFileName, string testNum1, string testNum2, string yearIn) : inputTree(0), inputFile(0)
{
  // Make the vector with the text file listed in the input
  cout << "Input File Path: " << inputFileName << endl;
  ifstream listFiles;
  listFiles.open(inputFileName);

  string file = "";
  int i = 0;
  int start = atoi(testNum1.c_str());
  int end = atoi(testNum2.c_str());
  cout << "TestNum 1: " << start << " and TestNum 2: " << end << endl;

  if (listFiles.is_open())
  {
    while (listFiles >> file)
    {
      if (i >= start && i <= end) {
        files.push_back(file);
        cout << "files: " << file << endl;

      }
      i++;
    }
  }
  cout << "Number of Entries: " << files.size() << endl;

  TString sampleName = file;
  cout << "Sample Name: " << sampleName << endl;

  // Parse the incoming file names to assign labels  
  isSig = (sampleName.Contains("Bprime"));
  isMadgraphBkg = (sampleName.Contains("QCD") || sampleName.Contains("madgraphMLM"));
  isTOP = (sampleName.Contains("Mtt") || sampleName.Contains("ST") || sampleName.Contains("ttZ") || sampleName.Contains("ttW") || sampleName.Contains("ttH") || sampleName.Contains("TTTo"));
  isTT = (sampleName.Contains("TT_Tune") || sampleName.Contains("Mtt") || sampleName.Contains("TTTo"));
  isVV = (sampleName.Contains("WW_") || sampleName.Contains("WZ_") || sampleName.Contains("ZZ_"));
  isMC = !(sampleName.Contains("Single") || sampleName.Contains("Data18") || sampleName.Contains("EGamma"));
  isSM = sampleName.Contains("SingleMuon");
  isSE = (sampleName.Contains("SingleElectron") || sampleName.Contains("EGamma"));

  if(isTT) samplebin = 0;
  else if(sampleName.Contains("ST_")) samplebin = 1;
  else if(sampleName.Contains("TTW") || sampleName.Contains("TTZ") || sampleName.Contains("ttH")) samplebin = 2;
  else if(sampleName.Contains("WJetsToLNu")) samplebin = 3;
  else if(sampleName.Contains("DYJets")) samplebin = 4;
  else if(isVV) samplebin = 5;
  else if(sampleName.Contains("QCD")) samplebin = 6;
  else if(isSig){
    if(sampleName.Contains("M-800")) samplebin = 7;
    if(sampleName.Contains("M-1000")) samplebin = 8;
    if(sampleName.Contains("M-1200")) samplebin = 9;
    if(sampleName.Contains("M-1300")) samplebin = 10;
    if(sampleName.Contains("M-1400")) samplebin = 11;
    if(sampleName.Contains("M-1500")) samplebin = 12;
    if(sampleName.Contains("M-1600")) samplebin = 13;
    if(sampleName.Contains("M-1700")) samplebin = 14;
    if(sampleName.Contains("M-1800")) samplebin = 15;
    if(sampleName.Contains("M-2000")) samplebin = 16;
    if(sampleName.Contains("M-2200")) samplebin = 17;
  }

  year = yearIn; // May need to change this line to get things to work

  TObjArray *tokens = sampleName.Tokenize("/");
  sample = ((TObjString *)(tokens->At(5)))->String();
  if(!isMC){ 
    string runera = (((TObjString *)(tokens->At(4)))->String()).Data();
    string process = (((TObjString *)(tokens->At(7)))->String()).Data();
    era = runera.back();
    if(year == "2016APV" && era == "B" && process.find("ver1") != std::string::npos) era = "A";
  } 
  delete tokens;

  jecera = era;
  if(year == "2016" && !isMC) jecera = "FGH";
  if(year == "2016APV" && !isMC){
    if(era == "A" or era == "B" or era == "C" or era == "D") jecera = "BCD";
    else jecera = "EF";
  }

  if(isMC){
    if(sampleName.Contains("_ext1")) era = "ext1";
    if(sampleName.Contains("_ext2")) era = "ext2";
    if(sampleName.Contains("_ext3")) era = "ext3";
  }

  btagpts = {15.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1200.0, 99999.0};

  if(year == "2016APV"){

    btageffs = { // B, C, Light
      {0.829835, 0.534866, 0.510229},
      {0.852295, 0.487347, 0.334173},
      {0.872671, 0.440800, 0.186784},
      {0.881861, 0.411432, 0.131229},
      {0.892275, 0.405335, 0.119147},
      {0.901108, 0.416303, 0.118356},
      {0.911451, 0.447708, 0.129164},
      {0.921408, 0.488303, 0.153077},
      {0.929370, 0.587864, 0.195291},
      {0.933337, 0.639332, 0.265921},
      {0.933098, 0.630793, 0.306818},
      {0.932616, 0.630037, 0.332766},
      {0.936644, 0.667409, 0.424512},
      {0.931860, 0.734737, 0.535264},
      {0.948849, 0.788618, 0.625000},
    };

    elIDSF = {  // AN-2019-199 v11 Fig 6(7,8,9 for other eras) (SUSY Double Disco)
      {1.005, 1.001, 0.988, 0.988, 0.968, 0.961, 0.977, 0.977, 0.980, 0.999},
      {1.024, 1.018, 1.082, 0.993, 0.977, 0.987, 0.975, 1.103, 1.001, 1.005},
      {0.941, 1.028, 0.976, 1.007, 0.988, 0.978, 1.012, 1.084, 0.959, 1.032}
    };
    elIDSFUnc = {
      {0.036, 0.017, 0.041, 0.026, 0.019, 0.019, 0.026, 0.041, 0.017, 0.036},
      {0.016, 0.013, 0.072, 0.011, 0.007, 0.007, 0.011, 0.072, 0.012, 0.016},
      {0.064, 0.033, 0.122, 0.028, 0.019, 0.019, 0.028, 0.119, 0.033, 0.062},
    };
    elISOSF = {  // same AN, Fig 10(11,12,13 for other eras)
      {1.001, 1.000, 1.004, 1.001, 1.000, 1.000, 0.999, 0.999, 1.001, 1.002},
      {0.999, 0.999, 0.999, 1.000, 1.000, 1.000, 1.001, 0.999, 1.000, 1.002},
      {0.994, 0.998, 0.988, 1.001, 1.000, 1.000, 1.001, 0.992, 0.998, 0.995}
    };
    elISOSFUnc = 0.002; //largest of any bin, not worth a histogram
    muISOSF = { // same AN, Fig 14(15,16,17 for other eras)
      {0.999, 1.000, 1.000, 1.000},
      {1.000, 1.000, 0.999, 1.001}
    };
    muHLTSF = { // rows are etabin, columns ptbin
      {0.97774, 0.9726, 0.96488, 0.9417, 0.87642},
      {0.97579, 0.95524, 0.96572, 0.96525, 0.75831},
      {0.96482, 0.98496, 0.95403, 0.98674, 0.96618},
      {0.96386, 0.97067, 0.98424, 0.76465, 1.19231}
    };
    muHLTSFUnc = {
      {0.00858, 0.00354, 0.00686, 0.01936, 0.04207},
      {0.01471, 0.00748, 0.01295, 0.02976, 0.10028},
      {0.0165, 0.00618, 0.01334, 0.0294, 0.06186},
      {0.04455, 0.0211, 0.04766, 0.13205, 0.09391}
    };
    elHLTSF = {
	       {1.06797, 0.85151, 1.04762},
	       {1.06713, 0.97127, 1.01831},
	       {1.06616, 1.0303, 1.08065},
	       {0.9936, 1.01903, 0.92554},
	       {1.00941, 0.98559, 0.99891},
	       {1.01336, 0.98826, 1.04159},
	       {0.9998, 1.02166, 0.94219},
	       {0.99527, 0.74098, 1.16667},
	       {1.00724, 1.00424, 1.09774},
	       {0.96812, 1.02572, 1.11392},
    };
    elHLTSFUnc = {
		  {0.03836, 0.12645, 0.02437},
		  {0.02399, 0.05044, 0.07431},
		  {0.04751, 0.07618, 0.02651},
		  {0.01814, 0.02484, 0.04601},
		  {0.01313, 0.01998, 0.03108},
		  {0.01321, 0.0198, 0.02091},
		  {0.01696, 0.02353, 0.04886},
		  {0.05214, 0.09298, 0.04366},
		  {0.0288, 0.04914, 0.02008},
		  {0.04551, 0.06189, 0.04008},
    };
  }else if(year == "2016"){

    btageffs = { // B, C, Light
      {0.849801, 0.565951, 0.519292},
      {0.871238, 0.515046, 0.347213},
      {0.884010, 0.461947, 0.199234},
      {0.895665, 0.436844, 0.140126},
      {0.902270, 0.425167, 0.125008},
      {0.910727, 0.441395, 0.123484},
      {0.922480, 0.454206, 0.132626},
      {0.928428, 0.506281, 0.155521},
      {0.936479, 0.586857, 0.194188},
      {0.938037, 0.644332, 0.258005},
      {0.936037, 0.638778, 0.297550},
      {0.938835, 0.636197, 0.319172},
      {0.939513, 0.661658, 0.393939},
      {0.931919, 0.708171, 0.461538},
      {0.940379, 0.786585, 0.595000},
    };

    elIDSF = {
      {0.999, 1.001, 0.982, 0.985, 0.978, 0.994, 0.988, 0.970, 0.985, 0.991},
      {1.039, 1.026, 1.060, 1.017, 0.998, 1.009, 0.995, 1.083, 0.994, 1.012},
      {1.037, 1.012, 0.888, 0.997, 1.010, 1.016, 0.992, 0.922, 1.047, 0.935}
    };
    elIDSFUnc = {
      {0.037, 0.026, 0.046, 0.026, 0.016, 0.016, 0.026, 0.046, 0.026, 0.037},
      {0.022, 0.011, 0.072, 0.010, 0.006, 0.006, 0.010, 0.072, 0.011, 0.022},
      {0.073, 0.029, 0.110, 0.022, 0.016, 0.016, 0.022, 0.108, 0.031, 0.072},
    };
    elISOSF = {
      {1.001, 1.000, 1.001, 0.999, 1.000, 0.999, 0.999, 1.001, 1.000, 1.001},
      {0.999, 1.001, 1.004, 1.000, 1.001, 1.000, 1.000, 1.002, 1.000, 1.000},
      {0.993, 0.998, 0.989, 0.999, 0.999, 1.001, 0.999, 0.989, 0.998, 0.994}
    };
    elISOSFUnc = 0.002;
    muISOSF = {
      {1.000, 0.999, 0.999, 1.001},
      {1.000, 1.001, 1.000, 1.002}
    };
    muHLTSF = {
      {0.98362, 0.98623, 0.99124, 0.95174, 0.99199},
      {0.98779, 0.97787, 0.96513, 0.93283, 0.91462},
      {1.01839, 1.0065, 0.99663, 0.94006, 0.9382},
      {1.00977, 1.0079, 1.02854, 1.08213, 0.71795}
    };
    muHLTSFUnc = {
      {0.00815, 0.00328, 0.00586, 0.01699, 0.02888},
      {0.01371, 0.0062, 0.0124, 0.03308, 0.06896},
      {0.01255, 0.00556, 0.01195, 0.0324, 0.09489},
      {0.05077, 0.01773, 0.03382, 0.02072, 0.18982},
    };
    elHLTSF = {
	{1.03014, 0.83014, 1.06122},
	{1.02285, 0.9592, 0.94865},
	{1.08569, 1.12368, 0.81944},
	{1.00538, 0.96318, 0.95721},
	{1.00818, 0.9856, 0.98429},
	{0.97618, 0.96249, 1.02682},
	{0.98282, 0.97902, 0.88658},
	{0.94034, 0.99777, 0.80556},
	{1.03174, 1.03791, 0.88994},
	{1.02163, 1.02687, 1.12222},
    };
    elHLTSFUnc = {
	{0.03885, 0.11161, 0.02575},
	{0.02963, 0.05519, 0.08352},
	{0.04337, 0.01912, 0.18768},
	{0.01819, 0.02902, 0.04385},
	{0.01345, 0.02124, 0.03216},
	{0.01389, 0.02354, 0.02766},
	{0.01966, 0.02683, 0.05123},
	{0.05395, 0.09239, 0.18087},
	{0.02901, 0.04605, 0.14282},
	{0.04571, 0.06245, 0.03904},
    };
  }else if(year == "2017"){

    btageffs = {
      {0.859705, 0.555959, 0.492094},
      {0.888146, 0.536376, 0.326982},
      {0.910808, 0.502367, 0.175392},
      {0.921982, 0.478112, 0.121900},
      {0.926561, 0.465780, 0.111225},
      {0.933731, 0.470531, 0.110502},
      {0.940000, 0.481269, 0.118697},
      {0.946382, 0.526891, 0.139537},
      {0.950598, 0.605253, 0.175671},
      {0.950137, 0.645084, 0.227488},
      {0.948931, 0.635245, 0.250078},
      {0.947138, 0.632706, 0.264945},
      {0.949487, 0.662092, 0.336838},
      {0.938316, 0.704684, 0.421262},
      {0.946124, 0.732601, 0.516129}
    };

    elIDSF = {
      {1.013, 0.994, 1.001, 0.983, 0.979, 0.980, 0.983, 0.967, 0.994, 1.023},
      {1.019, 1.001, 1.026, 1.011, 0.992, 0.992, 0.990, 0.990, 1.023, 1.028},
      {0.977, 0.972, 1.009, 0.996, 0.995, 1.008, 1.013, 1.087, 1.023, 1.035}
    };
    elIDSFUnc = {
      {0.033, 0.018, 0.047, 0.025, 0.018, 0.018, 0.025, 0.047, 0.018, 0.033},
      {0.016, 0.009, 0.036, 0.010, 0.006, 0.006, 0.009, 0.036, 0.009, 0.015},
      {0.055, 0.029, 0.099, 0.017, 0.016, 0.016, 0.017, 0.097, 0.029, 0.054}
    };
    elISOSF = {
      {1.000, 1.002, 1.004, 1.000, 1.000, 1.000, 1.000, 1.003, 1.001, 1.000},
      {1.002, 1.000, 1.000, 1.000, 1.001, 1.000, 1.001, 0.999, 1.001, 0.999},
      {0.998, 0.999, 0.996, 0.999, 1.003, 1.000, 1.000, 0.997, 0.999, 1.008}
    };
    elISOSFUnc = 0.008;
    muISOSF = {
      {1.000, 1.000, 1.000, 1.000},
      {1.000, 0.999, 1.000, 1.000}
    };
    muHLTSF = {
      {0.9801, 0.98764, 0.98809, 0.98698, 0.98259},
      {0.9757, 0.96532, 0.97838, 0.97019, 0.91306},
      {0.98607, 1.0102, 1.01146, 0.97565, 0.99721},
      {0.99818, 1.01533, 1.03488, 1.10558, 1.04225}
    };
    muHLTSFUnc = {
      {0.00567, 0.00212, 0.00382, 0.00868, 0.0197},
      {0.01082, 0.00466, 0.00732, 0.01755, 0.04674},
      {0.01071, 0.00346, 0.00561, 0.01727, 0.038},
      {0.03958, 0.01413, 0.0239, 0.04597, 0.02491}
    };
    elHLTSF = {
	       {0.98337, 1.01174, 0.88939},
	       {0.9556, 1.03358, 0.97555},
	       {0.92756, 1.07689, 0.84828},
	       {0.94214, 0.94608, 1.00729},
	       {0.9444, 0.95752, 0.98874},
	       {0.96579, 0.98604, 1.00515},
	       {0.96451, 1.01399, 0.98523},
	       {1.01327, 1.10539, 0.93941},
	       {0.93989, 0.94243, 0.96526},
	       {0.89329, 0.99857, 0.96711},
    };
    elHLTSFUnc = {
		  {0.03809, 0.07219, 0.10022},
		  {0.02614, 0.03252, 0.05911},
		  {0.04183, 0.04526, 0.15793},
		  {0.01487, 0.0215, 0.02741},
		  {0.01111, 0.01632, 0.0215},
		  {0.01076, 0.01535, 0.01959},
		  {0.01455, 0.02028, 0.02924},
		  {0.04134, 0.04604, 0.135},
		  {0.02548, 0.03985, 0.07436},
		  {0.04121, 0.06019, 0.10878},
    };
  }else{

    btageffs = {
      {0.869006, 0.569043, 0.495051},
      {0.892377, 0.544933, 0.325319},
      {0.911736, 0.507194, 0.178191},
      {0.919209, 0.482783, 0.123977},
      {0.925770, 0.467890, 0.110386},
      {0.932663, 0.463212, 0.110406},
      {0.938477, 0.488811, 0.118632},
      {0.945820, 0.530118, 0.140454},
      {0.950766, 0.614258, 0.179932},
      {0.951791, 0.649368, 0.233308},
      {0.951093, 0.643851, 0.258580},
      {0.950231, 0.638892, 0.279099},
      {0.947389, 0.670874, 0.347167},
      {0.951818, 0.700409, 0.449732},
      {0.941429, 0.785455, 0.561611},
    };

    elIDSF = {
      {0.986, 0.980, 1.021, 0.964, 0.973, 0.973, 0.973, 0.980, 0.976, 0.969},
      {0.984, 0.987, 1.083, 0.992, 0.977, 0.991, 1.001, 1.076, 0.991, 0.986},
      {1.014, 1.041, 0.902, 0.987, 0.997, 1.011, 1.010, 1.011, 1.022, 0.963}
    };
    elIDSFUnc = {
      {0.035, 0.021, 0.024, 0.028, 0.019, 0.019, 0.028, 0.024, 0.021, 0.035},
      {0.013, 0.009, 0.048, 0.011, 0.005, 0.005, 0.011, 0.048, 0.009, 0.013},
      {0.043, 0.029, 0.183, 0.022, 0.012, 0.012, 0.022, 0.183, 0.030, 0.042}
    };
    elISOSF = {
      {1.001, 1.000, 1.002, 0.999, 0.999, 0.999, 1.000, 1.001, 1.000, 1.001},
      {1.001, 1.000, 1.000, 1.001, 1.000, 1.000, 1.000, 0.998, 1.001, 0.999},
      {0.999, 1.000, 0.998, 1.001, 1.000, 1.000, 1.001, 0.998, 1.000, 0.999}
    };
    elISOSFUnc = 0.002;
    muISOSF = {
      {1.000, 1.000, 1.000, 1.000},
      {1.000, 0.999, 1.000, 1.001}
    };
    muHLTSF = {
      {0.98675, 0.988, 0.98389, 0.98667, 0.99726},
      {0.9679, 0.98152, 0.97593, 1.00438, 0.97792},
      {0.99798, 1.0017, 1.00678, 1.01222, 1.02847},
      {0.95344, 1.01615, 0.97478, 1.04948, 0.91837}
    };
    muHLTSFUnc = {
      {0.00416, 0.0017, 0.00327, 0.00716, 0.01507},
      {0.0087, 0.00315, 0.006, 0.01015, 0.03443},
      {0.00668, 0.00253, 0.00463, 0.01087, 0.02492},
      {0.02957, 0.01041, 0.02395, 0.03289, 0.15075}
    };
    elHLTSF = {
	       {0.92956, 0.99004, 0.91985},
	       {0.98449, 0.92816, 0.99065},
	       {0.93523, 1.01502, 0.98342},
	       {0.98793, 0.99445, 1.02878},
	       {0.99355, 1.00935, 0.99988},
	       {0.98114, 1.00203, 0.99266},
	       {0.98862, 1.00698, 0.96965},
	       {0.96667, 0.87563, 0.97996},
	       {0.94752, 0.97961, 1.05568},
	       {0.97815, 0.97344, 1.01281},
    };
    elHLTSFUnc = {
		  {0.02893, 0.05119, 0.09858},
		  {0.01834, 0.03043, 0.05082},
		  {0.03328, 0.03938, 0.08435},
		  {0.01069, 0.01643, 0.0196},
		  {0.00789, 0.01093, 0.01485},
		  {0.00816, 0.01195, 0.01597},
		  {0.01064, 0.01574, 0.02273},
		  {0.03133, 0.05214, 0.07876},
		  {0.01918, 0.02997, 0.03604},
		  {0.02921, 0.05045, 0.06342},
    };
  }    

  pnetpts = {200,300,400,500,600,800,1000,1200,1500,2000};

  pnet_t_t = {
    {0.05033, 0.13082, 0.43044, 0.72995, 0.8321, 0.87784, 0.90759, 0.92422, 0.92662, 0.93072}, // ttbar
    {0.04751, 0.11688, 0.44316, 0.73689, 0.8329, 0.86828, 0.89489, 0.90323, 0.86555, 0.76667}, // single-t
    {0.06337, 0.14203, 0.44779, 0.73293, 0.8255, 0.87586, 0.89527, 0.91557, 0.93516, 0.93704}, // ttbarVH
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // Wjets
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // Zjets
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // VV
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // qcd
    {0.04378, 0.13127, 0.43255, 0.75929, 0.84744, 0.8879, 0.87632, 0.89286, 0.94118, 1.0}, // BprimeM800
    {0.08304, 0.14514, 0.43504, 0.72683, 0.84946, 0.87801, 0.89231, 0.91089, 0.93478, 0.83333}, // BprimeM1000
    {0.06522, 0.14837, 0.43998, 0.709, 0.82532, 0.88602, 0.8892, 0.87708, 0.93103, 0.91667}, // BprimeM1200
    {0.04094, 0.14751, 0.43537, 0.70658, 0.8165, 0.87963, 0.89155, 0.88071, 0.88987, 1.0}, // BprimeM1300
    {0.11024, 0.15295, 0.43624, 0.70662, 0.81354, 0.8794, 0.89127, 0.88397, 0.87432, 0.90909}, // BprimeM1400
    {0.0873, 0.1572, 0.43706, 0.7046, 0.81041, 0.87904, 0.9057, 0.89016, 0.87819, 0.81356}, // BprimeM1500
    {0.07759, 0.15661, 0.43508, 0.70376, 0.81355, 0.87821, 0.91069, 0.89035, 0.89715, 0.82796}, // BprimeM1600
    {0.04505, 0.15385, 0.43288, 0.70283, 0.81013, 0.87571, 0.91015, 0.90497, 0.89469, 0.87013}, // BprimeM1700
    {0.10619, 0.1583, 0.42802, 0.69989, 0.80799, 0.87374, 0.91031, 0.91453, 0.89697, 0.85281}, // BprimeM1800
    {0.05814, 0.15448, 0.42128, 0.69091, 0.80437, 0.8715, 0.90861, 0.92511, 0.9067, 0.90244}, // BprimeM2000
    {0.08, 0.15698, 0.42509, 0.69019, 0.79967, 0.86946, 0.90622, 0.92562, 0.92441, 0.90094}, // BprimeM2200
  };
  pnet_t_W = {
    {0.29069, 0.51768, 0.5187, 0.46474, 0.42692, 0.37788, 0.31712, 0.25896, 0.21443, 0.14689}, // ttbar
    {0.23175, 0.47566, 0.46442, 0.41178, 0.36661, 0.30275, 0.23379, 0.19713, 0.14286, 0.16667}, // single-t
    {0.28523, 0.48612, 0.48541, 0.42735, 0.38719, 0.33112, 0.26725, 0.21491, 0.15944, 0.13177}, // ttbarVH
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // Wjets
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // Zjets
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // VV
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // qcd
    {0.21891, 0.41622, 0.44196, 0.39807, 0.3539, 0.31121, 0.25, 0.10714, 0.29412, 0.0}, // BprimeM800
    {0.22491, 0.40149, 0.42781, 0.38597, 0.35577, 0.30183, 0.23504, 0.20792, 0.08696, 0.16667}, // BprimeM1000
    {0.19565, 0.39722, 0.41752, 0.37715, 0.34558, 0.31434, 0.23699, 0.22259, 0.14655, 0.16667}, // BprimeM1200
    {0.22222, 0.39226, 0.41893, 0.37342, 0.33715, 0.31306, 0.24961, 0.18465, 0.15419, 0.05}, // BprimeM1300
    {0.23622, 0.39608, 0.41684, 0.37087, 0.33372, 0.30464, 0.24917, 0.17628, 0.16667, 0.15909}, // BprimeM1400
    {0.16667, 0.39519, 0.40986, 0.37232, 0.3295, 0.29447, 0.25842, 0.18673, 0.14735, 0.11864}, // BprimeM1500
    {0.21552, 0.38077, 0.40486, 0.36882, 0.33152, 0.28696, 0.26067, 0.18069, 0.17224, 0.12903}, // BprimeM1600
    {0.15315, 0.38418, 0.41061, 0.36584, 0.32839, 0.27839, 0.25542, 0.19394, 0.14791, 0.11039}, // BprimeM1700
    {0.23009, 0.38182, 0.40079, 0.36623, 0.32949, 0.27454, 0.24811, 0.20147, 0.16318, 0.07359}, // BprimeM1800
    {0.32558, 0.35089, 0.39574, 0.35498, 0.32403, 0.27127, 0.2275, 0.21216, 0.16031, 0.1182}, // BprimeM2000
    {0.2, 0.34001, 0.40134, 0.35966, 0.31822, 0.27031, 0.2171, 0.19802, 0.16824, 0.09151}, // BprimeM2200
  };
  pnet_W_t = {
    {0.01277, 0.02969, 0.05524, 0.08818, 0.11627, 0.14131, 0.15569, 0.17035, 0.1814, 0.13889}, // ttbar
    {0.011, 0.03063, 0.05206, 0.08047, 0.09679, 0.10783, 0.10853, 0.12461, 0.07194, 0.13514}, // single-t
    {0.02553, 0.06311, 0.15874, 0.25389, 0.31448, 0.34642, 0.37025, 0.39559, 0.41667, 0.33333}, // ttbarVH
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // Wjets
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // Zjets
    {0.00964, 0.02071, 0.02326, 0.03458, 0.03421, 0.04085, 0.0437, 0.03356, 0.04225, 0.15385}, // VV
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // qcd
    {0.02469, 0.03117, 0.0407, 0.08473, 0.18171, 0.17982, 0.12987, 0.15584, 0.15789, 0.25}, // BprimeM800
    {0.01288, 0.0324, 0.03924, 0.04824, 0.08667, 0.16544, 0.165, 0.18889, 0.2, 0.0}, // BprimeM1000
    {0.01258, 0.03488, 0.04099, 0.04792, 0.05362, 0.09544, 0.17018, 0.14863, 0.13072, 0.3125}, // BprimeM1200
    {0.02985, 0.03457, 0.04274, 0.04695, 0.05217, 0.07502, 0.15969, 0.12362, 0.15668, 0.03571}, // BprimeM1300
    {0.0354, 0.03391, 0.04401, 0.04645, 0.051, 0.0662, 0.14788, 0.15303, 0.12973, 0.16279}, // BprimeM1400
    {0.02632, 0.03445, 0.04214, 0.04687, 0.05266, 0.06084, 0.1108, 0.16254, 0.13811, 0.09589}, // BprimeM1500
    {0.02632, 0.03303, 0.0428, 0.05064, 0.05182, 0.05825, 0.09118, 0.15584, 0.14455, 0.13281}, // BprimeM1600
    {0.0, 0.03384, 0.04407, 0.04951, 0.052, 0.057, 0.07371, 0.14328, 0.15856, 0.13298}, // BprimeM1700
    {0.0, 0.03487, 0.04464, 0.04867, 0.05257, 0.05653, 0.06717, 0.12802, 0.14691, 0.15909}, // BprimeM1800
    {0.02083, 0.03502, 0.04402, 0.05005, 0.05293, 0.05539, 0.05805, 0.08426, 0.13256, 0.12174}, // BprimeM2000
    {0.06897, 0.04233, 0.04892, 0.05077, 0.05882, 0.05676, 0.05693, 0.06809, 0.12259, 0.13186}, // BprimeM2200
  };
  pnet_W_W = {
    {0.43478, 0.7508, 0.81291, 0.78317, 0.75379, 0.70002, 0.59875, 0.5163, 0.33953, 0.41667}, // ttbar
    {0.41112, 0.75028, 0.83064, 0.82115, 0.81658, 0.81613, 0.8154, 0.84112, 0.8705, 0.89189}, // single-t
    {0.43062, 0.69105, 0.72211, 0.63681, 0.55644, 0.46717, 0.38036, 0.32206, 0.28431, 0.35294}, // ttbarVH
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // Wjets
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // Zjets
    {0.37952, 0.69324, 0.77458, 0.76625, 0.77577, 0.75156, 0.72773, 0.81879, 0.73239, 0.69231}, // VV
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // qcd
    {0.4716, 0.78158, 0.86932, 0.86233, 0.773, 0.71211, 0.67532, 0.7013, 0.63158, 0.75}, // BprimeM800
    {0.45064, 0.76248, 0.85794, 0.88362, 0.87777, 0.79247, 0.73017, 0.72778, 0.64444, 1.0}, // BprimeM1000
    {0.42767, 0.74288, 0.8398, 0.87113, 0.89505, 0.88171, 0.77164, 0.73829, 0.75817, 0.6875}, // BprimeM1200
    {0.37313, 0.72515, 0.83273, 0.86382, 0.89052, 0.89614, 0.80392, 0.75386, 0.78341, 0.78571}, // BprimeM1300
    {0.32743, 0.72169, 0.8227, 0.85863, 0.88377, 0.89937, 0.83486, 0.76048, 0.75405, 0.90698}, // BprimeM1400
    {0.35965, 0.71237, 0.81792, 0.85142, 0.87514, 0.90157, 0.87549, 0.78157, 0.75699, 0.72603}, // BprimeM1500
    {0.38158, 0.70743, 0.81154, 0.84476, 0.87065, 0.90034, 0.89655, 0.81403, 0.77133, 0.75781}, // BprimeM1600
    {0.28814, 0.69758, 0.80981, 0.83632, 0.86422, 0.89542, 0.90689, 0.84266, 0.79121, 0.73936}, // BprimeM1700
    {0.375, 0.68641, 0.7965, 0.83113, 0.85905, 0.89031, 0.91022, 0.86583, 0.78926, 0.82576}, // BprimeM1800
    {0.27083, 0.66367, 0.78858, 0.82644, 0.84983, 0.87903, 0.90961, 0.90794, 0.83934, 0.77217}, // BprimeM2000
    {0.34483, 0.65472, 0.77266, 0.8081, 0.83773, 0.87087, 0.90152, 0.91597, 0.8797, 0.79692}, // BprimeM2200
  };
  pnet_J_t = {
    {0.01688, 0.05067, 0.15491, 0.17732, 0.14971, 0.10297, 0.06159, 0.04593, 0.04153, 0.03788}, // ttbar
    {0.00462, 0.0128, 0.03845, 0.05412, 0.0573, 0.05657, 0.06196, 0.06284, 0.05905, 0.0963}, // single-t
    {0.0233, 0.0622, 0.1792, 0.20714, 0.17736, 0.13668, 0.09021, 0.06942, 0.05114, 0.04344}, // ttbarVH
    {0.00048, 0.00149, 0.00646, 0.01326, 0.01806, 0.02213, 0.02689, 0.03467, 0.04218, 0.04925}, // Wjets
    {0.00143, 0.00268, 0.00859, 0.0158, 0.02109, 0.02502, 0.02974, 0.03756, 0.04598, 0.05098}, // Zjets
    {0.00043, 0.00315, 0.00794, 0.01153, 0.01459, 0.01581, 0.01736, 0.03143, 0.0, 0.0}, // VV
    {0.00072, 0.00234, 0.00777, 0.01463, 0.02018, 0.02537, 0.02847, 0.03259, 0.03713, 0.04664}, // qcd
    {0.01518, 0.05193, 0.20695, 0.25971, 0.17495, 0.10512, 0.09286, 0.08333, 0.06618, 0.02857}, // BprimeM800
    {0.00834, 0.03562, 0.16382, 0.29473, 0.30285, 0.17755, 0.10894, 0.07916, 0.107, 0.08}, // BprimeM1000
    {0.01159, 0.02679, 0.12065, 0.24379, 0.33191, 0.30242, 0.11403, 0.09393, 0.10024, 0.04396}, // BprimeM1200
    {0.00834, 0.02415, 0.10782, 0.21511, 0.31623, 0.34812, 0.16031, 0.10033, 0.09639, 0.11321}, // BprimeM1300
    {0.00181, 0.02284, 0.09655, 0.19457, 0.28448, 0.36597, 0.20574, 0.10976, 0.10282, 0.09722}, // BprimeM1400
    {0.01016, 0.0198, 0.08669, 0.1731, 0.25669, 0.36152, 0.26962, 0.12182, 0.11127, 0.10383}, // BprimeM1500
    {0.00286, 0.01855, 0.0797, 0.16125, 0.23242, 0.34149, 0.32516, 0.13859, 0.1075, 0.1327}, // BprimeM1600
    {0.00969, 0.01689, 0.07385, 0.14659, 0.21566, 0.31865, 0.37874, 0.17483, 0.12526, 0.09312}, // BprimeM1700
    {0.0092, 0.01628, 0.06626, 0.13538, 0.19838, 0.28993, 0.391, 0.23846, 0.12441, 0.08602}, // BprimeM1800
    {0.00661, 0.01448, 0.05979, 0.11784, 0.17219, 0.24855, 0.37473, 0.3698, 0.17786, 0.10209}, // BprimeM2000
    {0.00794, 0.01342, 0.05283, 0.10619, 0.15294, 0.2148, 0.32735, 0.41912, 0.25192, 0.13276}, // BprimeM2200
  };
  pnet_J_W = {
    {0.10397, 0.15608, 0.1428, 0.09795, 0.07299, 0.05157, 0.03589, 0.02708, 0.02464, 0.02125}, // ttbar
    {0.05112, 0.06807, 0.06423, 0.05089, 0.04219, 0.03523, 0.02935, 0.03791, 0.02904, 0.00741}, // single-t
    {0.12954, 0.17277, 0.14507, 0.09496, 0.07001, 0.05221, 0.03765, 0.03359, 0.02714, 0.01998}, // ttbarVH
    {0.05828, 0.07422, 0.07245, 0.05858, 0.04886, 0.04028, 0.03224, 0.02761, 0.02552, 0.02135}, // Wjets
    {0.07574, 0.08568, 0.08106, 0.06481, 0.05317, 0.04346, 0.03423, 0.02845, 0.02653, 0.022}, // Zjets
    {0.05946, 0.11366, 0.14796, 0.12911, 0.10897, 0.09236, 0.08679, 0.06483, 0.06075, 0.06349}, // VV
    {0.05896, 0.07495, 0.06826, 0.05514, 0.0462, 0.03485, 0.02915, 0.02528, 0.02377, 0.0177}, // qcd
    {0.07063, 0.1266, 0.15031, 0.10409, 0.06125, 0.0389, 0.03429, 0.04167, 0.03676, 0.05714}, // BprimeM800
    {0.06589, 0.10312, 0.13669, 0.128, 0.09838, 0.05523, 0.03168, 0.02908, 0.03704, 0.0}, // BprimeM1000
    {0.0713, 0.0906, 0.10913, 0.11792, 0.12332, 0.0873, 0.03758, 0.03229, 0.02934, 0.01099}, // BprimeM1200
    {0.05746, 0.0847, 0.10174, 0.11057, 0.11962, 0.09656, 0.04054, 0.03915, 0.02008, 0.07547}, // BprimeM1300
    {0.05144, 0.08262, 0.09385, 0.09945, 0.11272, 0.10485, 0.05397, 0.033, 0.03151, 0.02778}, // BprimeM1400
    {0.05817, 0.07871, 0.08837, 0.09357, 0.10211, 0.10691, 0.06732, 0.03192, 0.03424, 0.01093}, // BprimeM1500
    {0.05052, 0.07798, 0.08769, 0.08495, 0.09833, 0.10315, 0.08371, 0.04225, 0.0325, 0.02844}, // BprimeM1600
    {0.05717, 0.07616, 0.08446, 0.07946, 0.08755, 0.10318, 0.09551, 0.04563, 0.04071, 0.02429}, // BprimeM1700
    {0.0644, 0.07444, 0.07948, 0.0757, 0.08083, 0.09571, 0.10085, 0.05002, 0.03016, 0.03584}, // BprimeM1800
    {0.06232, 0.07421, 0.07341, 0.07089, 0.07741, 0.0846, 0.09518, 0.08345, 0.03359, 0.02784}, // BprimeM2000
    {0.0576, 0.07046, 0.07036, 0.06666, 0.06679, 0.07624, 0.08908, 0.09301, 0.05214, 0.02141}, // BprimeM2200
  };
  
  // Build polynomials for corrections
  poly2U->SetParameter(0,    0.998164);  
  poly2U->SetParameter(1, 8.51652e-05); 
  poly2U->SetParameter(2,-6.62919e-07);
  poly2U->SetParameter(3,  4.1205e-10); 
  poly2U->SetParameter(4,-9.57795e-14); 
  poly2U->SetParameter(5, 7.67371e-18); 
  poly2U->SetParameter(6,0.454456);
  poly2D->SetParameter(0,    0.998183);  
  poly2D->SetParameter(1, 8.30069e-05); 
  poly2D->SetParameter(2,-6.63629e-07);
  poly2D->SetParameter(3, 4.06495e-10); 
  poly2D->SetParameter(4,-9.42671e-14); 
  poly2D->SetParameter(5, 7.51924e-18); 
  poly2D->SetParameter(6,0.351156);
  polyHT->SetParameter(0,     1.09245);
  polyHT->SetParameter(1,-0.000220375);
  polyHT->SetParameter(2,    0.607311);
  polyHTU->SetParameter(0,     1.09245);
  polyHTU->SetParameter(1,-0.000220375);
  polyHTU->SetParameter(2,    0.607311);
  polyHTU->SetParameter(3, 0.000531602);
  polyHTU->SetParameter(4,-3.99715e-07);
  polyHTU->SetParameter(5, 3.19837e-10);
  polyHTU->SetParameter(6, 0.000541832);
  polyHTD->SetParameter(0,     1.09245);
  polyHTD->SetParameter(1,-0.000220375);
  polyHTD->SetParameter(2,    0.607311);
  polyHTD->SetParameter(3, 0.000531602);
  polyHTD->SetParameter(4,-3.99715e-07);
  polyHTD->SetParameter(5, 3.19837e-10);
  polyHTD->SetParameter(6, 0.000541832);


}

rdf::~rdf()
{
  if (!inputTree)
    return;
  delete inputTree->GetCurrentFile();
}


#endif // #ifdef rdf_cxx
