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
#include "../../correctionlib/include/correction.h"

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

  std::vector<float> btagpts;
  std::vector<std::vector<float>> btageffs;

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
  virtual void analyzer_trigEff(TString testNum, TString jesvar);
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

    elIDSF = {
      {1.005, 1.001, 0.988, 0.988, 0.968, 0.961, 0.977, 0.977, 0.980, 0.999},
      {1.024, 1.018, 1.082, 0.993, 0.977, 0.987, 0.975, 1.103, 1.001, 1.005},
      {0.941, 1.028, 0.976, 1.007, 0.988, 0.978, 1.012, 1.084, 0.959, 1.032}
    };
    elIDSFUnc = {
      {0.036, 0.017, 0.041, 0.026, 0.019, 0.019, 0.026, 0.041, 0.017, 0.036},
      {0.016, 0.013, 0.072, 0.011, 0.007, 0.007, 0.011, 0.072, 0.012, 0.016},
      {0.064, 0.033, 0.122, 0.028, 0.019, 0.019, 0.028, 0.119, 0.033, 0.062},
    };
    elISOSF = {
      {1.001, 1.000, 1.004, 1.001, 1.000, 1.000, 0.999, 0.999, 1.001, 1.002},
      {0.999, 0.999, 0.999, 1.000, 1.000, 1.000, 1.001, 0.999, 1.000, 1.002},
      {0.994, 0.998, 0.988, 1.001, 1.000, 1.000, 1.001, 0.992, 0.998, 0.995}
    };
    elISOSFUnc = 0.002; //largest of any bin, not worth a histogram
    muISOSF = {
      {0.999, 1.000, 1.000, 1.000},
      {1.000, 1.000, 0.999, 1.001}
    };
    elHLTSF = {
      {0.969994,0.999272,1.000000,0.955299,0.972568},
      {0.979406,0.976827,0.948699,0.925570,0.904762},
      {0.998922,0.948318,1.086539,0.935962,0.990274},
      {0.981869,0.979395,1.035172,0.991033,1.042024},
      {0.981854,0.999758,1.023505,1.006556,0.980000},
      {0.966784,0.990014,1.000000,1.029556,1.000000},
      {0.986579,1.000000,1.000000,1.000000,1.000000}
    };
    elHLTSFUnc = {
      {0.113757,0.104129,0.000000,0.083539,0.074347},
      {0.087648,0.097016,0.157331,0.128106,0.095238},
      {0.078416,0.075684,0.094028,0.118257,0.080618},
      {0.078566,0.042963,0.036409,0.034765,0.043790},
      {0.060464,0.025592,0.024057,0.006599,0.020000},
      {0.071746,0.016706,0.000000,0.030430,0.000000},
      {0.043074,0.000000,1.000000,0.000000,1.000000}
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
    elHLTSF = {
      {0.969994,0.999272,1.000000,0.955299,0.972568},
      {0.979406,0.976827,0.948699,0.925570,0.904762},
      {0.998922,0.948318,1.086539,0.935962,0.990274},
      {0.981869,0.979395,1.035172,0.991033,1.042024},
      {0.981854,0.999758,1.023505,1.006556,0.980000},
      {0.966784,0.990014,1.000000,1.029556,1.000000},
      {0.986579,1.000000,1.000000,1.000000,1.000000}
    };
    elHLTSFUnc = {
      {0.113757,0.104129,0.000000,0.083539,0.074347},
      {0.087648,0.097016,0.157331,0.128106,0.095238},
      {0.078416,0.075684,0.094028,0.118257,0.080618},
      {0.078566,0.042963,0.036409,0.034765,0.043790},
      {0.060464,0.025592,0.024057,0.006599,0.020000},
      {0.071746,0.016706,0.000000,0.030430,0.000000},
      {0.043074,0.000000,1.000000,0.000000,1.000000}
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
    elHLTSF = {
      {0.958734,0.879279,0.800000,0.930872,0.819886},
      {0.972084,0.890878,1.091776,0.814388,0.934196},
      {0.936706,0.887334,1.062766,0.932065,0.935490},
      {0.939349,0.976644,0.868087,0.989084,1.057438},
      {0.993655,0.979765,0.940140,0.972614,0.955556},
      {0.985202,0.995961,1.000000,0.909091,1.000000},
      {0.977273,0.942952,1.000000,1.000000,1.000000}
    };
    elHLTSFUnc = {
      {0.131348,0.166195,0.687125,0.291188,0.370062},
      {0.098831,0.124641,0.772932,0.203191,0.347070},
      {0.101623,0.145139,0.693954,0.236831,0.382964},
      {0.083849,0.121101,0.499959,0.204401,0.338514},
      {0.060334,0.083162,0.354400,0.153739,0.239221},
      {0.126580,0.187983,0.764234,0.345859,0.712837},
      {0.251020,0.428866,1.000000,0.851152,1.000000}
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
    elHLTSF = {
      {0.948653,0.946580,1.030069,0.789884,0.914935},
      {0.959176,0.975594,0.982969,0.878223,0.978510},
      {0.914667,0.921914,0.859925,0.907098,0.856968},
      {0.936572,0.964878,0.790190,0.930776,1.000734},
      {0.977059,0.974735,0.939183,0.934640,0.870107},
      {0.988220,0.981175,1.000000,0.956091,1.011764},
      {1.002034,0.974624,0.000000,0.857143,1.000000}
    };
    elHLTSFUnc = {
      {0.089153,0.129760,0.493564,0.183087,0.280660},
      {0.068586,0.097296,0.411438,0.140571,0.206748},
      {0.071284,0.102117,0.422389,0.154512,0.228542},
      {0.059381,0.082177,0.313267,0.126428,0.191029},
      {0.042942,0.058278,0.269807,0.094799,0.145511},
      {0.090644,0.142771,1.000000,0.234939,0.530056},
      {0.187096,0.275891,0.000000,0.505971,1.000000}
    };
  }    


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
