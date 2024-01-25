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
#include "../correctionlib/include/correction.h"

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
  virtual void analyzer_taggerEff(TString testNum, TString jesvar);
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

}

rdf::~rdf()
{
  if (!inputTree)
    return;
  delete inputTree->GetCurrentFile();
}


#endif // #ifdef rdf_cxx
