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

// Header file for the classes stored in the TTree if any.
#include "vector"
#include <string>
#include "TLorentzVector.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

enum shift : char;

using namespace std;
using namespace ROOT::VecOps;

class rdf
{
private:
  TTree *inputTree; //! pointer to the analyzed TTree or TChain
  TFile *inputFile;
  TString psOutName, fsOutName;
  Int_t fCurrent; //! current Tree number in a TChain
  std::vector<std::string> files;

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
  Bool_t isNominal;
  Bool_t isBUp;
  Bool_t isBDn;
  Bool_t isLUp;
  Bool_t isLDn;

  TString sample;

  // Fixed size dimensions of array or collections stored in the TTree if any.
public:
  rdf(std::string inputFileName, TString preselFileName, TString finalselFileName);
  virtual ~rdf();
  virtual void analyzer_RDF(TString testNum);
  // virtual void     histoDrawRDF(vector<std::string> inputFile);

  //auto Bprime_gen_info(unsigned int nGenPart, ROOT::VecOps::RVec<int> &GenPart_pdgId, ROOT::VecOps::RVec<float> &GenPart_mass, ROOT::VecOps::RVec<float> &GenPart_pt, ROOT::VecOps::RVec<float> &GenPart_phi, ROOT::VecOps::RVec<float> &GenPart_eta, ROOT::VecOps::RVec<int> &GenPart_genPartIdxMother, ROOT::VecOps::RVec<int> &GenPart_status, ROOT::VecOps::RVec<int> &GenPart_statusFlags);
  //auto t_gen_info(unsigned int nGenPart, ROOT::VecOps::RVec<int> &GenPart_pdgId, ROOT::VecOps::RVec<float> &GenPart_mass, ROOT::VecOps::RVec<float> &GenPart_pt, ROOT::VecOps::RVec<float> &GenPart_phi, ROOT::VecOps::RVec<float> &GenPart_eta, ROOT::VecOps::RVec<int> &GenPart_genPartIdxMother, ROOT::VecOps::RVec<int> &GenPart_status);
};

#endif

#ifdef rdf_cxx

rdf::rdf(std::string inputFileName, TString preselFileName, TString finalselFileName) : inputTree(0), inputFile(0)
{
  // Make the vector with the text file listed in the input
  //inputFileName = "condor/BPrime.txt"; Delete after testing before submitting condor job

  std::cout << "Input File Path: " << inputFileName << std::endl;
  std::ifstream listFiles;
  listFiles.open(inputFileName);

  std::string file;
  if (listFiles.is_open())
  {
    while (listFiles >> file)
    {
      //std::cout << "Files: " << file << std::endl;
      files.push_back(file);
    }
  }
  std::cout << "Number of Entries: " << files.size() << std::endl;

  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  TString sampleName = file;

  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.

  psOutName = preselFileName;
  fsOutName = finalselFileName;

  isSig = (sampleName.Contains("Bprime"));
  isMadgraphBkg = (sampleName.Contains("QCD") || sampleName.Contains("madgraphMLM"));
  isTOP = (sampleName.Contains("Mtt") || sampleName.Contains("ST") || sampleName.Contains("ttZ") || sampleName.Contains("ttW") || sampleName.Contains("ttH") || sampleName.Contains("TTTo"));
  isTT = (sampleName.Contains("TT_Tune") || sampleName.Contains("Mtt") || sampleName.Contains("TTTo"));
  isVV = (sampleName.Contains("WW_") || sampleName.Contains("WZ_") || sampleName.Contains("ZZ_"));
  isMC = !(sampleName.Contains("Single") || sampleName.Contains("Data18"));
  isSM = sampleName.Contains("SingleMuon");
  isSE = (sampleName.Contains("SingleElectron") || sampleName.Contains("EGamma"));

  isBUp = false; // these will now get changed in makeRdfDnn.C
  isBDn = false;
  isLUp = false;
  isLDn = false;
  isNominal = true;
  isTTincMtt0to700 = preselFileName.Contains("Mtt0to700");
  isTTincMtt0to1000 = preselFileName.Contains("Mtt0to1000");
  isTTincMtt700to1000 = preselFileName.Contains("Mtt700to1000");
  isTTincMtt1000toInf = preselFileName.Contains("Mtt1000toInf");

  // Samples will be signal, ttbar (a background), and possibly others...
  sample = "singletop";
  if (isTT == true)
  {
    sample = "ttbar";
  }
  else if (isSig == true)
  {
    sample = "Bprime";
  }
  else if (isMadgraphBkg == true)
  {
    sample = "wjets";
  }
  std::cout << "Sample: " << sample << std::endl;
}

rdf::~rdf()
{
  if (!inputTree)
    return;
  delete inputTree->GetCurrentFile();
}

#endif // #ifdef rdf_cxx
