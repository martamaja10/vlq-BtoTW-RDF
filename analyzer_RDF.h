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
  vector<string> files;

  // Polynominals for LHE-based WJets HT scaling
  TF1 *poly2 = new TF1("poly2","max(0.402806,0.998174 + (8.40861e-05)*x + (-6.63274e-07)*x*x + (4.09272e-10)*x*x*x + (-9.50233e-14)*x*x*x*x + (7.59648e-18)*x*x*x*x*x)",100,5000);  
  TF1 *poly2U = new TF1("poly2U","max([6],[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)",100,5000);   
  TF1 *poly2D = new TF1("poly2D","max([6],[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)",100,5000);
  // Polynomials for data-derived TOP HT scaling
  TF1 *polyHT = new TF1("polyHT","min(1.0,max([0] + [1]*x,[2]))",700,5000);
  TF1 *polyHTU = new TF1("polyHTU","min(1.0,max([0] + [1]*x + sqrt([3] + 2*x*[4] + x*x*[5]),[2] + [6]))",700,5000);
  TF1 *polyHTD = new TF1("polyHTD","min(1.0,max([0] + [1]*x - sqrt([3] + 2*x*[4] + x*x*[5]),[2] - [6]))",700,5000);

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
  string jesvar;
  string sample;
  string year;

  // Fixed size dimensions of array or collections stored in the TTree if any.
public:
  // Main Methods
  rdf(string inputFileName, TString preselFileName, TString finalselFileName, string testNum1, string testNum2, string yearIn);
  virtual ~rdf();
  virtual void analyzer_RDF(TString testNum1);
};

#endif

#ifdef rdf_cxx

rdf::rdf(string inputFileName, TString preselFileName, TString finalselFileName, string testNum1, string testNum2, string yearIn) : inputTree(0), inputFile(0)
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

  psOutName = preselFileName;
  fsOutName = finalselFileName;

  jesvar = "Nominal"; // FIXME, swap Nominal, JECup, JECdn, JERup, JERdn somewhere;

  // Parse the incoming file names to assign labels
  isSig = (sampleName.Contains("Bprime"));
  isMadgraphBkg = (sampleName.Contains("QCD") || sampleName.Contains("madgraphMLM"));
  isTOP = (sampleName.Contains("Mtt") || sampleName.Contains("ST") || sampleName.Contains("ttZ") || sampleName.Contains("ttW") || sampleName.Contains("ttH") || sampleName.Contains("TTTo"));
  isTT = (sampleName.Contains("TT_Tune") || sampleName.Contains("Mtt") || sampleName.Contains("TTTo"));
  isVV = (sampleName.Contains("WW_") || sampleName.Contains("WZ_") || sampleName.Contains("ZZ_"));
  isMC = !(sampleName.Contains("Single") || sampleName.Contains("Data18") || sampleName.Contains("EGamma"));
  isSM = sampleName.Contains("SingleMuon");
  isSE = (sampleName.Contains("SingleElectron") || sampleName.Contains("EGamma"));


  TObjArray *tokens = sampleName.Tokenize("/");
  sample = ((TObjString *)(tokens->At(5)))->String();
  delete tokens;

  year = yearIn; // May need to change this line to get things to work

  isTTincMtt0to700 = preselFileName.Contains("Mtt0to700"); // FIXME -- not using this yet, but we need to for high-mass ttbar (or some more clever way to weight based on mass value).
  isTTincMtt0to1000 = preselFileName.Contains("Mtt0to1000");
  isTTincMtt700to1000 = preselFileName.Contains("Mtt700to1000");
  isTTincMtt1000toInf = preselFileName.Contains("Mtt1000toInf");

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
