//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 16 20:39:17 2015 by ROOT version 6.02/05
// from TTree ljmet/ljmet
// found on file: /eos/uscms/store/user/lpcljm/LJMet_1lep_110915/X53X53_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns/X53X53_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns_1.root
//////////////////////////////////////////////////////////

#ifndef rdf_h
#define rdf_h

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

#include "lwtnn/lwtnn/interface/parse_json.hh"
#include "lwtnn/lwtnn/interface/LightweightNeuralNetwork.hh"

enum shift:char;

using namespace std;
using namespace ROOT::VecOps;

class rdf {
public :
   TTree          *inputTree;   //!pointer to the analyzed TTree or TChain
   TFile          *inputFile;
   TString         psOutName, fsOutName;
   Int_t           fCurrent; //!current Tree number in a TChain

   Bool_t          isSig;
   Bool_t          isTpTp;
   Bool_t          isBpBp;
   Bool_t          isTOP;
   Bool_t          isMadgraphBkg; // W, Z, QCD
   Bool_t          isMC;
   Bool_t          isSM;
   Bool_t          isSE;
   Bool_t          isTT;
   Bool_t          isVV;
   Bool_t          isTTincMtt0to1000;
   Bool_t          isTTincMtt0to700;
   Bool_t          isTTincMtt700to1000;
   Bool_t          isTTincMtt1000toInf;
   Bool_t          outBWBW;
   Bool_t          outTZBW;
   Bool_t          outTHBW;
   Bool_t          outTZTH;
   Bool_t          outTZTZ;
   Bool_t          outTHTH;
   Bool_t          outTWTW;
   Bool_t          outBZTW;
   Bool_t          outBHTW;
   Bool_t          outBZBH;
   Bool_t          outBZBZ;
   Bool_t          outBHBH;
   Bool_t          isNominal;
   Bool_t          isBUp;
   Bool_t          isBDn;
   Bool_t          isLUp;
   Bool_t          isLDn;

   lwt::LightweightNeuralNetwork* lwtnnTT;    

   // Fixed size dimensions of array or collections stored in the TTree if any.
 
   rdf(TString inputFileName, TString preselFileName, TString finalselFileName);
   virtual ~rdf();
   virtual void     step1RDF_forLJMet(std::string sample, TString chan, TString testNum, int year);
   //virtual void     histoDrawRDF(vector<std::string> inputFile);
};

#endif

#ifdef rdf_cxx
rdf::rdf(TString inputFileName, TString preselFileName, TString finalselFileName) : inputTree(0), inputFile(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  psOutName = preselFileName;
  fsOutName = finalselFileName;

  string dnnFileTT = "vlq_mlp_June_08_20_TT.json";
  ifstream input_cfgTT( dnnFileTT );
  lwt::JSONConfig cfgTT = lwt::parse_json(input_cfgTT);
  lwtnnTT = new lwt::LightweightNeuralNetwork(cfgTT.inputs, cfgTT.layers, cfgTT.outputs);

  isSig  = (inputFileName.Contains("prime") || inputFileName.Contains("X53") || inputFileName.Contains("ChargedHiggs_Hplus"));
  if(isSig){
    if(inputFileName.Contains("Tprime")) isTpTp = true;
    else if(inputFileName.Contains("Bprime")) isBpBp = true;
  }

  isMadgraphBkg = (inputFileName.Contains("QCD") || inputFileName.Contains("madgraphMLM"));
  isTOP = (inputFileName.Contains("Mtt") || inputFileName.Contains("ST") || inputFileName.Contains("ttZ") || inputFileName.Contains("ttW") || inputFileName.Contains("ttH") || inputFileName.Contains("TTTo"));
  isTT = (inputFileName.Contains("TT_Tune") || inputFileName.Contains("Mtt") || inputFileName.Contains("TTTo"));
  isVV = (inputFileName.Contains("WW_") || inputFileName.Contains("WZ_") || inputFileName.Contains("ZZ_"));
  isMC = !(inputFileName.Contains("Single") || inputFileName.Contains("Data18"));
  isSM = inputFileName.Contains("SingleMuon");
  isSE = (inputFileName.Contains("SingleElectron") || inputFileName.Contains("EGamma"));

  isBUp = false; // these will now get changed in makeRdfDnn.C
  isBDn = false;
  isLUp = false;
  isLDn = false;
  isNominal = true;
  isTTincMtt0to700    = preselFileName.Contains("Mtt0to700");
  isTTincMtt0to1000   = preselFileName.Contains("Mtt0to1000");
  isTTincMtt700to1000 = preselFileName.Contains("Mtt700to1000");
  isTTincMtt1000toInf = preselFileName.Contains("Mtt1000toInf");
  outBWBW = preselFileName.Contains("BWBW");
  outTZBW = preselFileName.Contains("TZBW");
  outTHBW = preselFileName.Contains("THBW");
  outTZTH = preselFileName.Contains("TZTH");
  outTZTZ = preselFileName.Contains("TZTZ");
  outTHTH = preselFileName.Contains("THTH");
  outTWTW = preselFileName.Contains("TWTW");
  outBZTW = preselFileName.Contains("BZTW");
  outBHTW = preselFileName.Contains("BHTW");
  outBZBH = preselFileName.Contains("BZBH");
  outBZBZ = preselFileName.Contains("BZBZ");
  outBHBH = preselFileName.Contains("BHBH");

  /* std::cout<<"Opening file: "<<inputFileName<<std::endl; */
  /* if(!(inputFile=TFile::Open(inputFileName))){ */
  /*   std::cout<<"WARNING! File doesn't exist! Exiting" << std::endl; */
  /*   exit(1); */
  /* } */

  // Now done in the .cc
  /* inputTree=(TTree*)inputFile->Get("ljmet"); */
  /* if(inputTree->GetEntries()==0){ */
  /*   std::cout<<"WARNING! Found 0 events in the tree!!!!"<<std::endl;; */
  /*   exit(1); */
  /* } */

  //outputFile=new TFile(preselFileName,"RECREATE");   
  
}

rdf::~rdf()
{
  delete lwtnnTT;
  if (!inputTree) return;
  delete inputTree->GetCurrentFile();
}



#endif // #ifdef rdf_cxx
