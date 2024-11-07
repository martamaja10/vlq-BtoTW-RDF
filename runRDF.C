#include "analyzer_RDF.cc"
#include "BPrime.cc"
#include "cleanJet.cc"
#include "cut_ptrel.cc"
#include "dnnPrep.cc"
#include "generatorInfo.cc"
#include "utilities.cc"
#include "W_t_reco.cc"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

void runRDF()
{
  string testNum1 = "0";
  string testNum2 = "100";
  string inputFile = "samples_files.txt";
  string year = "2016";
  rdf t(inputFile, testNum1, testNum2, year); // names get set to class members, should be known w/o passing

  bool isData = false;
  if(inputFile.find("Single") != std::string::npos || inputFile.find("EGamma") != std::string::npos) isData = true;

  if(isData) t.analyzer_RDF(testNum1,"Nominal");
  else{
    //vector<TString> shifts = {"Nominal","JECup","JECdn","JERup","JERdn"};
    vector<TString> shifts = {"Nominal"};
    for(size_t i = 0; i < shifts.size(); i++){
      cout << "\nRunning shift " << shifts[i] << endl;

      t.analyzer_RDF(testNum1,shifts[i]);

      cout << "\nFinished shift " << shifts[i] << endl;
    }
  }
  cout << "\nFinished all analyzing" << endl;

};
