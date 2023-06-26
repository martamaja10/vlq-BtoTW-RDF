#include "analyzer_RDF.cc"
//#include "constructor.cc"
//#include "lambdaFxns.cc"
#include "cleanJet.cc"
#include "dnnPrep.cc"
#include "W_t_reco.cc"
#include "BPrime.cc"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

void runRDF(TString testNum, std::string inputFile)
{
	rdf t(inputFile, "preselTree_" + testNum, "finalselTree_" + testNum); // names get set to class members, should be known w/o passing
	t.analyzer_RDF(testNum);
};
