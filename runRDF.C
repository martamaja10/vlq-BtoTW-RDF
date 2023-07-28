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

void runRDF(string testNum1, string testNum2, string inputFile, string year)
{
	rdf t(inputFile, "preselTree_" + testNum1, "finalselTree_" + testNum1, testNum1, testNum2, year); // names get set to class members, should be known w/o passing
	t.analyzer_RDF(testNum1);
};
