#include "analyzer_RDF.cpp"
#include "cleanJet.cc"
#include "dnnPrep.cc"
#include "W_t_reco.cc"
#include "BPrime.cc"
// Going to need new root files

void runRDF(TString channel,TString testNum, std::string inputFile)
{

	rdf t(inputFile,"preselTree_"+channel+"_"+testNum,"finalselTree_"+channel+"_"+testNum); // names get set to class members, should be known w/o passing
	int year = 2017;

	t.step1RDF_forLJMet(inputFile, channel, testNum ,year);
};
