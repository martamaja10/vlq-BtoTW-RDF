// Run histogram Drawing Script
#include "histoDrawRDF.C"
//#include "step1RDF_forLJMet.cpp"
//#include "cleanJet.cc"
//#include "dnnPrep.cc"
//#include "W_t_reco.cc"
//#include "TBPrime.cc"

void runHistoDraw(TString trial)
{
        //      Tree Files to run from
	std::vector<std::string> inputFiles = {"RDF_Signal_seniorResearch_Muon__m800.root"};
        //      std::vector<std::string> inputFiles = {"RDF_Signal_seniorResearch_Muon__m1400.root"};
        //      std::vector<std::string> inputFiles = {"RDF_Signal_seniorResearch_Muon__m2000.root"};

        rdf t(TString(inputFiles[1]),"preselTree_"+trial,"finalselTree_"+trial);
        t.histoDrawRDF(inputFiles);
};
