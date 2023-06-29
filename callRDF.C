// callRDF.C with lwtnn includes
void callRDF(TString testNum, std::string inputfile)
{
	gROOT->ProcessLine(".x runRDF.C(\""+testNum+"\",\""+inputfile+"\")");
};
