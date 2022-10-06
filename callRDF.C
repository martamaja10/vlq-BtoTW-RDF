// callRDF.C with lwtnn includes
void callRDF(TString testNum, std::string inputfile)
{
	gSystem->Load("liblwtnnlwtnn.so");
	gROOT->ProcessLine(".x runRDF.C(\""+testNum+"\",\""+inputfile+"\")");
};
