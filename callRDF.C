// callRDF.C with lwtnn includes
void callRDF(TString channel, TString testNum, std::string inputfile)
{
	gSystem->Load("liblwtnnlwtnn.so");
	gROOT->ProcessLine(".x runRDF.C(\""+channel+"\",\""+testNum+"\",\""+inputfile+"\")");
};
