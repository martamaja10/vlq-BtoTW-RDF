// callRDF.C with lwtnn includes
void callRDF(TString testNum, TString inputfile)
{
	gSystem->Load("../lwtnn/build/lib/liblwtnn.so");
	gROOT->ProcessLine(".x runRDF.C(\""+testNum+"\",\""+inputfile+"\")");
};
