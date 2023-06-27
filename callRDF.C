// callRDF.C with lwtnn includes
void callRDF(TString testNum, TString inputfile)
{
	gSystem->Load("../../../LCGenv/lwtnn/build/lib/liblwtnn.so");
	gROOT->ProcessLine(".x runRDF.C(\""+testNum+"\",\""+inputfile+"\")");
};
