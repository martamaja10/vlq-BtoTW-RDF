// callRDF.C with lwtnn includes
void callHistDraw(std::vector<std::string> inputFile)
{
	gSystem->Load("~/nobackup/YOURWORKINGAREA/CMSSW_11_0_0/lib/slc7_amd64_gcc820/liblwtnnlwtnn.so");
	gROOT->ProcessLine(".x histoDrawRDF.C("/+inputFile+"/");
};

