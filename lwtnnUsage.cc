// -----------------------------------------------
//   LWTNN IMPLIMENTATION AND MYMAP CALCULATION
// -----------------------------------------------
auto varMap = [lwtnnBB](float corr_met_MultiLepCalc, float AK4HT, int NJets_JetSubCalc, int NJetsAK8_JetSubCalc, float AK4HTpMETpLepPt, float jetPt_1, float jetPt_2, float jetPt_3, float sdMass_1, float sdMass_2, float sdMass_3, float dnnJ_1, float dnnJ_2, float dnnJ_3, float dnnT_1, float dnnT_2, float dnnT_3, float dnnH_1, float dnnH_2, float dnnH_3, float dnnZ_1, float dnnZ_2, float dnnZ_3, float dnnW_1, float dnnW_2, float dnnW_3, float dnnB_1, float dnnB_2, float dnnB_3, int dnnLargest_1, int dnnLargest_2, int dnnLargest_3, int nJ_DeepAK8, int nT_DeepAK8, int nH_DeepAK8, int nZ_DeepAK8, int nW_DeepAK8, int nB_DeepAK8, float tau21_1, float tau21_2, float tau21_3, float minDR_leadAK8otherAK8)
{
	ROOT::VecOps::RVec<float> dnn_SigWjetTtbar (6,0);
	std::map<std::string,double> varMapBB;
	std::map<std::string,double> myMapBB;
	
	myMapBB = {
	   {"WjetsBB",-999},
	   {"ttbarBB",-999},
	   {"Bprime",-999},
	};
	
	varMapBB = {
	   {"corr_met_MultiLepCalc", corr_met_MultiLepCalc},
	   {"AK4HTpMETpLepPt", AK4HTpMETpLepPt},
	   {"AK4HT", AK4HT},
	   {"NJets_JetSubCalc", NJets_JetSubCalc},
	   {"NJetsAK8_JetSubCalc", NJetsAK8_JetSubCalc},
	   {"minDR_leadAK8otherAK8", minDR_leadAK8otherAK8},
	   {"nH_DeepAK8", nH_DeepAK8},
	   {"nT_DeepAK8", nT_DeepAK8},
	   {"jetPt_1", jetPt_1},
	   {"jetPt_2", jetPt_2},
	   {"jetPt_3", jetPt_3},
	   {"sdMass_1", sdMass_1},
	   {"sdMass_2", sdMass_2},
	   {"sdMass_3", sdMass_3},
	   {"dnnLargest_1", dnnLargest_1},
	   {"dnnLargest_2", dnnLargest_2},
	   {"dnnLargest_3", dnnLargest_3},
	   {"dnnJ_1", dnnJ_1},
	   {"dnnJ_2", dnnJ_2},
	   {"dnnJ_3", dnnJ_3},
	   {"dnnH_2", dnnH_2},
	   {"dnnH_3", dnnH_3},
	   {"dnnT_1", dnnT_1},
	};
	
	myMapBB = lwtnnBB->compute(varMapBB);
	dnn_SigWjetTtbar[0] = myMapBB["Bprime"];
	dnn_SigWjetTtbar[1] = myMapBB["WjetsBB"];
	dnn_SigWjetTtbar[2] = myMapBB["ttbarBB"];

	return dnn_SigWjetTtbar;
};
