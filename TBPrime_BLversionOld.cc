// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
// Fxn to return any and all float TPrime and BPrime variables needed for plotting
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
// This depends on decay products, so there are going to need to be modifications to those. 
// SHould be able to be shortened a great deal, since we don't need that many decay products. The structure of the decay IDs will need to be explored when the time comes to trim this down.
RVec<float> three_jet(TLorentzVector top_lv, TLorentzVector Wlv, int isLeptonic, RVec<float>& ak8_pt, RVec<float>& ak8_eta, RVec<float>& ak8_phi, RVec<float>& ak8_mass, RVec<float>& dnn_T_DeepAK8Calc_PtOrdered, RVec<float>& dnn_H_DeepAK8Calc_PtOrdered, RVec<float>& dnn_Z_DeepAK8Calc_PtOrdered, RVec<float>& dnn_W_DeepAK8Calc_PtOrdered, RVec<float>& dnn_B_DeepAK8Calc_PtOrdered, RVec<int>& dnn_largest_DeepAK8Calc_PtOrdered, RVec<float>& theJetAK8SoftDropCorr_PtOrdered)
{
        RVec<int> validBTagged (2,0);
        int val = -99;
        int tag = -99;
        TLorentzVector jet_lv;
        std::vector<pair<TLorentzVector,int>> jets_lv;
        float deltaR_leptonicparticle_AK8_PtOrdered = 0;
        bool isLeptonic_W = false;
        bool isLeptonic_t = false;

	for(unsigned int ijet=0; ijet < ak8_pt.size(); ijet++)
	{
		jet_lv.SetPtEtaPhiM(ak8_pt.at(ijet),ak8_eta.at(ijet),ak8_phi.at(ijet),ak8_mass.at(ijet));
		if(isLeptonic == 0)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(Wlv);
			isLeptonic_W = true;
		}
		if(isLeptonic == 1)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(top_lv);
			isLeptonic_t = true;
		}
		// Get 3 highest-pT jets that are not close to t/W (deltaR > .8) and store AK8 index and 4-vector
		if(jets_lv.size() > 1)
		{
			continue;
		}
		if(jet_lv.DeltaR(top_lv) > 0.8 and isLeptonic_t) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
		if(jet_lv.DeltaR(Wlv) > 0.8 and isLeptonic_W) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
	}

	float highPtAK8Jet1_SoftDropCorrectedMass = -999;
	float highPtAK8Jet2_SoftDropCorrectedMass = -999;
	float Bprime1_DeepAK8_Mass = -999;
	float Bprime1_DeepAK8_Pt = -9999;
	float Bprime1_DeepAK8_Eta = 9;
	float Bprime1_DeepAK8_Phi = 9;
	float Bprime1_DeepAK8_deltaR = -9;
	float probSum_DeepAK8_decay = -999;
	float probSum_DeepAK8_four = -999;
	float jet2_DeepAK8 = -999;
	float leptonicBprimeJetIDs_DeepAK8 = -1;
	bool validBDecay_DeepAK8 = false;
	bool taggedTjet = false;
	bool taggedWbjetJet = false;
	bool taggedWjet = false;
	bool is2Jet = false;
	RVec<float> hadronicBprimeJetIDs_DeepAK8 (2,0);
	
	// ----------------------------------------------------------------------------
	// VLQ Decay -- 3 AK8 jets away from leptonic particle
	// ----------------------------------------------------------------------------
	if(jets_lv.size() > 2) {std::cout << "Problem: > 2 AK8s for Tprime reco" << std::endl;}
	else if(jets_lv.size() == 0){std::cout << "Empty jets_lv" << std::endl; }
	int npass_ThreeJets = 0;
	if(jets_lv.size() == 1 && jets_lv.size() == 2)
	{
//		std::cout << "Calculating..." << std::endl;
		npass_ThreeJets = npass_ThreeJets + 1;
		probSum_DeepAK8_decay = 0;
		probSum_DeepAK8_four = 0;
		for (unsigned int i = 0; i < jets_lv.size(); i++)
		{
			// "four" raw sum of heavy particle probabilities
			probSum_DeepAK8_four += dnn_W_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second) + dnn_T_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second) + dnn_Z_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second) + dnn_H_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second);
		}
		
		// Save masses of the 3 jets for plotting
		highPtAK8Jet1_SoftDropCorrectedMass = theJetAK8SoftDropCorr_PtOrdered.at(jets_lv.at(0).second);
		if(jets_lv.size() == 2)
		{
			is2Jet = true;
			highPtAK8Jet2_SoftDropCorrectedMass = theJetAK8SoftDropCorr_PtOrdered.at(jets_lv.at(1).second);
		}
		// get the tags
		float jet1_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(0).second);
		if(jets_lv.size() == 2)
		{
			float jet2_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(1).second);
		}
		// pair up the jet tag with the pT index 0,1,2 and sort by tag (orders J, T, H, Z, W, B)
		std::vector <pair<float,float>> decayJets_DeepAK8;
		decayJets_DeepAK8.push_back(std::make_pair(jet1_DeepAK8,0));
		if(jets_lv.size() == 2)
		{
			decayJets_DeepAK8.push_back(std::make_pair(jet2_DeepAK8,1));
		}
		std::sort(decayJets_DeepAK8.begin(),decayJets_DeepAK8.end());
		TLorentzVector Bprime1_DeepAK8_lv;
		validBDecay_DeepAK8 = false;
	
		// Only 1 Decay
		if(isLeptonic_W)
		{
			if(decayJets_DeepAK8.at(0).first == 1)// && decayJets_DeepAK8.at(1).first==4)
			{
				validBDecay_DeepAK8 = true; val = 1;
				tag = 8;
				taggedTjet = true;
				Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
                                Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);

			}
			// if jet1 ==4 and is2Jets and jet2 == 5
			else if(decayJets_DeepAK8.at(0).first==4 && is2Jet == true && decayJets_DeepAK8.at(1).first==5)
			{
                               validBDecay_DeepAK8 = true; val = 1;
                               tag = 8;
                               taggedWbjetJet = true;
                               Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first; // decayJets.second gives the jets_lv index to get 4-vec
                               Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first);
                        }
			else
			{
				if(jets_lv.size() == 1)
				{
	                                validBDecay_DeepAK8 = false; val = 1;
                 	                Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
	                    	        Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
				}
				else
				{
	                                validBDecay_DeepAK8 = false; val = 1;
                 	                Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
                        	        Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first);
				}
			} 
		}
		else //isLeptonic_t
		{
			if(decayJets_DeepAK8.at(0).first==4)// && decayJets_DeepAK8.at(1).first==4)
                        {
                                validBDecay_DeepAK8 = true; val = 1;
                                tag = 8;
                                taggedWjet = true;
                                Bprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(0).second).first; // decayJets.second gives the jets_lv index to get 4-vec
                                Bprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
                        }
			else
			{
	                 	validBDecay_DeepAK8 = false; val = 1;                                                      			
                       	tag = 8;
	                     	Bprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
	                   	Bprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);

			}
		}	
		// Edited here 
//		if(Bprime1_DeepAK8_lv.M() != -999)
//		{
		Bprime1_DeepAK8_Mass = Bprime1_DeepAK8_lv.M();
		Bprime1_DeepAK8_Pt = Bprime1_DeepAK8_lv.Pt();
		Bprime1_DeepAK8_Eta = Bprime1_DeepAK8_lv.Eta();
		Bprime1_DeepAK8_Phi = Bprime1_DeepAK8_lv.Phi();
//		}
	}
	// Include validBDecy_DeepAK8 variable and what it's decayed as 
	//leptonic and hadronic jetIDs can be slimmed out, no longer exist
	RVec<float> TandBPrimeVec = {Bprime1_DeepAK8_Mass,Bprime1_DeepAK8_Pt,Bprime1_DeepAK8_Eta,Bprime1_DeepAK8_Phi,Bprime1_DeepAK8_deltaR,highPtAK8Jet1_SoftDropCorrectedMass,highPtAK8Jet2_SoftDropCorrectedMass,leptonicBprimeJetIDs_DeepAK8,hadronicBprimeJetIDs_DeepAK8[0],hadronicBprimeJetIDs_DeepAK8[1]};
	return TandBPrimeVec;
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
// 									TTAGGING FXN
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
RVec<int> three_jet_TTag(TLorentzVector top_lv, TLorentzVector Wlv, int isLeptonic, RVec<float>& ak8_pt, RVec<float>& ak8_eta, RVec<float>& ak8_phi, RVec<float>& ak8_mass, RVec<int>& dnn_largest_DeepAK8Calc_PtOrdered, int hadronicTprimeJetIDs1_DeepAK8, int hadronicTprimeJetIDs2_DeepAK8)
{
	RVec<int> validTTagged (2,0);
	int val = -99;
	int tag = -99;
	TString tagCut = "";
	TLorentzVector jet_lv;
	std::vector<pair<TLorentzVector,int>> jets_lv;
	float deltaR_leptonicparticle_AK8_PtOrdered = 0;
	bool isLeptonic_W = false;
	bool isLeptonic_t = false;
	for(unsigned int ijet=0; ijet < ak8_pt.size(); ijet++)
	{
		jet_lv.SetPtEtaPhiM(ak8_pt.at(ijet),ak8_eta.at(ijet),ak8_phi.at(ijet),ak8_mass.at(ijet));
		if(isLeptonic == 0)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(Wlv);
			isLeptonic_W = true;
		}
		if(isLeptonic == 1)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(top_lv);
			isLeptonic_t = true;
		}
		// Get 3 highest-pT jets that are not close to t/W (deltaR > .8) and store AK8 index and 4-vector
		if(jets_lv.size() >= 2){continue;}
		if(jet_lv.DeltaR(top_lv) > 0.8 and isLeptonic_t) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
		if(jet_lv.DeltaR(Wlv) > 0.8 and isLeptonic_W) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
	}
	
	bool validTDecay_DeepAK8 = false; // [0,y]
	bool taggedBWBW_DeepAK8 = false; // [x,0]
	bool taggedTHBW_DeepAK8 = false; // [x,1]
	bool taggedTHTH_DeepAK8 = false; // [x,2]
	bool taggedTZBW_DeepAK8 = false; // [x,3]
	bool taggedTZTH_DeepAK8 = false; // [x,4]
	bool taggedTZTZ_DeepAK8 = false; // [x,5]
	
	// ----------------------------------------------------------------------------
	// VLQ Decay -- 3 AK8 jets away from leptonic particle
	// ----------------------------------------------------------------------------
	if(jets_lv.size() > 3) {std::cout << "Problem: > 3 AK8s for Tprime reco" << std::endl;}
	if(jets_lv.size() == 3)
	{
		// ----------------------------------------------------------------------------
		// DeepAK8 SECTION -- TT
		// ----------------------------------------------------------------------------
		
		// get the tags
		int jet1_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(0).second);
		int jet2_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(1).second);
		int jet3_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(2).second);
		// pair up the jet tag with the pT index 0,1,2 and sort by tag (orders J, T, H, Z, W, B)
		std::vector <pair<int,int>> decayJets_DeepAK8;
		decayJets_DeepAK8.push_back(std::make_pair(jet1_DeepAK8,0));
		decayJets_DeepAK8.push_back(std::make_pair(jet2_DeepAK8,1));
		decayJets_DeepAK8.push_back(std::make_pair(jet3_DeepAK8,2));
		std::sort(decayJets_DeepAK8.begin(),decayJets_DeepAK8.end());
		
		// Start forming 4 particle groups
		validTDecay_DeepAK8 = false;
		if(isLeptonic_t)
		{
			if(decayJets_DeepAK8.at(0).first==2 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> tH bW, BB -> tW bH
				validTDecay_DeepAK8 = true; val = 0;
				taggedTHBW_DeepAK8 = true; tag = 1;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==2 && decayJets_DeepAK8.at(2).first==2)
			{ // TTbar --> tH and tH
				validTDecay_DeepAK8 = true; val = 0;
				taggedTHTH_DeepAK8 = true; tag = 2;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==2 && decayJets_DeepAK8.at(2).first==3)
			{ // TTbar --> tH and tZ
				validTDecay_DeepAK8 = true; val = 0;
				taggedTZTH_DeepAK8 = true; tag = 4;
			}
			else if(decayJets_DeepAK8.at(0).first==3 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> tZ bW, BB -> tW bZ
				validTDecay_DeepAK8 = true; val = 0;
				taggedTZBW_DeepAK8 = true; tag = 3;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==3 && decayJets_DeepAK8.at(2).first==3)
			{ // TTbar --> tZ tZ
				validTDecay_DeepAK8 = true; val = 0;
				taggedTZTZ_DeepAK8 = true; tag = 5;
			}
		}
		else
		{ // isLeptonic_W
			if(decayJets_DeepAK8.at(0).first==4 && decayJets_DeepAK8.at(1).first==5 && decayJets_DeepAK8.at(2).first==5)
			{ // bW bW
				validTDecay_DeepAK8 = true; val = 0;
				taggedBWBW_DeepAK8 = true; tag = 0;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==3 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> bW tZ, BB -> tW bZ
				validTDecay_DeepAK8 = true; val = 0;
				taggedTZBW_DeepAK8 = true; tag = 3;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==2 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> bW tH, BB -> tW bH
				validTDecay_DeepAK8 = true; val = 0;
				taggedTHBW_DeepAK8 = true; tag = 1;
			}
		}
	}
	if(!validTDecay_DeepAK8)
	{
		// signal categories for only hadronic VLQ valid
		if(hadronicTprimeJetIDs1_DeepAK8 == 1 && hadronicTprimeJetIDs2_DeepAK8 == 2) {tag = 9;} // 'tH'
		else if(hadronicTprimeJetIDs1_DeepAK8 == 1 && hadronicTprimeJetIDs2_DeepAK8 == 3) {tag = 10;} // 'tZ'
		else if(hadronicTprimeJetIDs1_DeepAK8 == 4 && hadronicTprimeJetIDs2_DeepAK8 == 5) {tag = 11;} // 'bW'
		// if whichSig == 'TT' and 'tH' not in tag and 'tZ' not in tag and 'bW' not in tag:
		else{tag = -5;}
	}
	validTTagged = {val,tag};
	return validTTagged;
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
//  									 BTAGGING FXN 
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
RVec<int> three_jet_BTag(TLorentzVector top_lv, TLorentzVector Wlv, int isLeptonic, RVec<float>& ak8_pt, RVec<float>& ak8_eta, RVec<float>& ak8_phi, RVec<float>& ak8_mass, RVec<int>& dnn_largest_DeepAK8Calc_PtOrdered, int hadronicBprimeJetIDs1_DeepAK8, int hadronicBprimeJetIDs2_DeepAK8)
{
	RVec<int> validBTagged (2,0);
	int val = -99;
	int tag = -99;
	TLorentzVector jet_lv;
	std::vector<pair<TLorentzVector,int>> jets_lv;
	float deltaR_leptonicparticle_AK8_PtOrdered = 0;
	bool isLeptonic_W = false;
	bool isLeptonic_t = false;
	for(unsigned int ijet=0; ijet < ak8_pt.size(); ijet++)
	{
		jet_lv.SetPtEtaPhiM(ak8_pt.at(ijet),ak8_eta.at(ijet),ak8_phi.at(ijet),ak8_mass.at(ijet));
		if(isLeptonic == 0)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(Wlv);
			isLeptonic_W = true;
		}
		if(isLeptonic == 1)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(top_lv);
			isLeptonic_t = true;
		}
		// Get 3 highest-pT jets that are not close to t/W (deltaR > .8) and store AK8 index and 4-vector
		if(jets_lv.size() >= 3){continue;}
		if(jet_lv.DeltaR(top_lv) > 0.8 and isLeptonic_t) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
		if(jet_lv.DeltaR(Wlv) > 0.8 and isLeptonic_W) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
	}
	
	bool validBDecay_DeepAK8 = false; // [1,y]
	bool taggedTZBW_DeepAK8 = false; // [x,3] T and B
	bool taggedTWTW_DeepAK8 = false; // [x,6] B
	bool taggedBZTW_DeepAK8 = false; // [x,7] B
	bool taggedBHTW_DeepAK8 = false; // [x,8] B
	
	// ----------------------------------------------------------------------------
	// VLQ Decay -- 3 AK8 jets away from leptonic particle
	// ----------------------------------------------------------------------------
	if(jets_lv.size() > 3) {std::cout << "Problem: > 3 AK8s for Tprime reco" << std::endl;}
	if(jets_lv.size() == 3)
	{
		// ----------------------------------------------------------------------------
		// DeepAK8 SECTION -- TT
		// ----------------------------------------------------------------------------
		
		// get the tags
		int jet1_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(0).second);
		int jet2_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(1).second);
		int jet3_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(2).second);
		// pair up the jet tag with the pT index 0,1,2 and sort by tag (orders J, T, H, Z, W, B)
		std::vector <pair<int,int>> decayJets_DeepAK8;
		decayJets_DeepAK8.push_back(std::make_pair(jet1_DeepAK8,0));
		decayJets_DeepAK8.push_back(std::make_pair(jet2_DeepAK8,1));
		decayJets_DeepAK8.push_back(std::make_pair(jet3_DeepAK8,2));
		std::sort(decayJets_DeepAK8.begin(),decayJets_DeepAK8.end());
		
		// Start forming 4 particle groups
		validBDecay_DeepAK8 = false;
		if(isLeptonic_t)
		{
			if(decayJets_DeepAK8.at(0).first==2 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> tH bW, BB -> tW bH
				validBDecay_DeepAK8 = true; val = 1;
				taggedBHTW_DeepAK8 = true; tag = 8;
			}
			else if(decayJets_DeepAK8.at(0).first==3 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> tZ bW, BB -> tW bZ
				validBDecay_DeepAK8 = true; val = 1;
				taggedBZTW_DeepAK8 = true; tag = 7;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==4)
			{ // BB/XX -> tW tW, jets t W W
				validBDecay_DeepAK8 = true; val = 1;
				taggedTWTW_DeepAK8 = true; tag = 6;
			}
		}
		else
		{ // isLeptonic_W
			if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==3 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> bW tZ, BB -> tW bZ
				validBDecay_DeepAK8 = true; val = 1;
				taggedBZTW_DeepAK8 = true; tag = 7;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==2 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> bW tH, BB -> tW bH
				validBDecay_DeepAK8 = true; val = 1;
				taggedBHTW_DeepAK8 = true; tag = 8;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==1 && decayJets_DeepAK8.at(2).first==4)
			{ // BB -> tW tW, jets t t W
				validBDecay_DeepAK8 = true; val = 1;
				taggedTWTW_DeepAK8 = true; tag = 6;
			}
		}
	}
	if(!validBDecay_DeepAK8)
	{
		// signal categories for only hadronic VLQ valid
		if(hadronicBprimeJetIDs1_DeepAK8 == 1 && hadronicBprimeJetIDs2_DeepAK8 == 4) {tag = 12;} // 'tW'
		else if(hadronicBprimeJetIDs1_DeepAK8 == 3 && hadronicBprimeJetIDs2_DeepAK8 == 5) {tag = 13;} // 'bZ'
		else if(hadronicBprimeJetIDs1_DeepAK8 == 2 && hadronicBprimeJetIDs2_DeepAK8 == 5){tag = 14;} // 'bH'
		// if(whichSig == "BB" and 'bH' not in tag and 'bZ' not in tag and 'tW' not in tag)
		else{tag = -5;}
	}
	validBTagged = {val,tag};
	return validBTagged;
};
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------



