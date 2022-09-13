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
  
  for(unsigned int ijet=0; ijet < ak8_pt.size(); ijet++) {
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
      if(jets_lv.size() > 1) continue;

      if(jet_lv.DeltaR(top_lv) > 0.8 and isLeptonic_t) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
      if(jet_lv.DeltaR(Wlv) > 0.8 and isLeptonic_W) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
  }
//	std::cout << "jets_lv after the loop is now " << jets_lv.size() << std::endl;

  float Bprime1_DeepAK8_Mass = -999;
  float Bprime1_DeepAK8_Pt = -9999;
  float Bprime1_DeepAK8_Eta = 9;
  float Bprime1_DeepAK8_Phi = 9;
  float Bprime1_DeepAK8_deltaR = -9;
  
  TLorentzVector Bprime1_DeepAK8_lv;
  float validBDecay_DeepAK8 = -1;
  
  float taggedTjet = 0;
  float taggedWbjetJet = 0;
  float taggedWjet = 0;
  bool is2Jet = false;
  
  // ----------------------------------------------------------------------------
  // VLQ Decay -- 3 AK8 jets away from leptonic particle
  // ----------------------------------------------------------------------------
  // Change print stmnt
  // Get rid of jet 3, make structure of always jet1, if jet2 exists, then ....
  if(jets_lv.size() > 2) {std::cout << "Problem: > 2 AK8s for Tprime reco" << std::endl;}
  else if(jets_lv.size() == 0){std::cout << "Problem: Empty jets_lv" << std::endl; }
  int npass_ThreeJets = 0;

  if(jets_lv.size() ==1 || jets_lv.size() ==2){
    npass_ThreeJets = npass_ThreeJets + 1;
      
    // Save masses of the 3 jets for plotting
    // Just do jet 1
    if(jets_lv.size() == 2){
      is2Jet = true;
    }
    // get the tags
    float jet1_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(0).second);
    float jet2_DeepAK8 = -1;
    if(jets_lv.size() == 2) jet2_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(1).second);
    
    // pair up the jet tag with the pT index 0,1,2 and sort by tag (orders J, T, H, Z, W, B)
    std::vector <pair<float,float>> decayJets_DeepAK8;
    decayJets_DeepAK8.push_back(std::make_pair(jet1_DeepAK8,0));
    if(jets_lv.size() == 2) decayJets_DeepAK8.push_back(std::make_pair(jet2_DeepAK8,1));
    
    std::sort(decayJets_DeepAK8.begin(),decayJets_DeepAK8.end());
    
    // Start forming 4 particle groups
    validBDecay_DeepAK8 = -1;
    
    // Only 1 Decay
    // Derived from Previous logic: TT -> tZ bW, BB -> tW bZ
    if(isLeptonic_W){
      if(decayJets_DeepAK8.at(0).first == 1){ // leading jet is top tagged
	validBDecay_DeepAK8 = 1;
	taggedTjet = 1;
	Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first; // decayJets.second gives the jets_lv index to get 4-vec
	Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);	    
      }
      else{ // leading jet is not top tagged
	if(!is2Jet){ // all other 1-jet events are invalid
	  validBDecay_DeepAK8 = 0;
	  Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
	  Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
	}
	else{ // it's a 2-jet event
	  if(decayJets_DeepAK8.size() < 2) std::cout << "PROBLEM, need 2 elements in decayJets and have " << decayJets_DeepAK8.size() << std::endl;
	  
	  if(decayJets_DeepAK8.at(0).first==4 && decayJets_DeepAK8.at(1).first==5){ // W and b tags
	    validBDecay_DeepAK8 = 1;
	    taggedWbjetJet = 1;
	    Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first; // decayJets.second gives the jets_lv index to get 4-vec
	    Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first);
	  }
	  else{ // any other set of tags
	    validBDecay_DeepAK8 = 0;
	    Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
	    Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first);
	  }
	}
      }
    }else{ //isLeptonic_t
      if(decayJets_DeepAK8.at(0).first==4){ // leading jet is W-tagged	    
	validBDecay_DeepAK8 = 1;
	taggedWjet = 1;
	Bprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(0).second).first; // decayJets.second gives the jets_lv index to get 4-vec
	Bprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
      }else{
	validBDecay_DeepAK8 = 0;
	Bprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
	Bprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
      }
    }

    Bprime1_DeepAK8_Mass = Bprime1_DeepAK8_lv.M();
    Bprime1_DeepAK8_Pt = Bprime1_DeepAK8_lv.Pt();
    Bprime1_DeepAK8_Eta = Bprime1_DeepAK8_lv.Eta();
    Bprime1_DeepAK8_Phi = Bprime1_DeepAK8_lv.Phi();
  }

  // Include validBDecy_DeepAK8 variable and what it's decayed as 
  //leptonic and hadronic jetIDs can be slimmed out, no longer exist
  RVec<float> TandBPrimeVec = {Bprime1_DeepAK8_Mass,Bprime1_DeepAK8_Pt,Bprime1_DeepAK8_Eta,Bprime1_DeepAK8_Phi,Bprime1_DeepAK8_deltaR,validBDecay_DeepAK8,taggedWbjetJet,taggedTjet,taggedWjet};
  return TandBPrimeVec;
};

