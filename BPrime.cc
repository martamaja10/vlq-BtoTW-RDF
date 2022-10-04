// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
// Fxn to return any and all float TPrime and BPrime variables needed for plotting
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------

RVec<float> BPrime_reco(TLorentzVector top_lv, TLorentzVector Wlv, int leptonicParticle, RVec<float>& ak8_pt, RVec<float>& ak8_eta, RVec<float>& ak8_phi, RVec<float>& ak8_mass, RVec<int>& dpak8_tag, RVec<float>& ak8_sdmass)
{
  RVec<int> validBTagged (2,0);
  int val = -99;
  int tag = -99;
  TLorentzVector jet_lv;
  std::vector<pair<TLorentzVector,int>> jets_lv;
  float DR_LP_AK8 = 0;
  bool isLeptonic_W = false;
  bool isLeptonic_t = false;
  
  for(unsigned int ijet=0; ijet < ak8_pt.size(); ijet++) {
      jet_lv.SetPtEtaPhiM(ak8_pt.at(ijet),ak8_eta.at(ijet),ak8_phi.at(ijet),ak8_mass.at(ijet));
      if(leptonicParticle == 0)
	{
	  DR_LP_AK8 = jet_lv.DeltaR(Wlv);
	  isLeptonic_W = true;
	}
      if(leptonicParticle == 1)
	{
	  DR_LP_AK8 = jet_lv.DeltaR(top_lv);
	  isLeptonic_t = true;
	}
      if(jets_lv.size() > 1) continue; // we only take 2 for single Bprime

      if(jet_lv.DeltaR(top_lv) > 0.8 and isLeptonic_t) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
      if(jet_lv.DeltaR(Wlv) > 0.8 and isLeptonic_W) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
  }

  float Bprime_mass = -999;
  float Bprime_pt = -9999;
  float Bprime_eta = 9;
  float Bprime_phi = 9;
  float Bprime_DR = -9;
  
  TLorentzVector Bprime_lv;
  float validBDecay = -1;
  
  float taggedTjet = 0;
  float taggedWbjetJet = 0;
  float taggedWjet = 0;
  bool is2Jet = false;
  
  // ----------------------------------------------------------------------------
  // VLQ Decay -- 3 AK8 jets away from leptonic particle
  // ----------------------------------------------------------------------------
  if(jets_lv.size() > 2) {std::cout << "Problem: > 2 AK8s for Tprime reco" << std::endl;}
  else if(jets_lv.size() == 0){std::cout << "Problem: Empty jets_lv" << std::endl; }
	
  int npass_reco = 0;

  if(jets_lv.size() ==1 || jets_lv.size() ==2){
    npass_reco = npass_reco + 1;
      
    if(jets_lv.size() == 2){
      is2Jet = true;
    }
	  
    // get the tags
    float jet1_tag = dpak8_tag.at(jets_lv.at(0).second);
    float jet2_tag = -1;
    if(is2Jet) jet2_tag = dpak8_tag.at(jets_lv.at(1).second);
    
    // pair up the jet tag with the pT index 0,1,2 and sort by tag (orders J, T, W,)
    std::vector <pair<float,float>> decayJets;
    decayJets.push_back(std::make_pair(jet1_tag,0));
    if(is2Jet) decayJets.push_back(std::make_pair(jet2_tag,1));
    std::sort(decayJets.begin(),decayJets.end());
    
    // Start forming 2 particle groups
    validBDecay = -1;
    
    // Only 1 Decay
    // Derived from Previous logic: TT -> tZ bW, BB -> tW bZ
    if(isLeptonic_W){
      if(decayJets.at(0).first == 1){ // leading jet is top tagged
	validBDecay = 1;
	taggedTjet = 1;
	Bprime_lv = Wlv+jets_lv.at(decayJets.at(0).second).first; // decayJets.second gives the jets_lv index to get 4-vec
	Bprime_deltaR = Wlv.DeltaR(jets_lv.at(decayJets.at(0).second).first);	    
      }
      else{ // leading jet is not top tagged
	if(!is2Jet){ // all other 1-jet events are invalid
	  validBDecay = 0;
	  Bprime_lv = Wlv+jets_lv.at(decayJets.at(0).second).first;
	  Bprime_deltaR = Wlv.DeltaR(jets_lv.at(decayJets.at(0).second).first);
	}
	else{ // it's a 2-jet event
	  if(decayJets.size() < 2) std::cout << "PROBLEM, need 2 elements in decayJets and have " << decayJets.size() << std::endl;
	  
	  if(decayJets.at(0).first==2){ // W tag found
	    validBDecay_DeepAK8 = 1;
	    taggedWbjetJet = 1;
	    Bprime_lv = Wlv+jets_lv.at(decayJets.at(0).second).first+jets_lv.at(decayJets.at(1).second).first; // decayJets.second gives the jets_lv index to get 4-vec
	    Bprime_deltaR = Wlv.DeltaR(jets_lv.at(decayJets.at(0).second).first+jets_lv.at(decayJets.at(1).second).first);
	  }
	  else{ // any other set of tags
	    validBDecay = 0;
	    Bprime_lv = Wlv+jets_lv.at(decayJets.at(0).second).first+jets_lv.at(decayJets.at(1).second).first;
	    Bprime_deltaR = Wlv.DeltaR(jets_lv.at(decayJets.at(0).second).first+jets_lv.at(decayJets.at(1).second).first);
	  }
	}
      }
    }else{ //isLeptonic_t
      if(decayJets.at(0).first==2){ // leading jet is W-tagged	    
	validBDecay = 1;
	taggedWjet = 1;
	Bprime_lv = top_lv+jets_lv.at(decayJets.at(0).second).first; // decayJets.second gives the jets_lv index to get 4-vec
	Bprime_deltaR = top_lv.DeltaR(jets_lv.at(decayJets.at(0).second).first);
      }else{
	validBDecay = 0;
	Bprime_lv = top_lv+jets_lv.at(decayJets.at(0).second).first;
	Bprime_deltaR = top_lv.DeltaR(jets_lv.at(decayJets.at(0).second).first);
      }
    }

    Bprime_mass = Bprime_lv.M();
    Bprime_pt = Bprime_lv.Pt();
    Bprime_eta = Bprime_lv.Eta();
    Bprime_phi = Bprime_lv.Phi();
  }

  // Include validBDecy_DeepAK8 variable and what it's decayed as 
  //leptonic and hadronic jetIDs can be slimmed out, no longer exist
  RVec<float> BPrimeVec = {Bprime_mass,Bprime_pt,Bprime_eta,Bprime_phi,Bprime_deltaR,validBDecay,taggedWbjetJet,taggedTjet,taggedWjet};
  return BPrimeVec;
};

