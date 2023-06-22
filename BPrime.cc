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
  float Bprime_pt = -999;
  float Bprime_eta = 9;
  float Bprime_phi = 9;
  float Bprime_DR = -9;
  float Bprime_ptbal = -999;
  float Bprime_chi2 = -999;
  
  TLorentzVector Wcand_lv;
  TLorentzVector Tcand_lv;
  TLorentzVector Bprime_lv;
  float validBDecay = -1;
  
  float taggedTjet = 0;
  float taggedTjet2nd = 0;
  float taggedWbjetJet = 0;
  float taggedWjet = 0;
  float taggedWjet2nd = 0;
  bool is2Jet = false;
  
  // ----------------------------------------------------------------------------
  // VLQ Decay -- 1 or 2 AK8 jets away from leptonic particle
  // ----------------------------------------------------------------------------
  if(jets_lv.size() > 2) {std::cout << "Problem: > 2 AK8s for Bprime reco" << std::endl;}
  else if(jets_lv.size() == 0){std::cout << "Failed reco event, no jets" << std::endl; }
	
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
    
    // Derived from Previous logic: TT -> tZ bW, BB -> tW bZ
    if(isLeptonic_W){

      Wcand_lv = Wlv; // W candidate is always the leptonic W here

      if(decayJets.at(0).first == 1){ // leading jet is top tagged
	validBDecay = 1;
	taggedTjet = 1;
	Tcand_lv = jets_lv.at(decayJets.at(0).second).first; // decayJets.second gives the jets_lv index to get 4-vec
      }
      else if(is2Jet && decayJets.at(1).first == 1){ // 2nd-leading jet is top tagged
	validBDecay = 0; // 0 for suboptimal validity
	taggedTjet = 1;
	Tcand_lv = jets_lv.at(decayJets.at(1).second).first;
      }
      else{ // leading jets are not top tagged
	if(!is2Jet){ 
	  validBDecay = -1; // all other 1-jet events are really mistagged
	  Tcand_lv = jets_lv.at(decayJets.at(0).second).first; //take the jet as the T candidate with a bad tag
	}
	else{ // it's a 2-jet event, will assume partially-merged top and add the two jets, looking for a W tag
	  if(decayJets.at(0).first == 2){ // W tag found 
	    validBDecay = 1;
	    taggedWbjetJet = 1;
	    Tcand_lv = jets_lv.at(decayJets.at(0).second).first+jets_lv.at(decayJets.at(1).second).first;
	  }
	  else if(decayJets.at(1).first == 2){ // W tag found 
	    validBDecay = 0; // sub-optimal
	    taggedWbjetJet = 1;
	    Tcand_lv = jets_lv.at(decayJets.at(0).second).first+jets_lv.at(decayJets.at(1).second).first;
	  }
	  else{ // any other set of tags
	    validBDecay = -1;
	    Tcand_lv = jets_lv.at(decayJets.at(0).second).first+jets_lv.at(decayJets.at(1).second).first;
	  }
	}
      }

      Bprime_DR = Wcand_lv.DeltaR(Tcand_lv);
      Bprime_ptbal = Wcand_lv.Pt()/Tcand_lv.Pt();
      Bprime_chi2 = pow(Tcand_lv.M() - 170.1,2)/(14.29*14.29) + pow(Bprime_DR - TMath::Pi(),2)/(0.2*0.2) + pow(Bprime_ptbal - 1,2)/(0.8*0.8);
      // hadronic t chi2 with values from 2016 B2G-17-018 (W isn't there because W_lv uses mass constraint)

    }else{ //isLeptonic_t

      Tcand_lv = top_lv; // T candidate is always the leptonic top here

      if(decayJets.at(0).first == 2){ // leading jet is W-tagged	    
	validBDecay = 1;
	taggedWjet = 1;
	Wcand_lv = jets_lv.at(decayJets.at(0).second).first;
      }else if(is2Jet && decayJets.at(1).first == 2){ // 2nd-leading jet is W-tagged
	validBDecay = 0; // sub-optimal
	taggedWjet = 1;
	Wcand_lv = jets_lv.at(decayJets.at(1).second).first;
      }else{
	validBDecay = -1;
	Wcand_lv = jets_lv.at(decayJets.at(0).second).first; // let's assume a 1-jet mistag scenario
      }

      Bprime_DR = Wcand_lv.DeltaR(Tcand_lv);
      Bprime_ptbal = Wcand_lv.Pt()/Tcand_lv.Pt();
      Bprime_chi2 = pow(Tcand_lv.M() - 172.6,2)/(19.1*19.1) + pow(Wcand_lv.M() - 85.5,2)/(8.7*8.7) + pow(Bprime_DR - TMath::Pi(),2)/(0.2*0.2) + pow(Bprime_ptbal - 1,2)/(0.8*0.8); // leptonic t, hadronic W chi2 with values from 2016 B2G-17-018

    }

    Bprime_lv = Wcand_lv + Tcand_lv;
    Bprime_mass = Bprime_lv.M();
    Bprime_pt = Bprime_lv.Pt();
    Bprime_eta = Bprime_lv.Eta();
    Bprime_phi = Bprime_lv.Phi();
  }

  // Include validBDecy_DeepAK8 variable and what it's decayed as 
  //leptonic and hadronic jetIDs can be slimmed out, no longer exist
  RVec<float> BPrimeVec = {Bprime_mass,Bprime_pt,Bprime_eta,Bprime_phi,Bprime_DR,Bprime_ptbal,Bprime_chi2,validBDecay,taggedWbjetJet,taggedTjet,taggedWjet};
  return BPrimeVec;
};


///////////////////////
//    Alternative    //
///////////////////////
RVec<float> BPrime_reco_alt(TLorentzVector lepton_lv, TLorentzVector top_lv, TLorentzVector Wlv, int leptonicParticle, RVec<float>& ak8_pt, RVec<float>& ak8_eta, RVec<float>& ak8_phi, RVec<float>& ak8_mass, RVec<int>& dpak8_tag, RVec<float>& ak8_sdmass)
{
  RVec<int> validBTagged (2,0);
  int val = -99;
  int tag = -99;
  TLorentzVector jet_lv;
  std::vector<pair<TLorentzVector,int>> jets_lv;
  float DR_LP_AK8 = 0;
  bool isLeptonic_W = false;
  bool isLeptonic_t = false;
  RVec<float> DR_lep(ak8_pt.size(), 9);
  int n_ak8 = ak8_pt.size();
     
  for(unsigned int ijet=0; ijet<n_ak8; ijet++) {
    jet_lv.SetPtEtaPhiM(ak8_pt.at(ijet),ak8_eta.at(ijet),ak8_phi.at(ijet),ak8_mass.at(ijet));
    DR_lep[ijet] = abs(jet_lv.DeltaR(lepton_lv) - 3);
  }
  DR_lep = Argsort(DR_lep);

  for(unsigned int i=0; i<n_ak8; i++){
    int ijet = DR_lep[i];
    
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
  float Bprime_pt = -999;
  float Bprime_eta = 9;
  float Bprime_phi = 9;
  float Bprime_DR = -9;
  float Bprime_ptbal = -999;
  float Bprime_chi2 = -999;
  
  TLorentzVector Wcand_lv;
  TLorentzVector Tcand_lv;
  TLorentzVector Bprime_lv;
  float validBDecay = -1;
  
  float taggedTjet = 0;
  float taggedTjet2nd = 0;
  float taggedWbjetJet = 0;
  float taggedWjet = 0;
  float taggedWjet2nd = 0;
  bool is2Jet = false;
  
  // ----------------------------------------------------------------------------
  // VLQ Decay -- 1 or 2 AK8 jets away from leptonic particle
  // ----------------------------------------------------------------------------
  if(jets_lv.size() > 2) {std::cout << "Problem: > 2 AK8s for Bprime reco" << std::endl;}
  else if(jets_lv.size() == 0){std::cout << "Failed reco event, no jets" << std::endl; }
  
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
    
    // Derived from Previous logic: TT -> tZ bW, BB -> tW bZ
    if(isLeptonic_W){

      Wcand_lv = Wlv; // W candidate is always the leptonic W here

      if(decayJets.at(0).first == 1){ // leading jet is top tagged
	validBDecay = 1;
	taggedTjet = 1;
	Tcand_lv = jets_lv.at(decayJets.at(0).second).first; // decayJets.second gives the jets_lv index to get 4-vec
      }
      else if(is2Jet && decayJets.at(1).first == 1){ // 2nd-leading jet is top tagged
	validBDecay = 0; // 0 for suboptimal validity
	taggedTjet = 1;
	Tcand_lv = jets_lv.at(decayJets.at(1).second).first;
      }
      else{ // leading jets are not top tagged
	if(!is2Jet){ 
	  validBDecay = -1; // all other 1-jet events are really mistagged
	  Tcand_lv = jets_lv.at(decayJets.at(0).second).first; //take the jet as the T candidate with a bad tag
	}
	else{ // it's a 2-jet event, will assume partially-merged top and add the two jets, looking for a W tag
	  if(decayJets.at(0).first == 2){ // W tag found 
	    validBDecay = 1;
	    taggedWbjetJet = 1;
	    Tcand_lv = jets_lv.at(decayJets.at(0).second).first+jets_lv.at(decayJets.at(1).second).first;
	  }
	  else if(decayJets.at(1).first == 2){ // W tag found 
	    validBDecay = 0; // sub-optimal
	    taggedWbjetJet = 1;
	    Tcand_lv = jets_lv.at(decayJets.at(0).second).first+jets_lv.at(decayJets.at(1).second).first;
	  }
	  else{ // any other set of tags
	    validBDecay = -1;
	    Tcand_lv = jets_lv.at(decayJets.at(0).second).first+jets_lv.at(decayJets.at(1).second).first;
	  }
	}
      }

      Bprime_DR = Wcand_lv.DeltaR(Tcand_lv);
      Bprime_ptbal = Wcand_lv.Pt()/Tcand_lv.Pt();
      Bprime_chi2 = pow(Tcand_lv.M() - 170.1,2)/(14.29*14.29) + pow(Bprime_DR - TMath::Pi(),2)/(0.2*0.2) + pow(Bprime_ptbal - 1,2)/(0.8*0.8);
      // hadronic t chi2 with values from 2016 B2G-17-018 (W isn't there because W_lv uses mass constraint)

    }else{ //isLeptonic_t

      Tcand_lv = top_lv; // T candidate is always the leptonic top here

      if(decayJets.at(0).first == 2){ // leading jet is W-tagged    
	validBDecay = 1;
	taggedWjet = 1;
	Wcand_lv = jets_lv.at(decayJets.at(0).second).first;
      }else if(is2Jet && decayJets.at(1).first == 2){ // 2nd-leading jet is W-tagged
	validBDecay = 0; // sub-optimal
	taggedWjet = 1;
	Wcand_lv = jets_lv.at(decayJets.at(1).second).first;
      }else{
	validBDecay = -1;
	Wcand_lv = jets_lv.at(decayJets.at(0).second).first; // let's assume a 1-jet mistag scenario
      }

      Bprime_DR = Wcand_lv.DeltaR(Tcand_lv);
      Bprime_ptbal = Wcand_lv.Pt()/Tcand_lv.Pt();
      Bprime_chi2 = pow(Tcand_lv.M() - 172.6,2)/(19.1*19.1) + pow(Wcand_lv.M() - 85.5,2)/(8.7*8.7) + pow(Bprime_DR - TMath::Pi(),2)/(0.2*0.2) + pow(Bprime_ptbal - 1,2)/(0.8*0.8); // leptonic t, hadronic W chi2 with values from 2016 B2G-17-018

    }

    Bprime_lv = Wcand_lv + Tcand_lv;
    Bprime_mass = Bprime_lv.M();
    Bprime_pt = Bprime_lv.Pt();
    Bprime_eta = Bprime_lv.Eta();
    Bprime_phi = Bprime_lv.Phi();
  }

  // Include validBDecy_DeepAK8 variable and what it's decayed as 
  //leptonic and hadronic jetIDs can be slimmed out, no longer exist
  RVec<float> BPrimeVec = {Bprime_mass,Bprime_pt,Bprime_eta,Bprime_phi,Bprime_DR,Bprime_ptbal,Bprime_chi2,validBDecay,taggedWbjetJet,taggedTjet,taggedWjet};
    return BPrimeVec;
};
