// Methods in this file ---- all the methods in this file are in the commented section of analyzer_RDF.cc
// maxFxn(), JetDiscriminator(), DR_calc(), minDR_ptRel_lead_calc()

// --------------------------------------------------------
// 		     MAX VARIABLE FXN
// 	(Fxn to find max variable for dnnLargest)
// --------------------------------------------------------

// Commented fxn
ROOT::VecOps::RVec<int> maxFxn(ROOT::VecOps::RVec<float>& dnnJ, ROOT::VecOps::RVec<float>& dnnT, ROOT::VecOps::RVec<float>& dnnW)
{
	ROOT::VecOps::RVec<int> maxInt (dnnJ.size(),-1);
	for(int i = 0; i < dnnJ.size(); i++)
	{
		float maxVar = -999;
		if(maxVar < dnnJ[i]){maxVar = dnnJ[i]; maxInt[i] = 0;}
		if(maxVar < dnnT[i]){maxVar = dnnT[i]; maxInt[i] = 1;}
		if(maxVar < dnnW[i]){maxVar = dnnW[i]; maxInt[i] = 2;}
	}
	return maxInt;
};

// Commented fxn
ROOT::VecOps::RVec<int> JetDiscriminator(ROOT::VecOps::RVec<float>& dnnT, ROOT::VecOps::RVec<float>& dnnW){
  int nJets = dnnT.size();
  ROOT::VecOps::RVec<int> tag (nJets, -1);
  
  for(int i=0; i<nJets; i++){
    if(dnnT[i] > 0.58){tag[i] = 1;} // 1% bkg ala https://twiki.cern.ch/twiki/bin/view/CMS/ParticleNetSFs
    else if(dnnW[i] > 0.70){tag[i] = 2;} // TOO TIGHT at 0.94 1%, reducing to the 0.70 5% bkgd
    else{tag[i] = 0;}
  }
  return tag;
};

// --------------------------------------------
// 		 DR CALCULATOR
// --------------------------------------------
// Commented fxn
ROOT::VecOps::RVec<float> DR_calc(ROOT::VecOps::RVec<float>& jet_pt, ROOT::VecOps::RVec<float>& jet_eta, ROOT::VecOps::RVec<float>& jet_phi, ROOT::VecOps::RVec<float>& jet_mass, float lep_pt, float lep_eta, float lep_phi, float lep_mass)
{
	ROOT::VecOps::RVec<float> DR (jet_pt.size(),0);
	for(int i = 0; i < jet_pt.size(); i++) {DR[i] = DeltaR(jet_eta[i],lep_eta,jet_phi[i],lep_phi);}
	return DR;
};

// --------------------------------------------
// 	   MINDR & PTREL CALCULATOR
// --------------------------------------------
// Commented fxn
ROOT::VecOps::RVec<float> minDR_ptRel_lead_calc(ROOT::VecOps::RVec<float>& jet_pt, ROOT::VecOps::RVec<float>& jet_eta, ROOT::VecOps::RVec<float>& jet_phi, ROOT::VecOps::RVec<float>& jet_mass, TLorentzVector lepton_lv)
{
	TLorentzVector jet_TLV, leadJet;
	float deltaR_lepJets = 0;
	float ptRel_lepJets = 0;
	float minDR_lepJets = 1000;
	float minDR_leadJetotherJet = 1000;
	leadJet.SetPtEtaPhiM(jet_pt[0],jet_eta[0],jet_phi[0],jet_mass[0]);
	if(jet_pt.size() < 1) {minDR_lepJets = -99.0;}
	if(jet_pt.size() < 2) {minDR_leadJetotherJet = -99.0;}
	for(int i = 0; i < jet_pt.size(); i++)
	{
		jet_TLV.SetPtEtaPhiM(jet_pt[i],jet_eta[i],jet_phi[i],jet_mass[i]);
		deltaR_lepJets = lepton_lv.DeltaR(jet_TLV);
		if(deltaR_lepJets < minDR_lepJets)
		{
			minDR_lepJets = lepton_lv.DeltaR(jet_TLV);
			ptRel_lepJets = lepton_lv.P()*(jet_TLV.Vect().Cross(lepton_lv.Vect()).Mag()/jet_TLV.P()/lepton_lv.P());
		}
		if(i > 0)
		{
			float tempDR = leadJet.DeltaR(jet_TLV);
			if(tempDR < minDR_leadJetotherJet){minDR_leadJetotherJet = tempDR;}
		}
	}
	ROOT::VecOps::RVec<float> returnVec = {minDR_lepJets,ptRel_lepJets,minDR_leadJetotherJet};
	return returnVec;
};
