// --------------------------------------------------------
// 		     MAX VARIABLE FXN
// 	(Fxn to find max variable for dnnLargest)
// --------------------------------------------------------

ROOT::VecOps::RVec<int> maxFxn(ROOT::VecOps::RVec<float>& dnnJ, ROOT::VecOps::RVec<float>& dnnT, ROOT::VecOps::RVec<float>& dnnH, ROOT::VecOps::RVec<float>& dnnZ, ROOT::VecOps::RVec<float>& dnnW, ROOT::VecOps::RVec<float>& dnnB)
{
	ROOT::VecOps::RVec<int> maxInt (dnnJ.size(),-1);
	for(int i = 0; i < dnnJ.size(); i++)
	{
		float maxVar = -999;
		if(maxVar < dnnJ[i]){maxVar = dnnJ[i]; maxInt[i] = 0;}
		if(maxVar < dnnT[i]){maxVar = dnnT[i]; maxInt[i] = 1;}
		if(maxVar < dnnH[i]){maxVar = dnnH[i]; maxInt[i] = 2;}
		if(maxVar < dnnZ[i]){maxVar = dnnZ[i]; maxInt[i] = 3;}
		if(maxVar < dnnW[i]){maxVar = dnnW[i]; maxInt[i] = 4;}
		if(maxVar < dnnB[i]){maxVar = dnnB[i]; maxInt[i] = 5;}
	}
	return maxInt;
};

// -------------------------------------------
// 	  TLORENTZVECTOR CONSTRUCTOR
// -------------------------------------------
TLorentzVector fVectorConstructor(ROOT::VecOps::RVec<float>& pt, ROOT::VecOps::RVec<float>& eta, ROOT::VecOps::RVec<float>& phi, ROOT::VecOps::RVec<float>& mass)
{
	TLorentzVector lv;
	for(int i = 0; i < pt.size(); i++){lv.SetPtEtaPhiM(pt[i],eta[i],phi[i],mass[i]);}
	return lv;
};

// --------------------------------------------
// 		 DR CALCULATOR
// --------------------------------------------
ROOT::VecOps::RVec<float> DR_calc(ROOT::VecOps::RVec<float>& jet_pt, ROOT::VecOps::RVec<float>& jet_eta, ROOT::VecOps::RVec<float>& jet_phi, ROOT::VecOps::RVec<float>& jet_mass, ROOT::VecOps::RVec<float>& lep_pt, ROOT::VecOps::RVec<float>& lep_eta, ROOT::VecOps::RVec<float>& lep_phi, ROOT::VecOps::RVec<float>& lep_mass)
{
	ROOT::VecOps::RVec<float> DR (jet_pt.size(),0);
	for(int i = 0; i < jet_pt.size(); i++) {DR[i] = DeltaR(jet_eta[i],lep_eta[0],jet_phi[i],lep_phi[0]);}
	return DR;
};

// --------------------------------------------
// 	   MINDR & PTREL CALCULATOR
// --------------------------------------------
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
