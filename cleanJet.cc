// --------------------------------------------------------
// 		    JET CLEANING FXN
// --------------------------------------------------------
//auto cleanJets = [](ROOT::VecOps::RVec<float>& jt_pt, ROOT::VecOps::RVec<float>& jt_mass, ROOT::VecOps::RVec<int>& jt_id, ROOT::VecOps::RVec<float>& jt_eta, ROOT::VecOps::RVec<float>& jt_phi, ROOT::VecOps::RVec<float>& lep_pt, ROOT::VecOps::RVec<float>& lep_mass, ROOT::VecOps::RVec<float>& lep_eta, ROOT::VecOps::RVec<float>& lep_phi, float dR_LIM)
ROOT::VecOps::RVec<float> cleanJets(ROOT::VecOps::RVec<float>& jt_pt, ROOT::VecOps::RVec<float>& jt_mass, ROOT::VecOps::RVec<int>& jt_id, ROOT::VecOps::RVec<float>& jt_eta, ROOT::VecOps::RVec<float>& jt_phi, ROOT::VecOps::RVec<float>& lep_pt, ROOT::VecOps::RVec<float>& lep_mass, ROOT::VecOps::RVec<float>& lep_eta, ROOT::VecOps::RVec<float>& lep_phi, float dR_LIM)
{
	ROOT::VecOps::RVec<float> cleanJets_ (jt_id.size(),0);
	ROOT::VecOps::RVec<float> dR (jt_id.size(),0);
	ROOT::VecOps::RVec<float> pt_rel (jt_id.size(),0);
	auto isClean = true;
	int j = 0;
	for(auto &i: Nonzero(jt_id))
	{
		TLorentzVector Jets, Leptons;
		Jets.SetPtEtaPhiM(jt_pt[i],jt_eta[i],jt_phi[i],jt_mass[i]);
		Leptons.SetPtEtaPhiM(lep_pt[0],lep_eta[0],lep_phi[0],lep_mass[0]);
		dR[i] = DeltaR(jt_eta[i],lep_eta[0],jt_phi[i],lep_phi[0]);
		pt_rel[i] = (Jets.Vect().Cross(Leptons.Vect())).Mag()/Jets.P();
		isClean = true;
		if(dR[i] < dR_LIM && pt_rel[i] < 20){isClean = false;}
		if(isClean == false){continue;}
		if(isClean == true){cleanJets_[i] = 1;}
	}
	return cleanJets_;
};
