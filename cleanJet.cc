// --------------------------------------------------------
// 		    JET CLEANING FXN
// --------------------------------------------------------
//auto cleanJets = [](ROOT::VecOps::RVec<float>& jt_pt, ROOT::VecOps::RVec<float>& jt_mass, ROOT::VecOps::RVec<int>& jt_id, ROOT::VecOps::RVec<float>& jt_eta, ROOT::VecOps::RVec<float>& jt_phi, ROOT::VecOps::RVec<float>& lep_pt, ROOT::VecOps::RVec<float>& lep_mass, ROOT::VecOps::RVec<float>& lep_eta, ROOT::VecOps::RVec<float>& lep_phi, float dR_LIM)

ROOT::VecOps::RVec<float> assign_leps(bool isMu, bool isEl, ROOT::VecOps::RVec<int>& TPassMu,ROOT::VecOps::RVec<int>& TPassEl,ROOT::VecOps::RVec<float>& Muon_pt,ROOT::VecOps::RVec<float>& Muon_eta,ROOT::VecOps::RVec<float>& Muon_phi,ROOT::VecOps::RVec<float>& Muon_mass,ROOT::VecOps::RVec<float>& Muon_miniPFRelIso_all,ROOT::VecOps::RVec<float>& Electron_pt,ROOT::VecOps::RVec<float>& Electron_eta,ROOT::VecOps::RVec<float>& Electron_phi,ROOT::VecOps::RVec<float>& Electron_mass,ROOT::VecOps::RVec<float>& Electron_miniPFRelIso_all){

  float lep_pt = -9;
  float lep_eta = -9;
  float lep_phi = -9; 
  float lep_mass = -9;
  float lep_miniIso = -9;

  if(isMu){
    for(unsigned int imu=0; imu < Muon_pt.size(); imu++) {
      if(TPassMu.at(imu) == 1){
	if (lep_pt > -1) std::cout << "Problem: found two muons with TPassMu = 1" << std::endl;
	lep_pt = Muon_pt.at(imu);
	lep_eta = Muon_eta.at(imu);
	lep_phi = Muon_phi.at(imu);
	lep_mass = Muon_mass.at(imu);
	lep_miniIso = Muon_miniPFRelIso_all.at(imu);
      }
    }
  }else if(isEl){
    for(unsigned int iel=0; iel < Electron_pt.size(); iel++) {
      if(TPassEl.at(iel) == 1){
	if (lep_pt > -1) std::cout << "Problem: found two electrons with TPassEl = 1" << std::endl;
	lep_pt = Electron_pt.at(iel);
	lep_eta = Electron_eta.at(iel);
	lep_phi = Electron_phi.at(iel);
	lep_mass = Electron_mass.at(iel);
	lep_miniIso = Electron_miniPFRelIso_all.at(iel);
      }
    }
  }
  
  RVec<float> lepVec = {lep_pt,lep_eta,lep_phi,lep_mass,lep_miniIso};
  return lepVec;

}


ROOT::VecOps::RVec<float> cleanJets(ROOT::VecOps::RVec<float>& jt_pt, ROOT::VecOps::RVec<float>& jt_mass, ROOT::VecOps::RVec<int>& jt_id, ROOT::VecOps::RVec<float>& jt_eta, ROOT::VecOps::RVec<float>& jt_phi, float lep_pt, float lep_mass, float lep_eta, float lep_phi, float dR_LIM)
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
		Leptons.SetPtEtaPhiM(lep_pt,lep_eta,lep_phi,lep_mass);
		dR[i] = DeltaR(jt_eta[i],lep_eta,jt_phi[i],lep_phi);
		pt_rel[i] = (Jets.Vect().Cross(Leptons.Vect())).Mag()/Jets.P();
		isClean = true;
		if(dR[i] < dR_LIM && pt_rel[i] < 25){isClean = false;}
		if(isClean == false){continue;}
		if(isClean == true){cleanJets_[i] = 1;}
	}
	return cleanJets_;

	// FIXME: this is just the 2D cut. Try Electron/Muon/Jet_cleanmask, try Electron/Muon_jetidx, try Electron/Muon_jetPtRelv2, 
	// Can I implement the adjustment and re-JEC of jets? 
	// Better to just add it to CRAB job when getting the JEC unc?
	// Leptonic W events...this should be fine
	// Leptonic t events...is the ptRel enough to keep a b-jet?
};
