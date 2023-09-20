// Methods in this file:
// assign_leps(C), cleanJets()

// --------------------------------------------------------
// 		    JET CLEANING FXN
// --------------------------------------------------------

using namespace ROOT::VecOps;
using namespace std;

//  correction::CompoundCorrection::Ref ak4corr;
//  correction::Correction::Ref ak4corrUnc;
//  correction::Correction::Ref ak4ptres;
//  correction::Correction::Ref ak4jer;


// Commented Method Only
RVec<float> assign_leps(bool isMu, bool isEl, RVec<int> &TPassMu, RVec<int> &TPassEl, RVec<float> &Muon_pt, RVec<float> &Muon_eta, RVec<float> &Muon_phi, RVec<float> &Muon_mass, RVec<float> &Muon_miniPFRelIso_all, RVec<float> &Electron_pt, RVec<float> &Electron_eta, RVec<float> &Electron_phi, RVec<float> &Electron_mass, RVec<float> &Electron_miniPFRelIso_all)
{

  float lep_pt = -9;
  float lep_eta = -9;
  float lep_phi = -9; 
  float lep_mass = -9;
  float lep_miniIso = -9;

  if(isMu){
    for(unsigned int imu=0; imu < Muon_pt.size(); imu++) {
      if(TPassMu.at(imu) == 1){
	      if (lep_pt > -1) cout << "Problem: found two muons with TPassMu = 1" << endl;
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
	      if (lep_pt > -1) cout << "Problem: found two electrons with TPassEl = 1" << endl;
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

