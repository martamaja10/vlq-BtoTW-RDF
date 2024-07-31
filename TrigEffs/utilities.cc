// Methods in this file:
// assign_leps(C), cleanJets()

using namespace ROOT::VecOps;
using namespace std;

RVec<float> DeltaR_VecAndFloat(RVec<float>& jet_eta, RVec<float>& jet_phi, const float& lep_eta, const float& lep_phi)
{
  RVec<float> DR (jet_eta.size(),999);
  for(int i = 0; i < jet_eta.size(); i++) { DR[i] = DeltaR(jet_eta[i],lep_eta,jet_phi[i],lep_phi); }
  return DR;
};

RVec<float> ptRel(RVec<float>& jet_pt, RVec<float>& jet_eta, RVec<float>& jet_phi, RVec<float>& jet_mass, const float& lepton_pt, const float& lepton_eta, const float& lepton_phi, const float& lepton_mass)
{
  RVec<float> ptrel (jet_pt.size(),-1);
  TLorentzVector jet;
  TLorentzVector lepton;
  lepton.SetPtEtaPhiM(lepton_pt, lepton_eta, lepton_phi, lepton_mass);
  for(int i = 0; i < jet_pt.size(); i++) {
      jet.SetPtEtaPhiM(jet_pt[i], jet_eta[i], jet_phi[i], jet_mass[i]);
      ptrel[i] = (jet.Vect().Cross(lepton.Vect())).Mag() / jet.P();
  }
  return ptrel;
};

RVec<float> assign_leps(bool isMu, bool isEl, RVec<int> &TPassMu, RVec<int> &TPassEl, RVec<float> &Muon_pt, RVec<float> &Muon_eta, RVec<float> &Muon_phi, RVec<float> &Muon_mass, RVec<int> &Muon_jetIdx, RVec<float> &Electron_pt, RVec<float> &Electron_eta, RVec<float> &Electron_phi, RVec<float> &Electron_mass, RVec<int> &Electron_jetIdx)
{

  float lep_pt = -9;
  float lep_eta = -9;
  float lep_phi = -9;
  float lep_mass = -9;
  float lep_jetidx = -9;

  if(isMu){
    for(unsigned int imu=0; imu < Muon_pt.size(); imu++) {
      if(TPassMu.at(imu) == 1){
	if (lep_pt > -1) cout << "Problem: found two muons with TPassMu = 1" << endl;
	lep_pt = Muon_pt.at(imu);
	lep_eta = Muon_eta.at(imu);
	lep_phi = Muon_phi.at(imu);
	lep_mass = Muon_mass.at(imu);
	lep_jetidx = Muon_jetIdx.at(imu);
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
	lep_jetidx = Electron_jetIdx.at(iel);
      }
    }
  }
  
  RVec<float> lepVec = {lep_pt,lep_eta,lep_phi,lep_mass,lep_jetidx};
  return lepVec;
};

RVec<int> drop_jets(const RVec<int> &goodjets, const float &probejet, const float &tagjet)
{
  RVec<int> passfail = goodjets;
  if(probejet > -1) passfail.at(probejet) = 0;
  if(tagjet > -1) passfail.at(tagjet) = 0;
  return passfail;
};

auto Electron_cutBasedIdNoIso_tight(unsigned int nElectron, RVec<int> &Electron_vidNestedWPBitmap)
{
  RVec<int> noIso_tight(nElectron, 0);
  for (unsigned int i = 0; i < nElectron; i++)
    {
      noIso_tight[i] = 1;    
      list<int> vars{0, 1, 2, 3, 4, 5, 6, 8, 9}; // checking this
      for (int x : vars)
        {
	  if (((Electron_vidNestedWPBitmap[i] >> (x * 3)) & 0x7) < 4)
            {
	      noIso_tight[i] = 0;
            }
        }
    }
  return noIso_tight;
};
