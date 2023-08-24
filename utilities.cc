// Methods in this file:
// leptonicCheck(), Electron_cutBasedIdNoIso_tight(), genttbarMassCalc(C), fVectorConstructor(), lvConstructor(C), DeltaR_VecAndFloat(), ptRel()
// I put some of the methods in here becasue I didn't know where else to put them.

using namespace std;
using namespace ROOT::VecOps;


// The following functions could probably all go to the plotting marco
auto leptonicCheck(string sample, int trueLeptonicT, int trueLeptonicW)
{
  if (sample.find("Bprime") == std::string::npos)
    {
      return -9;
    } // not sure if this line is needed. check.

  int trueLeptonicMode = -9;

  if ((trueLeptonicT != 1) && (trueLeptonicW == 1))
    {
      trueLeptonicMode = 0;
    } // leptonic W
  else if ((trueLeptonicT == 1) && (trueLeptonicW != 1))
    {
      trueLeptonicMode = 1;
    } // leptonic T
  else if ((trueLeptonicT == 1) && (trueLeptonicW == 1))
    {
      trueLeptonicMode = 2;
    } // dileptonic
  else if ((trueLeptonicT == 0) && (trueLeptonicW == 0))
    {
      trueLeptonicMode = -1;
    } // hadronic
  
  return trueLeptonicMode;
};

auto Electron_cutBasedIdNoIso_tight(string sample, unsigned int nElectron, RVec<int> &Electron_vidNestedWPBitmap)
{
    RVec<int> noIso_tight(nElectron, 0);
    for (unsigned int i = 0; i < nElectron; i++)
    {
        list<int> vars{0, 1, 2, 3, 4, 5, 6, 8, 9}; // checking this
        for (int x : vars)
        {
            if (((Electron_vidNestedWPBitmap[i] >> (x * 3)) & 0x7) >= 4)
            {
                noIso_tight[i] = 1;
            }
        }
    }
    return noIso_tight;
};


// -------------------------------------------
// 	  TLORENTZVECTOR CONSTRUCTOR
// -------------------------------------------
RVec<TLorentzVector> fVectorConstructor(RVec<float> &pt, RVec<float> &eta, RVec<float> &phi, RVec<float> &mass)
{
    RVec<TLorentzVector> lv;
    TLorentzVector tlv;
    for (int i = 0; i < pt.size(); i++)
    {
        tlv.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
        lv.push_back(tlv);
    }
    return lv;
};

// -------------------------------------------------
//    TLORENTZVECTOR CONSTRUCTOR FOR FLOAT ONLYS
// -------------------------------------------------
// Commented Method Only
TLorentzVector lvConstructor(float pt, float eta, float phi, float mass)
{
	TLorentzVector lv;
	lv.SetPtEtaPhiM(pt,eta,phi,mass);
	return lv;
};

auto getHighestPt(RVec<float> &VLepton_pt, RVec<float> &VLepton_eta, RVec<float> &VLepton_phi){
  RVec<float> Lepton0_PtEtaPhi(3,-1);

  if(VLepton_pt.size()==0){return Lepton0_PtEtaPhi;}

  auto maxIdx = ArgMax(VLepton_pt);

  Lepton0_PtEtaPhi[0] = VLepton_pt[maxIdx];
  Lepton0_PtEtaPhi[1] = VLepton_eta[maxIdx];
  Lepton0_PtEtaPhi[2] = VLepton_phi[maxIdx];

  return Lepton0_PtEtaPhi;
};

// --------------------------------------------
//               DR CALCULATOR
// --------------------------------------------
ROOT::VecOps::RVec<float> DeltaR_VecAndFloat(ROOT::VecOps::RVec<float>& jet_eta, ROOT::VecOps::RVec<float>& jet_phi, float& lep_eta, float& lep_phi)
{
  ROOT::VecOps::RVec<float> DR (jet_eta.size(),999);
  for(int i = 0; i < jet_eta.size(); i++) {DR[i] = DeltaR(jet_eta[i],lep_eta,jet_phi[i],lep_phi);}
  return DR;
};

ROOT::VecOps::RVec<float> ptRel(ROOT::VecOps::RVec<float>& jet_pt, ROOT::VecOps::RVec<float>& jet_eta, ROOT::VecOps::RVec<float>& jet_phi, ROOT::VecOps::RVec<float>& jet_mass, float& lepton_pt, float& lepton_eta, float& lepton_phi, float& lepton_mass)
{
  ROOT::VecOps::RVec<float> ptrel (jet_pt.size(),-1);
  TLorentzVector jet;
  TLorentzVector lepton;
  lepton.SetPtEtaPhiM(lepton_pt, lepton_eta, lepton_phi, lepton_mass);
  for(int i = 0; i < jet_pt.size(); i++) {
      jet.SetPtEtaPhiM(jet_pt[i], jet_eta[i], jet_phi[i], jet_mass[i]);
      ptrel[i] = (jet.Vect().Cross(lepton.Vect())).Mag() / jet.P();
  }
  return ptrel;
}
