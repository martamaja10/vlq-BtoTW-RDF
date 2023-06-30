// Methods in this file:
// cut_ptrel()

// ----------------------------------------------------------------------------------------------
//     Calculate pass/fail for (minDRAK4 > 0.4 OR ptRelAK4 > 25) on lepton candidates
// ----------------------------------------------------------------------------------------------

using namespace std;
using namespace ROOT::VecOps;

auto cut_ptrel(double dR_LIM_AK4, double ptrel_LIM, RVec<TLorentzVector> leptons, RVec<double> gcJet_eta, RVec<double> gcJet_phi, RVec<double> gcJet_pt, RVec<double> gcJet_mass)
{
    // First, loop through all the muons.  We do each one seperately
  int Nlep = leptons.size();
  RVec<int> passOrFail(Nlep,1);

  int NgcJets = gcJet_eta.size();
  if(NgcJets==0){return passOrFail;}

  for (unsigned int i=0; i<Nlep; i++)
    {
      // Now, call the built in DeltaR function with a for loop to calculate it individually for each jet
      RVec<double> dr;
      RVec<double> leptonEta(NgcJets, leptons[i].Eta());
      RVec<double> leptonPhi(NgcJets, leptons[i].Phi());

      dr = DeltaR(leptonEta, gcJet_eta, leptonPhi, gcJet_phi);
/*
        double dr_temp;
        for (int i = 0; i < gcforwJet_eta.size(); i++)
        {
            // Calculate the dr for the jet compared to lepton and save to RVec<float> dr
            dr_temp = DeltaR(lepton.Eta(), gcforwJet_eta[i], lepton.Phi(), gcforwJet_phi[i]);
            dr.push_back(dr_temp);
        }
*/
      // Find the index of the smallest value in dr, and use it to save the minimum value to minDR
      auto minIndex = ArgMin(dr);
      auto minDR = dr[minIndex];
		
      TLorentzVector jet;
      jet.SetPtEtaPhiM(gcJet_pt[minIndex], gcJet_eta[minIndex], gcJet_phi[minIndex], gcJet_mass[minIndex]);
      auto ptRel = (jet.Vect().Cross(leptons[i].Vect())).Mag() / jet.P();
      
      if (minDR > dR_LIM_AK4 || ptRel > ptrel_LIM)
        {
	  passOrFail[i] = 1;
        }
      else
        {
	  passOrFail[i] = 0;
        }
    }
  return passOrFail;
}
