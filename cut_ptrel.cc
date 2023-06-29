// Methods in this file:
// cut_ptrel()

// ----------------------------------------------------------------------------------------------
//     Calculate pass/fail for (minDRAK4 > 0.4 OR ptRelAK4 > 25) on lepton candidates
// ----------------------------------------------------------------------------------------------

using namespace std;
using namespace ROOT::VecOps;

auto cut_ptrel(double dR_LIM_AK4, double ptrel_LIM, RVec<TLorentzVector> leptons, double NJets_forward, RVec<double> gcforwJet_eta, RVec<double> gcforwJet_phi, RVec<double> gcforwJet_pt, RVec<double> gcforwJet_mass)
{
    // First, loop through all the muons.  We do each one seperately
    RVec<int> passOrFail;

    for (const TLorentzVector &lepton : leptons)
    {
        // Now, call the built in DeltaR function with a for loop to calculate it individually for each jet
        RVec<double> dr;
        RVec<double> leptonEta(gcforwJet_eta.size(), lepton.Eta());
        RVec<double> leptonPhi(gcforwJet_phi.size(), lepton.Phi());

        dr = DeltaR(leptonEta, gcforwJet_eta, leptonPhi, gcforwJet_phi);
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
        jet.SetPtEtaPhiM(gcforwJet_pt[minIndex], gcforwJet_eta[minIndex], gcforwJet_phi[minIndex], gcforwJet_mass[minIndex]);
        auto ptRel = (jet.Vect().Cross(lepton.Vect())).Mag() / jet.P();
        if (minDR > dR_LIM_AK4 || ptRel > ptrel_LIM)
        {
            passOrFail.push_back(1);
        }
        else
        {
            passOrFail.push_back(0);
        }
    }
    return passOrFail;
}
