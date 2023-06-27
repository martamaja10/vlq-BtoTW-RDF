// ----------------------------------------------------------------------------------------------
//     Calculate pass/fail for (minDRAK4 > 0.4 OR ptRelAK4 > 25) on lepton candidates
// ----------------------------------------------------------------------------------------------

auto cut_ptrel(double dR_LIM_AK4, double ptrel, RVec<TLorentzVector> leptons, double NJets_forward, RVec<double> gcforwJet_eta, RVec<double> gcforwJet_phi, RVec<double> gcforwJet_pt, RVec<double> gcforwJet_mass)
{
    // First, loop through all the muons.  We do each one seperately
    ROOT::RVec<int> passOrFail;
    //std::cout << "Number of Jets: " << NJets_forward << std::endl;

    for (const TLorentzVector &lepton : leptons)
    {
        // Now, call the built in DeltaR function with a for loop to calculate it individually for each jet
        ROOT::RVec<double> dr;
        double dr_temp;
        //std::cout << "-------New Lepton-------" << std::endl;
        
        for (int i = 0; i < gcforwJet_eta.size(); i++)
        {
            // Calculate the dr for the jet compared to lepton and save to ROOT::RVec<float> dr
            dr_temp = ROOT::VecOps::DeltaR(lepton.Eta(), gcforwJet_eta[i], lepton.Phi(), gcforwJet_phi[i]);
            //std::cout << "dr_temp: " << dr_temp << std::endl;
            dr.push_back(dr_temp);
        }
        
        // Find the index of the smallest value in dr, and use it to save the minimum value to minDR
        auto minIndex = ROOT::VecOps::ArgMin(dr);
        //std::cout << "Min Index: " << minIndex << std::endl;
        auto minDR = dr[0];
        //std::cout << "Min DR: " << minDR << std::endl;
        
        TLorentzVector jet;
        jet.SetPtEtaPhiM(gcforwJet_pt[minIndex], gcforwJet_eta[minIndex], gcforwJet_phi[minIndex], gcforwJet_mass[minIndex]); 
        auto ptRel = (jet.Vect().Cross(lepton.Vect())).Mag() / jet.P();
        if (minDR > dR_LIM_AK4 || ptRel > ptrel)
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
