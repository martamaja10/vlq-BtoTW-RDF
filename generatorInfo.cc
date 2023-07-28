// Methods in this file:
// Bprime_gen_info(), t_gen_info(), W_gen_info(), t_bkg_idx(), W_bkg_idx(), FatJet_matching_sig(C), FatJet_matching_bkg(C)

// ----------------------------------------------------
//           Bprime truth extraction:
// ----------------------------------------------------

using namespace std;
using namespace ROOT::VecOps;

auto Bprime_gen_info(string sample, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_mass, RVec<float> &GenPart_pt, RVec<float> &GenPart_phi, RVec<float> &GenPart_eta, RVec<int> &GenPart_genPartIdxMother, RVec<int> &GenPart_status, RVec<int> &GenPart_statusFlags)
{
  RVec<float> BPrimeInfo(6, -999);
  if (sample.find("Bprime") == std::string::npos)
  {
    return BPrimeInfo;
  }

  for (unsigned int i = 0; i < nGenPart; i++)
  {
    int id = GenPart_pdgId[i];
    if (abs(id) != 6000007)
    {
      continue;
    }

    bitset<15> statusFlags(GenPart_statusFlags[i]);
    if (statusFlags.to_string()[1] == '0')
    {
      continue;
    } // takes the last B'
    BPrimeInfo[0] = GenPart_pt[i];
    BPrimeInfo[1] = GenPart_eta[i];
    BPrimeInfo[2] = GenPart_phi[i];
    BPrimeInfo[3] = GenPart_mass[i];
    BPrimeInfo[4] = GenPart_pdgId[i];
    BPrimeInfo[5] = GenPart_status[i];
  }

  return BPrimeInfo; // if entries -999, then no Bprime was found
};

// ----------------------------------------------------
//           t truth extraction:
// ----------------------------------------------------
auto t_gen_info(string sample, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_mass, RVec<float> &GenPart_pt, RVec<float> &GenPart_phi, RVec<float> &GenPart_eta, RVec<int> &GenPart_genPartIdxMother, RVec<int> &GenPart_status)
{
  RVec<float> t_gen_info(30, -999);
  if (sample.find("Bprime") == std::string::npos)
  {
    return t_gen_info;
  }

  int trueLeptonicT = -1;
  for (unsigned int i = 0; i < nGenPart; i++)
  {
    int id = GenPart_pdgId[i];
    int motherIdx = GenPart_genPartIdxMother[i];

    if (abs(GenPart_pdgId[motherIdx]) != 6)
    {
      continue;
    } // find t daughters
    if (abs(id) != 24 && abs(id) != 5)
    {
      continue;
    }

    // store t info
    t_gen_info[0] = GenPart_pt[motherIdx];
    t_gen_info[1] = GenPart_eta[motherIdx];
    t_gen_info[2] = GenPart_phi[motherIdx];
    t_gen_info[3] = GenPart_mass[motherIdx];
    t_gen_info[4] = GenPart_pdgId[motherIdx];
    t_gen_info[5] = GenPart_status[motherIdx];

    int igen = i;
    for (unsigned int j = i; j < nGenPart; j++)
    {
      if (GenPart_pdgId[j] != id)
      {
        continue;
      }
      if (GenPart_genPartIdxMother[j] != igen)
      {
        continue;
      }
      igen = j; // take the last copy of t daughter
    }

    if (abs(id) == 5)
    { // store b info
      t_gen_info[6] = GenPart_pt[igen];
      t_gen_info[7] = GenPart_eta[igen];
      t_gen_info[8] = GenPart_phi[igen]; // did not record gen mass, because =0 for all b
      t_gen_info[9] = GenPart_pdgId[igen];
      t_gen_info[10] = GenPart_status[igen];
    }
    else
    { // store W info
      t_gen_info[11] = GenPart_pt[igen];
      t_gen_info[12] = GenPart_eta[igen];
      t_gen_info[13] = GenPart_phi[igen];
      t_gen_info[14] = GenPart_mass[igen];
      t_gen_info[15] = GenPart_pdgId[igen];
      t_gen_info[16] = GenPart_status[igen];
      for (unsigned int j = igen; j < nGenPart; j++)
      {
        if (GenPart_genPartIdxMother[j] != igen)
        {
          continue;
        } // look for W daughters
        int j_id = GenPart_pdgId[j];

        int jgen = j;
        for (unsigned int k = j; k < nGenPart; k++)
        {
          if (GenPart_pdgId[k] != j_id)
          {
            continue;
          }
          if (GenPart_genPartIdxMother[k] != j_id)
          {
            continue;
          }
          jgen = k; // take the last copy of W daughter
        }

        int n = 0;
        if (abs(j_id) == 11 || abs(j_id) == 13 || abs(j_id) == 15)
        {
          trueLeptonicT = 1;
        } // store e/mu/tau first
        else if (abs(j_id) == 12 || abs(j_id) == 14 || abs(j_id) == 16)
        {
          trueLeptonicT = 1;
          n = 6;
        } // then neutrinos
        else if (trueLeptonicT == -1)
        {
          trueLeptonicT = 0;
        } // quark 1
        else if (trueLeptonicT == 0)
        {
          n = 6;
        } // quark 2
        else
        {
          cout << "Error in t_gen_info: not  e/mu/tau, neutrino, or quark" << endl;
        }

        t_gen_info[17 + n] = GenPart_pt[jgen];
        t_gen_info[18 + n] = GenPart_eta[jgen];
        t_gen_info[19 + n] = GenPart_phi[jgen];
        t_gen_info[20 + n] = GenPart_mass[jgen];
        t_gen_info[21 + n] = GenPart_pdgId[jgen];
        t_gen_info[22 + n] = GenPart_status[jgen];
      }
    }
  }
  t_gen_info[29] = trueLeptonicT;

  return t_gen_info;
};

// ----------------------------------------------------
//           W truth extraction:
// ----------------------------------------------------
auto W_gen_info(string sample, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_mass, RVec<float> &GenPart_pt, RVec<float> &GenPart_phi, RVec<float> &GenPart_eta, RVec<int> &GenPart_genPartIdxMother, RVec<int> &GenPart_status, int daughterW_gen_pdgId)
{
  RVec<float> W_gen_info(19, -999);
  if (sample.find("Bprime") == std::string::npos)
  {
    return W_gen_info;
  }
  int trueLeptonicW = -1;

  for (unsigned int i = 0; i < nGenPart; i++)
  {
    int id = GenPart_pdgId[i];
    int motherIdx = GenPart_genPartIdxMother[i];

    if (abs(id) > 17 || GenPart_pdgId[motherIdx] != (-daughterW_gen_pdgId))
    {
      continue;
    } // look for daughters of W

    if (trueLeptonicW == -1)
    {
      W_gen_info[0] = GenPart_pt[motherIdx];
      W_gen_info[1] = GenPart_eta[motherIdx];
      W_gen_info[2] = GenPart_phi[motherIdx];
      W_gen_info[3] = GenPart_mass[motherIdx];
      W_gen_info[4] = GenPart_pdgId[motherIdx];
      W_gen_info[5] = GenPart_status[motherIdx];
    }

    int igen = i;
    for (unsigned int j = igen; j < nGenPart; j++)
    {
      if (GenPart_pdgId[j] != id)
      {
        continue;
      }
      if (GenPart_genPartIdxMother[j] != igen)
      {
        continue;
      }
      igen = j; // take the last copy of W daughter
    }

    int n = 0;
    if (abs(id) == 11 || abs(id) == 13 || abs(id) == 15)
    {
      trueLeptonicW = 1;
    } // store e/mu/tau first
    else if (abs(id) == 12 || abs(id) == 14 || abs(id) == 16)
    {
      trueLeptonicW = 1;
      n = 6;
    } // then neutrinos
    else if (trueLeptonicW == -1)
    {
      trueLeptonicW = 0;
    } // quark 1
    else if (trueLeptonicW == 0)
    {
      n = 6;
    } // quark 2
    else
    {
      cout << "Error in W_gen_info: not  e/mu/tau, neutrino, or quark" << endl;
    }

    W_gen_info[6 + n] = GenPart_pt[i];
    W_gen_info[7 + n] = GenPart_eta[i];
    W_gen_info[8 + n] = GenPart_phi[i];
    W_gen_info[9 + n] = GenPart_mass[i];
    W_gen_info[10 + n] = GenPart_pdgId[i];
    W_gen_info[11 + n] = GenPart_status[i];
  }
  W_gen_info[18] = trueLeptonicW;

  return W_gen_info;
};

// ----------------------------------------------------
//           W,t truth extraction for bkg:
// ----------------------------------------------------

auto t_bkg_idx(string sample, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<int> &GenPart_genPartIdxMother, RVec<int> &GenPart_statusFlags)
{
  if (sample.find("Bprime") != std::string::npos)
  {
    RVec<int> t_daughter_idx;
    return t_daughter_idx;
  }

  RVec<int> t_idx;
  for (unsigned int i = 0; i < nGenPart; i++)
  {
    if (abs(GenPart_pdgId[i]) != 6)
    {
      continue;
    }
    bitset<15> statusFlags(GenPart_statusFlags[i]);

    if (statusFlags.to_string()[1] == '0')
    {
      continue;
    } // take last copy of t
    t_idx.push_back(i);
  }

  int Nt = t_idx.size();
  RVec<int> t_daughter_idx(Nt * 3, -99);

  for (unsigned int i = 0; i < Nt; i++)
  {
    for (unsigned int j = t_idx[i]; j < nGenPart; j++)
    {
      if (GenPart_genPartIdxMother[j] != t_idx[i])
      {
        continue;
      } // pick out daughters of t

      int id = GenPart_pdgId[j];
      if (abs(id) != 5 && abs(id) != 24)
      {
        continue;
      } // pick out daughter b, W

      if (abs(id) == 5)
      {
        t_daughter_idx[i * 3] = j;
      } // record the first copy of b
      else
      {
        int jgen = j;
        for (unsigned int k = j; k < nGenPart; k++)
        {
          if (GenPart_pdgId[k] != id)
          {
            continue;
          }
          if (GenPart_genPartIdxMother[k] != jgen)
          {
            continue;
          }
          jgen = k; // take the last copy of W
        }

        int n = 1;
        for (unsigned int k = j; k < nGenPart; k++)
        {
          if (GenPart_genPartIdxMother[k] != jgen)
          {
            continue;
          } // pick out daughters of W
          if (abs(GenPart_pdgId[k]) > 17)
          {
            continue;
          }                              // to exclude 24->22,24
          t_daughter_idx[i * 3 + n] = k; // record the first copy of W daughter
          n += 1;
        }
      }
    }
  }
  return t_daughter_idx;
};

auto W_bkg_idx(string sample, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<int> &GenPart_genPartIdxMother, RVec<int> &GenPart_statusFlags, RVec<int> &t_bkg_idx)
{
  if (sample.find("Bprime") != std::string::npos)
  {
    RVec<int> W_daughter_idx;
    return W_daughter_idx;
  }
  //cout << "Event" << endl;
  RVec<int> W_idx;
  for (unsigned int i = 0; i < nGenPart; i++)
  {
    if (abs(GenPart_pdgId[i]) != 24)
    {
      continue;
    }

    bitset<15> statusFlags(GenPart_statusFlags[i]);
    if (statusFlags.to_string()[1] == '0')
    {
      continue;
    } // take last copy of W

    bool exclude = false;
    for (unsigned int j = 0; j < t_bkg_idx.size(); j += 3)
    {
      if (i == GenPart_genPartIdxMother[t_bkg_idx[j + 1]])
      {
        exclude = true;
        break;
      } // exclude W's from t
    }

    if (exclude)
    {
      continue;
    }
    W_idx.push_back(i);
  }

  int nW = W_idx.size();
  RVec<int> W_daughter_idx(nW * 2, -99);

  for (unsigned int i = 0; i < nW; i++)
  {
    int n = 0;
    for (unsigned int j = 0; j < nGenPart; j++)
    {
      if (GenPart_genPartIdxMother[j] != W_idx[i])
      {
        continue;
      } // pick out daughters of W
      if (abs(GenPart_pdgId[j]) > 17)
      {
        continue;
      } // to exclude 24->22,24

      W_daughter_idx[i * 2 + n] = j; // record the first copy of W daughter
      n += 1;
    }
  }

  return W_daughter_idx;
};

// The following functions could probably all go to the plotting marco
// Commented Method Only
auto FatJet_matching_sig(string sample, RVec<float> &goodcleanFatJets, RVec<float> &gcFatJet_eta, RVec<float> &gcFatJet_phi, int NFatJets, RVec<int> &FatJet_subJetIdx1, unsigned int nSubJet, RVec<int> &SubJet_hadronFlavour, RVec<int> &GenPart_pdgId, double daughterb_gen_eta, double daughterb_gen_phi, double tDaughter1_gen_eta, double tDaughter1_gen_phi, int tDaughter1_gen_pdgId, double tDaughter2_gen_eta, double tDaughter2_gen_phi, int tDaughter2_gen_pdgId, double WDaughter1_gen_eta, double WDaughter1_gen_phi, int WDaughter1_gen_pdgId, double WDaughter2_gen_eta, double WDaughter2_gen_phi, int WDaughter2_gen_pdgId)
{
  RVec<int> matched_GenPart(NFatJets, -9);
  if (sample.find("Bprime") == std::string::npos)
  {
    return matched_GenPart;
  }

  RVec<int> gcFatJet_subJetIdx1 = FatJet_subJetIdx1[goodcleanFatJets];

  // cout << "Event: " << endl;
  for (unsigned int i = 0; i < NFatJets; i++)
  {
    // cout << "\n" << "Fatjet: " << endl;
    double fatjet_eta = gcFatJet_eta[i];
    double fatjet_phi = gcFatJet_phi[i];

    double dR_b = DeltaR(fatjet_eta, daughterb_gen_eta, fatjet_phi, daughterb_gen_phi);
    double dR_q1 = DeltaR(fatjet_eta, tDaughter1_gen_eta, fatjet_phi, tDaughter1_gen_phi);
    double dR_q2 = DeltaR(fatjet_eta, tDaughter2_gen_eta, fatjet_phi, tDaughter2_gen_phi);

    double dR_q3 = DeltaR(fatjet_eta, WDaughter1_gen_eta, fatjet_phi, WDaughter1_gen_phi);
    double dR_q4 = DeltaR(fatjet_eta, WDaughter2_gen_eta, fatjet_phi, WDaughter2_gen_phi);

    if (dR_b < 0.8 && dR_q1 < 0.8 && dR_q2 < 0.8)
    {
      if (abs(tDaughter1_gen_pdgId) < 6)
      {
        matched_GenPart[i] = 6;
      } // pos stands for hadronic t
      else
      {
        matched_GenPart[i] = -6;
      } // neg stands for leptonic t
    }
    else if (dR_q1 < 0.8 && dR_q2 < 0.8)
    {
      if (abs(tDaughter1_gen_pdgId) < 6)
      {
        matched_GenPart[i] = 24;
      }
      else
      {
        matched_GenPart[i] = -24;
      }
    }

    if (dR_q3 < 0.8 && dR_q4 < 0.8)
    {
      if (abs(WDaughter1_gen_pdgId) < 6)
      {
        matched_GenPart[i] = 24;
      }
      else
      {
        matched_GenPart[i] = -24;
      }
    }

    if (matched_GenPart[i] != -9)
    {
      continue;
    }

    int firstsub = FatJet_subJetIdx1[i];
    for (int isub = firstsub; isub < nSubJet; isub++)
    {
      if (SubJet_hadronFlavour[isub] != 0)
      {
        matched_GenPart[i] = SubJet_hadronFlavour[isub];
      }
      else
      {
        matched_GenPart[i] = 0;
      }
    }
  }
  return matched_GenPart;
};

// Commented Method Only
auto FatJet_matching_bkg(string sample, RVec<float> &goodcleanFatJets, RVec<float> &gcFatJet_eta, RVec<float> &gcFatJet_phi, int NFatJets, RVec<int> &FatJet_subJetIdx1, unsigned int nSubJet, RVec<int> &SubJet_hadronFlavour, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_phi, RVec<float> &GenPart_eta, RVec<int> &GenPart_genPartIdxMother, RVec<int> &t_bkg_idx, RVec<int> &W_bkg_idx)
{
  RVec<int> matched_GenPart(NFatJets, -9);
  if (sample.find("Bprime") != std::string::npos)
  {
    return matched_GenPart;
  }

  RVec<int> gcFatJet_subJetIdx1 = FatJet_subJetIdx1[goodcleanFatJets];

  int ntD = t_bkg_idx.size();
  int nWD = W_bkg_idx.size();

  RVec<double> t_eta(ntD, -9);
  RVec<double> t_phi(ntD, -9);
  RVec<double> W_eta(nWD, -9);
  RVec<double> W_phi(nWD, -9);

  if (ntD != 0)
  {
    for (unsigned int i = 0; i < ntD; i++)
    {
      int igen = t_bkg_idx[i];
      int id = GenPart_pdgId[igen];
      for (unsigned int j = igen; j < nGenPart; j++)
      {
        if (GenPart_pdgId[j] != id)
        {
          continue;
        }
        if (GenPart_genPartIdxMother[j] != igen)
        {
          continue;
        }
        igen = j; // take the last copy of t daughter
      }
      t_eta[i] = GenPart_eta[igen];
      t_phi[i] = GenPart_phi[igen];
    }
  }

  if (nWD != 0)
  {
    for (unsigned int i = 0; i < nWD; i++)
    {
      int igen = W_bkg_idx[i];
      int id = GenPart_pdgId[igen];
      for (unsigned int j = igen; j < nGenPart; j++)
      {
        if (GenPart_pdgId[j] != id)
        {
          continue;
        }
        if (GenPart_genPartIdxMother[j] != igen)
        {
          continue;
        }
        igen = j; // take the last copy of W daughter
      }
      W_eta[i] = GenPart_eta[igen];
      W_phi[i] = GenPart_phi[igen];
    }
  }

  for (unsigned int i = 0; i < NFatJets; i++)
  {
    double fatjet_eta = gcFatJet_eta[i];
    double fatjet_phi = gcFatJet_phi[i];

    for (unsigned int j = 0; j < t_bkg_idx.size() / 3; j++)
    {
      double dR_b = DeltaR(fatjet_eta, t_eta[j * 3], fatjet_phi, t_phi[j * 3]);
      double dR_q1 = DeltaR(fatjet_eta, t_eta[j * 3 + 1], fatjet_phi, t_phi[j * 3 + 1]);
      double dR_q2 = DeltaR(fatjet_eta, t_eta[j * 3 + 2], fatjet_phi, t_phi[j * 3 + 2]);

      if (dR_b < 0.8 && dR_q1 < 0.8 && dR_q2 < 0.8)
      {
        if (abs(GenPart_pdgId[t_bkg_idx[j * 3 + 1]]) < 6)
        {
          matched_GenPart[i] = 6;
          break;
        } // pos stands for hadronic t
        else
        {
          matched_GenPart[i] = -6;
          break;
        } // neg stands for leptonic t
      }
      else if (dR_q1 < 0.8 && dR_q2 < 0.8)
      {
        if (abs(GenPart_pdgId[t_bkg_idx[j * 3 + 1]]) < 6)
        {
          matched_GenPart[i] = 24;
          break;
        }
        else
        {
          matched_GenPart[i] = -24;
          break;
        }
      }
    }

    if (matched_GenPart[i] != -9)
    {
      continue;
    }
    for (unsigned int j = 0; j < W_bkg_idx.size() / 2; j++)
    {
      double dR_q1 = DeltaR(fatjet_eta, W_eta[j * 2], fatjet_phi, W_phi[j * 2]);
      double dR_q2 = DeltaR(fatjet_eta, W_eta[j * 2 + 1], fatjet_phi, W_phi[j * 2 + 1]);

      if (dR_q1 < 0.8 && dR_q2 < 0.8)
      {
        if (abs(GenPart_pdgId[j * 2]) < 6)
        {
          matched_GenPart[i] = 24;
          break;
        }
        else
        {
          matched_GenPart[i] = -24;
          break;
        }
      }
    }

    if (matched_GenPart[i] != -9)
    {
      continue;
    }
    int firstsub = FatJet_subJetIdx1[i];
    for (int isub = firstsub; isub < nSubJet; isub++)
    {
      if (SubJet_hadronFlavour[isub] != 0)
      {
        matched_GenPart[i] = SubJet_hadronFlavour[isub];
      }
      else
      {
        matched_GenPart[i] = 0;
      }
    }
  }
  return matched_GenPart;
};
