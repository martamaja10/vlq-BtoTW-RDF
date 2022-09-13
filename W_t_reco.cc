// -------------------------------------------------
//         W-BOSON TLORENTZVECTOR CALCULATOR
// -------------------------------------------------
//This is needed for single search
TLorentzVector lpNu_WCalc(float corr_met, float corr_met_phi, TLorentzVector lepton_lv)
{
	const double MW = 80.4;
	double metpx = corr_met * cos(corr_met_phi);
	double metpy = corr_met * sin(corr_met_phi);
	double metpt = corr_met;
	TLorentzVector Wlv_1, Wlv_2, Wlv, Nulv_1, Nulv_2;
	double nuPz_1;
	double nuPz_2;
	
	double Dtmp = MW*MW - lepton_lv.M()*lepton_lv.M() + 2*(lepton_lv.Px()*metpx + lepton_lv.Py()*metpy);
	double Atmp = 4.0*(lepton_lv.Energy()*lepton_lv.Energy() - lepton_lv.Pz()*lepton_lv.Pz());
	double Btmp = -4.0*Dtmp*lepton_lv.Pz();
	double Ctmp = 4.0*lepton_lv.Energy()*lepton_lv.Energy()*metpt*metpt - Dtmp*Dtmp;
	
	double DETtmp = Btmp*Btmp - 4.0*Atmp*Ctmp;
	
	if(DETtmp >= 0)
	{
		nuPz_1 = (-Btmp + TMath::Sqrt(DETtmp))/(2.0*Atmp);
		nuPz_2 = (-Btmp - TMath::Sqrt(DETtmp))/(2.0*Atmp);
		Nulv_1.SetPxPyPzE(metpx,metpy,nuPz_1,TMath::Sqrt((metpt)*(metpt)+(nuPz_1)*(nuPz_1)));
		Nulv_2.SetPxPyPzE(metpx,metpy,nuPz_2,TMath::Sqrt((metpt)*(metpt)+(nuPz_2)*(nuPz_2)));
	}
	else
	{
		nuPz_1 = -Btmp/(2.0*Atmp);
		nuPz_2 = -Btmp/(2.0*Atmp);
		// does another quad solution for pT and scales pT in result. Reduces M, pT, DR.
		double alpha = lepton_lv.Px()*metpx/metpt + lepton_lv.Py()*metpy/metpt;
		double Delta = MW*MW - lepton_lv.M()*lepton_lv.M();
		Atmp = 4.0*(lepton_lv.Pz()*lepton_lv.Pz() - lepton_lv.Energy()*lepton_lv.Energy() + alpha*alpha);
		Btmp = 4.0*alpha*Delta;
		Ctmp = Delta*Delta;
		float DETtmp2 = Btmp*Btmp - 4.0*Atmp*Ctmp;
		double pTnu_1 = (-Btmp + TMath::Sqrt(DETtmp2))/(2.0*Atmp);
		double pTnu_2 = (-Btmp - TMath::Sqrt(DETtmp2))/(2.0*Atmp);
		Nulv_1.SetPxPyPzE(metpx*pTnu_1/metpt,metpy*pTnu_1/metpt,nuPz_1,TMath::Sqrt(pTnu_1*pTnu_1 + nuPz_1*nuPz_1));
		Nulv_2.SetPxPyPzE(metpx*pTnu_2/metpt,metpy*pTnu_2/metpt,nuPz_2,TMath::Sqrt(pTnu_2*pTnu_2 + nuPz_2*nuPz_2));
	}
	Wlv_1 = Nulv_1 + lepton_lv;
	Wlv_2 = Nulv_2 + lepton_lv;
	
	if(fabs(Wlv_1.M() - MW) < fabs(Wlv_2.M() - MW)) {Wlv = Wlv_1;}
	else {Wlv = Wlv_2;}
	return Wlv;
};

// ----------------------------------------------------
// 		  BABY dR_Wt CALCULATOR
// ----------------------------------------------------

double dR_Wt_Calc(TLorentzVector Wlv, TLorentzVector lepton_lv){return Wlv.DeltaR(lepton_lv);};

// -----------------------------------------------------
// 		    isLEPTONIC W or t
// -----------------------------------------------------
bool isLeptonic_X(float minMleppJet)
{
	int isLep_X = -1; // 0 for W --- 1 for t
	bool isLeptonic_W = false;
	bool isLeptonic_t = false;
	// best combo of W vs t truth match with this
	if(minMleppJet > 150){isLeptonic_W = true; isLep_X = 0;}
	else{isLeptonic_t = true; isLep_X = 1;}
	return isLep_X;
};

// -----------------------------------------------------
// 	  HOMEMADE TLORENTZVECTOR CONSTRUCTOR
// -----------------------------------------------------

ROOT::VecOps::RVec<float> lpNu_t_Calc(int isLeptonic, ROOT::VecOps::RVec<float>& jet_pt, ROOT::VecOps::RVec<float>& jet_eta, ROOT::VecOps::RVec<float>& jet_phi, ROOT::VecOps::RVec<float>& jet_mass, TLorentzVector Wlv, float minMleppJet, int ind_MinMlj)
{
	float t_mass = -999;
	float t_pt = -999;
	float t_eta = -999;
	float t_phi = -999;
	float t_dRWb = -999;
	float deltaRbW = 999;
	int bIndex = 999;
	TLorentzVector jet_lv, top_lv;
	for(unsigned int ijet=0; ijet < jet_pt.size(); ijet++)
	{
		jet_lv.SetPtEtaPhiM(jet_pt.at(ijet),jet_eta.at(ijet),jet_phi.at(ijet),jet_mass.at(ijet));
		if(jet_lv.DeltaR(Wlv) < deltaRbW)
		{
			deltaRbW = jet_lv.DeltaR(Wlv);
			bIndex = ijet;
		}
	}
	
	// Form a leptonic top candidate if the b is close enough
	if(isLeptonic == 1)
	{
		if(deltaRbW > 0.8) {bIndex = ind_MinMlj;} // use a close b unless it doesn't exist
		TLorentzVector bottom_lv;
		bottom_lv.SetPtEtaPhiM(jet_pt.at(bIndex),jet_eta.at(bIndex),jet_phi.at(bIndex),jet_mass.at(bIndex));
		top_lv = bottom_lv + Wlv;
		t_mass = top_lv.M();
		t_pt = top_lv.Pt();
		t_eta = top_lv.Eta();
		t_phi = top_lv.Phi();
		t_dRWb = bottom_lv.DeltaR(Wlv);
	}
	else
	{
		t_pt = 9999;
		t_eta = 9;
		t_phi = 9;
		t_mass = -999;
	}
	ROOT::VecOps::RVec<float> t_FiveVec = {t_pt,t_eta,t_phi,t_mass,t_dRWb};
	return t_FiveVec;
};

// -------------------------------------------------
//    TLORENTZVECTOR CONSTRUCTOR FOR FLOAT ONLYS
// -------------------------------------------------
TLorentzVector top_lvConstructor(float pt, float eta, float phi, float mass)
{
	TLorentzVector lv;
	lv.SetPtEtaPhiM(pt,eta,phi,mass);
	return lv;
};

