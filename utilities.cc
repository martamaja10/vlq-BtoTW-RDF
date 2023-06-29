// Methods in this file:
// leptonicCheck(), Electron_cutBasedIdNoIso_tight(), genttbarMassCalc(C), fVectorConstructor(), lvConstructor(C)
// I put some of the methods in here becasue I didn't know where else to put them.

using namespace std;
using namespace ROOT::VecOps;

// The following functions could probably all go to the plotting marco
auto leptonicCheck(string sample, int trueLeptonicT, int trueLeptonicW)
{
    if (sample != "Bprime")
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

// ----------------------------------------------------
//   		ttbar background mass CALCULATOR:
// ----------------------------------------------------

// Commented Method Only
auto genttbarMassCalc(string sample, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_mass, RVec<float> &GenPart_pt, RVec<float> &GenPart_phi, RVec<float> &GenPart_eta, RVec<int> &GenPart_genPartIdxMother, RVec<int> &GenPart_status)
{
    int returnVar = 0;
    if (sample == "ttbar")
    {
        int genTTbarMass = -999;
        double topPtWeight = 1.0;
        TLorentzVector top, antitop;
        bool gottop = false;
        bool gotantitop = false;
        bool gottoppt = false;
        bool gotantitoppt = false;
        float toppt, antitoppt;
        for (unsigned int p = 0; p < nGenPart; p++)
        {
            int id = GenPart_pdgId[p];
            if (abs(id) != 6)
            {
                continue;
            }
            if (GenPart_mass[p] < 10)
            {
                continue;
            }
            int motherid = GenPart_pdgId[GenPart_genPartIdxMother[p]];
            if (abs(motherid) != 6)
            {
                if (!gottop && id == 6)
                {
                    top.SetPtEtaPhiM(GenPart_pt[p], GenPart_eta[p], GenPart_phi[p], GenPart_mass[p]);
                    gottop = true;
                }
                if (!gotantitop && id == -6)
                {
                    antitop.SetPtEtaPhiM(GenPart_pt[p], GenPart_eta[p], GenPart_phi[p], GenPart_mass[p]);
                    gotantitop = true;
                }
            }
            if (GenPart_status[p] == 62)
            {
                if (!gottoppt && id == 6)
                {
                    toppt = GenPart_pt[p];
                    gottoppt = true;
                }
                if (!gotantitoppt && id == -6)
                {
                    antitoppt = GenPart_pt[p];
                    gotantitoppt = true;
                }
            }
        }
        if (gottop && gotantitop)
        {
            genTTbarMass = (top + antitop).M();
        }
        if (gottoppt && gotantitoppt)
        {
            float SFtop = TMath::Exp(0.0615 - 0.0005 * toppt);
            float SFantitop = TMath::Exp(0.0615 - 0.0005 * antitoppt);
            topPtWeight = TMath::Sqrt(SFtop * SFantitop);
        }
        returnVar = genTTbarMass;
    }
    return returnVar;
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