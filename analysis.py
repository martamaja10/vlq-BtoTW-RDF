import ROOT
from ROOT import RDataFrame, TFile, gInterpreter, TLorentzVector
import time
cpp_files = [
    "BPrime.cc",
    "cleanJet.cc",
    "cut_ptrel.cc",
    "dnnPrep.cc",
    "generatorInfo.cc",
    "utilities.cc",
    "W_t_reco.cc",
    "analyzer_RDF.h",
    "lumiMask.h"

]
for cpp_file in cpp_files:
    ROOT.gInterpreter.Declare(f'#include "{cpp_file}"')

poly2 = ROOT.TF1("poly2", "max(0.402806, 0.998174 + (8.40861e-05)*x + (-6.63274e-07)*x*x + (4.09272e-10)*x*x*x + (-9.50233e-14)*x*x*x*x + (7.59648e-18)*x*x*x*x*x)", 100, 5000)

poly2U = ROOT.TF1("poly2U", "max([6], [0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)", 100, 5000)
poly2U.SetParameters(0.998174, 8.40861e-05, -6.63274e-07, 4.09272e-10, -9.50233e-14, 7.59648e-18, 0.402806)

poly2D = ROOT.TF1("poly2D", "max([6], [0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)", 100, 5000)
poly2D.SetParameters(0.998174, 8.40861e-05, -6.63274e-07, 4.09272e-10, -9.50233e-14, 7.59648e-18, 0.402806)

# Define the polynomials for TOP HT scaling
polyHT = ROOT.TF1("polyHT", "min(1.0, max([0] + [1]*x, [2]))", 700, 5000)
polyHT.SetParameters(1.0, -1.0e-04, 0.1)

polyHTU = ROOT.TF1("polyHTU", "min(1.0, max([0] + [1]*x + sqrt([3] + 2*x*[4] + x*x*[5]), [2] + [6]))", 700, 5000)
polyHTU.SetParameters(1.0, -1.0e-04, 0.1, 1.0, -1.0e-06, 1.0e-09, 0.1)

polyHTD = ROOT.TF1("polyHTD", "min(1.0, max([0] + [1]*x - sqrt([3] + 2*x*[4] + x*x*[5]), [2] - [6]))", 700, 5000)
polyHTD.SetParameters(1.0, -1.0e-04, 0.1, 1.0, -1.0e-06, 1.0e-09, 0.1)

class RDFAnalyzer:
    def __init__(self, inputFiles, testNum1, testNum2, year):
        self.sample = inputFiles
        self.samplebin = testNum1
        self.year = year
        self.isMC = not any(x in inputFiles for x in ["Single", "EGamma"])
        self.files = inputFiles  

        self.poly2 = poly2
        self.poly2U = poly2U
        self.poly2D = poly2D
        self.polyHT = polyHT
        self.polyHTU = polyHTU
        self.polyHTD = polyHTD

    def topHTpoly(self, AK4HT):
        return ROOT.std.vector('double')([
            self.polyHT.Eval(AK4HT),
            self.polyHTU.Eval(AK4HT),
            self.polyHTD.Eval(AK4HT)
        ])
      
       

    def analyzer_RDF(self, testNum, jesvar):
        ROOT.ROOT.EnableImplicitMT(4)
        print(f"Number of threads: {ROOT.ROOT.GetThreadPoolSize()}")

        start_time = time.time()

        sample = self.sample
        samplebin = self.samplebin
        year = self.year
        isMC = self.isMC

        print(f"Sample in cc: {sample}, bin # {samplebin}")
        print(f"Year in cc: {year}")
        print(f"isMC? {isMC}, jesvar = {jesvar}")
        if not isMC:
            print(f"Data era = {self.era}, for jec {self.jecera}")
        #else:
         #   print(f"MC extension tag (blank or ext) = {self.era}")

        # Here the initialization and correction variables would be defined
        # Similar to the C++ code, we would define functions and other variables
        # [placeholder]
        # Skipping the initialization of TF1 functions and vectors as they need to be adapted to Python.

        rdf_input = RDataFrame("Events", self.files)

        # Apply MET filters
        METgeneralFilters = rdf_input.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1", "MET Filters") \
                                       .Filter("nJet > 0 && nFatJet > 0", "Event has > 1 AK4 and > 1 AK8")

        truth = METgeneralFilters
        #print(truth)
        print("Flag, all good here, here")
        # Define Lepton definitions
        LepDefs = truth.Define("Electron_cutBasedIdNoIso_tight", "Electron_cutBasedIdNoIso_tight(nElectron, Electron_vidNestedWPBitmap)") \
                        .Define("TPassMu", "abs(Muon_eta)<2.4 && Muon_mediumId==1 && Muon_miniIsoId>=3 && abs(Muon_dz) < 0.5 && Muon_dxy < 0.2") \
                        .Define("TPassEl", "abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)") \
                        .Define("VetoMu", "TPassMu && (Muon_pt>25)") \
                        .Define("VetoEl", "TPassEl && (Electron_pt>25)") \
                        .Define("SignalIsoMu", "TPassMu && (Muon_pt>=55)") \
                        .Define("SignalIsoEl", "TPassEl && (Electron_pt>=55)") \
                        #.Define("nVetoLep", "(int) (Sum(VetoMu)+Sum(VetoEl))")
        

        # Define Lepton selection
        tkmutrig = " || HLT_OldMu100 || HLT_TkMu100" if year != "2016" else " || HLT_TkMu50"
        eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165"
        LepSelect = LepDefs.Define("isMu", f"(nMuon>0) && (HLT_Mu50{tkmutrig}) && (nElectron == 0)") \
                            .Define("isEl", f"(nElectron>0) && ({eltrig}) && (nMuon == 0)") \
                            .Filter("isMu || isEl", "Event is either muon or electron")
        #tkmutrig = " || HLT_OldMu100 || HLT_TkMu100" if year != "2016" else " || HLT_TkMu50"
        #eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon200"
        #LepSelect = LepDefs.Define("isMu", f"(nMuon>0) && (HLT_Mu50{tkmutrig}) && (nSignalIsoMu==1) && (nVetoIsoLep==0) && (nElectron == 0 || nSignalIsoEl == 0)") \
                            #.Define("isEl", f"(nElectron>0) && ({eltrig}) && (nSignalIsoEl==1) && (nVetoIsoLep==0) && (nMuon == 0 || nSignalIsoMu == 0)") \
                            #.Filter("isMu || isEl", "Event is either muon or electron")

        # Lep Assign
        LepAssign = LepSelect.Define("assignleps", "assign_leps(isMu,isEl,SignalIsoMu,SignalIsoEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_miniPFRelIso_all,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_miniPFRelIso_all)") \
        .Define("lepton_pt", "assignleps[0]") \
        .Define("lepton_eta", "assignleps[1]") \
        .Define("lepton_phi", "assignleps[2]") \
        .Define("lepton_mass", "assignleps[3]") \
        .Define("lepton_miniIso", "assignleps[4]")
        # ---------------------------------------------------------
        #                    MET Selection
        # ---------------------------------------------------------
        METSelect = LepAssign.Filter("MET_pt > 60", "Pass corr MET > 60");

        # ---------------------------------------------------------
        #               HT Calculation and N Jets cuts
        # ---------------------------------------------------------
        JetSelect = (METSelect.Define("DR_lepJets", "DeltaR_VecAndFloat(Jet_eta, Jet_phi, lepton_eta, lepton_phi)")
                        .Define("ptrel_lepJets", "ptRel(Jet_pt, Jet_eta, Jet_phi, Jet_mass, lepton_pt, lepton_eta, lepton_phi, lepton_mass)")
                        .Define("goodJets", "Jet_pt > 30 && abs(Jet_eta) < 2.5 && Jet_jetId > 1 && (DR_lepJets > 0.4 || ptrel_lepJets > 20)")
                        .Define("gcJet_HT", "Sum(Jet_pt[goodJets == true])")
                        .Define("DR_lepFatJets", "DeltaR_VecAndFloat(FatJet_eta, FatJet_phi, lepton_eta, lepton_phi)")
                        .Define("goodFatJets", "FatJet_pt > 200 && abs(FatJet_eta) < 2.5 && FatJet_jetId > 1 && (DR_lepFatJets > 0.8)")
                        .Define("NFatJets", "(int) Sum(goodFatJets)")
                        .Define("NOS_gcFatJets", "(int) Sum(DR_lepFatJets[goodFatJets == true] > TMath::Pi()/2)")
                        .Filter("gcJet_HT > 250", "Pass HT > 250")
                        .Filter("NFatJets > 0", "Pass N good central AK8 > 0")
                        .Filter("NOS_gcFatJets > 0", "Pass N good central other side AK8 > 0"))
        # ---------------------------------------------------------
        #         Jet pt ordering, counting, lepton association
        # ---------------------------------------------------------
        JetVars = (JetSelect.Define("NJets_central", "(int) Sum(goodJets)")
                    .Define("gcJet_pt_unsort", "Jet_pt[goodJets == true]")
                    .Define("gcJet_ptargsort", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcJet_pt_unsort))")
                    .Define("gcJet_pt", "reorder(gcJet_pt_unsort, gcJet_ptargsort)")
                    .Define("gcJet_eta", "reorder(Jet_eta[goodJets == true], gcJet_ptargsort)")
                    .Define("gcJet_phi", "reorder(Jet_phi[goodJets == true], gcJet_ptargsort)")
                    .Define("gcJet_mass", "reorder(Jet_mass[goodJets == true], gcJet_ptargsort)")
                    .Define("gcJet_DeepFlav", "reorder(Jet_btagDeepFlavB[goodJets == true], gcJet_ptargsort)")
                    .Define("DR_gcJets_central", "reorder(DR_lepJets[goodJets == true], gcJet_ptargsort)")
                    .Define("minDR_lepJets", "ROOT::VecOps::Min(DR_gcJets_central)")
                    .Define("ptrel_atMinDR_lepJets", "reorder(ptrel_lepJets[goodJets == true], gcJet_ptargsort)[ROOT::VecOps::ArgMin(DR_gcJets_central)]")
                    .Define("OS_gcJets", "DR_gcJets_central > TMath::Pi()/2")
                    .Define("SS_gcJets", "DR_gcJets_central <= TMath::Pi()/2")
                    .Define("NOS_gcJets_central", "(int) Sum(OS_gcJets)")
                    .Define("NSS_gcJets_central", "(int) Sum(SS_gcJets)")
                    .Define("gcOSJet_pt", "gcJet_pt[OS_gcJets == true]")
                    .Define("gcOSJet_eta", "gcJet_eta[OS_gcJets == true]")
                    .Define("gcOSJet_phi", "gcJet_phi[OS_gcJets == true]")
                    .Define("gcOSJet_mass", "gcJet_mass[OS_gcJets == true]")
                    .Define("gcSSJet_pt", "gcJet_pt[SS_gcJets == true]")
                    .Define("gcSSJet_eta", "gcJet_eta[SS_gcJets == true]")
                    .Define("gcSSJet_phi", "gcJet_phi[SS_gcJets == true]")
                    .Define("gcSSJet_mass", "gcJet_mass[SS_gcJets == true]"))

        ForwardJetVars = JetVars.Define("goodcleanForwardJets", "Jet_pt > 30 and abs(Jet_eta) >= 2.5 and Jet_jetId > 1") \
        .Define("NJets_forward", "(int) Sum(goodcleanForwardJets)") \
        .Define("gcforwJet_pt_unsort", "Jet_pt[goodcleanForwardJets == true]") \
        .Define("gcforwJet_ptargsort", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcforwJet_pt_unsort))") \
        .Define("gcforwJet_pt", "reorder(gcforwJet_pt_unsort, gcforwJet_ptargsort)") \
        .Define("gcforwJet_eta", "reorder(Jet_eta[goodcleanForwardJets == true], gcforwJet_ptargsort)") \
        .Define("gcforwJet_phi", "reorder(Jet_phi[goodcleanForwardJets == true], gcforwJet_ptargsort)") \
        .Define("gcforwJet_mass", "reorder(Jet_mass[goodcleanForwardJets == true], gcforwJet_ptargsort)") \
        .Define("gcforwJet_DeepFlav", "reorder(Jet_btagDeepFlavB[goodcleanForwardJets == true], gcforwJet_ptargsort)")

        FatJetVars = ForwardJetVars.Define("gcFatJet_pt_unsort", "FatJet_pt[goodFatJets == true]") \
        .Define("gcFatJet_ptargsort", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcFatJet_pt_unsort))") \
        .Define("gcFatJet_pt", "reorder(gcFatJet_pt_unsort, gcFatJet_ptargsort)") \
        .Define("gcFatJet_eta", "reorder(FatJet_eta[goodFatJets == true], gcFatJet_ptargsort)") \
        .Define("gcFatJet_phi", "reorder(FatJet_phi[goodFatJets == true], gcFatJet_ptargsort)") \
        .Define("gcFatJet_mass", "reorder(FatJet_mass[goodFatJets == true], gcFatJet_ptargsort)") \
        .Define("gcFatJet_sdmass", "reorder(FatJet_msoftdrop[goodFatJets == true], gcFatJet_ptargsort)") \
        .Define("DR_gcFatJets", "reorder(DR_lepFatJets[goodFatJets == true], gcFatJet_ptargsort)") \
        .Define("minDR_lepFatJets", "ROOT::VecOps::Min(DR_gcFatJets)") \
        .Define("ptrel_atMinDR_lepFatJets", "ptRel(gcFatJet_pt, gcFatJet_eta, gcFatJet_phi, gcFatJet_mass, lepton_pt, lepton_eta, lepton_phi, lepton_mass)[ROOT::VecOps::ArgMin(DR_gcFatJets)]") \
        .Define("SS_gcFatJets", "DR_gcFatJets <= TMath::Pi()/2") \
        .Define("NSS_gcFatJets", "(int) Sum(SS_gcFatJets)") \
        .Define("OS_gcFatJets", "DR_gcFatJets > TMath::Pi()/2") \
        .Define("gcOSFatJet_pt", "gcFatJet_pt[OS_gcFatJets == true]") \
        .Define("gcOSFatJet_eta", "gcFatJet_eta[OS_gcFatJets == true]") \
        .Define("gcOSFatJet_phi", "gcFatJet_phi[OS_gcFatJets == true]") \
        .Define("gcOSFatJet_mass", "gcFatJet_mass[OS_gcFatJets == true]") \
        .Define("gcOSFatJet_sdmass", "gcFatJet_sdmass[OS_gcFatJets == true]")

        # ---------------------------------------------------------
        # 		Add scale factors and MC jet-based calcs
        # ---------------------------------------------------------
        scaleFactors = FatJetVars
        # ---------------------------------------------------------
        # 		  JET Tagging variables
        # ---------------------------------------------------------
        Taggers = (scaleFactors.Define("lepton_lv", "lvConstructor(lepton_pt, lepton_eta, lepton_phi, lepton_mass)")
        .Define("gcJet_ST", "gcJet_HT + lepton_pt + MET_pt")
        .Define("gcFatJet_pNetJ", "reorder(FatJet_particleNet_QCD[goodFatJets == true], gcFatJet_ptargsort)")
        .Define("gcFatJet_pNetTvsQCD", "reorder(FatJet_particleNet_TvsQCD[goodFatJets == true], gcFatJet_ptargsort)")
        .Define("gcFatJet_pNetWvsQCD", "reorder(FatJet_particleNet_WvsQCD[goodFatJets == true], gcFatJet_ptargsort)")
        .Define("gcOSFatJet_pNetJ", "gcFatJet_pNetJ[OS_gcFatJets == true]")
        .Define("gcOSFatJet_pNetTvsQCD", "gcFatJet_pNetTvsQCD[OS_gcFatJets == true]")
        .Define("gcOSFatJet_pNetWvsQCD", "gcFatJet_pNetWvsQCD[OS_gcFatJets == true]")
        .Define("gcFatJet_pNetT", "(gcFatJet_pNetTvsQCD * gcFatJet_pNetJ) / (1 - gcFatJet_pNetTvsQCD)")
        .Define("gcFatJet_pNetW", "(gcFatJet_pNetWvsQCD * gcFatJet_pNetJ) / (1 - gcFatJet_pNetWvsQCD)")
        .Define("gcOSFatJet_pNetT", "gcFatJet_pNetT[OS_gcFatJets == true]")
        .Define("gcOSFatJet_pNetW", "gcFatJet_pNetW[OS_gcFatJets == true]"))

        # Continue with other processing steps
        # [place holder] 
        # Reconstruction

        Reconstruction = Taggers.Define("W_lv", "W_reco(MET_pt,MET_phi,lepton_lv)") \
        .Define("W_pt", "W_lv.Pt()") \
        .Define("W_eta", "W_lv.Eta()") \
        .Define("W_phi", "W_lv.Phi()") \
        .Define("W_mass", "W_lv.M()") \
        .Define("W_MT", "sqrt(2*lepton_pt*MET_pt*(1-cos(lepton_phi - MET_phi)))") \
        .Define("minMlj_output", "minM_lep_jet_calc(gcJet_pt, gcJet_eta, gcJet_phi, gcJet_mass, lepton_lv)") \
        .Define("DR_W_lep", "W_lv.DeltaR(lepton_lv)") \
        .Define("minM_lep_Jet", "minMlj_output[0]") \
        .Define("minM_lep_Jet_jetID", "(int) minMlj_output[1]") \
        .Define("minM_lep_Jet_TorW", "isLeptonic_X(minM_lep_Jet)") \
        # Save Snapshot to file
        print(f"-------------------------------------------------")
        print(f">>> Saving {sample} Snapshot...")
        finalFile = "output_snapshot.root"
         #finalFile = f"RDF_{sample}_{year}_{testNum}.root"
         #finalFile = f"RDF_{sample}{self.era}_{year}_{testNum}.root"
        snapCol = [] 

        opts = ROOT.ROOT.RDF.RSnapshotOptions()
        if jesvar != "Nominal":
            opts.fMode = "UPDATE"

        Reconstruction.Snapshot(f"Events_{jesvar}", finalFile, snapCol, opts)
        print(f"Output File: {finalFile}")
        print(f"-------------------------------------------------")

        end_time = time.time()
        print(f"Execution Time: {end_time - start_time} seconds")

        print("Cut statistics:")
        Reconstruction.Report().Print()

        if jesvar == "Nominal":
            print("Adding Counter tree to the file:")
            rdf_runs = RDataFrame("Runs", self.files)
            rdf_runs.Snapshot("Runs", finalFile, rdf_runs.GetColumnNames(), opts)

        print("Done!")

def runRDF(testNum1, testNum2, inputFile, year):
    # Check if the input is a .txt file
    if inputFile.endswith('.txt'):
        # Read the .txt file and extract ROOT file paths
        with open(inputFile, 'r') as f:
            root_files = [line.strip() for line in f if line.strip()]
    else:
        # If not a .txt file, assume it's a single ROOT file path or a TChain
        root_files = [inputFile]

    t = RDFAnalyzer(root_files, testNum1, testNum2, year)
    isData = any(x in inputFile for x in ["Single", "EGamma"])

    if isData:
        t.analyzer_RDF(testNum1, "Nominal")
    else:
        shifts = ["Nominal", "JECup", "JECdn", "JERup", "JERdn"]
        for shift in shifts:
            print(f"\nRunning shift {shift}")
            t.analyzer_RDF(testNum1, shift)
            print(f"\nFinished shift {shift}")

    print("\nFinished all analyzing")


runRDF("1", "3", "/Users/andreaolamejicanos/vlq-BtoTW-RDF/inputfile.txt", "2016")
