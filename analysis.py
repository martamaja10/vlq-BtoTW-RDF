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
    "lumiMask.h",

]

ROOT.gInterpreter.Load("functions_cpp.so")

for cpp_file in cpp_files:
    ROOT.gInterpreter.Declare(f'#include "{cpp_file}"')




class RDFAnalyzer:
    def __init__(self, inputFiles, testNum1, testNum2, year):
        self.sample = inputFiles
        self.samplebin = testNum1
        self.year = year
        self.isMC = not any(x in inputFiles for x in ["Single", "EGamma"])
        self.files = inputFiles  



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
        def get_year_settings(year):
            deepjetL = None
            mutrig = "TkMu50"
            yrstr = ""
            yr = ""
            jecyr = ""
            jeryr = ""
            jecver = ""

            if year == "2016APV":
                deepjetL = 0.0508
                yrstr = "2016preVFP"
                yr = "16"
                jecyr = "UL16APV"
                jeryr = "Summer20UL16APV_JRV3"
                jecver = "V7"
            elif year == "2016":
                deepjetL = 0.0480
                yrstr = "2016postVFP"
                yr = "16"
                jecyr = "UL16"
                jeryr = "Summer20UL16_JRV3"
                jecver = "V7"
            elif year == "2017":
                mutrig = "OldMu100_or_TkMu100"
                deepjetL = 0.0532
                yrstr = "2017"
                yr = "17"
                jecyr = "UL17"
                jeryr = "Summer19UL17_JRV2"
                jecver = "V5"
            elif year == "2018":
                mutrig = "OldMu100_or_TkMu100"
                deepjetL = 0.0490
                yrstr = "2018"
                yr = "18"
                jecyr = "UL18"
                jeryr = "Summer19UL18_JRV2"
                jecver = "V5"
            else:
                raise ValueError(f"ERROR: Can't parse the year to assign correctionLib json files. Expected 2016, 2016APV, 2017, or 2018. Got: {year}")

            return deepjetL, mutrig, yrstr, yr, jecyr, jeryr, jecver

        deepjetL, mutrig, yrstr, yr, jecyr, jeryr, jecver = get_year_settings("2016")
       
        def pnetWPs(year,dnnT, dnnW):
            # Set the working points based on the year
            year="2016"
            if year == "2016APV":
                wpT, wpW = 0.490, 0.677
            elif year == "2016":
                wpT, wpW = 0.495, 0.668
            elif year == "2017":
                wpT, wpW = 0.581, 0.709
            else:
                wpT, wpW = 0.580, 0.700

            # Initialize a vector for tags
            tag = RVec(int)()  # Use RVec for consistency with ROOT
            
            # Iterate over the input RVecs and apply the tag logic
            for i in range(len(dnnT)):
                tagval = 0  # Initialize tag value
                if dnnT[i] > wpT:
                    tagval += 1  # Add 1 if top score exceeds threshold
                if dnnW[i] > wpW:
                    tagval += 2  # Add 2 if W score exceeds threshold
                tag.push_back(tagval)  # Append the tag value to the tag vector

            return tag
           

        

        rdf_input = RDataFrame("Events", self.files)

        # Apply MET filters
        METgeneralFilters = rdf_input.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1", "MET Filters") \
                                       .Filter("nJet > 0 && nFatJet > 0", "Event has > 1 AK4 and > 1 AK8")

        truth = METgeneralFilters
        #print(truth)
        print("Flag, all good here, here")
        # Define Lepton definitions
        

        elHEMcut = ""
        if year == "2018":
            elHEMcut = " && (Electron_eta > -1.479 || (Electron_phi < -1.57 || Electron_phi > -0.87))"
        LepDefs = truth.Define("Electron_cutBasedIdNoIso_tight", "Electron_cutBasedIdNoIso_tight(nElectron, Electron_vidNestedWPBitmap)") \
                        .Define("TPassMu", "abs(Muon_eta)<2.4 && Muon_mediumId==1 && Muon_miniIsoId>=3 && abs(Muon_dz) < 0.5 && Muon_dxy < 0.2") \
                        .Define("TPassEl", "(abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_cutBasedIdNoIso_tight==1 && Electron_miniPFRelIso_all<0.1" + elHEMcut) \
                        .Define("VetoMu", "TPassMu && (Muon_pt>25)") \
                        .Define("VetoEl", "TPassEl && (Electron_pt>25)") \
                        .Define("SignalIsoMu", "TPassMu && (Muon_pt>=55)") \
                        .Define("SignalIsoEl", "TPassEl && (Electron_pt>=55)") \
                        .Define("nVetoLep", "(int) (Sum(VetoMu)+Sum(VetoEl))") \
                        .Define("SMuon_pt", "Muon_pt[SignalIsoMu == true]") \
                        .Define("SMuon_eta", "Muon_eta[SignalIsoMu == true]") \
                        .Define("SMuon_phi", "Muon_phi[SignalIsoMu == true]") \
                        .Define("SMuon_mass", "Muon_mass[SignalIsoMu == true]") \
                        .Define("SElectron_pt", "Electron_pt[SignalIsoEl == true]") \
                        .Define("SElectron_eta", "Electron_eta[SignalIsoEl == true]") \
                        .Define("SElectron_phi", "Electron_phi[SignalIsoEl == true]") \
                        .Define("SElectron_mass", "Electron_mass[SignalIsoEl == true]") \
                        .Define("Muon_P4", "fVectorConstructor(Muon_pt,Muon_eta,Muon_phi,Muon_mass)") \
                        .Define("SMuon_P4", "fVectorConstructor(SMuon_pt,SMuon_eta,SMuon_phi,SMuon_mass)") \
                        .Define("SElectron_P4", "fVectorConstructor(SElectron_pt,SElectron_eta,SElectron_phi,SElectron_mass)") \
                        .Define("SMuon_jetIdx", "Muon_jetIdx[SignalIsoMu == true]") \
                        .Define("SElectron_jetIdx", "Electron_jetIdx[SignalIsoEl == true]") \
                        .Define("nSignalIsoMu", "(int) Sum(SignalIsoMu)") \
                        .Define("nSignalIsoEl", "(int) Sum(SignalIsoEl)") \
                        .Define("VetoIsoMu", "(VetoMu == true && Muon_pt < 55)") \
                        .Define("VetoIsoEl", "(VetoEl == true && Electron_pt < 55)") \
                        .Define("nVetoIsoLep", "(int) (Sum(VetoIsoMu)+Sum(VetoIsoEl))") \
        
        tkmutrig = " || HLT_OldMu100 || HLT_TkMu100"
        eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon200"
        if year == "2017" and era == "B":
            tkmutrig = ""
            eltrig = "HLT_Ele35_WPTight_Gsf || HLT_Photon200"
        elif year == "2016" or year == "2016APV":
            tkmutrig = " || HLT_TkMu50"
            eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon175"
        if year == "2016APV" and (era == "A" or era == "B"):
            tkmutrig = ""
            eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon175"
        
        LepSelect = LepDefs.Define("isMu", f"(nMuon > 0) && (HLT_Mu50{tkmutrig}) && (nSignalIsoMu == 1) && (nVetoIsoLep == 0) && (nElectron == 0 or nSignalIsoEl == 0)") \
        .Define("isEl", f"(nElectron > 0) && ({eltrig}) && (nSignalIsoEl == 1) && (nVetoIsoLep == 0) && (nMuon == 0 or nSignalIsoMu == 0)") \
        .Filter("isMu || isEl", "Event is either muon or electron") \
        


        # Define the logic for different year and era conditions
        # Define Lepton selection
        #tkmutrig = " || HLT_OldMu100 || HLT_TkMu100" if year != "2016" else " || HLT_TkMu50"
        #eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165"
        #LepSelect = LepDefs.Define("isMu", f"(nMuon>0) && (HLT_Mu50{tkmutrig}) && (nElectron == 0)") \
        #                    .Define("isEl", f"(nElectron>0) && ({eltrig}) && (nMuon == 0)") \
         #                   .Filter("isMu || isEl", "Event is either muon or electron")
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
        METSelect = LepAssign.Filter("MET_pt > 60", "Pass corr MET > 60")

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
        #Define("gcHTCorr_top",ROOT.topHTpoly, ["gcJet_HT"])
        JetVars = (JetSelect.Define("NJets_central", "(int) Sum(goodJets)")
                    .Define("gcJet_pt_unsort", "Jet_pt[goodJets == true]")
                    .Define("gcJet_ptargsort", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcJet_pt_unsort))")
                    .Define("gcJet_pt", "reorder(gcJet_pt_unsort, gcJet_ptargsort)")
                    .Define("gcJet_eta", "reorder(Jet_eta[goodJets == true], gcJet_ptargsort)")
                    .Define("gcJet_phi", "reorder(Jet_phi[goodJets == true], gcJet_ptargsort)")
                    .Define("gcJet_mass", "reorder(Jet_mass[goodJets == true], gcJet_ptargsort)")
                    .Define("gcJet_DeepFlav", "reorder(Jet_btagDeepFlavB[goodJets == true], gcJet_ptargsort)")
                    .Define("gcJet_DeepFlavL", f"gcJet_DeepFlav > {deepjetL}") 
                    .Define("NJets_DeepFlavL", "(int) Sum(gcJet_DeepFlavL)")
                    .Define("DR_gcJets_central", "reorder(DR_lepJets[goodJets == true],gcJet_ptargsort)")
                    .Define("minDR_lepJets", "ROOT::VecOps::Min(DR_gcJets_central)")
                    .Define("ptrel_atMinDR_lepJets", "reorder(ptrel_lepJets[goodJets == true],gcJet_ptargsort)[ROOT::VecOps::ArgMin(DR_gcJets_central)]")
                    .Define("OS_gcJets", "DR_gcJets_central > TMath::Pi()/2")
                    .Define("SS_gcJets", "DR_gcJets_central <= TMath::Pi()/2")
                    .Define("NOS_gcJets_central", "(int) Sum(OS_gcJets)")
                    .Define("NSS_gcJets_central", "(int) Sum(SS_gcJets)")
                    .Define("gcOSJet_pt", "gcJet_pt[OS_gcJets == true]")
                    .Define("gcOSJet_eta", "gcJet_eta[OS_gcJets == true]")
                    .Define("gcOSJet_phi", "gcJet_phi[OS_gcJets == true]")
                    .Define("gcOSJet_mass", "gcJet_mass[OS_gcJets == true]")
                    .Define("gcOSJet_DeepFlavL", "gcJet_DeepFlavL[OS_gcJets == true]")
                    .Define("gcSSJet_pt", "gcJet_pt[SS_gcJets == true]")
                    .Define("gcSSJet_eta", "gcJet_eta[SS_gcJets == true]")
                    .Define("gcSSJet_phi", "gcJet_phi[SS_gcJets == true]")
                    .Define("gcSSJet_mass", "gcJet_mass[SS_gcJets == true]")
                    .Define("gcSSJet_DeepFlavL", "gcJet_DeepFlavL[SS_gcJets == true]")
                    .Define("NOS_gcJets_DeepFlavL", "(int) Sum(gcOSJet_DeepFlavL)")
                    .Define("NSS_gcJets_DeepFlavL", "(int) Sum(gcSSJet_DeepFlavL)")
                    )


                    #.Define("DR_gcJets_central", "reorder(DR_lepJets[goodJets == true], gcJet_ptargsort)")
                    #.Define("minDR_lepJets", "ROOT::VecOps::Min(DR_gcJets_central)")
                    #.Define("ptrel_atMinDR_lepJets", "reorder(ptrel_lepJets[goodJets == true], gcJet_ptargsort)[ROOT::VecOps::ArgMin(DR_gcJets_central)]")
                    #.Define("OS_gcJets", "DR_gcJets_central > TMath::Pi()/2")
                    #.Define("SS_gcJets", "DR_gcJets_central <= TMath::Pi()/2")
                    #.Define("NOS_gcJets_central", "(int) Sum(OS_gcJets)")
                    #.Define("NSS_gcJets_central", "(int) Sum(SS_gcJets)")
                    #.Define("gcOSJet_pt", "gcJet_pt[OS_gcJets == true]")
                    #.Define("gcOSJet_eta", "gcJet_eta[OS_gcJets == true]")
                    #.Define("gcOSJet_phi", "gcJet_phi[OS_gcJets == true]")
                    #.Define("gcSSJet_pt", "gcJet_pt[SS_gcJets == true]")
                    ##.Define("gcOSJet_mass", "gcJet_mass[OS_gcJets == true]")
                    #.Define("gcSSJet_eta", "gcJet_eta[SS_gcJets == true]")
                    #.Define("gcSSJet_phi", "gcJet_phi[SS_gcJets == true]")
                    #.Define("gcSSJet_mass", "gcJet_mass[SS_gcJets == true]")
                    #.Define("gcOSJet_DeepFlavL","gcJet_DeepFlavL[OS_gcJets == true]")
                    #.Define("gcSSJet_pt","gcJet_pt[SS_gcJets == true]")
                    #.Define("gcSSJet_eta","gcJet_eta[SS_gcJets == true]")
                    #.Define("gcSSJet_phi","gcJet_phi[SS_gcJets == true]")
                    #.Define("gcSSJet_mass","gcJet_mass[SS_gcJets == true]")
                    #.Define("gcSSJet_DeepFlavL","gcJet_DeepFlavL[SS_gcJets == true]")
                    #.Define("NOS_gcJets_DeepFlavL","(int) Sum(gcOSJet_DeepFlavL)")
                    #.Define("NSS_gcJets_DeepFlavL","(int) Sum(gcSSJet_DeepFlavL)"))

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
         #.Define("gcFatJet_pNetTag", 'Numba::pnetWPs(gcFatJet_pNetTvsQCD, gcFatJet_pNetWvsQCD)')
         #.Define("gcOSFatJet_pNetTag", "gcFatJet_pNetTag[OS_gcFatJets==true]")
         #.Define("gcFatJet_nJ", "Sum(gcFatJet_pNetTag == 0)")
         #.Define("gcFatJet_nT", "Sum(gcFatJet_pNetTag == 1)")
         #.Define("gcFatJet_nW", "Sum(gcFatJet_pNetTag == 2)")
         #.Define("gcFatJet_tau21", "reorder((FatJet_tau2 / FatJet_tau1)[goodFatJets == true],gcFatJet_ptargsort)")
         #.Define("gcOSFatJet_tau21", "gcFatJet_tau21[OS_gcFatJets == true]")
         #.Define("gcFatJet_tau32", "reorder((FatJet_tau3 / FatJet_tau2)[goodFatJets == true],gcFatJet_ptargsort)")
         #.Define("gcOSFatJet_tau32", "gcFatJet_tau32[OS_gcFatJets == true]")
         #.Define("minDR_leadAK8otherAK8", "minDR_leadJetOtherJet_calc(gcFatJet_eta,gcFatJet_phi)")
         #.Define("minDR_leadAK4otherAK4", "minDR_leadJetOtherJet_calc(gcJet_eta,gcJet_phi)")
         #.Define("minDR_AK8s_discrete","std::floor(minDR_leadAK8otherAK8/0.5)")
         #.Define("minDR_AK4s_discrete","std::floor(minDR_leadAK4otherAK4/0.5)"))


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
