import ROOT
import os 

# Question: A lot of the samples have _PSweights_.  The samples we were working with before didn't have this.  Is it good or bad?
# Sample Dictionaries: samples, samples_2016APVUL, samples_2016UL, samples_2017UL, samples_2018UL, samples_test, samples_QCD

class sample:
    def __init__(self, prefix, year, textlist, samplename): #, color
        self.prefix = prefix
        self.year = year
        self.textlist = textlist
        self.samplename = samplename
<<<<<<< HEAD
        #self.color = color
                        
Bprime_M1000_20UL16 = sample("Bprime_M1000_20UL16", "20UL16", "Bprime_M1000_20UL16NanoList.txt", "/BprimeBtoTW_M-1000_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
Bprime_M1000_20UL17 = sample("Bprime_M1000_20UL17", "20UL17", "Bprime_M1000_20UL17NanoList.txt", "/BprimeBtoTW_M-1000_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
Bprime_M1000_20UL18 = sample("Bprime_M1000_20UL18", "20UL18", "Bprime_M1000_20UL18NanoList.txt", "/BprimeBtoTW_M-1000_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
Bprime_M1200_20UL16 = sample("Bprime_M1200_20UL16", "20UL16", "Bprime_M1200_20UL16NanoList.txt", "/BprimeBtoTW_M-1200_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
Bprime_M1200_20UL17 = sample("Bprime_M1200_20UL17", "20UL17", "Bprime_M1200_20UL17NanoList.txt", "/BprimeBtoTW_M-1200_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
Bprime_M1200_20UL18 = sample("Bprime_M1200_20UL18", "20UL18", "Bprime_M1200_20UL18NanoList.txt", "/BprimeBtoTW_M-1200_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
Bprime_M1300_20UL16 = sample("Bprime_M1300_20UL16", "20UL16", "Bprime_M1300_20UL16NanoList.txt", "/BprimeBtoTW_M-1300_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
Bprime_M1300_20UL17 = sample("Bprime_M1300_20UL17", "20UL17", "Bprime_M1300_20UL17NanoList.txt", "/BprimeBtoTW_M-1300_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
Bprime_M1300_20UL18 = sample("Bprime_M1300_20UL18", "20UL18", "Bprime_M1300_20UL18NanoList.txt", "/BprimeBtoTW_M-1300_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
Bprime_M1400_20UL16 = sample("Bprime_M1400_20UL16", "20UL16", "Bprime_M1400_20UL16NanoList.txt", "/BprimeBtoTW_M-1400_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
Bprime_M1400_20UL17 = sample("Bprime_M1400_20UL17", "20UL17", "Bprime_M1400_20UL17NanoList.txt", "/BprimeBtoTW_M-1400_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
Bprime_M1400_20UL18 = sample("Bprime_M1400_20UL18", "20UL18", "Bprime_M1400_20UL18NanoList.txt", "/BprimeBtoTW_M-1400_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
Bprime_M1500_20UL16 = sample("Bprime_M1500_20UL16", "20UL16", "Bprime_M1500_20UL16NanoList.txt", "/BprimeBtoTW_M-1500_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
Bprime_M1500_20UL17 = sample("Bprime_M1500_20UL17", "20UL17", "Bprime_M1500_20UL17NanoList.txt", "/BprimeBtoTW_M-1500_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
Bprime_M1500_20UL18 = sample("Bprime_M1500_20UL18", "20UL18", "Bprime_M1500_20UL18NanoList.txt", "/BprimeBtoTW_M-1500_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
Bprime_M1600_20UL16 = sample("Bprime_M1600_20UL16", "20UL16", "Bprime_M1600_20UL16NanoList.txt", "/BprimeBtoTW_M-1600_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
Bprime_M1600_20UL17 = sample("Bprime_M1600_20UL17", "20UL17", "Bprime_M1600_20UL17NanoList.txt", "/BprimeBtoTW_M-1600_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
Bprime_M1600_20UL18 = sample("Bprime_M1600_20UL18", "20UL18", "Bprime_M1600_20UL18NanoList.txt", "/BprimeBtoTW_M-1600_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
Bprime_M1700_20UL16 = sample("Bprime_M1700_20UL16", "20UL16", "Bprime_M1700_20UL16NanoList.txt", "/BprimeBtoTW_M-1700_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
Bprime_M1700_20UL17 = sample("Bprime_M1700_20UL17", "20UL17", "Bprime_M1700_20UL17NanoList.txt", "/BprimeBtoTW_M-1700_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
Bprime_M1700_20UL18 = sample("Bprime_M1700_20UL18", "20UL18", "Bprime_M1700_20UL18NanoList.txt", "/BprimeBtoTW_M-1700_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
Bprime_M1800_20UL16 = sample("Bprime_M1800_20UL16", "20UL16", "Bprime_M1800_20UL16NanoList.txt", "/BprimeBtoTW_M-1800_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
Bprime_M1800_20UL17 = sample("Bprime_M1800_20UL17", "20UL17", "Bprime_M1800_20UL17NanoList.txt", "/BprimeBtoTW_M-1800_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
Bprime_M1800_20UL18 = sample("Bprime_M1800_20UL18", "20UL18", "Bprime_M1800_20UL18NanoList.txt", "/BprimeBtoTW_M-1800_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
Bprime_M2000_20UL16 = sample("Bprime_M2000_20UL16", "20UL16", "Bprime_M2000_20UL16NanoList.txt", "/BprimeBtoTW_M-2000_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
Bprime_M2000_20UL17 = sample("Bprime_M2000_20UL17", "20UL17", "Bprime_M2000_20UL17NanoList.txt", "/BprimeBtoTW_M-2000_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
Bprime_M2000_20UL18 = sample("Bprime_M2000_20UL18", "20UL18", "Bprime_M2000_20UL18NanoList.txt", "/BprimeBtoTW_M-2000_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
Bprime_M2200_20UL16 = sample("Bprime_M2200_20UL16", "20UL16", "Bprime_M2200_20UL16NanoList.txt", "/BprimeBtoTW_M-2200_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
Bprime_M2200_20UL17 = sample("Bprime_M2200_20UL17", "20UL17", "Bprime_M2200_20UL17NanoList.txt", "/BprimeBtoTW_M-2200_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
Bprime_M2200_20UL18 = sample("Bprime_M2200_20UL18", "20UL18", "Bprime_M2200_20UL18NanoList.txt", "/BprimeBtoTW_M-2200_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
Bprime_M800__20UL16 = sample("Bprime_M800__20UL16", "20UL16", "Bprime_M800__20UL16NanoList.txt", "/BprimeBtoTW_M-800_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
Bprime_M800__20UL17 = sample("Bprime_M800__20UL17", "20UL17", "Bprime_M800__20UL17NanoList.txt", "/BprimeBtoTW_M-800_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
Bprime_M800__20UL18 = sample("Bprime_M800__20UL18", "20UL18", "Bprime_M800__20UL18NanoList.txt", "/BprimeBtoTW_M-800_NWALO_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
DYMHT12002016APVUL = sample("DYMHT12002016APVUL", "2016APVUL", "DYMHT12002016APVULNanoList.txt", "/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
DYMHT12002016UL = sample("DYMHT12002016UL", "2016UL", "DYMHT12002016ULNanoList.txt", "/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
DYMHT12002017UL = sample("DYMHT12002017UL", "2017UL", "DYMHT12002017ULNanoList.txt", "/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
DYMHT12002018UL = sample("DYMHT12002018UL", "2018UL", "DYMHT12002018ULNanoList.txt", "/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
DYMHT2002016APVUL = sample("DYMHT2002016APVUL", "2016APVUL", "DYMHT2002016APVULNanoList.txt", "/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
DYMHT2002016UL = sample("DYMHT2002016UL", "2016UL", "DYMHT2002016ULNanoList.txt", "/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
DYMHT2002017UL = sample("DYMHT2002017UL", "2017UL", "DYMHT2002017ULNanoList.txt", "/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
DYMHT2002018UL = sample("DYMHT2002018UL", "2018UL", "DYMHT2002018ULNanoList.txt", "/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
DYMHT25002016APVUL = sample("DYMHT25002016APVUL", "2016APVUL", "DYMHT25002016APVULNanoList.txt", "/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
DYMHT25002016UL = sample("DYMHT25002016UL", "2016UL", "DYMHT25002016ULNanoList.txt", "/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
DYMHT25002017UL = sample("DYMHT25002017UL", "2017UL", "DYMHT25002017ULNanoList.txt", "/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
DYMHT25002018UL = sample("DYMHT25002018UL", "2018UL", "DYMHT25002018ULNanoList.txt", "/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
DYMHT4002016APVUL = sample("DYMHT4002016APVUL", "2016APVUL", "DYMHT4002016APVULNanoList.txt", "/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
DYMHT4002016UL = sample("DYMHT4002016UL", "2016UL", "DYMHT4002016ULNanoList.txt", "/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
DYMHT4002017UL = sample("DYMHT4002017UL", "2017UL", "DYMHT4002017ULNanoList.txt", "/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
DYMHT4002018UL = sample("DYMHT4002018UL", "2018UL", "DYMHT4002018ULNanoList.txt", "/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
DYMHT6002016APVUL = sample("DYMHT6002016APVUL", "2016APVUL", "DYMHT6002016APVULNanoList.txt", "/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
DYMHT6002016UL = sample("DYMHT6002016UL", "2016UL", "DYMHT6002016ULNanoList.txt", "/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
DYMHT6002017UL = sample("DYMHT6002017UL", "2017UL", "DYMHT6002017ULNanoList.txt", "/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
DYMHT6002018UL = sample("DYMHT6002018UL", "2018UL", "DYMHT6002018ULNanoList.txt", "/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
DYMHT8002016APVUL = sample("DYMHT8002016APVUL", "2016APVUL", "DYMHT8002016APVULNanoList.txt", "/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
DYMHT8002016UL = sample("DYMHT8002016UL", "2016UL", "DYMHT8002016ULNanoList.txt", "/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
DYMHT8002017UL = sample("DYMHT8002017UL", "2017UL", "DYMHT8002017ULNanoList.txt", "/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
DYMHT8002018UL = sample("DYMHT8002018UL", "2018UL", "DYMHT8002018ULNanoList.txt", "/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
JetHTRun2016APVB2016APVUL = sample("JetHTRun2016APVB2016APVUL", "2016APVUL", "JetHTRun2016APVB2016APVULNanoList.txt", "/JetHT/Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
JetHTRun2016APVC2016APVUL = sample("JetHTRun2016APVC2016APVUL", "2016APVUL", "JetHTRun2016APVC2016APVULNanoList.txt", "/JetHT/Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
JetHTRun2016APVD2016APVUL = sample("JetHTRun2016APVD2016APVUL", "2016APVUL", "JetHTRun2016APVD2016APVULNanoList.txt", "/JetHT/Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
JetHTRun2016APVE2016APVUL = sample("JetHTRun2016APVE2016APVUL", "2016APVUL", "JetHTRun2016APVE2016APVULNanoList.txt", "/JetHT/Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
JetHTRun2016APVF2016APVUL = sample("JetHTRun2016APVF2016APVUL", "2016APVUL", "JetHTRun2016APVF2016APVULNanoList.txt", "/JetHT/Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
JetHTRun2016F2016UL = sample("JetHTRun2016F2016UL", "2016UL", "JetHTRun2016F2016ULNanoList.txt", "/JetHT/Run2016F-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD")
JetHTRun2016G2016UL = sample("JetHTRun2016G2016UL", "2016UL", "JetHTRun2016G2016ULNanoList.txt", "/JetHT/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD")
JetHTRun2016H2016UL = sample("JetHTRun2016H2016UL", "2016UL", "JetHTRun2016H2016ULNanoList.txt", "/JetHT/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD")
JetHTRun2017B2017UL = sample("JetHTRun2017B2017UL", "2017UL", "JetHTRun2017B2017ULNanoList.txt", "/JetHT/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
JetHTRun2017C2017UL = sample("JetHTRun2017C2017UL", "2017UL", "JetHTRun2017C2017ULNanoList.txt", "/JetHT/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
JetHTRun2017D2017UL = sample("JetHTRun2017D2017UL", "2017UL", "JetHTRun2017D2017ULNanoList.txt", "/JetHT/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
JetHTRun2017E2017UL = sample("JetHTRun2017E2017UL", "2017UL", "JetHTRun2017E2017ULNanoList.txt", "/JetHT/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
JetHTRun2017F2017UL = sample("JetHTRun2017F2017UL", "2017UL", "JetHTRun2017F2017ULNanoList.txt", "/JetHT/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
JetHTRun2018A2018UL = sample("JetHTRun2018A2018UL", "2018UL", "JetHTRun2018A2018ULNanoList.txt", "/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD")
JetHTRun2018B2018UL = sample("JetHTRun2018B2018UL", "2018UL", "JetHTRun2018B2018ULNanoList.txt", "/JetHT/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD")
JetHTRun2018C2018UL = sample("JetHTRun2018C2018UL", "2018UL", "JetHTRun2018C2018ULNanoList.txt", "/JetHT/Run2018C-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD")
JetHTRun2018D2018UL = sample("JetHTRun2018D2018UL", "2018UL", "JetHTRun2018D2018ULNanoList.txt", "/JetHT/Run2018D-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD")
QCDHT10002016APVUL = sample("QCDHT10002016APVUL", "2016APVUL", "QCDHT10002016APVULNanoList.txt", "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
QCDHT10002016UL = sample("QCDHT10002016UL", "2016UL", "QCDHT10002016ULNanoList.txt", "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
QCDHT10002017UL = sample("QCDHT10002017UL", "2017UL", "QCDHT10002017ULNanoList.txt", "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
QCDHT10002018UL = sample("QCDHT10002018UL", "2018UL", "QCDHT10002018ULNanoList.txt", "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCDHT15002016APVUL = sample("QCDHT15002016APVUL", "2016APVUL", "QCDHT15002016APVULNanoList.txt", "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
QCDHT15002016UL = sample("QCDHT15002016UL", "2016UL", "QCDHT15002016ULNanoList.txt", "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
QCDHT15002017UL = sample("QCDHT15002017UL", "2017UL", "QCDHT15002017ULNanoList.txt", "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
QCDHT15002018UL = sample("QCDHT15002018UL", "2018UL", "QCDHT15002018ULNanoList.txt", "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCDHT20002016APVUL = sample("QCDHT20002016APVUL", "2016APVUL", "QCDHT20002016APVULNanoList.txt", "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
QCDHT20002016UL = sample("QCDHT20002016UL", "2016UL", "QCDHT20002016ULNanoList.txt", "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
QCDHT20002017UL = sample("QCDHT20002017UL", "2017UL", "QCDHT20002017ULNanoList.txt", "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
QCDHT20002018UL = sample("QCDHT20002018UL", "2018UL", "QCDHT20002018ULNanoList.txt", "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCDHT2002016APVUL = sample("QCDHT2002016APVUL", "2016APVUL", "QCDHT2002016APVULNanoList.txt", "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
QCDHT2002016UL = sample("QCDHT2002016UL", "2016UL", "QCDHT2002016ULNanoList.txt", "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
QCDHT2002017UL = sample("QCDHT2002017UL", "2017UL", "QCDHT2002017ULNanoList.txt", "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
QCDHT2002018UL = sample("QCDHT2002018UL", "2018UL", "QCDHT2002018ULNanoList.txt", "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCDHT3002016APVUL = sample("QCDHT3002016APVUL", "2016APVUL", "QCDHT3002016APVULNanoList.txt", "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
QCDHT3002016UL = sample("QCDHT3002016UL", "2016UL", "QCDHT3002016ULNanoList.txt", "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
QCDHT3002017UL = sample("QCDHT3002017UL", "2017UL", "QCDHT3002017ULNanoList.txt", "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
QCDHT3002018UL = sample("QCDHT3002018UL", "2018UL", "QCDHT3002018ULNanoList.txt", "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCDHT5002016APVUL = sample("QCDHT5002016APVUL", "2016APVUL", "QCDHT5002016APVULNanoList.txt", "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
QCDHT5002016UL = sample("QCDHT5002016UL", "2016UL", "QCDHT5002016ULNanoList.txt", "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
QCDHT5002018UL = sample("QCDHT5002018UL", "2018UL", "QCDHT5002018ULNanoList.txt", "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCDHT7002016APVUL = sample("QCDHT7002016APVUL", "2016APVUL", "QCDHT7002016APVULNanoList.txt", "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
QCDHT7002016UL = sample("QCDHT7002016UL", "2016UL", "QCDHT7002016ULNanoList.txt", "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
QCDHT7002017UL = sample("QCDHT7002017UL", "2017UL", "QCDHT7002017ULNanoList.txt", "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
QCDHT7002018UL = sample("QCDHT7002018UL", "2018UL", "QCDHT7002018ULNanoList.txt", "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
SingleElectronRun2016APVB2016APVUL = sample("SingleElectronRun2016APVB2016APVUL", "2016APVUL", "SingleElectronRun2016APVB2016APVULNanoList.txt", "/SingleElectron/Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleElectronRun2016APVC2016APVUL = sample("SingleElectronRun2016APVC2016APVUL", "2016APVUL", "SingleElectronRun2016APVC2016APVULNanoList.txt", "/SingleElectron/Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleElectronRun2016APVD2016APVUL = sample("SingleElectronRun2016APVD2016APVUL", "2016APVUL", "SingleElectronRun2016APVD2016APVULNanoList.txt", "/SingleElectron/Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleElectronRun2016APVE2016APVUL = sample("SingleElectronRun2016APVE2016APVUL", "2016APVUL", "SingleElectronRun2016APVE2016APVULNanoList.txt", "/SingleElectron/Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleElectronRun2016APVF2016APVUL = sample("SingleElectronRun2016APVF2016APVUL", "2016APVUL", "SingleElectronRun2016APVF2016APVULNanoList.txt", "/SingleElectron/Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleElectronRun2016F2016UL = sample("SingleElectronRun2016F2016UL", "2016UL", "SingleElectronRun2016F2016ULNanoList.txt", "/SingleElectron/Run2016F-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleElectronRun2016G2016UL = sample("SingleElectronRun2016G2016UL", "2016UL", "SingleElectronRun2016G2016ULNanoList.txt", "/SingleElectron/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleElectronRun2016H2016UL = sample("SingleElectronRun2016H2016UL", "2016UL", "SingleElectronRun2016H2016ULNanoList.txt", "/SingleElectron/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleElectronRun2017B2017UL = sample("SingleElectronRun2017B2017UL", "2017UL", "SingleElectronRun2017B2017ULNanoList.txt", "/SingleElectron/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleElectronRun2017C2017UL = sample("SingleElectronRun2017C2017UL", "2017UL", "SingleElectronRun2017C2017ULNanoList.txt", "/SingleElectron/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleElectronRun2017D2017UL = sample("SingleElectronRun2017D2017UL", "2017UL", "SingleElectronRun2017D2017ULNanoList.txt", "/SingleElectron/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleElectronRun2017E2017UL = sample("SingleElectronRun2017E2017UL", "2017UL", "SingleElectronRun2017E2017ULNanoList.txt", "/SingleElectron/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleElectronRun2017F2017UL = sample("SingleElectronRun2017F2017UL", "2017UL", "SingleElectronRun2017F2017ULNanoList.txt", "/SingleElectron/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleElectronRun2018A2018UL = sample("SingleElectronRun2018A2018UL", "2018UL", "SingleElectronRun2018A2018ULNanoList.txt", "/EGamma/Run2018A-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleElectronRun2018B2018UL = sample("SingleElectronRun2018B2018UL", "2018UL", "SingleElectronRun2018B2018ULNanoList.txt", "/EGamma/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleElectronRun2018C2018UL = sample("SingleElectronRun2018C2018UL", "2018UL", "SingleElectronRun2018C2018ULNanoList.txt", "/EGamma/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleElectronRun2018D2018UL = sample("SingleElectronRun2018D2018UL", "2018UL", "SingleElectronRun2018D2018ULNanoList.txt", "/EGamma/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleMuonRun2016APVB2016APVUL = sample("SingleMuonRun2016APVB2016APVUL", "2016APVUL", "SingleMuonRun2016APVB2016APVULNanoList.txt", "/SingleMuon/Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleMuonRun2016APVC2016APVUL = sample("SingleMuonRun2016APVC2016APVUL", "2016APVUL", "SingleMuonRun2016APVC2016APVULNanoList.txt", "/SingleMuon/Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleMuonRun2016APVD2016APVUL = sample("SingleMuonRun2016APVD2016APVUL", "2016APVUL", "SingleMuonRun2016APVD2016APVULNanoList.txt", "/SingleMuon/Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleMuonRun2016APVE2016APVUL = sample("SingleMuonRun2016APVE2016APVUL", "2016APVUL", "SingleMuonRun2016APVE2016APVULNanoList.txt", "/SingleMuon/Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleMuonRun2016APVF2016APVUL = sample("SingleMuonRun2016APVF2016APVUL", "2016APVUL", "SingleMuonRun2016APVF2016APVULNanoList.txt", "/SingleMuon/Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleMuonRun2016F2016UL = sample("SingleMuonRun2016F2016UL", "2016UL", "SingleMuonRun2016F2016ULNanoList.txt", "/SingleMuon/Run2016F-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleMuonRun2016G2016UL = sample("SingleMuonRun2016G2016UL", "2016UL", "SingleMuonRun2016G2016ULNanoList.txt", "/SingleMuon/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleMuonRun2016H2016UL = sample("SingleMuonRun2016H2016UL", "2016UL", "SingleMuonRun2016H2016ULNanoList.txt", "/SingleMuon/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleMuonRun2017B2017UL = sample("SingleMuonRun2017B2017UL", "2017UL", "SingleMuonRun2017B2017ULNanoList.txt", "/SingleMuon/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleMuonRun2017C2017UL = sample("SingleMuonRun2017C2017UL", "2017UL", "SingleMuonRun2017C2017ULNanoList.txt", "/SingleMuon/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleMuonRun2017D2017UL = sample("SingleMuonRun2017D2017UL", "2017UL", "SingleMuonRun2017D2017ULNanoList.txt", "/SingleMuon/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleMuonRun2017E2017UL = sample("SingleMuonRun2017E2017UL", "2017UL", "SingleMuonRun2017E2017ULNanoList.txt", "/SingleMuon/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleMuonRun2017F2017UL = sample("SingleMuonRun2017F2017UL", "2017UL", "SingleMuonRun2017F2017ULNanoList.txt", "/SingleMuon/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD")
SingleMuonRun2018A2018UL = sample("SingleMuonRun2018A2018UL", "2018UL", "SingleMuonRun2018A2018ULNanoList.txt", "/SingleMuon/Run2018A-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleMuonRun2018B2018UL = sample("SingleMuonRun2018B2018UL", "2018UL", "SingleMuonRun2018B2018ULNanoList.txt", "/SingleMuon/Run2018B-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleMuonRun2018C2018UL = sample("SingleMuonRun2018C2018UL", "2018UL", "SingleMuonRun2018C2018ULNanoList.txt", "/SingleMuon/Run2018C-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD")
SingleMuonRun2018D2018UL = sample("SingleMuonRun2018D2018UL", "2018UL", "SingleMuonRun2018D2018ULNanoList.txt", "/SingleMuon/Run2018D-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD")
STs2016APVUL = sample("STs2016APVUL", "2016APVUL", "STs2016APVULNanoList.txt", "/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
STs2016UL = sample("STs2016UL", "2016UL", "STs2016ULNanoList.txt", "/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
STs2017UL = sample("STs2017UL", "2017UL", "STs2017ULNanoList.txt", "/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
STs2018UL = sample("STs2018UL", "2018UL", "STs2018ULNanoList.txt", "/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
STt2016APVUL = sample("STt2016APVUL", "2016APVUL", "STt2016APVULNanoList.txt", "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
STt2016UL = sample("STt2016UL", "2016UL", "STt2016ULNanoList.txt", "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
STt2017UL = sample("STt2017UL", "2017UL", "STt2017ULNanoList.txt", "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
STt2018UL = sample("STt2018UL", "2018UL", "STt2018ULNanoList.txt", "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
STtb2016APVUL = sample("STtb2016APVUL", "2016APVUL", "STtb2016APVULNanoList.txt", "/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
STtb2016UL = sample("STtb2016UL", "2016UL", "STtb2016ULNanoList.txt", "/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
STtb2017UL = sample("STtb2017UL", "2017UL", "STtb2017ULNanoList.txt", "/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
STtb2018UL = sample("STtb2018UL", "2018UL", "STtb2018ULNanoList.txt", "/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
STtW2016APVUL = sample("STtW2016APVUL", "2016APVUL", "STtW2016APVULNanoList.txt", "/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
STtW2016UL = sample("STtW2016UL", "2016UL", "STtW2016ULNanoList.txt", "/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
STtW2017UL = sample("STtW2017UL", "2017UL", "STtW2017ULNanoList.txt", "/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
STtW2018UL = sample("STtW2018UL", "2018UL", "STtW2018ULNanoList.txt", "/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
STtWb2016APVUL = sample("STtWb2016APVUL", "2016APVUL", "STtWb2016APVULNanoList.txt", "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
STtWb2016UL = sample("STtWb2016UL", "2016UL", "STtWb2016ULNanoList.txt", "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
STtWb2017UL = sample("STtWb2017UL", "2017UL", "STtWb2017ULNanoList.txt", "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
STtWb2018UL = sample("STtWb2018UL", "2018UL", "STtWb2018ULNanoList.txt", "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
TTHB2016APVUL = sample("TTHB2016APVUL", "2016APVUL", "TTHB2016APVULNanoList.txt", "/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
TTHB2016UL = sample("TTHB2016UL", "2016UL", "TTHB2016ULNanoList.txt", "/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
TTHB2017UL = sample("TTHB2017UL", "2017UL", "TTHB2017ULNanoList.txt", "/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
TTHB2018UL = sample("TTHB2018UL", "2018UL", "TTHB2018ULNanoList.txt", "/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
TTHnonB2016APVUL = sample("TTHnonB2016APVUL", "2016APVUL", "TTHnonB2016APVULNanoList.txt", "/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
TTHnonB2016UL = sample("TTHnonB2016UL", "2016UL", "TTHnonB2016ULNanoList.txt", "/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
TTHnonB2017UL = sample("TTHnonB2017UL", "2017UL", "TTHnonB2017ULNanoList.txt", "/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
TTHnonB2018UL = sample("TTHnonB2018UL", "2018UL", "TTHnonB2018ULNanoList.txt", "/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
TTMT10002016APVUL = sample("TTMT10002016APVUL", "2016APVUL", "TTMT10002016APVULNanoList.txt", "/TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
TTMT10002016UL = sample("TTMT10002016UL", "2016UL", "TTMT10002016ULNanoList.txt", "/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
TTMT10002017UL = sample("TTMT10002017UL", "2017UL", "TTMT10002017ULNanoList.txt", "/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
TTMT10002018UL = sample("TTMT10002018UL", "2018UL", "TTMT10002018ULNanoList.txt", "/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
TTMT7002016APVUL = sample("TTMT7002016APVUL", "2016APVUL", "TTMT7002016APVULNanoList.txt", "/TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
TTMT7002016UL = sample("TTMT7002016UL", "2016UL", "TTMT7002016ULNanoList.txt", "/TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
TTMT7002017UL = sample("TTMT7002017UL", "2017UL", "TTMT7002017ULNanoList.txt", "/TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
TTMT7002018UL = sample("TTMT7002018UL", "2018UL", "TTMT7002018ULNanoList.txt", "/TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
TTTo2L2Nu2016APVUL = sample("TTTo2L2Nu2016APVUL", "2016APVUL", "TTTo2L2Nu2016APVULNanoList.txt", "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
TTTo2L2Nu2016UL = sample("TTTo2L2Nu2016UL", "2016UL", "TTTo2L2Nu2016ULNanoList.txt", "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
TTTo2L2Nu2017UL = sample("TTTo2L2Nu2017UL", "2017UL", "TTTo2L2Nu2017ULNanoList.txt", "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
TTTo2L2Nu2018UL = sample("TTTo2L2Nu2018UL", "2018UL", "TTTo2L2Nu2018ULNanoList.txt", "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
TTToHadronic2016APVUL = sample("TTToHadronic2016APVUL", "2016APVUL", "TTToHadronic2016APVULNanoList.txt", "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
TTToHadronic2016UL = sample("TTToHadronic2016UL", "2016UL", "TTToHadronic2016ULNanoList.txt", "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
TTToHadronic2017UL = sample("TTToHadronic2017UL", "2017UL", "TTToHadronic2017ULNanoList.txt", "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
TTToHadronic2018UL = sample("TTToHadronic2018UL", "2018UL", "TTToHadronic2018ULNanoList.txt", "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
TTToSemiLeptonic2016APVUL = sample("TTToSemiLeptonic2016APVUL", "2016APVUL", "TTToSemiLeptonic2016APVULNanoList.txt", "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
TTToSemiLeptonic2016UL = sample("TTToSemiLeptonic2016UL", "2016UL", "TTToSemiLeptonic2016ULNanoList.txt", "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
TTToSemiLeptonic2017UL = sample("TTToSemiLeptonic2017UL", "2017UL", "TTToSemiLeptonic2017ULNanoList.txt", "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
TTToSemiLeptonic2018UL = sample("TTToSemiLeptonic2018UL", "2018UL", "TTToSemiLeptonic2018ULNanoList.txt", "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
TTWl2016APVUL = sample("TTWl2016APVUL", "2016APVUL", "TTWl2016APVULNanoList.txt", "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
TTWl2016UL = sample("TTWl2016UL", "2016UL", "TTWl2016ULNanoList.txt", "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
TTWl2017UL = sample("TTWl2017UL", "2017UL", "TTWl2017ULNanoList.txt", "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
TTWl2018UL = sample("TTWl2018UL", "2018UL", "TTWl2018ULNanoList.txt", "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
TTWq2016APVUL = sample("TTWq2016APVUL", "2016APVUL", "TTWq2016APVULNanoList.txt", "/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
TTWq2016UL = sample("TTWq2016UL", "2016UL", "TTWq2016ULNanoList.txt", "/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
TTWq2017UL = sample("TTWq2017UL", "2017UL", "TTWq2017ULNanoList.txt", "/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
TTWq2018UL = sample("TTWq2018UL", "2018UL", "TTWq2018ULNanoList.txt", "/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
TTZM102016APVUL = sample("TTZM102016APVUL", "2016APVUL", "TTZM102016APVULNanoList.txt", "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
TTZM102016UL = sample("TTZM102016UL", "2016UL", "TTZM102016ULNanoList.txt", "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
TTZM102017UL = sample("TTZM102017UL", "2017UL", "TTZM102017ULNanoList.txt", "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
TTZM102018UL = sample("TTZM102018UL", "2018UL", "TTZM102018ULNanoList.txt", "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
TTZM1to102016APVUL = sample("TTZM1to102016APVUL", "2016APVUL", "TTZM1to102016APVULNanoList.txt", "/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
TTZM1to102016UL = sample("TTZM1to102016UL", "2016UL", "TTZM1to102016ULNanoList.txt", "/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
TTZM1to102017UL = sample("TTZM1to102017UL", "2017UL", "TTZM1to102017ULNanoList.txt", "/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
TTZM1to102018UL = sample("TTZM1to102018UL", "2018UL", "TTZM1to102018ULNanoList.txt", "/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
WJets2016APVUL = sample("WJets2016APVUL", "2016APVUL", "WJets2016APVULNanoList.txt", "/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
WJets2016UL = sample("WJets2016UL", "2016UL", "WJets2016ULNanoList.txt", "/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
WJets2017UL = sample("WJets2017UL", "2017UL", "WJets2017ULNanoList.txt", "/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
WJets2018UL = sample("WJets2018UL", "2018UL", "WJets2018ULNanoList.txt", "/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
WJetsHT12002016APVUL = sample("WJetsHT12002016APVUL", "2016APVUL", "WJetsHT12002016APVULNanoList.txt", "/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
WJetsHT12002016UL = sample("WJetsHT12002016UL", "2016UL", "WJetsHT12002016ULNanoList.txt", "/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
WJetsHT12002017UL = sample("WJetsHT12002017UL", "2017UL", "WJetsHT12002017ULNanoList.txt", "/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
WJetsHT12002018UL = sample("WJetsHT12002018UL", "2018UL", "WJetsHT12002018ULNanoList.txt", "/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
WJetsHT2002016APVUL = sample("WJetsHT2002016APVUL", "2016APVUL", "WJetsHT2002016APVULNanoList.txt", "/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
WJetsHT2002016UL = sample("WJetsHT2002016UL", "2016UL", "WJetsHT2002016ULNanoList.txt", "/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
WJetsHT2002017UL = sample("WJetsHT2002017UL", "2017UL", "WJetsHT2002017ULNanoList.txt", "/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
WJetsHT2002018UL = sample("WJetsHT2002018UL", "2018UL", "WJetsHT2002018ULNanoList.txt", "/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
WJetsHT25002016APVUL = sample("WJetsHT25002016APVUL", "2016APVUL", "WJetsHT25002016APVULNanoList.txt", "/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM")
WJetsHT25002016UL = sample("WJetsHT25002016UL", "2016UL", "WJetsHT25002016ULNanoList.txt", "/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM")
WJetsHT25002017UL = sample("WJetsHT25002017UL", "2017UL", "WJetsHT25002017ULNanoList.txt", "/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
WJetsHT25002018UL = sample("WJetsHT25002018UL", "2018UL", "WJetsHT25002018ULNanoList.txt", "/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM")
WJetsHT4002016APVUL = sample("WJetsHT4002016APVUL", "2016APVUL", "WJetsHT4002016APVULNanoList.txt", "/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
WJetsHT4002016UL = sample("WJetsHT4002016UL", "2016UL", "WJetsHT4002016ULNanoList.txt", "/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
WJetsHT4002017UL = sample("WJetsHT4002017UL", "2017UL", "WJetsHT4002017ULNanoList.txt", "/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
WJetsHT4002018UL = sample("WJetsHT4002018UL", "2018UL", "WJetsHT4002018ULNanoList.txt", "/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
WJetsHT6002016APVUL = sample("WJetsHT6002016APVUL", "2016APVUL", "WJetsHT6002016APVULNanoList.txt", "/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
WJetsHT6002016UL = sample("WJetsHT6002016UL", "2016UL", "WJetsHT6002016ULNanoList.txt", "/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
WJetsHT6002017UL = sample("WJetsHT6002017UL", "2017UL", "WJetsHT6002017ULNanoList.txt", "/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
WJetsHT6002018UL = sample("WJetsHT6002018UL", "2018UL", "WJetsHT6002018ULNanoList.txt", "/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
WJetsHT8002016APVUL = sample("WJetsHT8002016APVUL", "2016APVUL", "WJetsHT8002016APVULNanoList.txt", "/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
WJetsHT8002016UL = sample("WJetsHT8002016UL", "2016UL", "WJetsHT8002016ULNanoList.txt", "/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
WJetsHT8002017UL = sample("WJetsHT8002017UL", "2017UL", "WJetsHT8002017ULNanoList.txt", "/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
WJetsHT8002018UL = sample("WJetsHT8002018UL", "2018UL", "WJetsHT8002018ULNanoList.txt", "/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
WW2016APVUL = sample("WW2016APVUL", "2016APVUL", "WW2016APVULNanoList.txt", "/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
WW2016UL = sample("WW2016UL", "2016UL", "WW2016ULNanoList.txt", "/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
WW2017UL = sample("WW2017UL", "2017UL", "WW2017ULNanoList.txt", "/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
WW2018UL = sample("WW2018UL", "2018UL", "WW2018ULNanoList.txt", "/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
WZ2016APVUL = sample("WZ2016APVUL", "2016APVUL", "WZ2016APVULNanoList.txt", "/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
WZ2016UL = sample("WZ2016UL", "2016UL", "WZ2016ULNanoList.txt", "/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
WZ2017UL = sample("WZ2017UL", "2017UL", "WZ2017ULNanoList.txt", "/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
WZ2018UL = sample("WZ2018UL", "2018UL", "WZ2018ULNanoList.txt", "/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
ZZ2016APVUL = sample("ZZ2016APVUL", "2016APVUL", "ZZ2016APVULNanoList.txt", "/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM")
ZZ2016UL = sample("ZZ2016UL", "2016UL", "ZZ2016ULNanoList.txt", "/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM")
ZZ2017UL = sample("ZZ2017UL", "2017UL", "ZZ2017ULNanoList.txt", "/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM")
ZZ2018UL = sample("ZZ2018UL", "2018UL", "ZZ2018ULNanoList.txt", "/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCD200 = sample("QCD200", "2018", "QCD200NanoList.txt", "QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCD300 = sample("QCD300", "2018", "QCD300NanoList.txt", "QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCD500 = sample("QCD500", "2018", "QCD500NanoList.txt", "QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCD700 = sample("QCD700", "2018", "QCD700NanoList.txt", "QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCD1000 = sample("QCD1000", "2018", "QCD1000NanoList.txt", "QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCD1500 = sample("QCD1500", "2018", "QCD1500NanoList.txt", "QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")
QCD2000 = sample("QCD2000", "2018", "QCD2000NanoList.txt", "QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM")

samples={
    "Bprime_M1000_20UL16":Bprime_M1000_20UL16,
    "Bprime_M1000_20UL17":Bprime_M1000_20UL17,
    "Bprime_M1000_20UL18":Bprime_M1000_20UL18,
    "Bprime_M1200_20UL16":Bprime_M1200_20UL16,
    "Bprime_M1200_20UL17":Bprime_M1200_20UL17,
    "Bprime_M1200_20UL18":Bprime_M1200_20UL18,
    "Bprime_M1300_20UL16":Bprime_M1300_20UL16,
    "Bprime_M1300_20UL17":Bprime_M1300_20UL17,
    "Bprime_M1300_20UL18":Bprime_M1300_20UL18,
    "Bprime_M1400_20UL16":Bprime_M1400_20UL16,
    "Bprime_M1400_20UL17":Bprime_M1400_20UL17,
    "Bprime_M1400_20UL18":Bprime_M1400_20UL18,
    "Bprime_M1500_20UL16":Bprime_M1500_20UL16,
    "Bprime_M1500_20UL17":Bprime_M1500_20UL17,
    "Bprime_M1500_20UL18":Bprime_M1500_20UL18,
    "Bprime_M1600_20UL16":Bprime_M1600_20UL16,
    "Bprime_M1600_20UL17":Bprime_M1600_20UL17,
    "Bprime_M1600_20UL18":Bprime_M1600_20UL18,
    "Bprime_M1700_20UL16":Bprime_M1700_20UL16,
    "Bprime_M1700_20UL17":Bprime_M1700_20UL17,
    "Bprime_M1700_20UL18":Bprime_M1700_20UL18,
    "Bprime_M1800_20UL16":Bprime_M1800_20UL16,
    "Bprime_M1800_20UL17":Bprime_M1800_20UL17,
    "Bprime_M1800_20UL18":Bprime_M1800_20UL18,
    "Bprime_M2000_20UL16":Bprime_M2000_20UL16,
    "Bprime_M2000_20UL17":Bprime_M2000_20UL17,
    "Bprime_M2000_20UL18":Bprime_M2000_20UL18,
    "Bprime_M2200_20UL16":Bprime_M2200_20UL16,
    "Bprime_M2200_20UL17":Bprime_M2200_20UL17,
    "Bprime_M2200_20UL18":Bprime_M2200_20UL18,
    "Bprime_M800__20UL16":Bprime_M800__20UL16,
    "Bprime_M800__20UL17":Bprime_M800__20UL17,
    "Bprime_M800__20UL18":Bprime_M800__20UL18,
    "DYMHT12002016APVUL":DYMHT12002016APVUL,
    "DYMHT12002016UL":DYMHT12002016UL,
    "DYMHT12002017UL":DYMHT12002017UL,
    "DYMHT12002018UL":DYMHT12002018UL,
    "DYMHT2002016APVUL":DYMHT2002016APVUL,
    "DYMHT2002016UL":DYMHT2002016UL,
    "DYMHT2002017UL":DYMHT2002017UL,
    "DYMHT2002018UL":DYMHT2002018UL,
    "DYMHT25002016APVUL":DYMHT25002016APVUL,
    "DYMHT25002016UL":DYMHT25002016UL,
    "DYMHT25002017UL":DYMHT25002017UL,
    "DYMHT25002018UL":DYMHT25002018UL,
    "DYMHT4002016APVUL":DYMHT4002016APVUL,
    "DYMHT4002016UL":DYMHT4002016UL,
    "DYMHT4002017UL":DYMHT4002017UL,
    "DYMHT4002018UL":DYMHT4002018UL,
    "DYMHT6002016APVUL":DYMHT6002016APVUL,
    "DYMHT6002016UL":DYMHT6002016UL,
    "DYMHT6002017UL":DYMHT6002017UL,
    "DYMHT6002018UL":DYMHT6002018UL,
    "DYMHT8002016APVUL":DYMHT8002016APVUL,
    "DYMHT8002016UL":DYMHT8002016UL,
    "DYMHT8002017UL":DYMHT8002017UL,
    "DYMHT8002018UL":DYMHT8002018UL,
    "JetHTRun2016APVB2016APVUL":JetHTRun2016APVB2016APVUL,
    "JetHTRun2016APVC2016APVUL":JetHTRun2016APVC2016APVUL,
    "JetHTRun2016APVD2016APVUL":JetHTRun2016APVD2016APVUL,
    "JetHTRun2016APVE2016APVUL":JetHTRun2016APVE2016APVUL,
    "JetHTRun2016APVF2016APVUL":JetHTRun2016APVF2016APVUL,
    "JetHTRun2016F2016UL":JetHTRun2016F2016UL,
    "JetHTRun2016G2016UL":JetHTRun2016G2016UL,
    "JetHTRun2016H2016UL":JetHTRun2016H2016UL,
    "JetHTRun2017B2017UL":JetHTRun2017B2017UL,
    "JetHTRun2017C2017UL":JetHTRun2017C2017UL,
    "JetHTRun2017D2017UL":JetHTRun2017D2017UL,
    "JetHTRun2017E2017UL":JetHTRun2017E2017UL,
    "JetHTRun2017F2017UL":JetHTRun2017F2017UL,
    "JetHTRun2018A2018UL":JetHTRun2018A2018UL,
    "JetHTRun2018B2018UL":JetHTRun2018B2018UL,
    "JetHTRun2018C2018UL":JetHTRun2018C2018UL,
    "JetHTRun2018D2018UL":JetHTRun2018D2018UL,
    "QCDHT10002016APVUL":QCDHT10002016APVUL,
    "QCDHT10002016UL":QCDHT10002016UL,
    "QCDHT10002017UL":QCDHT10002017UL,
    "QCDHT10002018UL":QCDHT10002018UL,
    "QCDHT15002016APVUL":QCDHT15002016APVUL,
    "QCDHT15002016UL":QCDHT15002016UL,
    "QCDHT15002017UL":QCDHT15002017UL,
    "QCDHT15002018UL":QCDHT15002018UL,
    "QCDHT20002016APVUL":QCDHT20002016APVUL,
    "QCDHT20002016UL":QCDHT20002016UL,
    "QCDHT20002017UL":QCDHT20002017UL,
    "QCDHT20002018UL":QCDHT20002018UL,
    "QCDHT2002016APVUL":QCDHT2002016APVUL,
    "QCDHT2002016UL":QCDHT2002016UL,
    "QCDHT2002017UL":QCDHT2002017UL,
    "QCDHT2002018UL":QCDHT2002018UL,
    "QCDHT3002016APVUL":QCDHT3002016APVUL,
    "QCDHT3002016UL":QCDHT3002016UL,
    "QCDHT3002017UL":QCDHT3002017UL,
    "QCDHT3002018UL":QCDHT3002018UL,
    "QCDHT5002016APVUL":QCDHT5002016APVUL,
    "QCDHT5002016UL":QCDHT5002016UL,
    "QCDHT5002018UL":QCDHT5002018UL,
    "QCDHT7002016APVUL":QCDHT7002016APVUL,
    "QCDHT7002016UL":QCDHT7002016UL,
    "QCDHT7002017UL":QCDHT7002017UL,
    "QCDHT7002018UL":QCDHT7002018UL,
    "SingleElectronRun2016APVB2016APVUL":SingleElectronRun2016APVB2016APVUL,
    "SingleElectronRun2016APVC2016APVUL":SingleElectronRun2016APVC2016APVUL,
    "SingleElectronRun2016APVD2016APVUL":SingleElectronRun2016APVD2016APVUL,
    "SingleElectronRun2016APVE2016APVUL":SingleElectronRun2016APVE2016APVUL,
    "SingleElectronRun2016APVF2016APVUL":SingleElectronRun2016APVF2016APVUL,
    "SingleElectronRun2016F2016UL":SingleElectronRun2016F2016UL,
    "SingleElectronRun2016G2016UL":SingleElectronRun2016G2016UL,
    "SingleElectronRun2016H2016UL":SingleElectronRun2016H2016UL,
    "SingleElectronRun2017B2017UL":SingleElectronRun2017B2017UL,
    "SingleElectronRun2017C2017UL":SingleElectronRun2017C2017UL,
    "SingleElectronRun2017D2017UL":SingleElectronRun2017D2017UL,
    "SingleElectronRun2017E2017UL":SingleElectronRun2017E2017UL,
    "SingleElectronRun2017F2017UL":SingleElectronRun2017F2017UL,
    "SingleElectronRun2018A2018UL":SingleElectronRun2018A2018UL,
    "SingleElectronRun2018B2018UL":SingleElectronRun2018B2018UL,
    "SingleElectronRun2018C2018UL":SingleElectronRun2018C2018UL,
    "SingleElectronRun2018D2018UL":SingleElectronRun2018D2018UL,
    "SingleMuonRun2016APVB2016APVUL":SingleMuonRun2016APVB2016APVUL,
    "SingleMuonRun2016APVC2016APVUL":SingleMuonRun2016APVC2016APVUL,
    "SingleMuonRun2016APVD2016APVUL":SingleMuonRun2016APVD2016APVUL,
    "SingleMuonRun2016APVE2016APVUL":SingleMuonRun2016APVE2016APVUL,
    "SingleMuonRun2016APVF2016APVUL":SingleMuonRun2016APVF2016APVUL,
    "SingleMuonRun2016F2016UL":SingleMuonRun2016F2016UL,
    "SingleMuonRun2016G2016UL":SingleMuonRun2016G2016UL,
    "SingleMuonRun2016H2016UL":SingleMuonRun2016H2016UL,
    "SingleMuonRun2017B2017UL":SingleMuonRun2017B2017UL,
    "SingleMuonRun2017C2017UL":SingleMuonRun2017C2017UL,
    "SingleMuonRun2017D2017UL":SingleMuonRun2017D2017UL,
    "SingleMuonRun2017E2017UL":SingleMuonRun2017E2017UL,
    "SingleMuonRun2017F2017UL":SingleMuonRun2017F2017UL,
    "SingleMuonRun2018A2018UL":SingleMuonRun2018A2018UL,
    "SingleMuonRun2018B2018UL":SingleMuonRun2018B2018UL,
    "SingleMuonRun2018C2018UL":SingleMuonRun2018C2018UL,
    "SingleMuonRun2018D2018UL":SingleMuonRun2018D2018UL,
    "STs2016APVUL":STs2016APVUL,
    "STs2016UL":STs2016UL,
    "STs2017UL":STs2017UL,
    "STs2018UL":STs2018UL,
    "STt2016APVUL":STt2016APVUL,
    "STt2016UL":STt2016UL,
    "STt2017UL":STt2017UL,
    "STt2018UL":STt2018UL,
    "STtb2016APVUL":STtb2016APVUL,
    "STtb2016UL":STtb2016UL,
    "STtb2017UL":STtb2017UL,
    "STtb2018UL":STtb2018UL,
    "STtW2016APVUL":STtW2016APVUL,
    "STtW2016UL":STtW2016UL,
    "STtW2017UL":STtW2017UL,
    "STtW2018UL":STtW2018UL,
    "STtWb2016APVUL":STtWb2016APVUL,
    "STtWb2016UL":STtWb2016UL,
    "STtWb2017UL":STtWb2017UL,
    "STtWb2018UL":STtWb2018UL,
    "TTHB2016APVUL":TTHB2016APVUL,
    "TTHB2016UL":TTHB2016UL,
    "TTHB2017UL":TTHB2017UL,
    "TTHB2018UL":TTHB2018UL,
    "TTHnonB2016APVUL":TTHnonB2016APVUL,
    "TTHnonB2016UL":TTHnonB2016UL,
    "TTHnonB2017UL":TTHnonB2017UL,
    "TTHnonB2018UL":TTHnonB2018UL,
    "TTMT10002016APVUL":TTMT10002016APVUL,
    "TTMT10002016UL":TTMT10002016UL,
    "TTMT10002017UL":TTMT10002017UL,
    "TTMT10002018UL":TTMT10002018UL,
    "TTMT7002016APVUL":TTMT7002016APVUL,
    "TTMT7002016UL":TTMT7002016UL,
    "TTMT7002017UL":TTMT7002017UL,
    "TTMT7002018UL":TTMT7002018UL,
    "TTTo2L2Nu2016APVUL":TTTo2L2Nu2016APVUL,
    "TTTo2L2Nu2016UL":TTTo2L2Nu2016UL,
    "TTTo2L2Nu2017UL":TTTo2L2Nu2017UL,
    "TTTo2L2Nu2018UL":TTTo2L2Nu2018UL,
    "TTToHadronic2016APVUL":TTToHadronic2016APVUL,
    "TTToHadronic2016UL":TTToHadronic2016UL,
    "TTToHadronic2017UL":TTToHadronic2017UL,
    "TTToHadronic2018UL":TTToHadronic2018UL,
    "TTToSemiLeptonic2016APVUL":TTToSemiLeptonic2016APVUL,
    "TTToSemiLeptonic2016UL":TTToSemiLeptonic2016UL,
    "TTToSemiLeptonic2017UL":TTToSemiLeptonic2017UL,
    "TTToSemiLeptonic2018UL":TTToSemiLeptonic2018UL,
    "TTWl2016APVUL":TTWl2016APVUL,
    "TTWl2016UL":TTWl2016UL,
    "TTWl2017UL":TTWl2017UL,
    "TTWl2018UL":TTWl2018UL,
    "TTWq2016APVUL":TTWq2016APVUL,
    "TTWq2016UL":TTWq2016UL,
    "TTWq2017UL":TTWq2017UL,
    "TTWq2018UL":TTWq2018UL,
    "TTZM102016APVUL":TTZM102016APVUL,
    "TTZM102016UL":TTZM102016UL,
    "TTZM102017UL":TTZM102017UL,
    "TTZM102018UL":TTZM102018UL,
    "TTZM1to102016APVUL":TTZM1to102016APVUL,
    "TTZM1to102016UL":TTZM1to102016UL,
    "TTZM1to102017UL":TTZM1to102017UL,
    "TTZM1to102018UL":TTZM1to102018UL,
    "WJetsHT12002016APVUL":WJetsHT12002016APVUL,
    "WJetsHT12002016UL":WJetsHT12002016UL,
    "WJetsHT12002017UL":WJetsHT12002017UL,
    "WJetsHT12002018UL":WJetsHT12002018UL,
    "WJetsHT2002016APVUL":WJetsHT2002016APVUL,
    "WJetsHT2002016UL":WJetsHT2002016UL,
    "WJetsHT2002017UL":WJetsHT2002017UL,
    "WJetsHT2002018UL":WJetsHT2002018UL,
    "WJetsHT25002016APVUL":WJetsHT25002016APVUL,
    "WJetsHT25002016UL":WJetsHT25002016UL,
    "WJetsHT25002017UL":WJetsHT25002017UL,
    "WJetsHT25002018UL":WJetsHT25002018UL,
    "WJetsHT4002016APVUL":WJetsHT4002016APVUL,
    "WJetsHT4002016UL":WJetsHT4002016UL,
    "WJetsHT4002017UL":WJetsHT4002017UL,
    "WJetsHT4002018UL":WJetsHT4002018UL,
    "WJetsHT6002016APVUL":WJetsHT6002016APVUL,
    "WJetsHT6002016UL":WJetsHT6002016UL,
    "WJetsHT6002017UL":WJetsHT6002017UL,
    "WJetsHT6002018UL":WJetsHT6002018UL,
    "WJetsHT8002016APVUL":WJetsHT8002016APVUL,
    "WJetsHT8002016UL":WJetsHT8002016UL,
    "WJetsHT8002017UL":WJetsHT8002017UL,
    "WJetsHT8002018UL":WJetsHT8002018UL,
    "WW2016APVUL":WW2016APVUL,
    "WW2016UL":WW2016UL,
    "WW2017UL":WW2017UL,
    "WW2018UL":WW2018UL,
    "WZ2016APVUL":WZ2016APVUL,
    "WZ2016UL":WZ2016UL,
    "WZ2017UL":WZ2017UL,
    "WZ2018UL":WZ2018UL,
    "ZZ2016APVUL":ZZ2016APVUL,
    "ZZ2016UL":ZZ2016UL,
    "ZZ2017UL":ZZ2017UL,
    "ZZ2018UL":ZZ2018UL,
    "QCD200":QCD200,
    "QCD300":QCD300,
    "QCD500":QCD500,
    "QCD700":QCD700,
    "QCD1000":QCD1000,
    "QCD1500":QCD1500,
    "QCD2000":QCD2000,
}

samples_2016APVUL = {
    
    "DYMHT12002016APVUL":DYMHT12002016APVUL,
    "DYMHT2002016APVUL":DYMHT2002016APVUL,
    "DYMHT25002016APVUL":DYMHT25002016APVUL,
    "DYMHT4002016APVUL":DYMHT4002016APVUL,
    "DYMHT6002016APVUL":DYMHT6002016APVUL,
    "DYMHT8002016APVUL":DYMHT8002016APVUL,
    "JetHTRun2016APVB2016APVUL":JetHTRun2016APVB2016APVUL,
    "JetHTRun2016APVC2016APVUL":JetHTRun2016APVC2016APVUL,
    "JetHTRun2016APVD2016APVUL":JetHTRun2016APVD2016APVUL,
    "JetHTRun2016APVE2016APVUL":JetHTRun2016APVE2016APVUL,
    "JetHTRun2016APVF2016APVUL":JetHTRun2016APVF2016APVUL,
    "QCDHT10002016APVUL":QCDHT10002016APVUL,
    "QCDHT15002016APVUL":QCDHT15002016APVUL,
    "QCDHT20002016APVUL":QCDHT20002016APVUL,
    "QCDHT2002016APVUL":QCDHT2002016APVUL,
    "QCDHT3002016APVUL":QCDHT3002016APVUL,
    "QCDHT5002016APVUL":QCDHT5002016APVUL,
    "QCDHT7002016APVUL":QCDHT7002016APVUL,
    "SingleElectronRun2016APVB2016APVUL":SingleElectronRun2016APVB2016APVUL,
    "SingleElectronRun2016APVC2016APVUL":SingleElectronRun2016APVC2016APVUL,
    "SingleElectronRun2016APVD2016APVUL":SingleElectronRun2016APVD2016APVUL,
    "SingleElectronRun2016APVE2016APVUL":SingleElectronRun2016APVE2016APVUL,
    "SingleElectronRun2016APVF2016APVUL":SingleElectronRun2016APVF2016APVUL,
    "SingleMuonRun2016APVB2016APVUL":SingleMuonRun2016APVB2016APVUL,
    "SingleMuonRun2016APVC2016APVUL":SingleMuonRun2016APVC2016APVUL,
    "SingleMuonRun2016APVD2016APVUL":SingleMuonRun2016APVD2016APVUL,
    "SingleMuonRun2016APVE2016APVUL":SingleMuonRun2016APVE2016APVUL,
    "SingleMuonRun2016APVF2016APVUL":SingleMuonRun2016APVF2016APVUL,
    "STs2016APVUL":STs2016APVUL,
    "STt2016APVUL":STt2016APVUL,
    "STtb2016APVUL":STtb2016APVUL,
    "STtW2016APVUL":STtW2016APVUL,
    "STtWb2016APVUL":STtWb2016APVUL,
    "TTHB2016APVUL":TTHB2016APVUL,
    "TTHnonB2016APVUL":TTHnonB2016APVUL,
    "TTMT10002016APVUL":TTMT10002016APVUL,
    "TTMT7002016APVUL":TTMT7002016APVUL,
    "TTTo2L2Nu2016APVUL":TTTo2L2Nu2016APVUL,
    "TTToHadronic2016APVUL":TTToHadronic2016APVUL,
    "TTToSemiLeptonic2016APVUL":TTToSemiLeptonic2016APVUL,
    "TTWl2016APVUL":TTWl2016APVUL,
    "TTWq2016APVUL":TTWq2016APVUL,
    "TTZM102016APVUL":TTZM102016APVUL,
    "TTZM1to102016APVUL":TTZM1to102016APVUL,
    "WJetsHT12002016APVUL":WJetsHT12002016APVUL,
    "WJetsHT2002016APVUL":WJetsHT2002016APVUL,
    "WJetsHT25002016APVUL":WJetsHT25002016APVUL,
    "WJetsHT4002016APVUL":WJetsHT4002016APVUL,
    "WJetsHT6002016APVUL":WJetsHT6002016APVUL,
    "WJetsHT8002016APVUL":WJetsHT8002016APVUL,
    "WW2016APVUL":WW2016APVUL,
    "WZ2016APVUL":WZ2016APVUL,
    "ZZ2016APVUL":ZZ2016APVUL,

}

samples_2016UL = {
    "Bprime_M1000_20UL16":Bprime_M1000_20UL16,
    "Bprime_M1200_20UL16":Bprime_M1200_20UL16,
    "Bprime_M1300_20UL16":Bprime_M1300_20UL16,
    "Bprime_M1400_20UL16":Bprime_M1400_20UL16,
    "Bprime_M1500_20UL16":Bprime_M1500_20UL16,
    "Bprime_M1600_20UL16":Bprime_M1600_20UL16,
    "Bprime_M1700_20UL16":Bprime_M1700_20UL16,
    "Bprime_M1800_20UL16":Bprime_M1800_20UL16,
    "Bprime_M2000_20UL16":Bprime_M2000_20UL16,
    "Bprime_M2200_20UL16":Bprime_M2200_20UL16,
    "Bprime_M800__20UL16":Bprime_M800__20UL16,
    "DYMHT12002016UL":DYMHT12002016UL,
    "DYMHT2002016UL":DYMHT2002016UL,
    "DYMHT25002016UL":DYMHT25002016UL,
    "DYMHT4002016UL":DYMHT4002016UL,
    "DYMHT6002016UL":DYMHT6002016UL,
    "DYMHT8002016UL":DYMHT8002016UL,
    "JetHTRun2016F2016UL":JetHTRun2016F2016UL,
    "JetHTRun2016G2016UL":JetHTRun2016G2016UL,
    "JetHTRun2016H2016UL":JetHTRun2016H2016UL,
    "QCDHT10002016UL":QCDHT10002016UL,
    "QCDHT15002016UL":QCDHT15002016UL,
    "QCDHT20002016UL":QCDHT20002016UL,
    "QCDHT2002016UL":QCDHT2002016UL,
    "QCDHT3002016UL":QCDHT3002016UL,
    "QCDHT5002016UL":QCDHT5002016UL,
    "QCDHT7002016UL":QCDHT7002016UL,
    "SingleElectronRun2016F2016UL":SingleElectronRun2016F2016UL,
    "SingleElectronRun2016G2016UL":SingleElectronRun2016G2016UL,
    "SingleElectronRun2016H2016UL":SingleElectronRun2016H2016UL,
    "SingleMuonRun2016F2016UL":SingleMuonRun2016F2016UL,
    "SingleMuonRun2016G2016UL":SingleMuonRun2016G2016UL,
    "SingleMuonRun2016H2016UL":SingleMuonRun2016H2016UL,
    "STs2016UL":STs2016UL,
    "STt2016UL":STt2016UL,
    "STtb2016UL":STtb2016UL,
    "STtW2016UL":STtW2016UL,
    "STtWb2016UL":STtWb2016UL,
    "TTHB2016UL":TTHB2016UL,
    "TTHnonB2016UL":TTHnonB2016UL,
    "TTMT10002016UL":TTMT10002016UL,
    "TTMT7002016UL":TTMT7002016UL,
    "TTTo2L2Nu2016UL":TTTo2L2Nu2016UL,
    "TTToHadronic2016UL":TTToHadronic2016UL,
    "TTToSemiLeptonic2016UL":TTToSemiLeptonic2016UL,
    "TTWl2016UL":TTWl2016UL,
    "TTWq2016UL":TTWq2016UL,
    "TTZM102016UL":TTZM102016UL,
    "TTZM1to102016UL":TTZM1to102016UL,
    "WJetsHT12002016UL":WJetsHT12002016UL,
    "WJetsHT2002016UL":WJetsHT2002016UL,
    "WJetsHT25002016UL":WJetsHT25002016UL,
    "WJetsHT4002016UL":WJetsHT4002016UL,
    "WJetsHT6002016UL":WJetsHT6002016UL,
    "WJetsHT8002016UL":WJetsHT8002016UL,
    "WW2016UL":WW2016UL,
    "WZ2016UL":WZ2016UL,
    "ZZ2016UL":ZZ2016UL,
}

samples_2017UL = {
    "Bprime_M1000_20UL17":Bprime_M1000_20UL17,
    "Bprime_M1200_20UL17":Bprime_M1200_20UL17,
    "Bprime_M1300_20UL17":Bprime_M1300_20UL17,
    "Bprime_M1400_20UL17":Bprime_M1400_20UL17,
    "Bprime_M1500_20UL17":Bprime_M1500_20UL17,
    "Bprime_M1600_20UL17":Bprime_M1600_20UL17,
    "Bprime_M1700_20UL17":Bprime_M1700_20UL17,
    "Bprime_M1800_20UL17":Bprime_M1800_20UL17,
    "Bprime_M2000_20UL17":Bprime_M2000_20UL17,
    "Bprime_M2200_20UL17":Bprime_M2200_20UL17,
    "Bprime_M800__20UL17":Bprime_M800__20UL17,
    "DYMHT12002017UL":DYMHT12002017UL,
    "DYMHT2002017UL":DYMHT2002017UL,
    "DYMHT25002017UL":DYMHT25002017UL,
    "DYMHT4002017UL":DYMHT4002017UL,
    "DYMHT6002017UL":DYMHT6002017UL,
    "DYMHT8002017UL":DYMHT8002017UL,
    "JetHTRun2017B2017UL":JetHTRun2017B2017UL,
    "JetHTRun2017C2017UL":JetHTRun2017C2017UL,
    "JetHTRun2017D2017UL":JetHTRun2017D2017UL,
    "JetHTRun2017E2017UL":JetHTRun2017E2017UL,
    "JetHTRun2017F2017UL":JetHTRun2017F2017UL,
    "QCDHT10002017UL":QCDHT10002017UL,
    "QCDHT15002017UL":QCDHT15002017UL,
    "QCDHT20002017UL":QCDHT20002017UL,
    "QCDHT2002017UL":QCDHT2002017UL,
    "QCDHT3002017UL":QCDHT3002017UL,
    "QCDHT7002017UL":QCDHT7002017UL,
    "SingleElectronRun2017B2017UL":SingleElectronRun2017B2017UL,
    "SingleElectronRun2017C2017UL":SingleElectronRun2017C2017UL,
    "SingleElectronRun2017D2017UL":SingleElectronRun2017D2017UL,
    "SingleElectronRun2017E2017UL":SingleElectronRun2017E2017UL,
    "SingleElectronRun2017F2017UL":SingleElectronRun2017F2017UL,
    "SingleMuonRun2017B2017UL":SingleMuonRun2017B2017UL,
    "SingleMuonRun2017C2017UL":SingleMuonRun2017C2017UL,
    "SingleMuonRun2017D2017UL":SingleMuonRun2017D2017UL,
    "SingleMuonRun2017E2017UL":SingleMuonRun2017E2017UL,
    "SingleMuonRun2017F2017UL":SingleMuonRun2017F2017UL,
    "STs2017UL":STs2017UL,
    "STt2017UL":STt2017UL,
    "STtb2017UL":STtb2017UL,
    "STtW2017UL":STtW2017UL,
    "STtWb2017UL":STtWb2017UL,
    "TTHB2017UL":TTHB2017UL,
    "TTHnonB2017UL":TTHnonB2017UL,
    "TTMT10002017UL":TTMT10002017UL,
    "TTMT7002017UL":TTMT7002017UL,
    "TTTo2L2Nu2017UL":TTTo2L2Nu2017UL,
    "TTToHadronic2017UL":TTToHadronic2017UL,
    "TTToSemiLeptonic2017UL":TTToSemiLeptonic2017UL,
    "TTWl2017UL":TTWl2017UL,
    "TTWq2017UL":TTWq2017UL,
    "TTZM102017UL":TTZM102017UL,
    "TTZM1to102017UL":TTZM1to102017UL,
    "WJetsHT12002017UL":WJetsHT12002017UL,
    "WJetsHT2002017UL":WJetsHT2002017UL,
    "WJetsHT25002017UL":WJetsHT25002017UL,
    "WJetsHT4002017UL":WJetsHT4002017UL,
    "WJetsHT6002017UL":WJetsHT6002017UL,
    "WJetsHT8002017UL":WJetsHT8002017UL,
    "WW2017UL":WW2017UL,
    "WZ2017UL":WZ2017UL,
    "ZZ2017UL":ZZ2017UL,
}

samples_2018UL = {
    "Bprime_M1000_20UL18":Bprime_M1000_20UL18,
    "Bprime_M1200_20UL18":Bprime_M1200_20UL18,
    "Bprime_M1300_20UL18":Bprime_M1300_20UL18,
    "Bprime_M1400_20UL18":Bprime_M1400_20UL18,
    "Bprime_M1500_20UL18":Bprime_M1500_20UL18,
    "Bprime_M1600_20UL18":Bprime_M1600_20UL18,
    "Bprime_M1700_20UL18":Bprime_M1700_20UL18,
    "Bprime_M1800_20UL18":Bprime_M1800_20UL18,
    "Bprime_M2000_20UL18":Bprime_M2000_20UL18,
    "Bprime_M2200_20UL18":Bprime_M2200_20UL18,
    "Bprime_M800__20UL18":Bprime_M800__20UL18,
    "DYMHT12002018UL":DYMHT12002018UL,
    "DYMHT2002018UL":DYMHT2002018UL,
    "DYMHT25002018UL":DYMHT25002018UL,
    "DYMHT4002018UL":DYMHT4002018UL,
    "DYMHT6002018UL":DYMHT6002018UL,
    "DYMHT8002018UL":DYMHT8002018UL,
    # "JetHTRun2018A2018UL":JetHTRun2018A2018UL,
    # "JetHTRun2018B2018UL":JetHTRun2018B2018UL,
    # "JetHTRun2018C2018UL":JetHTRun2018C2018UL,
    # "JetHTRun2018D2018UL":JetHTRun2018D2018UL,
    "QCDHT10002018UL":QCDHT10002018UL,
    "QCDHT15002018UL":QCDHT15002018UL,
    "QCDHT20002018UL":QCDHT20002018UL,
    "QCDHT2002018UL":QCDHT2002018UL,
    "QCDHT3002018UL":QCDHT3002018UL,
    "QCDHT5002018UL":QCDHT5002018UL,
    "QCDHT7002018UL":QCDHT7002018UL,
    "SingleElectronRun2018A2018UL":SingleElectronRun2018A2018UL,
    "SingleElectronRun2018B2018UL":SingleElectronRun2018B2018UL,
    "SingleElectronRun2018C2018UL":SingleElectronRun2018C2018UL,
    "SingleElectronRun2018D2018UL":SingleElectronRun2018D2018UL,
    "SingleMuonRun2018A2018UL":SingleMuonRun2018A2018UL,
    "SingleMuonRun2018B2018UL":SingleMuonRun2018B2018UL,
    "SingleMuonRun2018C2018UL":SingleMuonRun2018C2018UL,
    "SingleMuonRun2018D2018UL":SingleMuonRun2018D2018UL,
    "STs2018UL":STs2018UL,
    "STt2018UL":STt2018UL,
    "STtb2018UL":STtb2018UL,
    "STtW2018UL":STtW2018UL,
    "STtWb2018UL":STtWb2018UL,
    "TTHB2018UL":TTHB2018UL,
    "TTHnonB2018UL":TTHnonB2018UL,
    "TTMT10002018UL":TTMT10002018UL,
    "TTMT7002018UL":TTMT7002018UL,
    "TTTo2L2Nu2018UL":TTTo2L2Nu2018UL,
    "TTToHadronic2018UL":TTToHadronic2018UL,
    "TTToSemiLeptonic2018UL":TTToSemiLeptonic2018UL,
    "TTWl2018UL":TTWl2018UL,
    "TTWq2018UL":TTWq2018UL,
    "TTZM102018UL":TTZM102018UL,
    "TTZM1to102018UL":TTZM1to102018UL,
    "WJetsHT12002018UL":WJetsHT12002018UL,
    "WJetsHT2002018UL":WJetsHT2002018UL,
    "WJetsHT25002018UL":WJetsHT25002018UL,
    "WJetsHT4002018UL":WJetsHT4002018UL,
    "WJetsHT6002018UL":WJetsHT6002018UL,
    "WJetsHT8002018UL":WJetsHT8002018UL,
    "WW2018UL":WW2018UL,
    "WZ2018UL":WZ2018UL,
    "ZZ2018UL":ZZ2018UL,
}

samples_test = {
    "SingleMuonRun2018A2018UL":SingleMuonRun2018A2018UL,
}

samples_BPrime = {
    "Bprime_M1600_20UL18":Bprime_M1600_20UL18,
}

samples_WJets = {
    "WJetsHT12002018UL":WJetsHT12002018UL,
    "WJetsHT2002018UL":WJetsHT2002018UL,
    "WJetsHT25002018UL":WJetsHT25002018UL,
    "WJetsHT4002018UL":WJetsHT4002018UL,
    "WJetsHT6002018UL":WJetsHT6002018UL,
    "WJetsHT8002018UL":WJetsHT8002018UL,
}

samples_QCD = {
=======
            
singleTb = sample("singleTb", "singleTbNanoList.txt", 
                  "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8") 
singleT = sample("singleT", "singleTNanoList.txt", 
                  "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8") 
ttbar = sample("ttbar", "ttbarNanoList.txt", 
                  "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8") 
ttjetsTb = sample("ttjetsTb", "ttjetsTbNanoList.txt", 
                  "TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8") 
ttjetsT = sample("ttjetsT", "ttjetsTNanoList.txt", 
                  "TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets200 = sample("WJets200", "WJets200NanoList.txt", 
                  "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets400 = sample("WJets400", "WJets400NanoList.txt", 
                  "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets600 = sample("WJets600", "WJets600NanoList.txt", 
                  "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets800 = sample("WJets800", "WJets800NanoList.txt", 
                  "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets1200 = sample("WJets1200", "WJets1200NanoList.txt", 
                  "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets2500 = sample("WJets2500", "WJets2500NanoList.txt", 
                  "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8") 
QCD200 = sample("QCD200", "QCD200NanoList.txt", 
                "QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD300 = sample("QCD300", "QCD300NanoList.txt", 
                "QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD500 = sample("QCD500", "QCD500NanoList.txt", 
                "QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD700 = sample("QCD700", "QCD700NanoList.txt", 
                "QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD1000 = sample("QCD1000", "QCD1000NanoList.txt", 
                "QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD1500 = sample("QCD1500", "QCD1500NanoList.txt", 
                "QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD2000 = sample("QCD2000", "QCD2000NanoList.txt", 
                "QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8")
Bprime800 = sample("Bprime800", "Bprime800NanoList.txt", 
                "BprimeBToTW_M-800.txt")
Bprime1400 = sample("Bprime1400", "Bprime1400NanoList.txt", 
                "BprimeBToTW_M-1400.txt")
Bprime2000 = sample("Bprime2000", "Bprime2000NanoList.txt", 
                "BprimeBToTW_M-2000.txt")

samples={
    #"singleTb":singleTb,
    #"singleT":singleT,
    #"ttbar":ttbar,
    #"ttjetsTb":ttjetsTb,
    #"ttjetsT":ttjetsT,
    #"WJets200":WJets200,
    #"WJets400":WJets400,
    #"WJets600":WJets600,
    #"WJets800":WJets800,
    #"WJets1200":WJets1200,
    #"WJets2500":WJets2500,
>>>>>>> becca63a8b6c385fbdbc9267bf9faf41de37fe8c
    "QCD200":QCD200,
    "QCD300":QCD300,
    "QCD500":QCD500,
    "QCD700":QCD700,
    "QCD1000":QCD1000,
    "QCD1500":QCD1500,
    "QCD2000":QCD2000,
<<<<<<< HEAD
}
=======
    "Bprime800":Bprime800,
    "Bprime1400":Bprime1400,
    "Bprime2000":Bprime2000,
}
>>>>>>> becca63a8b6c385fbdbc9267bf9faf41de37fe8c
