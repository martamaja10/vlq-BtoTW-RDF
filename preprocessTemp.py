import time
import sys
import os
import numpy as np
from ROOT import TTree, TH1D, TFile, RDataFrame
from root_numpy import tree2array

# %%
### Reading in basic parameters
start_time = time.time() # collected just for benchmarking
outdirName = sys.argv[1] # user can define this at runtime
testnum = 0 # currently hard-coded to 0 until it is needed


# %%
### Logging 

# Setting up directory
outdirName = outdirName + '/'
if not os.path.exists(outdirName):
    os.system('mkdir ' + outdirName)

# Creating file
logfile = open(outdirName + 'SingleBLog.txt', 'a+')


# %% 
### Signal Selection

# Defining weight classes
Bprime = 0.8
Bprime2 = 2.0

# %%
### Defining variables to be used with model (defined in RDataframe script)
vars = ['pNet_J_1','pNet_J_2',
        'pNet_T_1','pNet_T_2',
        'pNet_W_1','pNet_W_2',
        'dpak8_J_1','dpak8_J_2',
        'dpak8_T_1','dpak8_T_2',
        'dpak8_W_1','dpak8_W_2',
        'FatJet_pt_1','FatJet_pt_2',
        'FatJet_sdMass_1','FatJet_sdMass_2',
        'tau21_1','tau21_2',
        'nJ_dpak8','nT_dpak8','nW_dpak8',
        'nJ_pNet','nT_pNet','nW_pNet',
        'Jet_HT','Jet_ST','MET_pt',
        't_pt','t_mass',
        #'t_dRWb', # t_dRWb does not exist, should check RDF script
        'NJets_central', 'NJets_DeepFlavM','NFatJets','NJets_forward',
        'Bprime_DR','Bprime_ptbal','Bprime_chi2',
        'minDR_leadAK8otherAK8'] 

# %%
### Getting data from ROOT files
print('Opening files...')
eosdir = "root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/"

# Defining selection criteria for the events
seltrain = "NJets_forward == 0"
seltest = "NJets_forward != 0"

treeVars = vars

## Performing selections and reads

# Selections on ttbarT(b) and singlet(b)
fileTTbarT  = TFile.Open(eosdir + "ttbarT_hadd.root", "READ")
treeTTbarT  = fileTTbarT.Get("Events")
trainTTbarT = tree2array(treeTTbarT, treeVars, seltrain)
testTTbarT  = tree2array(treeTTbarT, treeVars, seltest)
print(trainTTbarT[0])
for i, event in enumerate(trainTTbarT):
    event.insert(0, 1)
    trainTTbarT[i] = event
for i, event in enumerate(testTTbarT):
    event.insert(0, 1)
    testTTbarT[i] = event
print(trainTTbarT[0])

fileTTbarTb  = TFile.Open(eosdir + "ttbarTb_hadd.root", "READ")
treeTTbarTb  = fileTTbarTb.Get("Events")
trainTTbarTb = tree2array(treeTTbarTb, treeVars, seltrain)
testTTbarTb  = tree2array(treeTTbarTb, treeVars, seltest)

fileSingleT  = TFile.Open(eosdir + "singleT_hadd.root", "READ")
treeSingleT  = fileSingleT.Get("Events")
trainSingleT = tree2array(treeSingleT, treeVars, seltrain)
testSingleT  = tree2array(treeSingleT, treeVars, seltest)

fileSingleTb  = TFile.Open(eosdir + "singleTb_hadd.root", "READ")
treeSingleTb  = fileSingleTb.Get("Events")
trainSingleTb = tree2array(treeSingleTb, treeVars, seltrain)
testSingleTb  = tree2array(treeSingleTb, treeVars, seltest)

# WJet selection
fileWJets2500  = TFile.Open(eosdir + "WJets2500_hadd.root", "READ")
treeWJets2500  = fileWJets2500.Get("Events")
trainWJets2500 = tree2array(treeWJets2500, treeVars, seltrain)
testWJets2500  = tree2array(treeWJets2500, treeVars, seltest)

fileWJets1200  = TFile.Open(eosdir + "WJets1200_hadd.root", "READ")
treeWJets1200  = fileWJets1200.Get("Events")
trainWJets1200 = tree2array(treeWJets1200, treeVars, seltrain)
testWJets1200  = tree2array(treeWJets1200, treeVars, seltest)

fileWJets800  = TFile.Open(eosdir + "WJets800_hadd.root", "READ")
treeWJets800  = fileWJets800.Get("Events")
trainWJets800 = tree2array(treeWJets800, treeVars, seltrain)
testWJets800  = tree2array(treeWJets800, treeVars, seltest)

fileWJets600  = TFile.Open(eosdir + "WJets600_hadd.root", "READ")
treeWJets600  = fileWJets600.Get("Events")
trainWJets600 = tree2array(treeWJets600, treeVars, seltrain)
testWJets600  = tree2array(treeWJets600, treeVars, seltest)

fileWJets400  = TFile.Open(eosdir + "WJets400_hadd.root", "READ")
treeWJets400  = fileWJets400.Get("Events")
trainWJets400 = tree2array(treeWJets400, treeVars, seltrain)
testWJets400  = tree2array(treeWJets400, treeVars, seltest)

fileWJets200  = TFile.Open(eosdir + "WJets200_hadd.root", "READ")
treeWJets200  = fileWJets200.Get("Events")
trainWJets200 = tree2array(treeWJets200, treeVars, seltrain)
testWJets200  = tree2array(treeWJets200, treeVars, seltest)

# Selection with signals
fileBp1  = TFile.Open(eosdir + "Bp800_hadd.root", "READ")
fileBp2 = TFile.Open(eosdir + "Bp1200_hadd.root", "READ")
# Swapping order depending on defined parameters
if Bprime == 2.0:
    temp = fileBp1
    fileBp1 = fileBp2
    fileBp2 = temp
treeBprime = fileBp1.Get("Events")
treeBprime2 = fileBp2.Get("Events")
trainBprime= tree2array(treeBprime, treeVars, seltrain)
testBprime= tree2array(treeBprime, treeVars, seltest)
testBprime2= tree2array(treeBprime2, treeVars, seltest)



# %%
