#!/usr/bin/env python
# coding: utf-8

# # Dense neural network with Keras
# Author: Javier Duarte
# Adapted by: Sam Carlson

# ## Loading `pandas` DataFrames
# Now we load two different `NumPy` arrays. One corresponding to the VV signal and one corresponding to the background.

# In[1]:


import uproot
import numpy as np
import pandas as pd
import h5py
import pickle

# fix random seed for reproducibility
seed = 7
np.random.seed(seed)

treename = 'Events'
filename = {}
weights = {}
upfile = {}
params = {}
df = {}

print('Performing initializations...')

#These two dictionaries make things easier

weights['ttbarT'] = 1
weights['ttbarTb'] = 1
#weights['singleT'] = 0.456
#weights['singleTb'] = 0.0506
weights['WJets2500'] = 0.0011
weights['WJets1200'] = 0.0148
weights['WJets800'] = 0.0544
weights['WJets600'] = 0.1128
weights['WJets400'] = 0.4749
weights['WJets200'] = 0.4466
weights['BpM2000'] = 1
weights['BpM1400'] = 1
weights['BpM800'] = 1

filename['BpM2000'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM2000_hadd.root'
filename['BpM1400'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM1400_hadd.root'
filename['BpM800'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM800_hadd.root'
filename['WJets1200'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets1200_hadd.root'
filename['WJets2500'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets2500_hadd.root'
filename['WJets800'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets800_hadd.root'
filename['WJets600'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets600_hadd.root'
filename['WJets400'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets400_hadd.root'
filename['WJets200'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets200_hadd.root'
filename['ttbarT'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/ttbarInc_hadd.root'
filename['ttbarTb'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/ttbarInc_hadd.root'

VARS = ['pNet_J_1',#'pNet_J_2',
        'pNet_T_1',#'pNet_T_2',
        'pNet_W_1',#'pNet_W_2',
        'FatJet_pt_1',#'FatJet_pt_2',
        'FatJet_sdMass_1',#'FatJet_sdMass_2',
        'tau21_1',#'tau21_2',
        'nJ_pNet','nT_pNet','nW_pNet',
        'Jet_HT','Jet_ST','MET_pt',
        't_pt','t_mass', 'Bprime_mass',
        #'t_dRWb', # t_dRWb does not exist, should check RDF script
        'NJets_central', 'NJets_DeepFlavM','NFatJets','NJets_forward',
        'Bprime_DR','Bprime_ptbal','Bprime_chi2'] # choose which vars to use (2d)

for i,key in enumerate(filename.keys()): ## Need .keys() here, or different list
    print('Now processing file ' + key + '...')
    upfile[key] = uproot.open(filename[key])
    params[key] = upfile[key][treename].arrays(VARS)
    df[key] = pd.DataFrame(params[key],columns=VARS)
    df[key] = df[key].loc[(df[key]['Bprime_mass'] > 0) & (df[key]['NJets_forward'] > 1)] # Selection criteria

    df[key]['Weight'] = weights[key]



## cut out undefined variables VARS[0] and VARS[1] > -999
#df[key]= df[key][(df[key][VARS[0]] > -999) & (df[key][VARS[1]] > -999)]

    # add isSignal variable
    if 'Bp' in key: ## Only label Bprime as signal
        df[key]['isSignal'] = np.full(len(df[key]), 2) 
    elif 'tt' in key:
        df[key]['isSignal'] = np.ones(len(df[key]))
    else:
        df[key]['isSignal'] = np.zeros(len(df[key]))
    print(df[key].shape)
    print(df[key]['Weight'])

print('Completed loading files, pickling input data...')
df_all = pd.concat([df['BpM2000'], df['BpM1400'], df['BpM800'], df['ttbarT'],df['ttbarTb'], df['WJets200'], 
                    df['WJets400'], df['WJets600'], df['WJets800'], df['WJets1200'], df['WJets2500']])
np.savez('NewAnalysisArrays', df_all.values)

