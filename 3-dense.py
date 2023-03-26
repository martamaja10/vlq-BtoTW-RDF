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
upfile = {}
params = {}
df = {}

print('Performing initializations...')
# These two lists make things easier later
filenames = ['ttbarT', 'ttbarTb', #'singleT', 'singleTb', 
             'WJets2500','WJets1200', 'WJets800', 'WJets600', 'WJets400', 'WJets200',
             'BpM2000', 'BpM14000', 'BpM800']

weights = [1, #0.456, 0.0506, 
           0.0011, 0.0148, 0.0544, 0.1128, 0.4749, 0.4466, 1, 1, 1]

filename['BpM2000'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM2000_hadd.root'
filename['BpM14000'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM14000_hadd.root'
filename['BpM800'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM800_hadd.root'
filename['WJets1200'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets1200_hadd.root'
filename['WJets2500'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets2500_hadd.root'
filename['WJets800'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets800_hadd.root'
filename['WJets600'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets600_hadd.root'
filename['WJets400'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets400_hadd.root'
filename['WJets200'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets200_hadd.root'
filename['ttbarT'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/ttbarT_hadd.root'
filename['ttbarTb'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/ttbarTb_hadd.root'

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

for i,key in enumerate(filenames): ## Need .keys() here, or different list
    print('Now open file ' + key + '...')
    upfile[key] = uproot.open(filename[key])
    print('Converting to arrays...')
    params[key] = upfile[key][treename].arrays(VARS)

    print('Generating dataframe...')
    df[key] = pd.DataFrame(params[key],columns=VARS)
    df[key] = df[key].loc[(df[key]['Bprime_mass'] > 0) & (df[key]['NJets_forward'] > 1)] # Selection criteria

    df[key]['Weight'] = weights[i]


## cut out undefined variables VARS[0] and VARS[1] > -999
#df[key]= df[key][(df[key][VARS[0]] > -999) & (df[key][VARS[1]] > -999)]

    # add isSignal variable
    if 'Bp' in key: ## Only label Bprime as signal
        df[key]['isSignal'] = np.full(len(df[key]), 2) 
    elif 'tt' in key:
        df[key]['isSignal'] = np.ones(len(df[key]))
    else:
        df[key]['isSignal'] = np.zeros(len(df[key]))
    df[key].to_pickle(key + '.pkl')

print('Completed loading files, pickling input data...')
