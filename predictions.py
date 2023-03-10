#To be run after 3-dense.py, imports 'model' and 'scaler'
#Opens all the .root files in the /eos/ storage area and applies prediction to all events

# In[1]:
import uproot
import numpy as np
from root_numpy import root2array, array2root
from tensorflow import keras
import pickle
from sklearn.preprocessing import StandardScaler
filename = {}

#import scaler and model
scaler = pickle.load(open('/uscms/home/khowey/nobackup/BtoTW/CMSSW_11_0_0/src/vlq-BtoTW-RDF/NewAnalysisModels/dnn_scaler.pkl'))
model = keras.models.load_model('/uscms/home/khowey/nobackup/BtoTW/CMSSW_11_0_0/src/vlq-BtoTW-RDF/NewAnalysisModels/MLP.h5')

#define dictionary with all events
filename['BpM2000'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM2000_hadd.root'
#filename['BpM1400'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM14000_hadd.root'
#filename['BpM800'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM800_hadd.root'
#filename['WJets1200'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets1200_hadd.root'
#filename['WJets2500'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets2500_hadd.root'
#filename['WJets800'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets800_hadd.root'
#filename['WJets600'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets600_hadd.root'
#filename['WJets400'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets400_hadd.root'
#filename['WJets200'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets200_hadd.root'
#filename['WJetsInc'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJetsInc_hadd.root'
#filename['singleT'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/singleT_hadd.root'
#filename['singleTb'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/singleTb_hadd.root'
filename['ttbarT'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/ttbarT_hadd.root'
#filename['ttbarTb'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/ttbarTb_hadd.root'

#this list of columns must match 3-dense.py
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
        'Bprime_DR','Bprime_ptbal','Bprime_chi2']

def get_features_from_file(filename='', treename='', branches=[]):
    t = root2array(filename, treename=treename, branches=branches) # structured numpy array 
    #print t.shape 
    t = t.view(np.float32).reshape(t.shape + (-1,)) # normal numpy array (trick from https://stackoverflow.com/questions/5957380/convert-structured-array-to-regular-numpy-array)
    t = t[~np.isinf(t).any(axis=1)] # Filters out infinite values
    #print t.shape
    return t

def write_prediction_to_file(features, model, filename='',treename='',branch=''):
    y_predict_all = model.predict(features) # normal numpy array
    #print y_predict_all.shape
    y_predict_all = np.array(y_predict_all, dtype=[(branch, np.float32)]) # structured numpy array
    #print y_predict_all.shape
    array2root(y_predict_all, filename, treename=treename, mode='recreate')

#running the get_features, scaler.transform, and write functions over all .root files
for key in filename.keys():
    print(key)
    X_all = get_features_from_file(filename[key], 
                                   treename='Events', 
                                   branches=VARS)
    X_all = scaler.transform(X_all)
    write_prediction_to_file(X_all, 
                             model, 
                             filename[key].replace('hadd','predict').replace('/jmanagan/BtoTW_RDF','/samuelca'), 
                             treename='Events', 
                             branch='dense')
print('Complete!')
