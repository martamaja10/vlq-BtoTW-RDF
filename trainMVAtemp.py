# %%
### imports

# external modules
import os
import sys
import time
import numpy as np
import math
import matplotlib.pyplot as plt
import importlib
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow import keras
from keras import backend as K
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.layers import Input, Dense, Concatenate
from tensorflow.keras.models import Model, Sequential, load_model
import importlib
from sklearn.preprocessing import StandardScaler
rom sklearn import svm, metrics, preprocessing
from plot_confusion_matrix import plot_confusion_matrix
import copy, time, os, sys, math, random, itertools
from sklearn import neural_network
from sklearn.linear_model import SGDClassifier
from sklearn import tree
from sklearn.metrics import roc_curve
from sklearn.externals import joblib
from sklearn.metrics import classification_report, f1_score, recall_score, precision_score, accuracy_score

# %% 
### Define helper functions
def millify(n):
   n = float(n)
   millnames = ['','k','M','G','T']
   millidx = max(0,min(len(millnames)-1,
                       int(math.floor(0 if n == 0 else math.log10(abs(n))/3))))

   return '{:.0f}{}'.format(n / 10**(3 * millidx), millnames[millidx])

def arch2tuple(n):
   layers, nodes = n.split('x',2)
   tup = (int(nodes),)
   out =()
   for x in range(int(layers)-1):
      out = out + tup
   return (out)

# %%
### User parameters
start_time = time.time()
arch = '3x10'
maxtest = 15000
outdir = sys.argv[1]
vararray = int(sys.argv[2])
testnum = int(sys.argv[3])
year = str(sys.argv[4])
if year == 'all': maxtest = 30000

# %%
### Configure logs
outdir = outdir + '/'

# Check if log file already exists
if testnum == 1:
   logfile = open(outdir + "NN_logs" + year + ".txt", "a+")
   logfile.write('\ntest, vararray, Testing Score (Accuracy), tt-as-BB, BB-as-BB, Precision, Recall, F-Score \n')
else:
   time.sleep(2)
   logfile = open(outdir+"BB_output_Jan21_"+year+".txt","a+")
   logfile.write('\n')

logfile.write(str(testnum)+", ")
logfile.write(str(vararray)+", ")

# %%
### Signal Selection
Bprime = 0.8
Bprime2 = 2.0
test2000 = False #use if Bprime = 2000

# %%
### Configure output
outStr = '_'+year+'BB_'+str(arch)+'_' + str(millify(maxtest)) +'test_vars'+str(vararray)+'_Test'+str(testnum)
print('Outstr:',outStr,'Outdir:',outdir)
if not os.path.exists(outdir): os.system('mkdir '+outdir)

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
### Perform permutations calculations
# TODO - Figure out what items must be present
combos = []
for size in range(7,len(varList)):
   thissize = list(itertools.combinations(varList,size))
   for item in thissize:
      ## Some sanity checks for vars we know for sure we want to include
      if 'dnnJ_1' not in item: continue
      if 'dnnJ_2' not in item and 'dnnJ_3' not in item: continue
      if 'AK4HT' not in item: continue
      if 'corr_met_MultiLepCalc' not in item and 'AK4HTpMETpLepPt' not in item: continue
      if 'NJetsDeepFlavwithSF_JetSubCalc' not in item and 'NJets_JetSubCalc' not in item: continue
      if 't_mass' not in item and 't_pt' not in item: continue
      if 't_dRWb' not in item and 'minDR_leadAK8otherAK8' not in item: continue
      combos.append(item)
combos.append(varList)

# %%
### Configure inputs
print('Reading in input data...')
indexKill = range(0,len(varList))
for item in vars:
   indexKill.remove(varList.index(item))

inStr = '_'+year+'BB_'+str(arch)+'_' + str(millify(maxtest)) +'test'
allmystuff = np.load(outdir+'Arrays'+inStr+'.npz')

trainData = (allmystuff['trainData']).tolist()
trainLabel = (allmystuff['trainLabel']).tolist()
testData = (allmystuff['testData']).tolist()
testLabel = (allmystuff['testLabel']).tolist()
testWJets = (allmystuff['testWJets']).tolist()
testTTbarT = (allmystuff['testTTbarT']).tolist()
testSingleT = (allmystuff['testSingleT']).tolist()
testBprime = (allmystuff['testBprime']).tolist()
testBprime2 = (allmystuff['testBprime2']).tolist()

## Get rid of the variables we aren't using this time
for i in range(0,len(trainData)):
   for j in sorted(indexKill, reverse = True):
      del trainData[i][j]
for i in range(0,len(testData)):
   for j in sorted(indexKill, reverse = True):
      del testData[i][j]
for i in range(0,len(testWJets)):
   for j in sorted(indexKill, reverse = True):
      del testWJets[i][j]
      del testBprime[i][j]
      del testBprime2[i][j]
      del testTTbarT[i][j]
      del SingleT[i][j]

# %% 
### Perform scaling
print('Building the scaler...')
scaler = preprocessing.StandardScaler().fit(trainData)
print('Transforming...')
trainData = scaler.transform(trainData)
testData = scaler.transform(testData)
testBprime2 = scaler.transform(testBprime2)
testBprime = scaler.transform(testBprime)
testTTbarT = scaler.transform(testTTbarT)
testSingleT = scaler.transform(testSingleT)
testWJets = scaler.transform(testWJets)

# %% 
### Neural Network Training
print('Training...')
mlp = neural_network.MLPClassifier(hidden_layer_sizes=arch2tuple(arch), activation='relu',early_stopping=True)
mlp.fit(trainData, trainLabel)

# %%
### Scores
print('Test data score =',mlp.score(testData, testLabel))
print('Train data score =',mlp.score(trainData, trainLabel))
logfile.write(str(round(mlp.score(testData, testLabel),5)) + ", ")

# %%
### Plot a confusion matrix of the scores
cm = metrics.confusion_matrix(mlp.predict(testData), testLabel)
plt.figure()
targetNames = [r'$\mathrm{W+jets}$',r'$\mathrm{Bprime}$',r'$\mathrm{T\overline{T}}$',r'$\mathrm{SingleT}$']
plot_confusion_matrix(cm.T, targetNames, normalize=True)

cm = (cm.T).astype('float') / (cm.T).sum(axis=1)[:, np.newaxis]

ofile.write(str(round(cm[1][2],5)) + ", " + str(round(cm[2][2],5)) + ", ")

