#%%
### imports
import time
import sys
import os
import copy
import numpy as np
import matplotlib.pyplot as plt
import random
import math
from ROOT import TTree, TH1D, TFile, RDataFrame
from root_numpy import tree2array
import itertools
from collections import Counter
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn import ensemble, svm
from sklearn.metrics import f1_score, recall_score, precision_score, accuracy_score
from sklearn.calibration import CalibratedClassifierCV
from sklearn.feature_selection import SelectKBest, f_regression
import pickle
import csv

# %%
### Defining assisting functions
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

def resample_with_replacement(X_train, sample_weight):
   # normalize sample_weights if not already
   sample_weight = sample_weight / sample_weight.sum(dtype=np.float64)
   
   X_train_resampled = np.zeros((len(X_train), len(X_train[0])), dtype=np.float32)
   print('\t Resampling:')
   for i in range(len(X_train)):
      if i%10000 == 0: print('\t\t ...'+ str(i) +'...')
      # draw a number from 0 to len(X_train)-1
      draw = np.random.choice(np.arange(len(X_train)), p=sample_weight)

      # place the X at the drawn number into the resampled X
      X_train_resampled[i] = X_train[draw]

   return X_train_resampled

# This function takes in a 1D array of lists (brokenArray) and a weight. 
# The function then adds the weight to the front of the list and creates
#       a 2D array
def addWeight(brokenArray, weight):
    reshapeArray = np.zeros((brokenArray.shape[0], len(vars) + 1))
    for i, event in enumerate(brokenArray):
        eventList = list(event)
        eventList.insert(0, weight) # (0 - index (first in this case), 1 - weight being assigned)
        reshapeArray[i] = np.array(eventList)
    return reshapeArray

def plot_confusion(actual_class, pred_class, title = 'Confusion Matrix'):
   confusion = np.zeros((3, 3))
   counts = Counter(actual_class)

   for i in range(len(pred_class)):
      confusion[actual_class[i]][pred_class[i]] += 1

   for i in counts.keys():
      confusion[i][:] /= counts[i]

   fig, ax = plt.subplots()
   ax.matshow(confusion)

   labels = ['WJet', 'TTbarT', 'Bprime']
   for (i, j), z in np.ndenumerate(confusion):
      ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
   
   plt.title(title)
   plt.xticks(range(3), labels[:3])
   plt.xlabel('Predicted label')
   plt.yticks(range(3), labels[:3])
   plt.ylabel('Actual label')
   plt.savefig(outdir+'CM_' + title + outStr+'.png')
   plt.show()
   return confusion


# %%
### Reading in basic parameters
outdir = 'LargeTestOutput'
inputdir = './Output/'
arch = '3x10'
testnum = 0 # currently hard-coded to 0 until it is needed
maxtest = 300000
year = '2018'
vararray = 'shortened'
outStr = '_'+year+'BB_'+str(arch)+'_' + str(millify(maxtest)) +'test_vars'+str(vararray)+'_Test'+str(testnum)


# %%
### Logging 

# Setting up directory
outdir = outdir + '/'
if not os.path.exists(outdir):
    os.system('mkdir ' + outdir)

if not os.path.exists(outdir): 
   os.system('mkdir '+ outdir)
if not os.path.exists(outdir + '/plots'): 
   os.system('mkdir "' + outdir + '/plots"')
if not os.path.exists(outdir + '/models'): 
   os.system('mkdir "' + outdir + '/models"')
outdir = './' + outdir + '/'


# Creating file
logfile = open(outdir + 'NewTestLog.txt', 'a+')


# %% 
### Signal Selection
use2000 = True # NOTE DIFFERENT MEANING - select if trying to test 2000 vs 800 data
Bprime = 2.0 if use2000 else 0.8
Bprime2 = 0.8 if use2000 else 2.0


#%%
### Importing variables
varFile = open(inputdir + 'MLPvars.txt')
varList = list(csv.reader(varFile, delimiter=','))[0]
varList = varList[:len(varList) - 1] # importing tacks on an unneeded tail


# %%
### Getting data from ROOT files
print('Opening files...')
eosdir = "root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/"
# Defining selection criteria for the events
seltest = "NJets_forward > 0 && Bprime_mass > 0 && isValidBDecay == 1"

# Potential variable for future loop controller
filenames = ['ttbarT', 'ttbarTb', #'singleT', 'singleTb', 
             'WJets2500','WJets1200', 'WJets800', 'WJets600', 'WJets400', 'WJets200']

weights = [1, 1, #0.456, 0.0506, 
           0.0011, 0.0148, 0.0544, 0.1128, 0.4749, 0.4466]

arraysTrain = []
arraysTest = []
for i,fname in enumerate(filenames):
    sys.stdout.write('\rNow processing file {}/'.format(i + 1) + str(len(filenames)) + ' - ' + fname + '...            ')
    sys.stdout.flush()
    weight = weights[i]
    fileOpener  = TFile.Open(eosdir + fname + "_hadd.root", "READ")
    treeMaker  = fileOpener.Get("Events")
    arraysTrain.append(addWeight(tree2array(treeMaker, treeVars, seltrain), weight))
    arraysTest.append(addWeight(tree2array(treeMaker, treeVars, seltest), weight))

trainTTbarT = arraysTrain.pop(0)
testTTbarT  = arraysTest.pop(0)
trainTTbarTb = arraysTrain.pop(0)
testTTbarTb  = arraysTest.pop(0)
#trainSingleT = arraysTrain.pop(0)
#testSingleT  = arraysTest.pop(0)
#trainSingleTb = arraysTrain.pop(0)
#testSingleTb  = arraysTest.pop(0)
trainWJets2500 = arraysTrain.pop(0)
testWJets2500  = arraysTest.pop(0)
trainWJets1200 = arraysTrain.pop(0)
testWJets1200  = arraysTest.pop(0)
trainWJets800 = arraysTrain.pop(0)
testWJets800  = arraysTest.pop(0)
trainWJets600 = arraysTrain.pop(0)
testWJets600  = arraysTest.pop(0)
trainWJets400 = arraysTrain.pop(0)
testWJets400  = arraysTest.pop(0)
trainWJets200 = arraysTrain.pop(0)
testWJets200  = arraysTest.pop(0)

treeVars = varList
sys.stdout.write('\rNow processing signal file...                                                                ')
sys.stdout.flush()
weight = 1
fileBp1  = TFile.Open(eosdir + "BpM14000_hadd.root", "READ")
fileBp2 = TFile.Open(eosdir + "BpM2000_hadd.root", "READ")
# Swapping order depending on defined parameters
treeBprime = fileBp1.Get("Events") if Bprime != 2.0 else fileBp2.Get("Events")
treeBprime2 = fileBp2.Get("Events") if Bprime != 2.0 else fileBp1.Get("Events")
testBprime = addWeight(tree2array(treeBprime, treeVars, seltest), 1).tolist()
testBprime2 = addWeight(tree2array(treeBprime2, treeVars, seltest), 1).tolist()
sys.stdout.write('\rDone\n')
sys.stdout.flush()


# %%
### Creaing new arrays for added data

print('Concatenating samples...')
## Add WJets together into a single sample and reshuffle
trainWJets = np.concatenate([trainWJets200, trainWJets400, trainWJets600, trainWJets800, trainWJets1200, trainWJets2500])
testWJets = np.concatenate([testWJets200, testWJets400, testWJets600, testWJets800, testWJets1200, testWJets2500])
np.random.shuffle(trainWJets)
np.random.shuffle(testWJets)

trainTTbarT = np.concatenate([trainTTbarT, trainTTbarTb])
testTTbarT = np.concatenate([testTTbarT, testTTbarTb])
np.random.shuffle(trainTTbarT)
np.random.shuffle(testTTbarT)

#trainSingleT = np.concatenate([trainSingleT, trainSingleTb])
#testSingleT = np.concatenate([testSingleT, testSingleTb])
#np.random.shuffle(trainSingleT)
#np.random.shuffle(testSingleT)

## Print initial information to the log file and the screen
logfile.write(str(len(trainTTbarT)) + ", " + str(len(trainBprime)) + ", " +str(len(trainWJets)) + ", " +str(len(testTTbarT)) + ", " +str(len(testBprime)) + ", " +str(len(testBprime2)) + ", " +str(len(testWJets))) # + ", " + str(len(testSingleT) + ", " + str(len(trainSingleT)) ))

np.random.shuffle(testTTbarT)
np.random.shuffle(testBprime)
np.random.shuffle(testWJets)
np.random.shuffle(testBprime2)

# %%
### Cleaning the imported arrays
cleanedList = []
for i, row in enumerate(testWJets):
   if np.inf in row or -np.inf in row or np.nan in row:
      cleanedList.pop(i)
testWJets = cleanedList

cleanedList = []
for i, row in enumerate(testTTbarT):
   if np.inf in row or -np.inf in row or np.nan in row:
      cleanedList.pop(i)
testTTbarT = cleanedList

cleanedList = []
for i, row in enumerate(testBprime):
   if np.inf in row or -np.inf in row or np.nan in row:
      cleanedList.pop(i)
testBprime = cleanedList

cleanedList = []
for i, row in enumerate(testBprime2):
   if np.inf in row or -np.inf in row or np.nan in row:
      cleanedList.pop(i)
testBprime2 = cleanedList


# %%
### Loading model and scaler from pickles
mlp = pickle.load(open(inputdir + 'models/MLP' +outStr+'.pkl'))
scaler = pickle.load(open(outdir+'models/Dnn_scaler_3bin'+outStr+'.pkl'))


# %%
### Scaling data
testBprime2 = scaler.transform(testBprime2)
testBprime = scaler.transform(testBprime)
testTTbarT = scaler.transform(testTTbarT)
testWJets = scaler.transform(testWJets)


# %%
### Getting probabilities from the classifier
print('\n--------------Evaluation of Model--------------')

# Get scores for non-training events on MLP
probs_WJetsMLP = mlp.predict_proba(testWJets)
probs_TTbarTMLP = mlp.predict_proba(testTTbarT)
probs_BprimeMLP = mlp.predict_proba(testBprime)
probs_Bprime2MLP = mlp.predict_proba(testBprime2)
# probs = [probsWJets, probsTTbarT, probsBprime]

# %%
### Plotting the comparison 
## WJets
# plt.close()
if not os.path.exists(outdir + 'plots'): os.system('mkdir '+outdir + 'plots')
plt.figure()
plt.xlabel('Predicted W boson score - MLP',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_WJetsMLP.T[0], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_TTbarTMLP.T[0], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_BprimeMLP.T[0], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.hist(probs_Bprime2MLP.T[0], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime2)+'\,TeV)}$', color='c', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig(outdir+'plots/score_WJetMLP'+outStr+'.png')

## TTbarT
# plt.close()
plt.figure()
plt.xlabel('Predicted top quark score - MLP',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_WJetsMLP.T[1], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_TTbarTMLP.T[1], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_BprimeMLP.T[1], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.hist(probs_Bprime2MLP.T[1], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime2)+'\,TeV)}$', color='c', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig(outdir+'plots/score_TTbarTMLP'+outStr+'.png')

## Signal
# plt.close()
plt.figure()
plt.xlabel('Predicted B quark score - MLP',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_WJetsMLP.T[2], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_TTbarTMLP.T[2], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_BprimeMLP.T[2], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.hist(probs_Bprime2MLP.T[2], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime2)+'\,TeV)}$', color='c', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig(outdir+'plots/score_BprimeMLP'+outStr+'.png')