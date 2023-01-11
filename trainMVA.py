# %%
### imports

# external modules
import os
import time
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from collections import Counter
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn import ensemble, svm
from sklearn.metrics import f1_score, recall_score, precision_score, accuracy_score
from sklearn.calibration import CalibratedClassifierCV
from sklearn.feature_selection import SelectKBest, f_regression
import pickle

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
### User parameters
start_time = time.time()
arch = '3x10'
maxtest = 300000

outdir = 'Output'
if not os.path.exists(outdir): 
   os.system('mkdir '+ outdir)
outdir = './' + outdir + '/'

inputdir = './Input'
inputdir += '/'

vararray = 'test'
testnum = 1
year = '2018'
feature_selection = True
dpak8 = True
cnn = False
if len(sys.argv) > 1:
    outdir = sys.argv[1]
    vararray = int(sys.argv[2])
    testnum = int(sys.argv[3])
    year = str(sys.argv[4])
if year == 'all': maxtest = 30000

# %%
### Configure logs
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
Bprime = 2.0
Bprime2 = 0.8
test2000 = True #use if Bprime = 2000

# %%
### Configure output
outStr = '_'+year+'BB_'+str(arch)+'_' + str(millify(maxtest)) +'test_vars'+str(vararray)+'_Test'+str(testnum)
print('Outstr:',outStr,'Outdir:',outdir)

# %%
### Defining variables to be used with model (defined in RDataframe script)
varList = ['pNet_J_1',#'pNet_J_2',
      'pNet_T_1',#'pNet_T_2',
      'pNet_W_1',#'pNet_W_2',
      'dpak8_J_1',#'dpak8_J_2',
      'dpak8_T_1',#'dpak8_T_2',
      'dpak8_W_1',#'dpak8_W_2',
      'FatJet_pt_1',#'FatJet_pt_2',
      'FatJet_sdMass_1',#'FatJet_sdMass_2',
      'tau21_1',#'tau21_2',
      'nJ_dpak8','nT_dpak8','nW_dpak8',
      'nJ_pNet','nT_pNet','nW_pNet',
      'Jet_HT','Jet_ST','MET_pt',
      't_pt','t_mass',
      #'t_dRWb', # t_dRWb does not exist, should check RDF script
      'NJets_central', 'NJets_DeepFlavM','NFatJets','NJets_forward',
      'Bprime_DR','Bprime_ptbal','Bprime_chi2',
      #'minDR_leadAK8otherAK8'
      ] 


# %%
### Importing data
inStr = '_'+year+'BB_'+str(arch)+'_' + str(millify(maxtest)) +'test'
print('Loading from file ' + inputdir + 'Arrays' + inStr + '.npz...')
allmystuff = np.load(inputdir+'Arrays'+inStr+'.npz')

trainData = (allmystuff['trainData']).tolist()
trainLabel = (allmystuff['trainLabel']).tolist()
testData = (allmystuff['testData']).tolist()
testLabel = (allmystuff['testLabel']).tolist()
testWJets = (allmystuff['testWJets']).tolist()
testTTbarT = (allmystuff['testTTbarT']).tolist()
testBprime = (allmystuff['testBprime']).tolist()
testBprime2 = (allmystuff['testBprime2']).tolist()

trainLabel2Txt = ['WJets', 'TTbarT', 'Signal']

# Fixes an old bug, so should be unnecessary post December 2022
if len(testLabel) == 0:
   testLabel = trainLabel[len(trainData):]
   trainLabel = trainLabel[:len(trainData)]

# Remove invalid rows
nInvalidRow = 0
for i,row in enumerate(trainData):
   if np.inf in row or -np.inf in row or np.nan in row:
      trainData.pop(i)
      trainLabel.pop(i)
      nInvalidRow += 1
if nInvalidRow > 0: print('Encountered and removed {} invalid train row(s).'.format(nInvalidRow))

nInvalidRow = 0
for i, row in enumerate(testData):
   if np.inf in row or -np.inf in row or np.nan in row:
      testData.pop(i)
      testLabel.pop(i)
      nInvalidRow += 1
if nInvalidRow > 0: print('Encountered and removed {} invalid test row(s).'.format(nInvalidRow))

for i, row in enumerate(testWJets):
   if np.inf in row or -np.inf in row or np.nan in row:
      testWJets.pop(i)

for i, row in enumerate(testTTbarT):
   if np.inf in row or -np.inf in row or np.nan in row:
      testTTbarT.pop(i)

print('Training on ' + str(len(trainData)) + ' events.')
print('Testing on {} events of the following makeup:'.format(len(testData)))
print('{} WJet events'.format(len(testWJets)))
print('{} TTbarT events'.format(len(testTTbarT)))
print('{} {} TeV B events'.format(len(testBprime), Bprime))
print('{} {} TeV B events'.format(len(testBprime2), Bprime2))



## Code for selecting out deepak8
if not dpak8:
   cols = []
   newVarList = []
   for i, var in enumerate(varList):
      if not 'dpak8' in var:
         cols.append(i)
         newVarList.append(var)

   varList = newVarList

   print('Selecting out deepak8 features from training data...')
   selectedTrain = []
   for event in trainData:
      newEvent = []
      for col in cols:
         newEvent.append(event[col])
      selectedTrain.append(newEvent)

   print('Selecting out deepak8 features from testing data...')
   selectedTest = []
   for event in testData:
      newEvent = []
      for col in cols:
         newEvent.append(event[col])
      selectedTest.append(newEvent)

   selectedBprime2 = []
   for event in testBprime2:
      newEvent = []
      for col in cols:
         newEvent.append(event[col])
      selectedBprime2.append(newEvent)

   selectedBprime = []
   for event in testBprime:
      newEvent = []
      for col in cols:
         newEvent.append(event[col])
      selectedBprime.append(newEvent)

   selectedTTbarT = []
   for event in testTTbarT:
      newEvent = []
      for col in cols:
         newEvent.append(event[col])
      selectedTTbarT.append(newEvent)

   selectedWJets = []
   for event in testWJets:
      newEvent = []
      for col in cols:
         newEvent.append(event[col])
      selectedWJets.append(newEvent)

   trainData = selectedTrain
   testData = selectedTest
   testBprime2 = selectedBprime2
   testBprime = selectedBprime
   testTTbarT = selectedTTbarT
   testWJets = selectedWJets

# %% 
### Perform scaling
print('\nBuilding the scaler...')
scaler = StandardScaler().fit(trainData)
print('Transforming...')
trainData = scaler.transform(trainData)
testData = scaler.transform(testData)
testBprime2 = scaler.transform(testBprime2)
testBprime = scaler.transform(testBprime)
testTTbarT = scaler.transform(testTTbarT)
testWJets = scaler.transform(testWJets)


# %%
### Training a basic decision tree with SKlearn
print('\n--------------Random Forest Feature Selection--------------')
tstart = time.time()
dtModel = ensemble.RandomForestClassifier(random_state = 42, n_estimators = 100)
dtModel.fit(trainData, trainLabel)
importance = dtModel.feature_importances_
for i, value in enumerate(importance):
   print('Feature: %0d, Score: %.5f' % (i, value))
plt.bar([x for x in range (len(importance))], importance)
plt.show()


# %%
### Training a basic MLP with SKlearn

print('\n--------------Training Multilayer Perceptron--------------')
tstart = time.time()
mlp = MLPClassifier(max_iter = 500, solver = 'adam', activation = 'relu', alpha = 1e-5, 
      hidden_layer_sizes = (25, 100), random_state = 42, shuffle = True, verbose = False,
      early_stopping = True, validation_fraction = 0.3)
mlp.fit(trainData, trainLabel)
mlpTime = time.time() - tstart
print(mlpTime)

# MLP loss curve
losscurve = mlp.loss_curve_
plt.figure()
plt.xlabel('iterations')
plt.ylabel('training loss')
plt.plot(losscurve)
plt.savefig(outdir+'MLP_trainloss'+outStr+'.png')

## Draw validation sample (10% of the training data) score
testscore = mlp.validation_scores_
plt.figure()
plt.xlabel('iterations')
plt.ylabel('validation score')
plt.plot(testscore)
plt.savefig(outdir+'MLP_valscore'+outStr+'.png')
plt.show()

preds = mlp.predict(testData)
plot_confusion(testLabel, preds, title = 'MLP')

# %%
### Getting probabilities from the classifiers
print('\n--------------Evaluation of Model--------------')

# Get scores for non-training events on MLP
probs_WJetsMLP = mlp.predict_proba(testWJets)
probs_TTbarTMLP = mlp.predict_proba(testTTbarT)
probs_BprimeMLP = mlp.predict_proba(testBprime)
probs_Bprime2MLP = mlp.predict_proba(testBprime2)
# probs = [probsWJets, probsTTbarT, probsBprime]
probsMLP = [probs_WJetsMLP, probs_TTbarTMLP, probs_Bprime2MLP]

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

# %%
### Getting metrics
accMLP = accuracy_score(testLabel, mlp.predict(testData))
precMLP = precision_score(testLabel, mlp.predict(testData), average='weighted')
recallMLP = recall_score(testLabel, mlp.predict(testData), average='weighted')
fscoreMLP = f1_score(testLabel, mlp.predict(testData), average='weighted')

accDT = accuracy_score(testLabel, dtModel.predict(testData))
precDT = precision_score(testLabel, dtModel.predict(testData), average='weighted')
recallDT = recall_score(testLabel, dtModel.predict(testData), average='weighted')
fscoreDT = f1_score(testLabel, dtModel.predict(testData), average='weighted')

print('------MLP------')
print('Accuracy: ' + str(accMLP))
print('Precision: ' + str(precMLP))
print('Recall: ' + str(recallMLP))
print('F-measure: ' + str(fscoreMLP))
print('Trained in ' + str(mlpTime) + ' s')

# %%
## Saving models to files
if not os.path.exists(outdir + 'models'): os.system('mkdir '+outdir + 'models')
pickle.dump(mlp, open(outdir+'models/MLP'+outStr+'.pkl', 'wb'))
pickle.dump(dtModel, open(outdir+'models/DT' + outStr +'.pkl', 'wb'))
pickle.dump(scaler, open(outdir+'models/Dnn_scaler_3bin'+outStr+'.pkl', 'wb'))