# %%
### imports

# external modules
import os
import time
import numpy as np
import math
import matplotlib.pyplot as plt
from collections import Counter
import tensorflow as tf
import tempfile
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.layers import Conv1D, Dense, MaxPooling1D, Flatten, Dropout
from keras.models import Sequential, Model, save_model, load_model
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.neural_network import MLPClassifier
from sklearn import ensemble
from sklearn.metrics import f1_score, recall_score, precision_score, accuracy_score
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

   for (i, j), z in np.ndenumerate(confusion):
      ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
   plt.title(title)
   plt.xlabel('Predicted label')
   plt.ylabel('Actual label')
   plt.show()
   return confusion

# %%
### User parameters
start_time = time.time()
arch = '3x10'
maxtest = 300000
outdir = './temp'
vararray = 'test'
testnum = 1
year = '2018'
# if len(sys.argv) > 1:
#    outdir = sys.argv[1]
#    vararray = int(sys.argv[2])
#    testnum = int(sys.argv[3])
#    year = str(sys.argv[4])
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
Bprime = 2.0
Bprime2 = 0.8
test2000 = False #use if Bprime = 2000

# %%
### Configure output
outStr = '_'+year+'BB_'+str(arch)+'_' + str(millify(maxtest)) +'test_vars'+str(vararray)+'_Test'+str(testnum)
print('Outstr:',outStr,'Outdir:',outdir)
if not os.path.exists(outdir): os.system('mkdir '+outdir)

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
print('Loading from file ' + outdir + 'Arrays' + inStr + '.npz...')
allmystuff = np.load(outdir+'Arrays'+inStr+'.npz')

trainData = (allmystuff['trainData']).tolist()
trainLabel = (allmystuff['trainLabel']).tolist()
testData = (allmystuff['testData']).tolist()
testLabel = (allmystuff['testLabel']).tolist()
testWJets = (allmystuff['testWJets']).tolist()
testTTbarT = (allmystuff['testTTbarT']).tolist()
testBprime = (allmystuff['testBprime']).tolist()
testBprime2 = (allmystuff['testBprime2']).tolist()

trainLabel2Txt = ['WJets', 'TTbarT', 'Signal']

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
      
# %%
### Evaluation of features
trainSelector = SelectKBest(f_regression, k=15).fit(trainData, trainLabel)
cols = trainSelector.get_support(indices = True).tolist()

selectedFeatures = []
for col in cols:
   selectedFeatures.append(varList[col])
print('\nSelected the following features for training:')
print(selectedFeatures)

# Eliminating unhelpful features
for i, col in enumerate(cols):
   if col <=5 and col > 2:
      cols.pop(i)
   if col <=11 and col > 8:
      cols.pop(i)

print('Selecting out unhelpful features from training data...')
selectedTrain = []
for event in trainData:
   newEvent = []
   for col in cols:
      newEvent.append(event[col])
   selectedTrain.append(newEvent)

print('Selecting out unhelpful features from testing data...')
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

# %% 
### Perform scaling
print('\nBuilding the scaler...')
scaler = StandardScaler().fit(selectedTrain)
print('Transforming...')
trainData = scaler.transform(selectedTrain)
testData = scaler.transform(selectedTest)
testBprime2 = scaler.transform(selectedBprime2)
testBprime = scaler.transform(selectedBprime)
testTTbarT = scaler.transform(selectedTTbarT)
testWJets = scaler.transform(selectedWJets)

# %%
### Training a basic MLP with SKlearn

print('\n--------------Training Multilayer Perceptron--------------')
tstart = time.time()
mlp = MLPClassifier(max_iter = 500, solver = 'adam', activation = 'relu', alpha = 1e-5, 
      hidden_layer_sizes = (25, 100), random_state = 1, shuffle = True, verbose = False,
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
plt.savefig(outdir+'trainloss'+outStr+'.png')
plt.close()

## Draw validation sample (10% of the training data) score
testscore = mlp.validation_scores_
plt.figure()
plt.xlabel('iterations')
plt.ylabel('validation score')
plt.plot(testscore)
plt.savefig(outdir+'valscore'+outStr+'.png')
plt.close()

ConfusionMatrixDisplay.from_estimator(mlp, testData, testLabel, normalize = 'true')
# %%
### Training a basic decision tree with SKlearn
print('\n--------------Random Forest Classifier--------------')
tstart = time.time()
dtModel = ensemble.RandomForestClassifier(random_state = 0, n_estimators = 100)
dtModel.fit(trainData, trainLabel)
dtTime = time.time() - tstart
print(dtTime)
ConfusionMatrixDisplay.from_estimator(dtModel, testData, testLabel, normalize = 'true')

# %%
### Train CNN with Tensorflow
print('\n--------------Training Convolutional Neural Network--------------')
tstart = time.time()

trainAr = np.array(trainData)
labelAr = np.array(trainLabel)
# Reshaping train data for input to CNN
sample_size = trainAr.shape[0]
time_steps = trainAr.shape[1]
input_dimension = 1
callbackList = [EarlyStopping(monitor = 'val_loss', mode = 'min', patience = 25, verbose = 1),
               ModelCheckpoint('./ModelChkpt/checkpoint', save_weights_only=True, monitor='val_accuracy', mode = 'max', save_best_only = True)]
trainCNN = trainAr.reshape(sample_size, time_steps, input_dimension)
tf.random.set_seed(42) # ensures models are somewhat consistent

# Building the model
cnn = Sequential(name="model_conv1D")
cnn.add(Conv1D(filters=128, kernel_size = 8, activation = 'relu', input_shape = (trainCNN.shape[1], trainCNN.shape[2])))
cnn.add(Dropout(0.5))
cnn.add(Conv1D(filters=256, kernel_size = 5, activation = 'relu'))
cnn.add(Dropout(0.2))
cnn.add(Dense(16, activation="relu"))
cnn.add(MaxPooling1D())
cnn.add(Flatten())
cnn.add(Dense(100, activation="relu"))
cnn.add(Dense(25, activation = 'softmax'))
cnn.compile(loss = 'sparse_categorical_crossentropy', optimizer = 'adam', metrics = ['accuracy'])

cnn.fit(trainCNN, labelAr, epochs = 500, validation_split = 0.1, verbose = 0,
   callbacks = callbackList)
cnn.load_weights('./ModelChkpt/checkpoint')
cnnTime = time.time() - tstart
print(cnnTime)
# %%
### Getting predictions from CNN
testAr = np.array(testData)
labelTAr = np.array(testLabel)

sample_test = testAr.shape[0]
time_test = trainAr.shape[1]

testCNN = testAr.reshape(sample_test, time_test, input_dimension)

pred_cnn = cnn.predict(testCNN).argmax(axis = -1)

confusion = plot_confusion(labelTAr, pred_cnn, title = 'CNN Confusion Matrix')
# %%
### Getting probabilities from the classifiers
print('\n--------------Evaluation of Models--------------')

# Get scores for non-training events on MLP
probs_WJetsMLP = mlp.predict_proba(testWJets)
probs_TTbarTMLP = mlp.predict_proba(testTTbarT)
probs_BprimeMLP = mlp.predict_proba(testBprime)
probs_Bprime2MLP = mlp.predict_proba(testBprime2)
# probs = [probsWJets, probsTTbarT, probsBprime]
probsMLP = [probs_WJetsMLP, probs_TTbarTMLP, probs_Bprime2MLP]
   

# Get scores for non-training events on DT
probs_WJetsDT = dtModel.predict_proba(testWJets)
probs_TTbarTDT = dtModel.predict_proba(testTTbarT)
probs_BprimeDT = dtModel.predict_proba(testBprime)
probs_Bprime2DT = dtModel.predict_proba(testBprime2)
# probs = [probsWJets, probsTTbarT, probsBprime]
probsDT = [probs_WJetsDT, probs_TTbarTDT, probs_Bprime2DT]


# Get scores for non-training events on CNN
probs_WJetsCNN = cnn.predict(testWJets.reshape(testWJets.shape[0], testWJets.shape[1], input_dimension))
probs_TTbarTCNN = cnn.predict(testTTbarT.reshape(testTTbarT.shape[0], testTTbarT.shape[1], input_dimension))
probs_BprimeCNN = cnn.predict(testBprime.reshape(testBprime.shape[0], testBprime.shape[1], input_dimension))
probs_Bprime2CNN = cnn.predict(testBprime2.reshape(testBprime2.shape[0], testBprime2.shape[1], input_dimension))
# probs = [probsWJets, probsTTbarT, probsBprime]
probsCNN = [probs_WJetsCNN, probs_TTbarTCNN, probs_Bprime2CNN]

# %%
### Plotting the comparison 
## WJets
plt.close()
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

plt.close()
plt.figure()
plt.xlabel('Predicted W boson score - DT',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_WJetsDT.T[0], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_TTbarTDT.T[0], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_BprimeDT.T[0], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.hist(probs_Bprime2DT.T[0], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime2)+'\,TeV)}$', color='c', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig(outdir+'plots/score_WJetDT'+outStr+'.png')

plt.close()
plt.figure()
plt.xlabel('Predicted W boson score - CNN',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_WJetsCNN.T[0], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_TTbarTCNN.T[0], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_BprimeCNN.T[0], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.hist(probs_Bprime2CNN.T[0], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime2)+'\,TeV)}$', color='c', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig(outdir+'plots/score_WJetCNN'+outStr+'.png')


## TTbarT
plt.close()
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

plt.close()
plt.figure()
plt.xlabel('Predicted top quark score - DT',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_WJetsDT.T[1], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_TTbarTDT.T[1], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_BprimeDT.T[1], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.hist(probs_Bprime2DT.T[1], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime2)+'\,TeV)}$', color='c', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig(outdir+'plots/score_TTbarTDT'+outStr+'.png')

plt.close()
plt.figure()
plt.xlabel('Predicted top quark score - CNN',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_WJetsCNN.T[1], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_TTbarTCNN.T[1], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_BprimeCNN.T[1], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.hist(probs_Bprime2CNN.T[1], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime2)+'\,TeV)}$', color='c', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig(outdir+'plots/score_TTbarTCNN'+outStr+'.png')


## Signal
plt.close()
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

plt.close()
plt.figure()
plt.xlabel('Predicted B quark score - DT',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_WJetsDT.T[2], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_TTbarTDT.T[2], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_BprimeDT.T[2], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.hist(probs_Bprime2DT.T[2], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime2)+'\,TeV)}$', color='c', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig(outdir+'plots/score_BprimeDT'+outStr+'.png')

plt.close()
plt.figure()
plt.xlabel('Predicted B quark score - CNN',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_WJetsCNN.T[2], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_TTbarTCNN.T[2], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_BprimeCNN.T[2], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.hist(probs_Bprime2CNN.T[2], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime2)+'\,TeV)}$', color='c', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig(outdir+'plots/score_BprimeCNN'+outStr+'.png')
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

accCNN = (confusion[0][0] + confusion[1][1] + confusion[2][2]) / np.sum(confusion)
precCNN = ((confusion[0][0] / (confusion[1][0] + confusion[2][0] + confusion[0][0])) + (confusion[1][1]/(confusion[0][1] + confusion[2][1] + confusion[1][1])) + (confusion[2][2]/(confusion[0][2] + confusion[1][2] + confusion[2][2]))) / 3
recallCNN = ((confusion[0][0] / (confusion[0][0] + confusion[0][1] + confusion[0][2])) + (confusion[1][1]/(confusion[1][0] + confusion[1][1] + confusion[1][2])) + (confusion[2][2]/(confusion[2][0] + confusion[2][1] + confusion[2][2]))) / 3
fscoreCNN = (precCNN + recallCNN) / 2

print('------MLP------')
print('Precision: ' + str(precMLP))
print('Recall: ' + str(recallMLP))
print('F-measure: ' + str(fscoreMLP))
print('Trained in ' + str(mlpTime) + ' s')

print('------DT------')
print('Precision: ' + str(precDT))
print('Recall: ' + str(recallDT))
print('F-measure: ' + str(fscoreDT))
print('Trained in ' + str(dtTime) + ' s')

print('------CNN------')
print('Precision: ' + str(precCNN))
print('Recall: ' + str(recallCNN))
print('F-measure: ' + str(fscoreCNN))
print('Trained in ' + str(cnnTime) + ' s')
# %%
## Saving models to files
def make_keras_picklable():
    def __getstate__(self):
        model_str = ""
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=True) as fd:
            save_model(self, fd.name, overwrite=True)
            model_str = fd.read()
        d = {'model_str': model_str}
        return d

    def __setstate__(self, state):
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=True) as fd:
            fd.write(state['model_str'])
            fd.flush()
            model = load_model(fd.name)
        self.__dict__ = model.__dict__

    cls = Model
    cls.__getstate__ = __getstate__
    cls.__setstate__ = __setstate__

# Run the function
make_keras_picklable()

pickle.dump(mlp, open(outdir+'models/Dnn_mlp_3bin'+outStr+'.pkl', 'wb'))
pickle.dump(dtModel, open(outdir+'models/dt_3bin' + outStr +'.pkl', 'wb'))
pickle.dump(cnn, open(outdir+'models/dt_3bin' + outStr +'.pkl', 'wb'))
pickle.dump(scaler, open(outdir+'Dnn_scaler_3bin'+outStr+'.pkl', 'wb'))
# %%
cnn.save(outdir+'models/cnn', 'wb')