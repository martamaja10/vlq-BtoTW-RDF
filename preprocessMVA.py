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
import subprocess

# %%
### Reading in basic parameters
start_time = time.time() # collected just for benchmarking
outdirName = sys.argv[1] # user can define this at runtime
arch = '3x10'
testnum = 0 # currently hard-coded to 0 until it is needed
maxtest = 300000

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
Bprime = 2.0
Bprime2 = 0.8
test2000 = True

# Defining plotting parameters
WithBprimeVars = False
outStr = '_2018BB_'+str(arch)+'_' + str(millify(maxtest)) +'test'

# %%
### Defining variables to be used with model (defined in RDataframe script)
vars = ['pNet_J_1',#'pNet_J_2',
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
        'Bprime_DR','Bprime_ptbal','Bprime_chi2']
        #'minDR_leadAK8otherAK8'] 

# %%
### Getting data from ROOT files
print('Opening files...')
eosdir = "root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/"

# Defining selection criteria for the events
seltrain = "Bprime_mass > 0 && (NJets_forward == 0 || (NJets_forward > 0 && FatJet_phi[0] <= 0))"
seltest = "Bprime_mass > 0 && NJets_forward > 0 && FatJet_phi[0] > 0"

treeVars = vars

# Potential variable for future loop controller
filenames = ['ttbarT', 'ttbarTb', #'singleT', 'singleTb', 
             'WJets2500','WJets1200', 'WJets800', 'WJets600', 'WJets400', 'WJets200']

weights = [1, 1, #0.456, 0.0506, 
           0.0011, 0.0148, 0.0544, 0.1128, 0.4749, 0.4466]

## Performing selections and reads

# Selections on ttbarT(b) and singlet(b)

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

# Selection with signals
fName = fname
print()
sys.stdout.write('\rNow processing signal file...                                                                ')
sys.stdout.flush()
weight = 1
fileBp1  = TFile.Open(eosdir + "BpM14000_hadd.root", "READ")
fileBp2 = TFile.Open(eosdir + "BpM2000_hadd.root", "READ")
# Swapping order depending on defined parameters
if Bprime == 2.0:
    temp = fileBp1
    fileBp1 = fileBp2
    fileBp2 = temp
treeBprime = fileBp1.Get("Events")
treeBprime2 = fileBp2.Get("Events")
trainBprime= addWeight(tree2array(treeBprime, treeVars, seltrain), weight)
testBprime= addWeight(tree2array(treeBprime, treeVars, seltest), weight)
testBprime2= addWeight(tree2array(treeBprime2, treeVars, seltest), weight)
waitDone = True
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
# %%
### Update the user
print('------------ Before Cuts -------------')
print('Training Events:')
print('Number of ttBarT: ' + str(len(trainTTbarT)))
print('Number of Bprime: ' + str(len(trainBprime)))
print('Number of WJets: ' + str(len(trainWJets)))
#print('Number of singleT: ' + str(len(trainSingleT)) + '\n')
print('Testing Events:')
print('Number of ttBarT: ' + str(len(testTTbarT)))
print('Number of Bprime: ' + str(len(testBprime)))
print('Number of Bprime2: ' + str(len(testBprime2)))
print('Number of WJets: ' + str(len(testWJets)))
#print('Number of singleT: ' + str(len(testSingleT)) + '\n')

# %%
### Post-processing to prepare date for plotting and export

## Shorten the testing arrays to the chosen length
testBprime = testBprime[:maxtest]
testBprime2 = testBprime2[:maxtest]
testWJets = testWJets[:maxtest]
testTTbarT = testTTbarT[:maxtest]
#testSingleT = testSingleT[:maxtest]

## Calculate the maximum number of allowed training events in each group
## We'll allow up to 30% imbalance between the samples. 
maxpersample = int(round(1.3*min(len(trainTTbarT), len(trainBprime), len(trainWJets)),0))#, len(trainSingleT)),0))

## Shorten the training arrays to the max allowed length
## These are not shuffled yet, so we will chop off fewer "good" testing events
trainBprime = trainBprime[:maxpersample]
trainWJets = trainWJets[:maxpersample]
trainTTbarT = trainTTbarT[:maxpersample]
#trainSingleT = trainSingleT[:maxpersample]

# Final shuffle to eliminate any unwanted patterns
np.random.shuffle(trainTTbarT)
np.random.shuffle(trainBprime)
np.random.shuffle(trainWJets)
#np.random.shuffle(trainSingleT)
np.random.shuffle(testTTbarT)
np.random.shuffle(testBprime)
np.random.shuffle(testWJets)
#np.random.shuffle(testSingleT)

# Print final size information to user
print('------------ Final Sizes ------------')
print('Training Events:')
print('Number of ttBarT: ' + str(len(trainTTbarT)))
print('Number of Bprime: ' + str(len(trainBprime)))
print('Number of WJets: ' + str(len(trainWJets)))
#print('Number of singleT: ' + str(len(trainSingleT)) + '\n')
print('Testing Events:')
print('Number of ttBarT: ' + str(len(testTTbarT)))
print('Number of Bprime: ' + str(len(testBprime)))
print('Number of Bprime2: ' + str(len(testBprime2)))
print('Number of WJets: ' + str(len(testWJets)))
#print('Number of singleT: ' + str(len(testSingleT)) + '\n')

logfile.write(str(len(trainTTbarT)) + ", " + str(len(trainBprime)) + ", " +str(len(trainWJets)) + ", " +str(len(testTTbarT)) + ", " +str(len(testBprime)) + ", " +str(len(testBprime2)) + ", " +str(len(testWJets)))# + ", " + str(len(testSingleT) + ", " + str(len(trainSingleT)) ))
logfile.close()

# %%
### Final steps before merge and plotting

## Resample the training data so that events with larger weights are more likely to be included
## We'll do this within each type of sample so that we don't break the equivalency
weightsTrainTTbarT = np.array([sub[0] for sub in trainTTbarT])
weightsTrainBprime = np.array([sub[0] for sub in trainBprime])
weightsTrainWJets = np.array([sub[0] for sub in trainWJets])
#weightsTrainSingleT = np.array([sub[0] for sub in trainSingleT])
weightsTestTTbarT = np.array([sub[0] for sub in testTTbarT])
weightsTestBprime = np.array([sub[0] for sub in testBprime])
weightsTestBprime2 = np.array([sub[0] for sub in testBprime2])
weightsTestWJets = np.array([sub[0] for sub in testWJets])
#weightsTestSingleT = np.array([sub[0] for sub in testSingleT])
## Cut off weights from front of array
trainTTbarT = [sub[1:] for sub in trainTTbarT]
trainBprime = [sub[1:] for sub in trainBprime]
trainWJets = [sub[1:] for sub in trainWJets]
#trainSingleT = [sub[1:] for sub in trainSingleT]
testTTbarT = [sub[1:] for sub in testTTbarT]
testBprime = [sub[1:] for sub in testBprime]
testBprime2 = [sub[1:] for sub in testBprime2]
testWJets = [sub[1:] for sub in testWJets]
#testSingleT = [sub[1:] for sub in testSingleT]

RStrainTTbarT = (resample_with_replacement(trainTTbarT, weightsTrainTTbarT)).tolist()
RStrainBprime = (resample_with_replacement(trainBprime, weightsTrainBprime)).tolist()
RStrainWJets = (resample_with_replacement(trainWJets, weightsTrainWJets)).tolist()
#RStrainSingleT = (resample_with_replacement(trainSingleT, weightsTrainSingleT)).tolist()
RStestTTbarT = (resample_with_replacement(testTTbarT, weightsTestTTbarT)).tolist()
RStestBprime = (resample_with_replacement(testBprime, weightsTestBprime)).tolist()
RStestBprime2 = (resample_with_replacement(testBprime2, weightsTestBprime2)).tolist()
RStestWJets = (resample_with_replacement(testWJets, weightsTestWJets)).tolist()
#RStestSingleT = (resample_with_replacement(testSingleT, weightsTestSingleT)).tolist()

## New versions are used for merging and copies are used for unaltered plotting
trainTTbarT = copy.copy(RStrainTTbarT)
trainBprime = copy.copy(RStrainBprime)
trainWJets = copy.copy(RStrainWJets)
#trainSingleT = copy.copy(RStrainSingleT)
testTTbarT = copy.copy(RStestTTbarT)
testBprime = copy.copy(RStestBprime)
testBprime2 = copy.copy(RStestBprime2)
testWJets = copy.copy(RStestWJets)
#testSingleT = copy.copy(RStestSingleT)

## Transpose these arrays to get arrays for plotting
## Each entry is one variable for all the events
## We will make sure all samples are the same size for plots
#numPerSample = min(len(trainTTbarT),len(trainBprime), len(trainWJets))#, len(trainSingleT))
#print('Setting numPerSample = '+str(numPerSample))

histsTTbarT = np.array(trainTTbarT).T#[:numPerSample]).T
histsBprime = np.array(trainBprime).T#[:numPerSample]).T
#histsTprime2 = np.array(trainTprime2[:numPerSample]).T
histsWJets = np.array(trainWJets).T#[:numPerSample]).T
#histsSingleT = np.array(trainSingleT[:numPerSample]).T


# Remove invalid rows
histsTTbarT = histsTTbarT[~np.isnan(histsTTbarT).any(axis=1)]
histsWJets = histsWJets[~np.isnan(histsWJets).any(axis=1)]
histsBprime= histsBprime[~np.isnan(histsBprime).any(axis=1)]
# %%
### Plotting input variables
print('Plotting input variables...')
for index, hist in enumerate(histsWJets):
   sys.stdout.write('\rNow processing plot {}/'.format(index + 1) + str(len(histsWJets)) + ' - ' + vars[index] + '...       ')
   sys.stdout.flush()
   plt.figure()
   plt.hist(hist, bins=50, color='g', label=r'$\mathrm{W+jets}$', histtype='step', normed=True)
   plt.hist(histsBprime[index], bins=50, color='y', label=r'$\mathrm{T\overline{T}\,('+str(Bprime)+'\,TeV)}$', histtype='step', normed=True)
   #plt.hist(histsTprime2[index], bins=50, color='c', label=r'$\mathrm{T\overline{T}\,('+str(Tprime2)+'\,TeV)}$', histtype='step', normed=True)
   plt.hist(histsTTbarT[index], bins=50, color='r', label=r'$\mathrm{t\bar{t}}$', histtype='step', normed=True)
   # plt.hist(histsSingleT[index], bins=50, color='k', label=r'$\mathrm{singleT}$', histtype='step', normed=True)
   plt.title('CMS Simulation',loc='left',size=18)
   plt.title('Work in progress',loc='right',size=14,style='italic')
   plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
   plt.xlabel(vars[index],horizontalalignment='right',x=1.0,size=14)
   plt.legend(loc='best',fontsize=14)
   if not WithBprimeVars: plt.savefig(outdirName+'plots_'+str(vars[index])+outStr)
   if WithBprimeVars: plt.savefig(outdirName+'plots_'+str(vars[index])+outStr)
   plt.close()
   
   # Logarithmic option
   plt.figure()
   plt.hist(hist, bins=50, color='g', label=r'$\mathrm{W+jets}$', histtype='step', normed=True)
   plt.hist(histsBprime[index], bins=50, color='y', label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', histtype='step', normed=True)
   #plt.hist(histsTprime2[index], bins=50, color='c', label=r'$\mathrm{T\overline{T}\,('+str(Tprime2)+'\,TeV)}$', histtype='step', normed=True)
   plt.hist(histsTTbarT[index], bins=50, color='r', label=r'$\mathrm{t\bar{t}}$', histtype='step', normed=True)
   # plt.hist(histsSingleT[index], bins=50, color='k', label=r'$\mathrm{singleT}$', histtype='step', normed=True)
   plt.title('CMS Simulation',loc='left',size=18)
   plt.title('Work in progress',loc='right',size=14,style='italic')
   plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
   plt.yscale('log', nonposy='clip')
   plt.xlabel(vars[index],horizontalalignment='right',x=1.0,size=14)
   plt.legend(loc='best',fontsize=14)
   if not WithBprimeVars: plt.savefig(outdirName+'logscale_plots_'+str(vars[index])+outStr)
   if WithBprimeVars: plt.savefig(outdirName+'logscale_plots_'+str(vars[index])+outStr)
   plt.close()

# %%
### Make arrays of training and testing data
trainData = []
trainLabel = []
testWeights = []
nEvents = len(RStrainTTbarT) + len(RStrainBprime) + len(RStrainWJets) #+ len(RStrainSingleT)
print('nEvents in the training set '+str(nEvents))

# Pull random data based on a selected random integer
while nEvents > 0:
    rng = random.randint(0,2)
    if(rng == 0 and len(RStrainWJets) > 0):
        trainData.append(RStrainWJets.pop())
        trainLabel.append(rng)
        nEvents -= 1

    elif(rng == 1 and len(RStrainTTbarT) > 0):
        trainData.append(RStrainTTbarT.pop())
        trainLabel.append(rng)
        nEvents -= 1

    elif(rng == 2 and len(RStrainBprime) > 0):
        trainData.append(RStrainBprime.pop())
        trainLabel.append(rng)
        nEvents -= 1

     #elif(rng == 3 and len(RStrainSingleT) > 0):
     #   trainData.append(RStrainSingleT.pop())
     #   trainLabel.append(rng)
     #   nEvents -= 1

    

# If using larger weight class, we use a different variable set
if test2000: nEventsTest = len(RStestTTbarT) + len(RStestBprime2) + len(RStestWJets) #+ len(RStestSingleT)
else: nEventsTest = len(RStestTTbarT) + len(RStestBprime) + len(RStestWJets)# + len(RStestSingleT)
print('nEventsTest = '+str(nEventsTest))

testData = []
testLabel = []
testWeight = [] 
# Once again else clauses kick in when selected list is empty
while(nEventsTest > 0):
   rng = random.randint(0, 2)
   if(rng == 0 and len(RStestWJets) > 0):
      testData.append(RStestWJets.pop())
      testLabel.append(rng)
      nEventsTest -= 1

   elif(rng == 1 and len(RStestTTbarT) > 0):
      testData.append(RStestTTbarT.pop())
      testLabel.append(rng)
      nEventsTest -= 1

   elif(rng == 2):
      if(test2000 and len(RStestBprime2) > 0): 
         testData.append(RStestBprime2.pop())
         testLabel.append(rng)
         nEventsTest -= 1
      elif(not test2000 and len(RStestBprime) > 0):
         testData.append(RStestBprime.pop())
         testLabel.append(rng)
         nEventsTest -= 1

   # elif(rng == 3 and len(RStestSingleT) > 0):
   #    testData.append(RStestSingleT.pop()) 
   #    trainLabel.append(rng)
   #    nEvents -= 1

# Save the output to a file
print('Saving file to ' + outdirName + 'Arrays' + outStr)
np.savez(outdirName + 'Arrays' + outStr, trainData=trainData, trainLabel = trainLabel, testData = testData, 
    testLabel = testLabel, testWJets = testWJets, testTTbarT = testTTbarT, testBprime = testBprime,
    testBprime2 = testBprime2)#, testSingleT = testSingleT)

print('Done')
print('Time Taken: %s minutes' % (round(time.time() - start_time, 2)/60)) 
