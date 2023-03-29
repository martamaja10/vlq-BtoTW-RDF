# baseline keras model
import numpy as np
import pandas as pd
import h5py
import pickle
import os
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten
from keras.utils import np_utils
import tensorflow as tf

# %%
### Loading data
npzPath = 'NewAnalysisArrays.npz'
outPath = 'NewAnalysisModels'
if not os.path.exists(outPath): 
   os.system('mkdir '+ outPath)
outPath = './' + outPath + '/'

print('Importing data from {}...'.format(npzPath))
archive = np.load(npzPath)
allmystuff = archive['arr_0']
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

NDIM = len(VARS)
inputs = Input(shape=(NDIM,), name = 'input')
dense1 = Dense(10, activation = 'relu')(inputs)
dense2 = Dense(10, activation = 'relu')(dense1)
dense3 = Dense(10, activation = 'relu')(dense2)   
outputs = Dense(3, name = 'output', kernel_initializer='normal', activation='softmax')(dense3)

#%%
### Creating the modeli
print('Generating model...')
model = Model(inputs = inputs, outputs = outputs)

# compile the model
model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
# print the model summary
model.summary()

# %%
### Processing datai
print('Scaling and separating data...')
eventData = allmystuff[:, :NDIM]
weights = allmystuff[:, NDIM]
labels = allmystuff[:, NDIM+1]

from sklearn.model_selection import train_test_split
X_train_val, X_test, Y_train_val, Y_test, weights_train, weights_test = train_test_split(eventData, labels, weights, test_size=0.2, random_state=7)

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler().fit(X_train_val)
pickle.dump(scaler, open(outPath + 'dnn_scaler.pkl', 'wb'))
X_train_val = scaler.transform(X_train_val)
X_test = scaler.transform(X_test)

#%%
### Preparing the model
# early stopping callback
from keras.callbacks import EarlyStopping
early_stopping = EarlyStopping(monitor='val_loss', patience=10)

from keras.callbacks import ModelCheckpoint
model_checkpoint = ModelCheckpoint(outPath + 'MLP.h5', monitor='val_loss', 
                                   verbose=0, save_best_only=True, 
                                   save_weights_only=False, mode='auto', 
                                   period=1)

# %%
### Training the model
# Train classifieri
print('Training Model...')
history = model.fit(X_train_val, 
                    Y_train_val, 
                    epochs=1000, 
                    batch_size=1024, 
                    verbose=0, # switch to 1 for more verbosity 
                    callbacks=[early_stopping, model_checkpoint], 
                    validation_split=0.1,
                    sample_weight = weights_train)


pred_mlp = cnn.predict(X_test).argmax(axis = -1)

from tf.math import confusion_matrix

confusion_matrix(Y_test, pred_mlp, num_classes = 3, weights = weights_test)

wjets_test = X_test[Y_test == 0]
ttbar_test = X_test[Y_test == 1]
bprime_test = X_test[Y_test == 2]

probs_wjets = model.predict(wjets_test)
probs_ttbar = model.predict(ttbar_test)
probs_bprime = model.predict(bprime_test)

plt.figure()
plt.xlabel('Predicted W boson score',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_wjets.T[0], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_ttbar.T[0], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_bprime.T[0], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig('plots/score_WJetC.png')

plt.figure()
plt.xlabel('Predicted top quark score',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_wjets.T[1], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_ttbar.T[1], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_bprime.T[1], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig(outdir+'plots/score_TTbarT.png')


plt.figure()
plt.xlabel('Predicted top quark score - CNN',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_wjets.T[2], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_ttbar.T[2], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_pprime.T[2], bins=20, range=(0,1), label=r'$\mathrm{Bprime\,('+str(Bprime)+'\,TeV)}$', color='m', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig(outdir+'plots/score_bprime.png')
