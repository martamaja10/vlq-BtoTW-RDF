# baseline keras model
import numpy as np
import pandas as pd
import h5py
import pickle
import os
from keras.models import Sequential, Model
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten
from keras.utils import np_utils
import tensorflow as tf
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

# %%
### Loading data
npzPath = 'NewAnalysisArrays.npz'
outPath = 'NewAnalysisModels'
if not os.path.exists(outPath): 
   os.system('mkdir '+ outPath)
outPath = './' + outPath + '/'

print('Importing data from {}...'.format(npzPath))
archive = np.load(npzPath)

# Three datasets
background = archive['background']
Bp800 = archive['lowMass']
Bp2000 = archive['largeMass']

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
dense2 = Dense(25, activation = 'relu')(dense1)
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
# Background Events
eventData = background[:, :NDIM]
weights = background[:, NDIM]
labels = background[:, NDIM+1]
X_Back_train_val, X_Back_test, Y_Back_train_val, Y_Back_test, weights_Back_train, weights_Back_test = train_test_split(eventData, labels, weights, test_size=0.9, random_state=7)

eventData = Bp800[:, :NDIM]
weights = Bp800[:, NDIM]
labels = Bp800[:, NDIM+1]
X_LM_train_val, X_LM_test, Y_LM_train_val, Y_LM_test, weights_LM_train, weights_LM_test = train_test_split(eventData, labels, weights, test_size=0.1, random_state=7)

eventData = Bp2000[:, :NDIM]
weights = Bp2000[:, NDIM]
labels = Bp2000[:, NDIM + 1]
X_HM_train_val, X_HM_test, Y_HM_train_val, Y_HM_test, weights_HM_train, weights_HM_test = train_test_split(eventData, labels, weights, test_size=0.9, random_state=7)


X_train_val = np.concatenate((X_Back_train_val[:9000], X_LM_train_val, X_HM_train_val[:1000])) # Take low mass BP and background for training
Y_train_val = np.concatenate((Y_Back_train_val[:9000], Y_LM_train_val, Y_HM_train_val[:1000])) 
weights_train = np.concatenate((weights_Back_train[:9000], weights_LM_train, weights_HM_train[:1000]))

# Randomly shuffling training set
indices = np.random.permutation(len(X_train_val))
X_train_val = X_train_val[indices]
Y_train_val = Y_train_val[indices]
weights_train = weights_train[indices]

# Building training set out of all events left
X_test = np.concatenate((X_Back_test, X_LM_test, Bp2000[:, :NDIM]))
Y_test = np.concatenate((Y_Back_test, Y_LM_test, Bp2000[:, NDIM+1]))
weights_test = np.concatenate((weights_Back_test, weights_LM_test, Bp2000[:, NDIM]))

indices = np.random.permutation(len(X_test))
X_test = X_test[indices]
Y_test = Y_test[indices]
weights_test = weights_test[indices]

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
with tf.device('/CPU:0'):
        history = model.fit(X_train_val, 
                    tf.keras.utils.to_categorical(Y_train_val), 
                    epochs=1000, 
                    batch_size=1024, 
                    verbose=0, # switch to 1 for more verbosity 
                    callbacks=[early_stopping, model_checkpoint], 
                    validation_split=0.1,
                    sample_weight = weights_train
                    )


pred_mlp = model.predict(X_test).argmax(axis = -1)

from sklearn.metrics import confusion_matrix

confusion_matrix(Y_test, pred_mlp)

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
plt.hist(probs_bprime.T[0], bins=20, range=(0,1), label=r'$\mathrm{Bprime}$', color='m', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig('plots/score_WJet.png')

plt.figure()
plt.xlabel('Predicted top quark score',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_wjets.T[1], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_ttbar.T[1], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_bprime.T[1], bins=20, range=(0,1), label=r'$\mathrm{Bprime}$', color='m', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig('plots/score_TTbarT.png')


plt.figure()
plt.xlabel('Predicted Bprime score',horizontalalignment='right',x=1.0,size=14)
plt.ylabel('Events per bin',horizontalalignment='right',y=1.0,size=14)
plt.title('CMS Simulation',loc='left',size=18)
plt.title('Work in progress',loc='right',size=14,style='italic')
plt.ylim([0.01,10.**4])
plt.hist(probs_wjets.T[2], bins=20, range=(0,1), label=r'$\mathrm{W+jets}$', color='g', histtype='step', log=True, density=True)
plt.hist(probs_ttbar.T[2], bins=20, range=(0,1), label=r'$\mathrm{t\bar{t}}$', color='y', histtype='step', log=True, density=True)
plt.hist(probs_bprime.T[2], bins=20, range=(0,1), label=r'$\mathrm{Bprime}$', color='m', histtype='step', log=True, density=True)
plt.legend(loc='best')
plt.savefig('plots/score_bprime.png')
