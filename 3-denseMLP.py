# baseline keras model
import numpy as np
import pandas as pd
import h5py
import pickle
<<<<<<< HEAD
import os
=======
>>>>>>> 92718ba1ccb8f729de7be5419a8c2c8fb4ca3680
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten
from keras.utils import np_utils

# %%
### Loading data
<<<<<<< HEAD
npzPath = 'NewAnalysisArrays.npz'
outPath = 'NewAnalysisModels'
if not os.path.exists(outPath): 
   os.system('mkdir '+ outPath)
outPath = './' + outPath + '/'

print('Importing data from {}...'.format(npzPath))
archive = np.load(npzPath)
allmystuff = archive['arr_0']
=======
allmystuff = np.load('NewAnalysisArrays.npz')
>>>>>>> 92718ba1ccb8f729de7be5419a8c2c8fb4ca3680
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
outputs = Dense(1, name = 'output', kernel_initializer='normal', activation='sigmoid')(inputs)

#%%
<<<<<<< HEAD
### Creating the modeli
print('Generating model...')
=======
### Creating the model
>>>>>>> 92718ba1ccb8f729de7be5419a8c2c8fb4ca3680
model = Model(inputs=inputs, outputs=outputs)
# compile the model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
# print the model summary
model.summary()

# %%
<<<<<<< HEAD
### Processing datai
print('Scaling and separating data...')
=======
### Processing data
>>>>>>> 92718ba1ccb8f729de7be5419a8c2c8fb4ca3680
eventData = allmystuff[:, :NDIM]
weights = allmystuff[:, NDIM]
labels = allmystuff[:, NDIM+1]

from sklearn.model_selection import train_test_split
X_train_val, X_test, Y_train_val, Y_test, weights_train, weights_test = train_test_split(eventData, labels, weights, test_size=0.2, random_state=7)

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler().fit(X_train_val)
<<<<<<< HEAD
pickle.dump(scaler, open(outPath + 'dnn_scaler.pkl', 'wb'))
=======
>>>>>>> 92718ba1ccb8f729de7be5419a8c2c8fb4ca3680
X_train_val = scaler.transform(X_train_val)
X_test = scaler.transform(X_test)

#%%
### Preparing the model
# early stopping callback
from keras.callbacks import EarlyStopping
early_stopping = EarlyStopping(monitor='val_loss', patience=10)

from keras.callbacks import ModelCheckpoint
<<<<<<< HEAD
model_checkpoint = ModelCheckpoint(outPath + 'MLP.h5', monitor='val_loss', 
=======
model_checkpoint = ModelCheckpoint('dense_model.h5', monitor='val_loss', 
>>>>>>> 92718ba1ccb8f729de7be5419a8c2c8fb4ca3680
                                   verbose=0, save_best_only=True, 
                                   save_weights_only=False, mode='auto', 
                                   period=1)

<<<<<<< HEAD
# %%
### Training the model
# Train classifieri
print('Training Model...')
=======

# %%
### Training the model
# Train classifier
>>>>>>> 92718ba1ccb8f729de7be5419a8c2c8fb4ca3680
history = model.fit(X_train_val, 
                    Y_train_val, 
                    epochs=1000, 
                    batch_size=1024, 
                    verbose=0, # switch to 1 for more verbosity 
                    callbacks=[early_stopping, model_checkpoint], 
<<<<<<< HEAD
                    validation_split=0.1,
                    sample_weight = weights_train)
=======
                    validation_split=0.25,
                    sample_weight = weights_train)
>>>>>>> 92718ba1ccb8f729de7be5419a8c2c8fb4ca3680
