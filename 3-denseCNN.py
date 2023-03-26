# baseline keras model
import numpy as np
import pandas as pd
import h5py
import pickle
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling1D, Dropout, Flatten, Conv1D
from keras.utils import np_utils
import tensorflow as tf
from keras.callbacks import ModelCheckpoint, EarlyStopping

# %%
### Loading data
allmystuff = np.load('NewAnalysisArrays.npz')
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

# %%
### Processing data
print('Preprocessing data...')
eventData = allmystuff[:, :NDIM]
weights = allmystuff[:, NDIM]
labels = allmystuff[:, NDIM+1]

from sklearn.model_selection import train_test_split
X_train_val, X_test, Y_train_val, Y_test, weights_train, weights_test = train_test_split(eventData, labels, weights, test_size=0.2, random_state=7)

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler().fit(X_train_val)
X_train_val = scaler.transform(X_train_val)
X_test = scaler.transform(X_test)

#%%
### Creating the model
sample_size = X_train_val.shape[0]
time_steps = X_train_val.shape[1]
input_dimension = 1
callbackList = [EarlyStopping(monitor = 'val_loss', mode = 'min', patience = 25, verbose = 1),
               ModelCheckpoint('./ModelChkpt/checkpoint', save_weights_only=True, monitor='val_accuracy', mode = 'max', save_best_only = True)]
trainCNN = X_train_val.reshape(sample_size, time_steps, input_dimension)
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

cnn.fit(trainCNN, Y_train_val, epochs = 500, validation_split = 0.1, verbose = 0,
   callbacks = callbackList)
cnn.load_weights('./ModelChkpt/checkpoint')
