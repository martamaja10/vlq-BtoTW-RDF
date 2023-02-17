#!/usr/bin/env python
# coding: utf-8

# # Dense neural network with Keras
# Author: Javier Duarte
# Adapted by: Sam Carlson

# ## Loading `pandas` DataFrames
# Now we load two different `NumPy` arrays. One corresponding to the VV signal and one corresponding to the background.

# In[1]:


import uproot
import numpy as np
import pandas as pd
import h5py

# fix random seed for reproducibility
seed = 7
np.random.seed(seed)

treename = 'Events'
filename = {}
upfile = {}
params = {}
df = {}

filename['BpM2000'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM2000_hadd.root'
filename['BpM1400'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM1400_hadd.root'
filename['BpM800'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM800_hadd.root'
filename['WJets1200'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets1200_hadd.root'
filename['WJets2500'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets2500_hadd.root'
filename['WJets800'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets800_hadd.root'
filename['WJets600'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets600_hadd.root'
filename['WJets400'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets400_hadd.root'
filename['WJets200'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/WJets200_hadd.root'
filename['ttbarInc'] = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/ttbarInc_hadd.root'

VARS = ['pNet_J_1',#'pNet_J_2',
        'pNet_T_1',#'pNet_T_2',
        'pNet_W_1',#'pNet_W_2',
        'FatJet_pt_1',#'FatJet_pt_2',
        'FatJet_sdMass_1',#'FatJet_sdMass_2',
        'tau21_1',#'tau21_2',
        'nJ_pNet','nT_pNet','nW_pNet',
        'Jet_HT','Jet_ST','MET_pt',
        't_pt','t_mass',
        #'t_dRWb', # t_dRWb does not exist, should check RDF script
        'NJets_central', 'NJets_DeepFlavM','NFatJets','NJets_forward',
        'Bprime_DR','Bprime_ptbal','Bprime_chi2'] # choose which vars to use (2d)

for key in filename:

    upfile[key] = uproot.open(filename[key])
    params[key] = upfile[key][treename].arrays(VARS)

    df[key] = pd.DataFrame(params[key],columns=VARS)


## cut out undefined variables VARS[0] and VARS[1] > -999
#df[key]= df[key][(df[key][VARS[0]] > -999) & (df[key][VARS[1]] > -999)]

# add isSignal variable
    df[key]['isSignal'] = np.ones(len(df[key])) 



# ## Define the model
# We'll start with a dense (fully-connected) NN layer.
# Our model will have a single fully-connected hidden layer with the same number of neurons as input variables. 
# The weights are initialized using a small Gaussian random number. 
# We will switch between linear and tanh activation functions for the hidden layer.
# The output layer contains a single neuron in order to make predictions. 
# It uses the sigmoid activation function in order to produce a probability output in the range of 0 to 1.
# 
# We are using the `binary_crossentropy` loss function during training, a standard loss function for binary classification problems. 
# We will optimize the model with the Adam algorithm for stochastic gradient descent and we will collect accuracy metrics while the model is trained.

# In[2]:


# baseline keras model
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten
from keras.utils import np_utils

NDIM = len(VARS)
inputs = Input(shape=(NDIM,), name = 'input')  
outputs = Dense(1, name = 'output', kernel_initializer='normal', activation='sigmoid')(inputs)

# creae the model
model = Model(inputs=inputs, outputs=outputs)
# compile the model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
# print the model summary
model.summary()


# ## Dividing the data into testing and training dataset
# 
# We will split the data into two parts (one for training+validation and one for testing). 
# We will also apply "standard scaling" preprocessing: http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html i.e. making the mean = 0 and the RMS = 1 for all input variables (based **only** on the training/validation dataset).
# We will also define our early stopping criteria to prevent over-fitting and we will save the model based on the best `val_loss`.

# In[3]:


df_all = pd.concat([df['VV'],df['bkg']])
dataset = df_all.values
X = dataset[:,0:NDIM]
Y = dataset[:,NDIM]

from sklearn.model_selection import train_test_split
X_train_val, X_test, Y_train_val, Y_test = train_test_split(X, Y, test_size=0.2, random_state=7)

# preprocessing: standard scalar
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler().fit(X_train_val)
X_train_val = scaler.transform(X_train_val)
X_test = scaler.transform(X_test)

# early stopping callback
from keras.callbacks import EarlyStopping
early_stopping = EarlyStopping(monitor='val_loss', patience=10)

# model checkpoint callback
# this saves our model architecture + parameters into dense_model.h5
from keras.callbacks import ModelCheckpoint
model_checkpoint = ModelCheckpoint('dense_model.h5', monitor='val_loss', 
                                   verbose=0, save_best_only=True, 
                                   save_weights_only=False, mode='auto', 
                                   period=1)


# ## Run training 
# Here, we run the training.

# In[4]:


# Train classifier
history = model.fit(X_train_val, 
                    Y_train_val, 
                    epochs=1000, 
                    batch_size=1024, 
                    verbose=0, # switch to 1 for more verbosity 
                    callbacks=[early_stopping, model_checkpoint], 
                    validation_split=0.25)


# ## Plot performance
# Here, we plot the history of the training and the performance in a ROC curve

# In[5]:


import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
# plot loss vs epoch
plt.figure(figsize=(15,10))
ax = plt.subplot(2, 2, 1)
ax.plot(history.history['loss'], label='loss')
ax.plot(history.history['val_loss'], label='val_loss')
ax.legend(loc="upper right")
ax.set_xlabel('epoch')
ax.set_ylabel('loss')

# plot accuracy vs epoch
ax = plt.subplot(2, 2, 2)
ax.plot(history.history['acc'], label='acc')
ax.plot(history.history['val_acc'], label='val_acc')
ax.legend(loc="upper left")
ax.set_xlabel('epoch')
ax.set_ylabel('acc')

# Plot ROC
Y_predict = model.predict(X_test)
from sklearn.metrics import roc_curve, auc
fpr, tpr, thresholds = roc_curve(Y_test, Y_predict)
roc_auc = auc(fpr, tpr)
ax = plt.subplot(2, 2, 3)
ax.plot(fpr, tpr, lw=2, color='cyan', label='auc = %.3f' % (roc_auc))
ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='k', label='random chance')
ax.set_xlim([0, 1.0])
ax.set_ylim([0, 1.0])
ax.set_xlabel('false positive rate')
ax.set_ylabel('true positive rate')
ax.set_title('receiver operating curve')
ax.legend(loc="lower right")
plt.show()


# In[6]:


df_all['dense'] = model.predict(X) # add prediction to array
print(df_all.iloc[:5])


# # Plot NN output vs input variables
# Here, we will plot the NN output and devision boundary as a function of the input variables.
# 
# **Question 1:** How can we fill the correct numpy arrays for plotting?

# In[7]:


# Hint: We want to make a three 2D numpy arrays: 
# x values at each (x, y) grid point
# y values at each (x, y) grid point
# z values (model prediction) at each (x, y) grid point

myXI, myYI = np.meshgrid(np.linspace(-2, 2, 200), np.linspace(-2, 2, 200))
# print shape
print(myXI.shape)

for i in range(0, len(myXI)):
    for j in range(0, len(myYI)):
        myXI[i,j] # x value of xi, yj point
        myYI[i,j] # y value of xi, yj point
        #myZI[i,j] = ??? # change this

myZI = model.predict(np.c_[myXI.ravel(), myYI.ravel()])
myZI = myZI.reshape(myXI.shape)


# **Question 2:** The code below shoes how to plot the NN output. How can we plot the NN decision boundary?

# In[8]:


from matplotlib.colors import ListedColormap
plt.figure(figsize=(20,7))

# plot contour map of NN output
# overlaid with test data points
ax = plt.subplot(1, 2, 1)
cm = plt.cm.RdBu
cm_bright = ListedColormap(['#FF0000', '#0000FF'])
cont_plot = ax.contourf(myXI, myYI, myZI, cmap=cm, alpha=.8)
ax.scatter(X_test[:, 0], X_test[:, 1], c=Y_test, cmap=cm_bright, edgecolors='k')
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)
ax.set_xlabel(VARS[0])
ax.set_ylabel(VARS[1])
plt.colorbar(cont_plot,ax=ax, boundaries=[0,1],label='NN output')

# plot decision boundary
# overlaid with test data points
ax = plt.subplot(1, 2, 2)
cm = plt.cm.RdBu
cm_bright = ListedColormap(['#FF0000', '#0000FF'])
cont_plot = ax.contourf(myXI, myYI, myZI>0.5, cmap=cm, alpha=.8)
ax.scatter(X_test[:, 0], X_test[:, 1], c=Y_test, cmap=cm_bright, edgecolors='k')
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)
ax.set_xlabel(VARS[0])
ax.set_ylabel(VARS[1])
plt.colorbar(cont_plot,ax=ax, boundaries=[0,1],label='NN output')


# **Question 3:** What happens if you increase/decrease the number of hidden layers?
# 
# **Question 4:** What happens if you increase/decrease the number of nodes per hidden layer?
# 
# **Question 5:** What happens if you add/remove dropout?
# 
# **Question 6:** What happens if you add/remove early stopping?

# ## Add prediction to `ROOT` trees
# Here we'll add the precition we've computed to `ROOT` trees.

# In[9]:


from root_numpy import root2array, array2root


def get_features_from_file(filename='', treename='', branches=[]):
    t = root2array(filename, treename=treename, branches=branches) # structured numpy array 
    #print t.shape 
    t = t.view(np.float32).reshape(t.shape + (-1,)) # normal numpy array (trick from https://stackoverflow.com/questions/5957380/convert-structured-array-to-regular-numpy-array)
    #print t.shape
    return t

def write_prediction_to_file(features, model, filename='',treename='',branch=''):
    y_predict_all = model.predict(features) # normal numpy array
    #print y_predict_all.shape
    y_predict_all = np.array(y_predict_all, dtype=[(branch, np.float32)]) # structured numpy array
    #print y_predict_all.shape
    array2root(y_predict_all, filename, treename=treename, mode='recreate')
    
X_all = get_features_from_file('data/ntuple_4mu_VV.root', 
                               treename='HZZ4LeptonsAnalysisReduced', 
                               branches=VARS)

X_all = scaler.transform(X_all)

write_prediction_to_file(X_all, 
                         model, 
                         filename='data/ntuple_4mu_VV_predict.root', 
                         treename='HZZ4LeptonsAnalysisReduced', 
                         branch='dense')

X_all = get_features_from_file('data/ntuple_4mu_bkg.root', 
                               treename='HZZ4LeptonsAnalysisReduced', 
                               branches=VARS)

X_all = scaler.transform(X_all)

write_prediction_to_file(X_all, 
                         model, 
                         filename='data/ntuple_4mu_bkg_predict.root', 
                         treename='HZZ4LeptonsAnalysisReduced', 
                         branch='dense')


# In[ ]:




