# %%
### Imports 

from ROOT import TTree, TH1D, TFile
from root_numpy import tree2array
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow import keras
from keras import backend as K
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.layers import Input, Dense, Concatenate
from tensorflow.keras.models import Model, Sequential, load_model
import importlib
from sklearn.preprocessing import StandardScaler

from sklearn import tree



# %%
### Defining assisting functions

# %%
### Read in arguments

# %%
### Setup logs

# %%
### Signal selection

# %%
### Setup output directory

# %%
### Defining variables for model input

# %%
### Open ROOT files

# %%
### Perform selections on data

# %%
### Further data preprocessing

# %% 
###Plot input data

# %%
### Make arrays of testing data