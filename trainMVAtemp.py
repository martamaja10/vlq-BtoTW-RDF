# %%
### Imports
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

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

# %%
### Read in arguments
start_time = time.time()
arch = '3x10'
maxtest = 15000
outdir = sys.argv[1]
vararray = int(sys.argv[2])
testnum = int(sys.argv[3])
year = str(sys.argv[4])
if year == 'all': maxtest = 30000

# %%
### Set up logs
outdir = outdir + '/'
if testnum == 1:
    logfile = open(outdir + 'NN_log_'+ year + '.txt', "a+")
    logfile.write('\ntest, vararray, Testing Score (Accuracy), tt-as-BB, BB-as-BB, Precision, Recall, F-Score \n')
else:
    time.sleep(2)
    logfile = open(outdir + 'NN_log_'+ year + '.txt', "a+")
    logfile.write('\n')
logfile.write(str(testnum) + ', ')
logfile.write(str(vararray) + ', ')

# %%
### Signal Selection
Bprime = 0.8
Bprime2 = 2.0
test2000 = True

# %%
### configure output
outStr = '_' + year + 'BB_' + str(arch) + '_' + str(millify(maxtest)) +'test_vars'+str(vararray)+'_Test'+str(testnum)
print('Output location: ' + outdir + outStr)
if not os.path.exists(outdir): os.system('mkdir ' + outdir)

# %%
### Defining input variables

varList = ['pNet_J_1','pNet_J_2',
        'pNet_T_1','pNet_T_2',
        'pNet_W_1','pNet_W_2',
        'dpak8_J_1','dpak8_J_2',
        'dpak8_T_1','dpak8_T_2',
        'dpak8_W_1','dpak8_W_2',
        'FatJet_pt_1','FatJet_pt_2',
        'FatJet_sdMass_1','FatJet_sdMass_2',
        'tau21_1','tau21_2',
        'nJ_dpak8','nT_dpak8','nW_dpak8',
        'nJ_pNet','nT_pNet','nW_pNet',
        'Jet_HT','Jet_ST','MET_pt',
        't_pt','t_mass',
        #'t_dRWb', # t_dRWb does not exist, should check RDF script
        'NJets_central', 'NJets_DeepFlavM','NFatJets','NJets_forward',
        'Bprime_DR','Bprime_ptbal','Bprime_chi2',
        'minDR_leadAK8otherAK8'] 

# %%
### Perform permuations calculations
## Make a big list of lists of different lengths and combinations
# TODO
combos = []
for size in range(7,len(varList)):
   thissize = list(itertools.combinations(varList,size))
   for item in thissize:
      ## Some sanity checks for vars we know for sure we want to include
    #   if 'pNet_J_1' not in item: continue
    #   if 'pNet_J_2' not in item and 'dnnJ_3' not in item: continue
    #   if 'AK4HT' not in item: continue
    #   if 'corr_met_MultiLepCalc' not in item and 'AK4HTpMETpLepPt' not in item: continue
    #   if 'NJetsDeepFlavwithSF_JetSubCalc' not in item and 'NJets_JetSubCalc' not in item: continue
    #   if 't_mass' not in item and 't_pt' not in item: continue
    #   if 't_dRWb' not in item and 'minDR_leadAK8otherAK8' not in item: continue
      combos.append(item)
combos.append(varList)

if year == 'all':
   combos = []
   for size in range(12,14):
      thissize = list(itertools.combinations(varList,size))
      for item in thissize:
         ## Some sanity checks for vars we know for sure we want to include
        #  if 'dnnJ_1' not in item: continue
        #  if 'dnnJ_2' not in item: continue
        #  if 'dnnJ_3' not in item: continue
        #  if 'AK4HT' not in item: continue
        #  if 'AK4HTpMETpLepPt' not in item: continue
        #  if 'NJetsDeepFlavwithSF_JetSubCalc' not in item: continue
        #  if 'NJets_JetSubCalc' not in item: continue
        #  if 't_mass' not in item and 't_pt' not in item: continue
        #  if 'minDR_leadAK8otherAK8' not in item: continue
         combos.append(item)

# %%
### Select input variables for training and testing
vars = list(combos[vararray])
print('Vars = ', vars)

indexKill = range(0, len(varList))
for item in vars:
   indexKill.remove(varList.index(item))

