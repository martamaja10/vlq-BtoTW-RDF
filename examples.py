>>> from ROOT import *
>>> tree = RDataFrame("Events","root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/Bp800_hadd.root");
>>> seltrain = "NJets_forward == 0"
>>> seltest = "NJets_forward > 0"

### Filter and define "weight"
>>> train = tree.Filter(seltrain)
>>> test = tree.Filter(seltest)
>>> train2 = train.Define("weight","42.0");
>>> test2 = test.Define("weight","41.0");
>>> weightmean = train2.Mean('weight');

### Check that weight exists sensibly -- yes
>>> print weightmean.GetValue()
42.0

### Do the "as numpy" -- I had forgotten to use "columns="!!!
>>> npytrain = train2.AsNumpy(columns=['weight','NJets_forward'])
>>> npytest = test2.AsNumpy(columns=['weight','NJets_forward'])

### Print the output of AsNumpy to see the structure
>>> print npytrain
{'NJets_forward': numpy.array([0, 0, 0, ..., 0, 0, 0], dtype=int32), 'weight': numpy.array([42., 42., 42., ..., 42., 42., 42.])}
>>> print npytest 
{'NJets_forward': numpy.array([1, 2, 1, ..., 1, 1, 2], dtype=int32), 'weight': numpy.array([41., 41., 41., ..., 41., 41., 41.])}
### Interesting...dictionaries: key = input variable, value = 1-D numpy array
### Filtered contents and the numeric weights look right

### Compare to the old method
>>> tfile = TFile.Open("root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/Bp800_hadd.root");
>>> ttree = tfile.Get("Events");
>>> from root_numpy import tree2array
>>> trainOldMethod = []
>>> testOldMethod = []
>>> vars = ['NJets_forward','Jet_HT']  ## "weight" doesn't exist here yet
>>> trainOldMethod = tree2array(ttree,vars,seltrain)
>>> testOldMethod = tree2array(ttree,vars,seltest)
>>> print trainOldMethod
[(0, 777.1875 ) (0, 778.     ) (0, 725.25   ) ... (0, 557.96875)]
>>> print testOldMethod
[(1, 880.75   ) (2, 571.25   ) (1, 512.75   ) ... (1, 745.875  )]
### As expected, this output is [ (event 1 stuff) (event 2 stuff) ... (event N stuff)]
### Here, stuff = (NJets_forward, Jet_HT), and teh values make sense. 
