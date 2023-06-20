# vlq-BtoTW-RDF
## RDataFrame code for B->tW search on NanoAOD

This fork contains some machine learning extensions to the preexisting R-Dataframe code being used in the B->tW search. 

## Structure of this repository

Most of the critical files are located in the root directory of the repository, namely the .cc and .cpp files. Each of these files plays a particular role in running the analysis.

Additionally, the python files in this repository are important for analyzing and training neural networks on the data produced by the RDF scripts. These are also located in the root directory. 

## Critical Files

A few key files in this repository are listed below along with a brief description of how they play into the analysis.

- `callRDF.C` - This file is the primary interface between the user/condor and the code in this repository. The file takes input in the form "[channel]" [testnum] [inputfile] where input file is a file path. The file then loads in the liblwtnnlwtnn.so library (for access to light-weight neural network access) and calls runRDF.C with the same inputs.

- `runRDF.C` - Similar to the above file, `runRDF.C` primarily calls another file to do its work. In more detail, the file creates an rdf object t with the user defined arguments in addition to some preset arguments (including analysis year) and calls the analyzer_RDF function from `analyzer_RDF.cpp` on this object. The goal of this file is to manage the RDF instance and include any extra functions also included in this repository.

- `analyzer_RDF.cpp` - This file is where the data begins being altered. The key function in this file is analyzer_RDF. This function takes in the user input and creates a dataframe. Initially, the function performs some preprocessing, setting up necessary variables for analysis. Then, the function implements the LWTNN library. Next, the function applies some flags and filtering to narrow the dataset down to only include necessary and "good" events. Finally, the function performs an analysis on the now cut data and creates a final snapshot file as an output. Note this cpp file has a corresponding header file.

So, in summary, the flow is `callRDF.C` -> `runRDF.C` -> `analyzer_RDF.cpp` with `callRDF.C` being the only file the user needs to interact with at runtime. 

In addition to these files, there are three primary python files.

- `root2NPZ.py` - This file takes the output of the RDF scripts and performs some further selections and cuts to make the data useful for training machine learning models. The output of this file is a numpy array (stored as a .npz) which stores an archive of the data. The archive contains a set of subarrays, one for each background type and signal mass. The second dimension is events and the third is features.

- `trainMLP.py` - This file takes the output of `root2NPZ` and trains a model on it. This model is a multilayer perceptron. This file trains this models, testsit on the testing data, and creates some plots analyzing the performance of the classifier. At the end, this models are saved to an h5 file along with a data scaler so they can be passed back into the RDF scripts to generate better plots. The input directory for this script should be pointed to the output directory of `root2NPZ.py`.

## Getting started

In order to run the code, the user must have access to cmslpc-sl7.fnal.gov and the appropriate setup that comes with that (as described in [this tutorial](https://fnallpc.github.io/cms-das-pre-exercises/01-CMSDataAnalysisSchoolPreExerciseFirstSet/index.html)). 

Once access is set up and the user has the github repository configured along with a working GRID certificate, the program should be runnable.

To actually run the analysis, use the following command: root -l -b -q callRDF.C\\(\\"`Channel`\\",\\"`testNum`\\",\\"`inputfile`\\"\\)

- `Channel` = “Muon” or “Electron”
- `testNum` = an integer
- `inputfile` = path to input file

It's important to note that every parenthesis and quotation must be escaped for proper syntax. 
