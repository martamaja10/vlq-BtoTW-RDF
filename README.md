# vlq-BtoTW-RDF
## RDataFrame code for B->tW search on NanoAOD

This fork contains some machine learning extensions to the preexisting R-Dataframe code being used in the B->tW search. 

## Structure of this repository

Most of the critical files are located in the root directory of the repository, namely the .cc and .cpp files. Each of these files plays a particular role in running the analysis.

## Critical Files

A few key files in this repository are listed below along with a brief description of how they play into the analysis.

- `callRDF.C` - This file is the primary interface between the user/condor and the code in this repository. The file takes input in the form "[channel]" [testnum] [inputfile] where input file is a file path. The file then loads in the liblwtnnlwtnn.so library (for access to light-weight neural network access) and calls runRDF.C with the same inputs.

- `runRDF.C` - Similar to the above file, `runRDF.C` primarily calls another file to do its work. In more detail, the file creates an rdf object t with the user defined arguments in addition to some preset arguments (including analysis year) and calls the analyzer_RDF function from `analyzer_RDF.cpp` on this object. The goal of this file is to manage the RDF instance and include any extra functions also included in this repository.

- `analyzer_RDF.cpp` - This file is where the data begins being altered. The key function in this file is analyzer_RDF. This function takes in the user input and creates a dataframe. Initially, the function performs some preprocessing, setting up necessary variables for analysis. Then, the function implements the LWTNN library. Next, the function applies some flags and filtering to narrow the dataset down to only include necessary and "good" events. Finally, the function performs an analysis on the now cut data and creates a final snapshot file as an output. Note this cpp file has a corresponding header file.

So, in summary, the flow is `callRDF.C` -> `runRDF.C` -> `analyzer_RDF.cpp` with `callRDF.C` being the only file the user needs to interact with at runtime. 

## Getting started

In order to run the code, the user must have access to cmslpc-sl7.fnal.gov and the appropriate setup that comes with that (as described in [this tutorial](https://fnallpc.github.io/cms-das-pre-exercises/01-CMSDataAnalysisSchoolPreExerciseFirstSet/index.html)). 

Once access is set up and the user has the github repository configured along with a working GRID certificate, the program should be runnable.

To actually run the analysis, use the following command: root -l -b -q callRDF.C\\(\\"`Channel`\\",\\"`testNum`\\",\\"`inputfile`\\"\\)

- `Channel` = “Muon” or “Electron”
- `testNum` = an integer
- `inputfile` = path to input file

It's important to note that every parenthesis and quotation must be escaped for proper syntax. 
