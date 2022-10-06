#!/bin/bash

hostname

infilename=${1}
outputDir=${2}
testnum=${3}

scratch=${PWD}
macroDir=${PWD}
source /cvmfs/cms.cern.ch/cmsset_default.sh
scramv1 project CMSSW CMSSW_11_0_0
cd CMSSW_11_0_0

echo "unpacking tar"
tar -xf ${scratch}/rdfjobs.tar
rm ${scratch}/rdfjobs.tar

eval `scramv1 runtime -sh`
cd src/vlq-BtoTW-RDF/

export PATH=$PATH:$macroDir

echo "Running RDF:"
root -l -b -q callRDF.C\(\"${testnum}\",\"${infilename}\"\)

echo "ROOT Files:"
ls -l *.root

echo "Removing presel files:"
rm *preselTree*.root

# copy output to eos

echo "xrdcp output for condor"
for FILE in *.root
do
  echo "xrdcp -f ${FILE} root://cmseos.fnal.gov/${outputDir}/${FILE}"
  xrdcp -f ${FILE} root://cmseos.fnal.gov/${outputDir}/${FILE} 2>&1
  XRDEXIT=$?
  if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
  fi
  rm ${FILE}
done

echo "done"
