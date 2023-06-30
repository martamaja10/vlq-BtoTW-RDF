# Arguments: 
# The first argument is True or False depending on whether you want to query the dasgoclient to make lists of root files in Nano List text files
# The second argument is True or False depending on whether you want to submit a condor job with the sample specified in the thrid argument
# The third argument is the prefix of the sample you want to use when you submit a condor job

import os,sys,shutil,datetime,time
from samples import *

execfile("/uscms_data/d3/kjohnso/EOSSafeUtils.py") # this is a python2 command, so ignore the error and run with python2
start_time = time.time()

makelists = False
runanalyzer = False
if len(sys.argv) >= 2: makelists = bool(eval(sys.argv[1]))
if len(sys.argv) >= 3: runanalyzer = bool(eval(sys.argv[2]))
if len(sys.argv) >= 4: 
    prefix = sys.argv[3] # 'singleTb'
    textlist = prefix + "NanoList.txt" # 'singleTbNanoList.txt'

relbase = '/uscms/home/jmanagan/nobackup/BtoTW/CMSSW_11_0_0/'
outDir='/store/user/jmanagan/BtoTW_Jul2023/LeptonChecks/'
condorDir='/uscms/home/jmanagan/nobackup/BtoTW/condor_July2023/' # recommend this be outside git area!
tarfile = '/uscms/home/jmanagan/nobackup/rdfjobs.tar' # outside the CMSSW

runDir=os.getcwd()

cTime=datetime.datetime.now()
date='%i_%i_%i_%i_%i_%i'%(cTime.year,cTime.month,cTime.day,cTime.hour,cTime.minute,cTime.second)

print ('---Starting Submission---')
count=0

if makelists:
    
    queryStatement = '/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /'
    additionalPath = '/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM'
    
    # Loops through all the samples listed in samples.py
    for k,v in samples.items():
        # BPrime given to me in correct format in a text file
        if (v.prefix != 'WJets2500' and v.prefix != 'BPrime800' and v.prefix != 'BPrime1400' and v.prefix != 'BPrime2000'):
            query = queryStatement + v.samplename + additionalPath + "\" > " + v.textlist 
            os.system(query)
        elif (v.prefix == 'WJets2500'):
            # 2500 is a special case because it is using v2 instead of v1
            query = queryStatement + v.samplename + '/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM' + "\" > " + v.textlist
            os.system(query)
        
if runanalyzer:    

    print ('Making tar:')
    if os.path.exists(tarfile): print ('*********** tar already exists! I ASSUME YOU WANT TO MAKE A NEW ONE! *************')

    os.chdir(relbase)
    print ('tar --exclude="src/.git" --exclude="tmp/" --exclude="src/VLQRDF" --exclude="src/vlq-singleProd-RDF" --exclude="src/vlq-BtoTW-RDF/*.root" --exclude=".SCRAM" -zcf '+tarfile+' ./*')
    os.system('tar --exclude="src/.git" --exclude="tmp/" --exclude="src/VLQRDF" --exclude="src/vlq-singleProd-RDF" --exclude="src/vlq-BtoTW-RDF/*.root" --exclude=".SCRAM" -zcf '+tarfile+' ./*')
    os.chdir(runDir)

    for k,v in samples.items():
        prefix = v.prefix;
        textlist = v.textlist;
        
        #Indent everything below this and comment the two lines at the top of this section
        #to submit condor jobs for all samples at at once

    
        rootfiles = ''
    
        # Delimiter is  "\n" so I can write and read by line
        with open(os.path.abspath(textlist),'r') as rootlist:
            for line in rootlist:
                if 'Bprime' not in textlist:
                    rootfiles += "root://cmsxrootd.fnal.gov/" + line.strip() + "\n"
                else:
                    rootfiles += "root://cmseos.fnal.gov/" + line.strip() + "\n"
            rootfiles = rootfiles[:-1]
        print ('Length of Rootfiles: ',len(rootfiles))
    
        # Write rootfiles to a textfile. Read in analyzer_RDF.h in the constructor
        fileName = prefix+".txt" # JH removing "condorDir" because I want condorDir to be outside git area
        listArgument = open(fileName, "w")
        listArgument.write(rootfiles)
        listArgument.close()
        print ('Wrote to File')
    
        # Redefining fileName so it is accessed from the directory above by analyzer_RDF.h
        fileName = "condor/"+prefix+".txt"

        os.system('eos root://cmseos.fnal.gov/ mkdir -p '+outDir+'/'+prefix)
        os.system('mkdir -p '+condorDir+'/'+prefix)
                
        dict={'RUNDIR':runDir, 'CONDORDIR':condorDir+'/'+prefix, 'CMSSWBASE':relbase, 'OUTPUTDIR':outDir+'/'+prefix, 'TARBALL':tarfile, 'TESTNUM':count, 'PREFIX':prefix, 'LIST':rootfiles, 'FILENAME':fileName}
        jdfName=condorDir+prefix+'/%(PREFIX)s_%(TESTNUM)s.job'%dict
        print ("jdfname: ",jdfName)
        jdf=open(jdfName,'w')
        jdf.write(
            """use_x509userproxy = true
universe = vanilla
Executable = %(RUNDIR)s/condorRDF.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = %(TARBALL)s
Output = %(PREFIX)s_RDF_%(TESTNUM)s.out
Error = %(PREFIX)s_RDF_%(TESTNUM)s.err
Log = %(PREFIX)s_RDF_%(TESTNUM)s.log
Notification = Never
Arguments = %(FILENAME)s %(OUTPUTDIR)s %(TESTNUM)s 

Queue 1"""%dict)
        jdf.close()
        os.chdir('%s/'%(condorDir+'/'+prefix))
        os.system('condor_submit %(PREFIX)s_%(TESTNUM)s.job'%dict)
        os.system('sleep 0.5')                                
        os.chdir('%s'%(runDir))
        print (count, " jobs submitted!!!")
        
        #Formatting line 101

print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))




