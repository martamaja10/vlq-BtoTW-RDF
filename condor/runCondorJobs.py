# Arguments: 
# The first argument is True or False depending on whether you want to query the dasgoclient to make lists of root files in Nano List text files
# The second argument is True or False depending on whether you want to submit a condor job with the sample specified in the thrid argument
# The third argument is the prefix of the sample you want to use when you submit a condor job

import os,sys,shutil,datetime,time
from samples import *

execfile("/uscms_data/d3/kjohnso/EOSSafeUtils.py") # this is a python2 command, so ignore the error and run with python2
start_time = time.time()

# --- Size of Condor Job ---
filesPerJob = 10 # Not used when jobsPerSample is used
jobsPerSample = 10

# --- Sample Dictionary ---
sample_dic = samples_WJets # This is the name of the dictionary we want to run over (the yellow underline doesn't prevent it from runniing)

# Take in Arguments from the command line
makelists = False
runanalyzer = False
if len(sys.argv) >= 2: makelists = bool(eval(sys.argv[1]))
if len(sys.argv) >= 3: runanalyzer = bool(eval(sys.argv[2]))
if len(sys.argv) >= 4: 
    prefix = sys.argv[3] # 'singleTb'
<<<<<<< HEAD
    textlist = prefix + "NanoList.txt"
    
relbase = '/uscms/home/kjohnso/nobackup/BtoTW/CMSSW_11_0_0/'
outDir='/store/user/kjohnso/BtoTW_Jul2023/WJets_MergeTest/'
condorDir='/uscms/home/kjohnso/nobackup/BtoTW/CMSSW_11_0_0/src/vlq-BtoTW-RDF/condor_2/' # recommend this be outside git area!
tarfile = '/uscms/home/kjohnso/nobackup/rdfjobs.tar' # outside the CMSSW
=======
    textlist = prefix + "NanoList.txt" # 'singleTbNanoList.txt'

relbase = '/uscms/home/jmanagan/nobackup/BtoTW/CMSSW_11_0_0/'
outDir='/store/user/jmanagan/BtoTW_Jul2023/LeptonChecks/'
condorDir='/uscms/home/jmanagan/nobackup/BtoTW/condor_July2023/' # recommend this be outside git area!
tarfile = '/uscms/home/jmanagan/nobackup/rdfjobs.tar' # outside the CMSSW
>>>>>>> becca63a8b6c385fbdbc9267bf9faf41de37fe8c

runDir=os.getcwd()
cTime=datetime.datetime.now()
date='%i_%i_%i_%i_%i_%i'%(cTime.year,cTime.month,cTime.day,cTime.hour,cTime.minute,cTime.second)

print ('--- Starting Submission ---')

# --- Make Lists ---
if makelists:
    queryStatement = '/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /'
        
    if (os.path.isdir('NanoList')): os.system("rm -r NanoList")
    os.system("mkdir NanoList")
    
    # Loops through all the samples listed in samples.py
    for k,v in samples.items():
        query = queryStatement + v.samplename + "\" > " + "NanoList/" + v.textlist
        os.system(query)
        
# --- Submit Condor Jobs ---
if runanalyzer:    
<<<<<<< HEAD
    
    if (not os.path.isdir('Rootfiles')): os.system("mkdir Rootfiles")
    if (not os.path.isdir('Output')): os.system("mkdir Output")
    
    # ------ Making the Tar File ------
    
    print ('--- Making tar ---')
=======

    print ('Making tar:')
>>>>>>> becca63a8b6c385fbdbc9267bf9faf41de37fe8c
    if os.path.exists(tarfile): print ('*********** tar already exists! I ASSUME YOU WANT TO MAKE A NEW ONE! *************')
    os.chdir(relbase)
    print ('tar --exclude="src/.git" --exclude="tmp/" --exclude="src/VLQRDF" --exclude="src/vlq-singleProd-RDF" --exclude="src/vlq-BtoTW-RDF/*.root" --exclude=".SCRAM" -zcf '+tarfile+' ./*')
    os.system('tar --exclude="src/.git" --exclude="tmp/" --exclude="src/VLQRDF" --exclude="src/vlq-singleProd-RDF" --exclude="src/vlq-BtoTW-RDF/*.root" --exclude=".SCRAM" -zcf '+tarfile+' ./*')
    os.chdir(runDir)

<<<<<<< HEAD
    for k,v in sample_dic.items(): # This is fine.  It is a dictionary that exists in samples.py that gets imported at runtime
        
        prefix = v.prefix;
        textlist = "NanoList/" + v.textlist;
        year = v.year;
        
        rootfiles = ''
        num = count = newNum = j = 0
        
        # --- Making string of root files --- (Delimiter is  "\n")
        with open(os.path.abspath(textlist),'r') as rootlist:
            for line in rootlist:
                num += 1
                rootfiles += "root://cmsxrootd.fnal.gov/" + line.strip() + "\n"
            rootfiles = rootfiles[:-1]
        print ('Number of Rootfiles: ' + str(num))
    
        if "TTToSemiLeptonic" in prefix: 
            filesPerJob = int(max(1,round(num/jobsPerSample)))
        else:
            filesPerJob = 999
        
        os.system('eos root://cmseos.fnal.gov/ mkdir -p '+outDir+'/'+prefix)
        os.system('mkdir -p '+condorDir+'/'+prefix)
        
        # --- Write rootfiles to a textfile --- Read in analyzer_RDF.h in the constructor
        fileName = "Rootfiles/" + prefix + ".txt" # JH removing "condorDir" because I want condorDir to be outside git area
        listArgument = open(fileName, "w")
        listArgument.write(rootfiles)
        listArgument.close()
        print ('--- Wrote to File ---')
    
        # Redefining fileName so it is accessed from the directory above for analyzer_RDF.h
        fileName = "condor/Rootfiles/"+prefix+".txt"
        
        print ('--- Submitting Condor Jobs ---')
        for i in range(filesPerJob,num+filesPerJob,filesPerJob):
            oldNum = newNum
            newNum = i
            print("Num test: " + str(oldNum) + " -> " + str(newNum))
            
            dict={'RUNDIR':runDir, 'CONDORDIR':condorDir+'/'+prefix, 'CMSSWBASE':relbase, 'OUTPUTDIR':outDir+prefix, 'TARBALL':tarfile, 'TESTNUM1':oldNum, 'TESTNUM2':newNum-1, 'PREFIX':prefix, 'FILENAME':fileName, 'YEAR':year}
            jdfName=condorDir+prefix+'/%(PREFIX)s_%(TESTNUM1)s.job'%dict
            print ("jdfname: ",jdfName)
            jdf=open(jdfName,'w')
            jdf.write(
=======
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
>>>>>>> becca63a8b6c385fbdbc9267bf9faf41de37fe8c
            """use_x509userproxy = true
universe = vanilla
Executable = %(RUNDIR)s/condorRDF.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = %(TARBALL)s
Output = %(PREFIX)s_RDF_%(TESTNUM1)s.out
Error = %(PREFIX)s_RDF_%(TESTNUM1)s.err
Log = %(PREFIX)s_RDF_%(TESTNUM1)s.log
Notification = Never
<<<<<<< HEAD
Arguments = %(FILENAME)s %(OUTPUTDIR)s %(TESTNUM1)s %(TESTNUM2)s %(YEAR)s 

Queue 1"""%dict)
            jdf.close()
            os.chdir('%s/'%(condorDir+'/'+prefix))
            os.system('condor_submit %(PREFIX)s_%(TESTNUM1)s.job'%dict)
            os.system('sleep 0.5')                                
            os.chdir('%s'%(runDir))
            print ( j, " jobs submitted!!!")
            j += 1
        
        #Formatting line 101
                
    #for k,v in sample_dic.items():
    #    os.system('mv ' + prefix + ' Output/')
=======
Arguments = %(FILENAME)s %(OUTPUTDIR)s %(TESTNUM)s 

Queue 1"""%dict)
        jdf.close()
        os.chdir('%s/'%(condorDir+'/'+prefix))
        os.system('condor_submit %(PREFIX)s_%(TESTNUM)s.job'%dict)
        os.system('sleep 0.5')                                
        os.chdir('%s'%(runDir))
        print (count, " jobs submitted!!!")
        
        #Formatting line 101
>>>>>>> becca63a8b6c385fbdbc9267bf9faf41de37fe8c

print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))




