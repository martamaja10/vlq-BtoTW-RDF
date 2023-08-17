# Arguments: 
# The first argument is True or False depending on whether you want to query the dasgoclient to make lists of root files in Nano List text files
# The second argument is True or False depending on whether you want to submit a condor job with the sample specified in the thrid argument
# The third argument is the prefix of the sample you want to use when you submit a condor job

import os,sys,shutil,datetime,time, subprocess
from samples import *

execfile("/uscms_data/d3/kjohnso/EOSSafeUtils.py") # this is a python2 command, so ignore the error and run with python2
start_time = time.time()

# --- Sample Dictionary ---
sample_dic = samples_2018UL # This is the name of the list (using list of class objects to keep ordering)

# --- Size of Condor Job ---
filesPerJob = 999
jobsPerSample = 10 # This is only used for the TTToSemiLeptonic samples

# Take in Arguments from the command line
makelists = False
runanalyzer = False
if len(sys.argv) >= 2: makelists = bool(eval(sys.argv[1]))
if len(sys.argv) >= 3: runanalyzer = bool(eval(sys.argv[2]))
if len(sys.argv) >= 4: 
    prefix = sys.argv[3] # 'singleTb'
    textlist = prefix + "NanoList.txt"
    
relbase = '/uscms/home/jmanagan/nobackup/BtoTW/CMSSW_11_0_0/'
outDir='/store/user/jmanagan/BtoTW_Aug2023_2018/'
condorDir='/uscms/home/jmanagan/nobackup/BtoTW/rdfjobs_2018/' # recommend this be outside git area!
tarfile = '/uscms/home/jmanagan/nobackup/rdfjobs.tar' # outside the CMSSW

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
    print ('--- DAS queries ---')
    for k,v in sample_dic.items():
        query = queryStatement + v.samplename + "\" > " + "NanoList/" + v.textlist
        print ('\t'+v.textlist)
        os.system(query)
        
# --- Submit Condor Jobs ---
if runanalyzer:    
    
    if (not os.path.isdir('Rootfiles')): os.system("mkdir Rootfiles")
    if (not os.path.isdir('Output')): os.system("mkdir Output")
    
    # ------ Making the Tar File ------
    
    print ('--- Making tar ---')
    if os.path.exists(tarfile): print ('*********** tar already exists! I ASSUME YOU WANT TO MAKE A NEW ONE! *************')
    os.chdir(relbase)
    print ('tar --exclude="src/.git" --exclude="src/*/.git" --exclude="tmp/" --exclude="src/vlq-BtoTW-SLA" --exclude="src/VLQRDF" --exclude="src/vlq-singleProd-RDF" --exclude="src/vlq-BtoTW-RDF/*.root" --exclude=".SCRAM" -zcf '+tarfile+' ./*')
    os.system('tar --exclude="src/.git" --exclude="src/*/.git" --exclude="tmp/" --exclude="src/vlq-BtoTW-SLA" --exclude="src/VLQRDF" --exclude="src/vlq-singleProd-RDF" --exclude="src/vlq-BtoTW-RDF/*.root" --exclude=".SCRAM" -zcf '+tarfile+' ./*')
    os.chdir(runDir)

    for k,v in sample_dic.items(): # This is fine.  It is a dictionary that exists in samples.py that gets imported at runtime
        
        prefix = v[0].prefix;
        textlist = "NanoList/" + v[0].textlist;
        year = v[0].year;

        # --- Check sample size and aim for <= 50GB per job
        command = '/cvmfs/cms.cern.ch/common/dasgoclient --query="dataset dataset=/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM | grep dataset.size" '
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        (out, err) = proc.communicate()
        try: 
            samplesize = int(out.split('\n')[0])
        except:    
            try: 
                samplesize = int(out.split('\n')[1])
            except:
                try: 
                    samplesize = int(out.split('\n')[2])
                except: 
                    print 'need more than 2 levels to get sample size'
                    exit(1)

        jobsPerSample = max(1,round(samplesize/50000000000.))
        
        rootfiles = ''
        num = count = newNum = j = 0
        
        # --- Making string of root files --- (Delimiter is  "\n")
        with open(os.path.abspath(textlist),'r') as rootlist:
            for line in rootlist:
                num += 1
                rootfiles += "root://cmsxrootd.fnal.gov/" + line.strip() + "\n"
            rootfiles = rootfiles[:-1]
        print ('Number of Rootfiles: ' + str(num))
    
        #if "TTToSemiLeptonic" in prefix: 
        filesPerJob = int(max(1,round(num/jobsPerSample)))
        
        os.system('eos root://cmseos.fnal.gov/ mkdir -p '+outDir+'/')
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
            
            dict={'RUNDIR':runDir, 'CONDORDIR':condorDir+'/'+prefix, 'CMSSWBASE':relbase, 'OUTPUTDIR':outDir, 'TARBALL':tarfile, 'TESTNUM1':oldNum, 'TESTNUM2':newNum-1, 'PREFIX':prefix, 'FILENAME':fileName, 'YEAR':year}
            jdfName=condorDir+prefix+'/%(PREFIX)s_%(TESTNUM1)s.job'%dict
            print ("jdfname: ",jdfName)
            jdf=open(jdfName,'w')
            jdf.write(
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

print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))




