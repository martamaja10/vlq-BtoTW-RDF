# Arguments: 
# The first argument is True or False depending on whether you want to query the dasgoclient to make lists of root files in Nano List text files
# The second argument is True or False depending on whether you want to submit a condor job with the sample specified in the thrid argument
# The third argument is the prefix of the sample you want to use when you submit a condor job

import os,sys,shutil,datetime,time, subprocess
parent = os.path.dirname(os.getcwd())
sys.path.append(parent+"/condor/")
from samples import *

exec(open("/uscms_data/d3/jmanagan/EOSSafeUtils.py").read()) # this is a python2 command, so ignore the error and run with python2
start_time = time.time()

# --- Sample Dictionary ---
sample_dic = samples # This is the name of the list (using list of class objects to keep ordering)

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
    
relbase = '/uscms/home/jmanagan/nobackup/BtoTW/CMSSW_12_4_8/'
outDir='/store/user/jmanagan/BtoTW_Apr2024_pNetEffs_tWseparate/'
condorDir='/uscms/home/jmanagan/nobackup/BtoTW/rdfjobs_Apr2024_pNetEffs_tWseparate/' # recommend this be outside git area!
tarfile = '/uscms/home/jmanagan/nobackup/rdfjobs.tar' # outside the CMSSW

runDir=os.getcwd()
cTime=datetime.datetime.now()
date='%i_%i_%i_%i_%i_%i'%(cTime.year,cTime.month,cTime.day,cTime.hour,cTime.minute,cTime.second)

print ('--- Starting Submission ---')

# --- Make Lists ---
if makelists:
    queryStatement = '/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /'
        
    #if (os.path.isdir('NanoList')): os.system("rm -r NanoList")
    #os.system("mkdir NanoList")
    
    # Loops through all the samples listed in samples.py
    print ('--- DAS queries ---')
    for k,v in sample_dic.items():
        if os.path.exists("NanoList/"+v.textlist): 
            print(v.textlist+' exists, skipping')
            continue
        query = queryStatement + v.samplename + "\" > " + "NanoList/" + v.textlist
        print ('\t'+v.textlist)
        os.system(query)

        # --- Making string of root files --- (Delimiter is  "\n")
        with open(os.path.abspath("NanoList/"+v.textlist),'r') as rootlist:
            rootfiles = [''.join(["root://cmsxrootd-site.fnal.gov/", line.strip(), "\n"]) for line in rootlist.readlines()]

        with open(os.path.abspath("NanoList/"+v.textlist),'w') as rootlist:
            rootlist.writelines(rootfiles)
        
# --- Submit Condor Jobs ---
if runanalyzer:    
    
    if (not os.path.isdir('Rootfiles')): os.system("mkdir Rootfiles")
    if (not os.path.isdir('Output')): os.system("mkdir Output")
    
    # ------ Making the Tar File ------
    
    print ('--- Making tar ---')    
    if os.path.exists(tarfile): 
        print ('*********** tar already exists! I ASSUME YOU WANT TO MAKE A NEW ONE! *************')
        os.system('rm '+tarfile)
    os.chdir(relbase)
    print ('tar --exclude="src/.git" --exclude="src/*/.git" --exclude="tmp/" --exclude="src/vlq-BtoTW-SLA" --exclude="src/VLQRDF" --exclude="src/vlq-singleProd-RDF" --exclude="src/vlq-BtoTW-RDF/*.root" --exclude="src/vlq-BtoTW-RDF/*.npz" --exclude="src/vlq-BtoTW-RDF/DeepJetEffs/" --exclude="src/vlq-BtoTW-RDF/PNetEffs/FirstPass/ --exclude=".SCRAM" -zcf '+tarfile+' ./*')
    os.system('tar --exclude="src/.git" --exclude="src/*/.git" --exclude="tmp/" --exclude="src/vlq-BtoTW-SLA" --exclude="src/VLQRDF" --exclude="src/vlq-singleProd-RDF" --exclude="src/vlq-BtoTW-RDF/*.root" --exclude="src/vlq-BtoTW-RDF/*.npz" --exclude="src/vlq-BtoTW-RDF/DeepJetEffs/" --exclude="src/vlq-BtoTW-RDF/PNetEffs/FirstPass/" --exclude=".SCRAM" -zcf '+tarfile+' ./*')
    os.chdir(runDir)

    count = 0

    for k,v in sample_dic.items(): # This is fine.  It is a dictionary that exists in samples.py that gets imported at runtime
        prefix = v.prefix;
        textlist = "NanoList/" + v.textlist;
        year = v.year;

        if 'TTWq2016APV' not in prefix: continue

        print('Submitting ' + prefix)
            
        os.system('eos root://cmseos.fnal.gov/ mkdir -p '+outDir+'/')
        os.system('mkdir -p '+condorDir+'/'+prefix)
        
        # Redefining fileName so it is accessed from the directory above for analyzer_RDF.h
        fileName = "../condor/"+textlist
        
        dict={'RUNDIR':runDir, 'CONDORDIR':condorDir+'/'+prefix, 'CMSSWBASE':relbase, 'OUTPUTDIR':outDir, 'TARBALL':tarfile, 'PREFIX':prefix, 'FILENAME':fileName, 'YEAR':year}
        jdfName=condorDir+prefix+'/%(PREFIX)s.job'%dict
        print ("jdfname: ",jdfName)
        jdf=open(jdfName,'w')
        jdf.write(
            """use_x509userproxy = true
universe = vanilla
Executable = %(RUNDIR)s/condorRDF.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = %(TARBALL)s
Output = %(PREFIX)s.out
Error = %(PREFIX)s.err
Log = %(PREFIX)s.log
Notification = Never
Arguments = %(FILENAME)s %(OUTPUTDIR)s %(YEAR)s 

Queue 1"""%dict)
        jdf.close()
        os.chdir('%s/'%(condorDir+'/'+prefix))
        os.system('condor_submit %(PREFIX)s.job'%dict)
        os.system('sleep 0.5')                                
        os.chdir('%s'%(runDir))
        print ( str(count) + " jobs submitted!!!")
        count += 1
        

print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))




