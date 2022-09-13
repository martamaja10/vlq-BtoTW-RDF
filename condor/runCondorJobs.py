import os,sys,shutil,datetime,time
from ROOT import *

execfile("/uscms_data/d3/jmanagan/EOSSafeUtils.py")

start_time = time.time()


relbase = '/uscms_data/d3/jmanagan/CMSSW_11_0_0/'
outDir='/store/user/jmanagan/SingleProd_RDF'
condorDir='/uscms_data/d3/jmanagan/rdfjobs/'
tarfile = '/uscms_data/d3/jmanagan/rdfjobs.tar'

runDir=os.getcwd()

cTime=datetime.datetime.now()
date='%i_%i_%i_%i_%i_%i'%(cTime.year,cTime.month,cTime.day,cTime.hour,cTime.minute,cTime.second)

print 'Making tar:'
if os.path.exists(tarfile): print '*********** tar already exists! I ASSUME YOU WANT TO MAKE A NEW ONE! *************'

os.chdir(relbase)
print 'tar --exclude="src/.git" --exclude="tmp/" --exclude="src/AOD2NanoAODOutreachTool" --exclude="src/practice" --exclude="src/testRDataFrame" --exclude="src/vlq-1lepDnn-RDF" --exclude=".SCRAM" -zcf '+tarfile+' ./*'
os.system('tar --exclude="src/.git" --exclude="tmp/" --exclude="src/AOD2NanoAODOutreachTool" --exclude="src/practice" --exclude="src/testRDataFrame" --exclude="src/vlq-1lepDnn-RDF" --exclude=".SCRAM" -zcf '+tarfile+' ./*')
os.chdir(runDir)

print 'Starting submission'
count=0

#if not os.path.exists('TTbarNanoList.txt'):
#    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > TTbarNanoList.txt')

rootfiles = []
with open(os.path.abspath('TTbarNanoList.txt'),'r') as rootlist:
    for line in rootlist:
        rootfiles.append('root://cmsxrootd.fnal.gov/'+line.strip())
print '\tTotal files:',len(rootfiles)

os.system('eos root://cmseos.fnal.gov/ mkdir -p '+outDir)
os.system('mkdir -p '+condorDir)

for ifile in rootfiles:

    count+=1
    #if count == 1: continue

    dict={'RUNDIR':runDir, 'CONDORDIR':condorDir, 'FILENAME':ifile, 'CMSSWBASE':relbase, 'OUTPUTDIR':outDir, 'TARBALL':tarfile, 'TESTNUM':count}
    jdfName=condorDir+'/TTbar_%(TESTNUM)s.job'%dict
    print "jdfname: ",jdfName
    jdf=open(jdfName,'w')
    jdf.write(
        """use_x509userproxy = true
universe = vanilla
Executable = %(RUNDIR)s/condorRDF.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = %(TARBALL)s
Output = TTbarRDF_%(TESTNUM)s.out
Error = TTbarRDF_%(TESTNUM)s.err
Log = TTbarRDF_%(TESTNUM)s.log
Notification = Never
Arguments = %(FILENAME)s %(OUTPUTDIR)s %(TESTNUM)s

Queue 1"""%dict)
    jdf.close()
    os.chdir('%s/'%(condorDir))
    os.system('condor_submit TTbar_%(TESTNUM)s.job'%dict)
    os.system('sleep 0.5')                                
    os.chdir('%s'%(runDir))
    print count, "jobs submitted!!!"


print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))





