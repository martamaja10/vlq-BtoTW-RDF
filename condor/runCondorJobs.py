import os,sys,shutil,datetime,time
from ROOT import *

execfile("/uscms_data/d3/jmanagan/EOSSafeUtils.py")

start_time = time.time()

makelists = False
if len(sys.argv) > 1: makelists = eval(bool(sys.argv[1]))

relbase = '/uscms_data/d3/jmanagan/BtoTW/CMSSW_11_0_0/'
outDir='/store/user/jmanagan/BtoTW_RDF'
condorDir='/uscms_data/d3/jmanagan/BtoTW/rdfjobs/'
tarfile = '/uscms_data/d3/jmanagan/BtoTW/rdfjobs.tar'

runDir=os.getcwd()

cTime=datetime.datetime.now()
date='%i_%i_%i_%i_%i_%i'%(cTime.year,cTime.month,cTime.day,cTime.hour,cTime.minute,cTime.second)

print 'Making tar:'
if os.path.exists(tarfile): print '*********** tar already exists! I ASSUME YOU WANT TO MAKE A NEW ONE! *************'

os.chdir(relbase)
print 'tar --exclude="src/.git" --exclude="tmp/" --exclude="src/VLQRDF" --exclude="src/vlq-singleProd-RDF" --exclude=".SCRAM" -zcf '+tarfile+' ./*'
os.system('tar --exclude="src/.git" --exclude="tmp/" --exclude="src/VLQRDF" --exclude="src/vlq-singleProd-RDF" --exclude=".SCRAM" -zcf '+tarfile+' ./*')
os.chdir(runDir)

print 'Starting submission'
count=0

if makelists:
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > WJetsNanoList.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > ttbarNanoList.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > ttjetsTNanoList.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > ttjetsTbNanoList.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > singleTbNanoList.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > singleTNanoList.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > WJets200NanoList.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > WJets400NanoList.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > WJets600NanoList.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > WJets800NanoList.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > WJets1200NanoList.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = /WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM" > WJets2500NanoList.txt')

else:
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





