import os,sys,datetime,time, subprocess, math
from ROOT import *
execfile("/uscms_data/d3/khowey/EOSSafeUtils.py")

start_time = time.time()

#IO directories must be full paths

inDir='/store/user/khowey/BtoTW_RDF_MLPs/'
outDir='/store/user/khowey/BtoTW_RDF_MLPs_hadds/'
scratchDir='/uscmst1b_scratch/lpc1/3DayLifetime/khowey/'

if not os.path.exists(scratchDir): os.system('mkdir '+scratchDir)
os.system('eos root://cmseos.fnal.gov/ mkdir -p '+outDir)

dirList = [
#    'BpM1400',
    'BpM2000',
    'BpM800',
#    'WJets1200',
#    'WJets200',
#    'WJets2500',
#    'WJets400',
#    'WJets600',
#    'WJets800',
#    'singleT',
#    'singleTb',
    'ttbarInc',
]

for sample in dirList:
    outList = ['none']
 
    for outlabel in outList:

        outsample = sample+'_'+outlabel
        if outlabel == 'none': outsample = sample

        rootfiles = EOSlist_root_files(inDir+'/'+outsample)

        print ("------------ hadding Sample:",outsample,"---------------")
        print ('N root files in',outsample,'=',len(rootfiles))


        nFilesPerHadd = 900

        onefile = ' root://cmseos.fnal.gov/'+inDir+'/'+outsample+'/'+rootfiles[-1]
        manyfiles = nFilesPerHadd*onefile
        lengthcheck = len('hadd -f root://cmseos.fnal.gov/'+outDir+'/'+outsample+'_hadd.root '+manyfiles)
        if lengthcheck > 126000.:
            toolong = lengthcheck - 126000.
            num2remove = math.ceil(toolong/len(onefile))
            nFilesPerHadd = int(nFilesPerHadd - num2remove)
            print ('Length estimate reduced from',lengthcheck,'by',toolong,'via removing',num2remove,'files for nFilesPerHadd of',nFilesPerHadd)

        if len(rootfiles) < nFilesPerHadd:
            haddcommand = 'hadd -f '+scratchDir+'/'+outsample+'_hadd.root '
            for file in rootfiles:
                haddcommand+=' root://cmseos.fnal.gov/'+inDir+'/'+outsample+'/'+file
            print ('Length of hadd command =',len(haddcommand))
            subprocess.call(haddcommand,shell=True)
            
            xrdcpcommand = 'xrdcp -f '+scratchDir+'/'+outsample+'_hadd.root root://cmseos.fnal.gov/'+outDir+'/'+outsample+'_hadd.root'
            subprocess.call(xrdcpcommand,shell=True)

            if bool(EOSisfile(outDir+'/'+outsample+'_hadd.root')) != True:
                print ('COPY ME LATER')        
        else:
            for i in range(int(math.ceil(len(rootfiles)/float(nFilesPerHadd)))):
                haddcommand = 'hadd -f '+scratchDir+'/'+outsample+'_'+str(i+1)+'_hadd.root '

                begin=i*nFilesPerHadd
                end=begin+nFilesPerHadd
                if end > len(rootfiles): end=len(rootfiles)
                print ('begin:',begin,'end:',end)

                for j in range(begin,end):
                    haddcommand+=' root://cmseos.fnal.gov/'+inDir+'/'+outsample+'/'+rootfiles[j]
                print ('Length of hadd command =',len(haddcommand))
                subprocess.call(haddcommand,shell=True)

                xrdcpcommand = 'xrdcp -f '+scratchDir+'/'+outsample+'_'+str(i+1)+'_hadd.root root://cmseos.fnal.gov/'+outDir+'/'+outsample+'_'+str(i+1)+'_hadd.root'
                subprocess.call(xrdcpcommand,shell=True)

                if bool(EOSisfile(outDir+'/'+outsample+'_'+str(i+1)+'_hadd.root')) != True:
                    print ('COPY ME LATER')

print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))



