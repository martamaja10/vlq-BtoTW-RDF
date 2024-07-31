import os,sys,subprocess

#dogrep = True

logdir = '/uscms/home/jmanagan/nobackup/BtoTW/rdfjobs_Apr2024_pNetEffs_tWseparate/'

reruns = [
'Bprime_M1200_2017UL/Bprime_M1200_2017UL.out',
'QCDHT3002016UL/QCDHT3002016UL.out',
'QCDHT5002016UL/QCDHT5002016UL.out',
'TTWq2016APVUL/TTWq2016APVUL.out',
'WJetsHT12002018UL/WJetsHT12002018UL.out',
]

resubs = {}

#if dogrep:
for rerun in reruns:
    os.chdir(logdir)
    rerun = rerun.replace('.out','.job')
    folder = rerun.split('/')[0]
    jobfile = rerun.split('/')[1]

    os.chdir(logdir+folder)
    os.system('condor_submit '+jobfile)
            
