import os,sys,subprocess

#dogrep = True

logdir = '/uscms/home/jmanagan/nobackup/BtoTW/rdfjobs_Jan2024_pNetEffs/'

reruns = [
    'Bprime_M1200_2017UL/Bprime_M1200_2017UL.out',
    'DYMHT4002017UL/DYMHT4002017UL.out',
    'DYMHT8002018UL/DYMHT8002018UL.out',
    'QCDHT10002018UL/QCDHT10002018UL.out',
    'QCDHT3002016UL/QCDHT3002016UL.out',
    'QCDHT3002017UL/QCDHT3002017UL.out',
    'QCDHT3002018UL/QCDHT3002018UL.out',
    'QCDHT7002017UL/QCDHT7002017UL.out',
    'STs2017UL/STs2017UL.out',
    'STt2018UL/STt2018UL.out',
    'STtWb2017UL/STtWb2017UL.out',
    'STtWb2018UL/STtWb2018UL.out',
    'STtb2016UL/STtb2016UL.out',
    'TTHB2016UL/TTHB2016UL.out',
    'TTHB2018UL/TTHB2018UL.out',
    'TTTo2L2Nu2018UL/TTTo2L2Nu2018UL.out',
    'TTToHadronic2016APVUL/TTToHadronic2016APVUL.out',
    'TTToHadronic2018UL/TTToHadronic2018UL.out',
    'TTToSemiLeptonic2016APVUL/TTToSemiLeptonic2016APVUL.out',
    'TTToSemiLeptonic2016UL/TTToSemiLeptonic2016UL.out',
    'WJetsHT2002017UL/WJetsHT2002017UL.out',
    'WJetsHT4002017UL/WJetsHT4002017UL.out',
    'WJetsHT6002016APVUL/WJetsHT6002016APVUL.out',
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
            
