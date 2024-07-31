import os,sys,subprocess,fnmatch

#dogrep = True
resubmit = True
findfails = True
jobsrunning = False
splitjobs = True
dumpbadfile = False
nsplit = 3

logdir = '/uscms/home/jmanagan/nobackup/BtoTW/rdfjobs_Jul2024_fullRun2/'

def findfiles(path, filtre):
    for root, dirs, files in os.walk(path):
        for f in fnmatch.filter(files, filtre):
            yield os.path.join(root, f)

reruns = [
]
# last was 7/30 00:47 TTsemi


if dumpbadfile:
    for rerun in reruns:
        print('RERUN: '+rerun)
        os.chdir(logdir)
        os.system("grep 'Number of Entries' "+rerun)
        rerun = rerun.replace('.out','.err')
        os.system("grep '/store/' "+rerun)
        os.system("grep 'status code' "+rerun)
    
## fill reruns via finding jobs that have no .out', (time, memory)
if findfails:
    for file in findfiles('/uscms_data/d3/jmanagan/BtoTW/rdfjobs_Jul2024_fullRun2/', '*.job'):
        outname = file.replace('.job','.out').replace('F_','F_RDF_').replace('UL_','UL_RDF_').replace('A_','A_RDF_').replace('B_','B_RDF_').replace('C_','C_RDF_').replace('D_','D_RDF_').replace('E_','E_RDF_').replace('G_','G_RDF_').replace('H_','H_RDF_')
        if not os.path.exists(outname):
            if not jobsrunning: reruns.append(outname[55:])
            else: os.system('ls -l '+file) # use time stamp to remove new running jobs

print(reruns[0:2])
resubs = {}

if not resubmit: exit()
for rerun in reruns:
    os.chdir(logdir)
    rerun = rerun.replace('_RDF','').replace('.out','.job')
    command = 'grep "Arguments" '+rerun
    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    out = out.decode('utf-8')
    jobnums = [int(out.split(" ")[4]),int(out.split(" ")[5])]
    
    print('\''+rerun+'\':'+str(jobnums))
    resubs[rerun] = jobnums
    #os.system('ls -l '+rerun.replace('_','_RDF_').replace('.job','.out'))
    #os.system('tail -20 '+rerun.replace('_','_RDF_').replace('.job','.err'))

for rerun in resubs.keys():
    folder = rerun.split('/')[0]
    jobfile = rerun.split('/')[1]

    if splitjobs:
        if resubs[rerun][0] == resubs[rerun][1]:
            print('NOT SPLITTING THIS ONE, SINGLE FILE '+folder+'/'+jobfile)
            os.chdir(logdir+folder)
            os.system('condor_submit '+jobfile)
            continue

        if nsplit == 2: 
            #half = int(resubs[rerun][0] + (resubs[rerun][1] - resubs[rerun][0])/2.0)
            newjobfile = jobfile.replace('_'+str(resubs[rerun][0]),'_'+str(half+1))
            os.chdir(logdir+folder)
            os.system('cp '+jobfile+' '+newjobfile)
            os.system('sed -i "s| '+str(resubs[rerun][1])+' | '+str(half)+' |" '+jobfile) # only change upper in original file
            os.system('sed -i "s| '+str(resubs[rerun][0])+' | '+str(half+1)+' |" '+newjobfile) # only change lower in new file
            os.system('sed -i "s|_'+str(resubs[rerun][0])+'\.|_'+str(half+1)+'\.|" '+newjobfile)
            
            os.system('condor_submit '+jobfile)
            os.system('condor_submit '+newjobfile) 

        if nsplit == 3:
            third = int(resubs[rerun][0] + (resubs[rerun][1] - resubs[rerun][0])/3.0)
            third2 = int(resubs[rerun][0] + 2*(resubs[rerun][1] - resubs[rerun][0])/3.0)
            print(folder+': new ranges '+str(resubs[rerun][0])+' - '+str(third)+' , '+str(third+1)+' - '+str(third2)+', and '+str(third2+1)+' - '+str(resubs[rerun][1]))
    
            newjobfile = jobfile.replace('_'+str(resubs[rerun][0]),'_'+str(third+1))
            newjobfile2 = jobfile.replace('_'+str(resubs[rerun][0]),'_'+str(third2+1))
        
            os.chdir(logdir+folder)
            os.system('cp '+jobfile+' '+newjobfile)
            os.system('cp '+jobfile+' '+newjobfile2)
            os.system('sed -i "s| '+str(resubs[rerun][1])+' | '+str(third)+' |" '+jobfile) # only change upper in original file

            os.system('sed -i "s| '+str(resubs[rerun][1])+' | '+str(third2)+' |" '+newjobfile) # change upper in middle file
            os.system('sed -i "s| '+str(resubs[rerun][0])+' | '+str(third+1)+' |" '+newjobfile) # change lower in middle file
            os.system('sed -i "s|_'+str(resubs[rerun][0])+'\.|_'+str(third+1)+'\.|" '+newjobfile)

            os.system('sed -i "s| '+str(resubs[rerun][0])+' | '+str(third2+1)+' |" '+newjobfile2) # only change lower in last file
            os.system('sed -i "s|_'+str(resubs[rerun][0])+'\.|_'+str(third2+1)+'\.|" '+newjobfile2)
            
            os.system('condor_submit '+jobfile)
            os.system('condor_submit '+newjobfile)
            os.system('condor_submit '+newjobfile2) 
    else:
        os.chdir(logdir+folder)
        os.system('condor_submit '+jobfile)
            
