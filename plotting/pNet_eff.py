import os, time
from ROOT import gStyle, TH2F, TFile, TTree, TCanvas

start = time.time()

indir = "/store/user/sxiaohe/vlq-BtoTW-RDF/cut_update1_170/Bprime_hadds/"
outdir = os.getcwd()+'/plots_confusionMatrix/'
if not os.path.exists(outdir): os.system('mkdir -p '+outdir)
samples = {'Bp800':'Bprime800_hadd.root', 
           'Bp1400':'Bprime1400_hadd.root',
           'Bp2000':'Bprime2000_hadd.root',
           #'qcd200':'QCD200_hadd.root',                                  
           #'qcd300':'QCD300_hadd.root',                                  
           #'qcd500':'QCD500_hadd.root',                                  
           #'qcd700':'QCD700_hadd.root',                                  
           #'qcd1000':'QCD1000_hadd.root',                                
           #'qcd1500':'QCD1500_hadd.root',                                
           #'qcd2000':'QCD2000_hadd.root',
           #'ttbar':'ttbar_hadd.root',
           #'wjets200':'WJets200_hadd.root',
           #'wjets400':'WJets400_hadd.root',
           #'wjets600':'WJets600_hadd.root',
           #'wjets800':'WJets800_hadd.root',
           #'wjets1200':'WJets1200_hadd.root',
           #'wjets2500':'WJets2500_hadd.root',
           #'singleT':'singleT_hadd.root',
           #'singleTb':'singleTb_hadd.root',
           #'data_obs':'ttbar_hadd.root' # data is a copy of ttbar for now, need histograms in the file for limits
}

c1 = TCanvas("c1", "c1", 600, 600)
h = TH2F("h", "Confusion matrix of pNet_tag", 3, 0, 3, 3, 0, 3)
gStyle.SetOptStat(11);

for sample in samples.keys():
    tfile = TFile.Open("root://cmseos.fnal.gov/"+indir+"/"+samples[sample], "READ")
    ftree = tfile.Get("Events")

    nEntries = ftree.GetEntries()
    for i in range(nEntries):
        if (ftree.GetEntry(i) > 0):
            genFlv = ftree.genFatJet_matching
            pNetTag = ftree.pNet_tag # 0=J, 1=t, 2=W
        
            njets = len(genFlv)

            for j in range(njets):
                trueTag = -1
                if(genFlv[j]==0 or genFlv[j]==4 or genFlv[j]==5):
                    trueTag = 0 # light QCD jet
                elif(genFlv[j]==6):
                    trueTag = 1 # top jet
                elif(genFlv[j]==24): trueTag = 2 # W jet
            
                if(trueTag!=-1):
                    h.Fill(trueTag, pNetTag[j])

h.Draw("colz")
outputName = outdir+"confusionMatrix"
c1.SaveAs(outputName+".png")
c1.SaveAs(outputName+".root")
end = time.time()

print("time elapsed: ", end-start)
raw_input("Press Enter to continue...")
