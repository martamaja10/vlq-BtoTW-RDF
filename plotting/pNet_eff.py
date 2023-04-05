import time
from ROOT import gStyle, TH2F, TFile, TTree, TCanvas

start = time.time()

fFile = TFile.Open("root://cmseos.fnal.gov//store/user/sxiaohe/vlq-BtoTW-RDF/cut_update1_170/Bprime_hadds/Bprime800_hadd.root", "READ")
fTree = fFile.Get("Events")

fCanvas = TCanvas("c", "c", 600, 600)
h = TH2F("h", "Confusion matrix of pNet_tag", 3, 0, 3, 3, 0, 3)
gStyle.SetOptStat(0);

n = fTree.GetEntries()
for i in range(100):
    if fTree.GetEntry(i) > 0:
        genFlv = fTree.genFatJet_matching
        pNetTag = fTree.pNet_tag # 0=J, 1=t, 2=W
        
        njets = len(genFlv)

        for j in range(njets):
            trueTag = -1
            if(genFlv[j]==0 or genFlv[j]==4 or genFlv[j]==5): trueTag = 0 # light QCD jet
            elif(genFlv[j]==6): trueTag = 1 # top jet
            elif(genFlv[j]==24): trueTag = 2 # W jet
            #print(genFlv[j], trueTag)
            
            if(trueTag!=-1):
                print(trueTag, pNetTag[j])
                h.Fill(trueTag, pNetTag[j])
h.Draw("colz")

end = time.time()

print("time elapsed: ", end-start)
raw_input("Press Enter to continue...")
