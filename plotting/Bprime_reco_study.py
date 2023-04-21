from ROOT import TFile, TTree, gStyle, TH1F, TCanvas

tfile = TFile.Open("root://cmseos.fnal.gov//store/user/sxiaohe/vlq-BtoTW-RDF/cut_update1_170/Bprime_hadds_alt/Bprime800_hadd.root")
ftree = tfile.Get("Events")

countT = 0
countW = 0
nEntries = ftree.GetEntries()
for i in range(nEntries):
    if(ftree.GetEntry(i)!=0):
        genFatJet_matching = ftree.genFatJet_matching
        n_gcFatJet = len(genFatJet_matching)
        if(genFatJet_matching[0]==6): countT+=1
        elif(genFatJet_matching[0]==24): countW+=1
        
        if(n_gcFatJet>1):
            if(genFatJet_matching[1]==6): countT+=1
            elif(genFatJet_matching[1]==24): countW+=1

counts = countT+countW
print(countT, countW, nEntries)
print("We get the correct jet {:.2f}% of the time.".format(100*counts/nEntries))
print("Among the correct jets, {:.2f}% are t jets.".format(100*countT/counts))
print("                      , {:.2f}% are W jets.".format(100*countW/counts))



        

        
