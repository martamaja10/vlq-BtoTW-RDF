from ROOT import TFile, TMath
from ROOT.VecOps import DeltaR, Argsort, RVec

tfile = TFile.Open("root://cmseos.fnal.gov//store/user/sxiaohe/vlq-BtoTW-RDF/cut_update1_170/Bprime_hadds_alt/Bprime800_hadd.root")
ftree = tfile.Get("Events")

countT = 0
countW = 0
nEntries = ftree.GetEntries()

for i in range(nEntries):
    if(ftree.GetEntry(i)!=0):
        genFatJet_matching = ftree.genFatJet_matching
        n_gcFatJet = len(genFatJet_matching)

        lep_eta = ftree.lepton_lv.Eta()
        lep_phi = ftree.lepton_lv.Phi()
        jet_eta = ftree.gcFatJet_eta
        jet_phi = ftree.gcFatJet_phi

        DR_lep = RVec('float')(n_gcFatJet)
        for ijet in range(n_gcFatJet):
            DR_lep[ijet] = abs(DeltaR(lep_eta, jet_eta[ijet], lep_phi, jet_phi[ijet]) - TMath.Pi())

        DR_lep = Argsort(DR_lep)
        
        idx1 = DR_lep[0]
        if(genFatJet_matching[idx1]==6):
            countT+=1
        elif(genFatJet_matching[idx1]==24): 
            countW+=1

        if(n_gcFatJet>1):
            idx2 = DR_lep[1] 
            if(genFatJet_matching[idx2]==6):
                countT+=1
            elif(genFatJet_matching[idx2]==24):
                countW+=1

counts = countT+countW
print(countT, countW, nEntries)
print("We get the correct jet {:.2f}% of the time.".format(100*counts/nEntries))
print("Among the correct jets, {:.2f}% are t jets ".format(100*countT/counts))
print("                      , {:.2f}% are W jets.".format(100*countW/counts))



        

        
