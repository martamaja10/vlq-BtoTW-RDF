from ROOT import TFile, TTree, gStyle, TH1F, TCanvas
from ROOT.VecOps import DeltaR
import numpy as np

tfile = TFile.Open("root://cmseos.fnal.gov//store/user/sxiaohe/vlq-BtoTW-RDF/cut_update1_170/Bprime_hadds_alt/Bprime2000_hadd.root")
ftree = tfile.Get("Events")

countT = 0
countW = 0
nEntries = ftree.GetEntries()
for i in range(nEntries):
    if(ftree.GetEntry(i)!=0):
        if(ftree.isValidBDecay==1 and ftree.Bprime_gen_exist==1):
            if(ftree.leptonicParticle==0):
                dR = DeltaR(ftree.W_eta, ftree.W_gen_eta, ftree.W_phi, ftree.W_gen_phi)
                if(dR < 0.4):
                #dR = ftree.W_lv.DeltaR(ftree.BPrime_lv-ftree.W_lv)
                    t_lv = ftree.BPrime_lv-ftree.W_lv
                    dR = DeltaR(t_lv.Eta(), ftree.t_gen_eta, t_lv.Phi(), ftree.t_gen_phi)
                    if(dR < 0.4):
                        countW+=1

            elif(ftree.leptonicParticle==1):
                dR = DeltaR(ftree.t_eta, ftree.t_gen_eta, ftree.t_phi, ftree.t_gen_phi)
                if(dR < 0.4):
                    W_lv = ftree.BPrime_lv-ftree.t_lv
                    dR = DeltaR(W_lv.Eta(), ftree.W_gen_eta, W_lv.Phi(), ftree.W_gen_phi)
                    if(dR < 0.4):
                        countT+=1
counts = countW+countT
print(countW, countT, counts)
#print("We get ", , "unambiguous events.")

W_mass, t_mass, dR_Wt, pt_bal = np.zeros(counts), np.zeros(counts), np.zeros(counts), np.zeros(counts)

idx=0
for i in range(nEntries):
    if(ftree.GetEntry(i)!=0):
        if(ftree.isValidBDecay==1 and ftree.Bprime_gen_exist==1):
            if(ftree.leptonicParticle==0):
                if(DeltaR(ftree.W_eta, ftree.W_gen_eta, ftree.W_phi, ftree.W_gen_phi)<0.4):
                    W_lv = ftree.W_lv
                    t_lv = ftree.BPrime_lv-W_lv
                    if(DeltaR(t_lv.Eta(), ftree.t_gen_eta, t_lv.Phi(), ftree.t_gen_phi)<0.4):
                        W_mass[idx] = W_lv.M()
                        t_mass[idx] = t_lv.M()
                        dR_Wt[idx] = W_lv.DeltaR(t_lv)
                        pt_bal[idx] = W_lv.Pt()/t_lv.Pt()
                        idx+=1
            elif(ftree.leptonicParticle==1):
                if(DeltaR(ftree.t_eta, ftree.t_gen_eta, ftree.t_phi, ftree.t_gen_phi)<0.4):
                    t_lv = ftree.t_lv    
                    W_lv = ftree.BPrime_lv-t_lv
                    if(DeltaR(W_lv.Eta(), ftree.W_gen_eta, W_lv.Phi(), ftree.W_gen_phi)<0.4):
                        t_mass[idx] = t_lv.M()
                        W_mass[idx] = W_lv.M()
                        dR_Wt[idx] =W_lv.DeltaR(t_lv)
                        pt_bal[idx] = W_lv.Pt()/t_lv.Pt()
                        idx+=1

print("t_mass_avg: {:.2f} ".format(np.mean(t_mass)))
print("t_mass_std: {:.2f} ".format(np.std(t_mass)))
print("W_mass_avg: {:.2f} ".format(np.mean(W_mass)))
print("W_mass_std: {:.2f} ".format(np.std(W_mass)))
print("dR_Wt_avg: {:.2f} ".format(np.mean(dR_Wt)))
print("dR_Wt_std: {:.2f} ".format(np.std(dR_Wt)))
print("pt_bal_avg: {:.2f} ".format(np.mean(pt_bal)))
print("pt_bal_std: {:.2f} ".format(np.std(pt_bal)))
        
