import os
import numpy as np
import ROOT
from ROOT import *
import matplotlib.pyplot as plt

indir = "root://cmseos.fnal.gov//store/user/kjohnso/BtoTW_Jul2023/LeptonChecks/QCDBp_scenarios/"
outdir = os.getcwd()+'/LepIso_histos/'
write = ""
if not os.path.exists(outdir): os.system('mkdir -p '+outdir)

samples = {"Bprime800":"Bprime800_scenarios.root",
           "Bprime1400":"Bprime1400_scenarios.root",
           "Bprime2000":"Bprime2000_scenarios.root",                    
           "QCD200":"QCD200_scenarios.root",   
           "QCD300":"QCD300_scenarios.root",
           "QCD500":"QCD500_scenarios.root",
           "QCD700":"QCD700_scenarios.root",
           "QCD1000":"QCD1000_scenarios.root",
           "QCD1500":"QCD1500_scenarios.root",
           "QCD2000":"QCD2000_scenarios.root"
}

lumi = 138000.0
nRun = {'Bprime800':99800.,
        'Bprime1400':99600.,
        'Bprime2000':99600.,  # fixme, need a sum-of-gen-weights calculated...
        'QCD200':61542214.,
        'QCD300':56214199.,
        'QCD500':61097673.,
        'QCD700':47314826.,
        'QCD1000':15230975.,
        'QCD1500':11887406.,
        'QCD2000':5710430.
}
xsec = {'Bprime800':1.0,'Bprime1400':1.0,'Bprime2000':1.0, #signals all the same at 1pb for now, predictions vary 
        'QCD200':1712000,
        'QCD300':347700,
        'QCD500':32100,
        'QCD700':6831,
        'QCD1000':1207,
        'QCD1500':119.9,
        'QCD2000':25.24
}

weights = {'Bprime800': 0.,
           'Bprime1400': 0.,
           'Bprime2000': 0.,
           'QCD200': 0.,
           'QCD300': 0.,
           'QCD500': 0.,
           'QCD700': 0.,
           'QCD1000': 0.,
           'QCD1500': 0.,
           'QCD2000': 0.
}

dphi ='''
using namespace ROOT::VecOps;
float dphi(const float& phi1, const float& phi2) {
    return DeltaPhi(phi1, phi2);
}
'''
ROOT.gInterpreter.Declare(dphi)

for sample in samples:
    filename = "root://cmseos.fnal.gov//store/user/kjohnso/BtoTW_Jul2023/LeptonChecks/QCDBp_short/"+sample+"_MET.root"
    tfile = ROOT.TFile.Open(filename)   
    ftree = tfile.Get("Events")
    ftree.GetEntry(0)
    Generator_weight = ftree.Generator_weight
    print("Generator_weight: ", Generator_weight)
    tfile.Close()

    filename = indir+samples[sample]
    LepIsoC = ROOT.RDataFrame("Events", filename)
    ElIsoC = LepIsoC.Filter("isEl").Define("MET_Lep_DeltaPhi","dphi(LepIsoC_phi, MET_phi)")
    #MuIsoC = LepIsoC.Filter("isMu").Define("MET_Lep_DeltaPhi","dphi(LepIsoC_phi, MET_phi)")
   
    # MET_phi = ftree.MET_phi
    # El_phi = ftree.lep_phi
    # El_DelPhi = VecOps.DeltaPhi(El_phi, MET_phi)
    
    outname = sample+"_triangular_El.root"
    ElIsoC.Snapshot("Events", outname)
    outfile = ROOT.TFile.Open(outname,"UPDATE")
    
    histo = ElIsoC.Histo2D(("PhiPt_histo","PhiPt_histo", 25,0,3.15,25,50,400),"MET_Lep_DeltaPhi","MET_pt")
    Generator_weight = 1.0
    weight = Generator_weight*lumi*xsec[sample]/(nRun[sample]*abs(Generator_weight))
    weights[sample] = weight
    
    histo.Draw("colz")
    histo.Scale(weight)
    histo.Write("PhiPt_histo")
    outfile.Close()

    
    # outname = sample+"_triangular_Mu.root"
    # MuIsoC.Snapshot("Events", outname)
    # outfile = ROOT.TFile.Open(outname,"UPDATE")

    # histo = MuIsoC.Histo2D(("PhiPt_histo","PhiPt_histo", 25,0,3.15,25,50,400),"MET_Lep_DeltaPhi","MET_pt")
    # Generator_weight = 1.0
    # weight = Generator_weight*lumi*xsec[sample]/(nRun[sample]*abs(Generator_weight))
    # weights[sample] = weight

    # histo.Draw("colz")
    # histo.Scale(weight)
    # histo.Write("PhiPt_histo")
    # outfile.Close()
    
    # canv1 = TCanvas("c1","c1",800,600)
    # histo.Draw("colz")
    # canv1.SaveAs("PhiPt_histo.png")
    

N = 400
slopeList = np.arange(N)
counts_Bp800 = np.zeros(N)
counts_Bp1400 = np.zeros(N)
counts_Bp2000 = np.zeros(N)
counts_QCD200 = np.zeros(N)
counts_QCD300 = np.zeros(N)
counts_QCD500 = np.zeros(N)
counts_QCD700 = np.zeros(N)
counts_QCD1000 = np.zeros(N)
counts_QCD1500 = np.zeros(N)
counts_QCD2000 = np.zeros(N)


noPass_Bp800 = np.zeros(N)
noPass_Bp1400 = np.zeros(N)
noPass_Bp2000 = np.zeros(N)
noPass_QCD200 = np.zeros(N)
noPass_QCD300 = np.zeros(N)
noPass_QCD500 = np.zeros(N)
noPass_QCD700 = np.zeros(N)
noPass_QCD1000 = np.zeros(N)
noPass_QCD1500 = np.zeros(N)
noPass_QCD2000 = np.zeros(N)


for sample in samples:
    print(sample)
    outname = sample+"_triangular.root"
    tfile = TFile.Open(outname)
    ftree = tfile.Get("Events") #ElIsoC_
    
    nEntries = ftree.GetEntries()
 
    for i in range(nEntries):
        if(ftree.GetEntry(i)>0):
            
            MET_phi = ftree.MET_phi
            El_phi = ftree.LepIsoC_phi
        
            MET_pt = ftree.MET_pt

            El_DelPhi = VecOps.DeltaPhi(El_phi, MET_phi)

            MET_ptThreshold = (slopeList/1.5) * El_DelPhi - slopeList # array

            for j in range(N):
                if(MET_pt>MET_ptThreshold[j]):
                    if(sample=="Bprime800"):
                        counts_Bp800[j]+=1
                    elif(sample=="Bprime1400"):
                        counts_Bp1400[j]+=1
                    elif(sample=="Bprime2000"):
                        counts_Bp2000[j]+=1
                    elif(sample=="QCD200"):
                        counts_QCD200[j]+=1
                    elif(sample=="QCD300"):
                        counts_QCD300[j]+=1
                    elif(sample=="QCD500"):
                        counts_QCD500[j]+=1
                    elif(sample=="QCD700"):
                        counts_QCD700[j]+=1
                    elif(sample=="QCD1000"):
                        counts_QCD1000[j]+=1
                    elif(sample=="QCD1500"):
                        counts_QCD1500[j]+=1
                    elif(sample=="QCD2000"):
                        counts_QCD2000[j]+=1
                else:
                    if(sample=="Bprime800"):
                        noPass_Bp800[j]+=1
                    elif(sample=="Bprime1400"):
                        noPass_Bp1400[j]+=1
                    elif(sample=="Bprime2000"):
                        noPass_Bp2000[j]+=1
                    elif(sample=="QCD200"):
                        noPass_QCD200[j]+=1
                    elif(sample=="QCD300"):
                        noPass_QCD300[j]+=1
                    elif(sample=="QCD500"):
                        noPass_QCD500[j]+=1
                    elif(sample=="QCD700"):
                        noPass_QCD700[j]+=1
                    elif(sample=="QCD1000"):
                        noPass_QCD1000[j]+=1
                    elif(sample=="QCD1500"):
                        noPass_QCD1500[j]+=1
                    elif(sample=="QCD2000"):
                        noPass_QCD2000[j]+=1                    

# true is significance and false is efficiency
SorE = True
percentIndex = .05
counts_bkg = counts_QCD200*weights["QCD200"] + counts_QCD300*weights["QCD300"] + counts_QCD500*weights["QCD500"] + counts_QCD700*weights["QCD700"] + counts_QCD1000*weights["QCD1000"] + counts_QCD1500*weights["QCD1500"] + counts_QCD2000*weights["QCD2000"]

if (SorE):
    print("Calculating significance...")
    
    significance_Bp800 = (counts_Bp800*weights["Bprime800"])/(counts_bkg**0.5)
    significance_Bp1400 = counts_Bp1400*weights["Bprime1400"]/(counts_bkg**0.5)
    significance_Bp2000 = counts_Bp2000*weights["Bprime2000"]/(counts_bkg**0.5)
    
    # Index of 800
    index800 = max(enumerate(significance_Bp800),key=lambda x: x[1])[0]
    write = "Max Index Bp800: " + str(index800) + "\n"
    #print("Max Index Bp800: " + str(index800))
    
    value800 = significance_Bp800[index800]
    value90P800 = value800 - percentIndex * (value800 - significance_Bp800[0])
    index90P800 = min(range(len(significance_Bp800)/2), key=lambda i: abs(significance_Bp800[i]-value90P800))
    write += "Max Index 95% Bp800: " + str(index90P800) + "\n"
    #print("Max Index 95% Bp800: " + str(index90P800))
    
    # Index of 1400
    index1400 = max(enumerate(significance_Bp1400),key=lambda x: x[1])[0]
    write += "Max Index Bp1400: " + str(index1400) + "\n"
    #print("Max Index Bp1400: ", index1400)
    
    value1400 = significance_Bp1400[index1400]
    value90P1400 = value1400 - percentIndex * (value1400 - significance_Bp1400[0])
    index90P1400 = min(range(len(significance_Bp1400)/2), key=lambda i: abs(significance_Bp1400[i]-value90P1400))
    #print("Max Index 95% Bp1400: ", index90P1400)
    write += "Max Index 95% Bp1400: " + str(index90P1400) + "\n"
    
    # Index of 2000
    index2000 = max(enumerate(significance_Bp2000),key=lambda x: x[1])[0]
    #print("Max Index Bp2000: ", index2000)
    write += "Max Index Bp2000: " + str(index2000) + "\n"
    
    value2000 = significance_Bp2000[index2000]
    value90P2000 = value2000 - percentIndex * (value2000 - significance_Bp2000[0])
    index90P2000 = min(range(len(significance_Bp2000)), key=lambda i: abs(significance_Bp2000[i]-value90P2000))
    #print("Max Index 95% Bp2000: ", index90P2000)
    write += "Max Index 95% Bp2000: " + str(index90P2000) + "\n"
    
    print("Making plots...")

    plt.plot(slopeList, significance_Bp800, marker='o', label="Bprime800", markersize=0.5)
    plt.plot(slopeList, significance_Bp1400, marker='o', label="Bprime1400", markersize=0.5)
    plt.plot(slopeList, significance_Bp2000, marker='o', label="Bprime2000", markersize=0.5)
    
    plt.ylabel("significance")
    plt.xlabel("slope")
    plt.legend()
    #plt.show()
    plt.savefig('significance.png')

if (SorE):
    print("Calculating efficiency...")
    
    noPass_bkg = noPass_QCD200*weights["QCD200"] + noPass_QCD300*weights["QCD300"] + noPass_QCD500*weights["QCD500"] + noPass_QCD700*weights["QCD700"] + noPass_QCD1000*weights["QCD1000"] + noPass_QCD1500*weights["QCD1500"] + noPass_QCD2000*weights["QCD2000"]      
    eff_Bp800 = counts_Bp800 / (counts_Bp800 + noPass_Bp800)
    eff_Bp1400 = counts_Bp1400 / (counts_Bp1400 + noPass_Bp1400)
    eff_Bp2000 = counts_Bp2000 / (counts_Bp2000 + noPass_Bp2000)
    eff_bkg = counts_bkg / (counts_bkg + noPass_bkg)
    
    # Efficient MET choice -> text
    # print("eff_Bp800: ", eff_Bp800[130])
    # print("eff_Bp1400: ", eff_Bp1400[130])
    # print("eff_Bp2000: ", eff_Bp2000[130])
    # print("eff_bkg: ", eff_bkg[130])
    
    write += "\n130: \n"
    write += "eff_Bp800: " + str(eff_Bp800[130]) + "\n"
    write += "eff_Bp1400: " + str(eff_Bp1400[130]) + "\n"
    write += "eff_Bp2000: " + str(eff_Bp2000[130]) + "\n"
    write += "eff_bkg: " + str(eff_bkg[130]) + "\n"
    
    write += "\n145: \n"
    write += "eff_Bp800: " + str(eff_Bp800[145]) + "\n"
    write += "eff_Bp1400: " + str(eff_Bp1400[145]) + "\n"
    write += "eff_Bp2000: " + str(eff_Bp2000[145]) + "\n"
    write += "eff_bkg: " + str(eff_bkg[145]) + "\n"
    
    print("Making plots...")

    # plt.plot(slopeList, eff_Bp800, marker='o', label="Bprime800", markersize=0.5)
    # plt.plot(slopeList, eff_Bp1400, marker='o', label="Bprime1400", markersize=0.5)                
    # plt.plot(slopeList, eff_Bp2000, marker='o', label="Bprime2000", markersize=0.5)
    # plt.plot(slopeList, eff_bkg, marker='o', label="Background", markersize=0.5)

    #plt.ylabel("efficiency")
    #plt.xlabel("slope")
    #plt.legend()
    #plt.show()
    #plt.savefig('efficiency.png')'
    
with open('EfficiencyMETchoice.txt', 'w') as f:
    f.write(write)