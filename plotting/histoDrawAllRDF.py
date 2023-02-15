import os,sys,time
from numpy import linspace
from array import array
from ROOT import *
#from weights import *  ## later, get weights.py working for nRun and xSec

# give a branch name to plot
#iPlot = sys.argv[1] 

# open an output file for the histos
#outputFile = TFile.Open("histos_"+iPlot+".root","RECREATE")

# Samples to process and categories to plot
indir = "/store/user/jmanagan/BtoTW_RDF"
outdir = os.getcwd()+'/plots_forwardjetAll/'
if not os.path.exists(outdir): os.system('mkdir -p '+outdir)
samples = {'Bp800':'Bp800_hadd.root', #'Bp1400':'Bp1400_hadd.root',
           'Bp2000':'Bp2000_hadd.root',
           'ttbar':'ttbarInc_hadd.root','wjets':'WJetsInc_hadd.root',
           'singleT':'singleT_hadd.root','singleTb':'singleTb_hadd.root',
           'data_obs':'ttbarInc_hadd.root' # data is a copy of ttbar for now, need histograms in the file for limits
}
tags = {'tjet':'taggedTjet == true',
        'Wjet':'taggedWjet == true',
        'Wbjet':'taggedWbjetJet == true',
        'other':'isValidBDecay < 1', # later likely some other selection like a network score
        'all':'isValidBDecay == 1', # combine all tagged events
    }

# Numbers to weight the samples relative to each other
lumi = 138000.0
nRun = {'Bp800':200000,'Bp1400':200000,'Bp2000':200000,  # fixme, need a sum-of-gen-weights calculated...
        'ttbar':476408000,'wjets':81051269,
        'singleT':178336000.0,'singleTb':95627000.0,
        'data_obs':1}
xsec = {'Bp800':1.0,'Bp1400':1.0,'Bp2000':1.0, #signals all the same at 1pb for now, predictions vary
        'ttbar':831.76*0.438,'wjets':61526.7,
        'singleT':136.02,'singleTb':80.95,
        'data_obs':1}

# Settings for drawing the graphs
sig1leg='B (0.8 TeV)'
sig2leg='B (2.0 TeV)'
sig3leg='B (1.4 TeV)'
sigScaleFact = 10 # zoom in/out on signal
print 'Scale factor = ',sigScaleFact
bkgProcList = ['wjets','ttbar','singleT','singleTb']
bkgHistColors = {'ttbar':kAzure+8,'wjets':kMagenta-2,'singleT':kGreen-3,'singleTb':kGreen-3,} #TT
yLog  = True
blind = True
lumiSys = 0.20 # 20% uncertainty on background
doPlotting = True

# Settings for different iPlots
plotList = {#discriminantName:(discriminantLJMETName, binning, xAxisLabel)
    'lepPt' :('lepton_pt',linspace(0, 1000, 51).tolist(),';Lepton p_{T} [GeV];'),
    'lepEta':('lepton_eta',linspace(-4, 4, 41).tolist(),';Lepton #eta;'),
    'lepPhi':('lepton_phi',linspace(-3.2,3.2,65).tolist(),';#phi(l)'),
    'lepIso':('lepton_miniIso',linspace(0,0.2,51).tolist(),';lepton mini isolation'),

    'JetEta':('gcJet_eta',linspace(-4, 4, 41).tolist(),';AK4 jet #eta;'),
    'JetPt' :('gcJet_pt',linspace(0, 1500, 51).tolist(),';AK4 jet p_{T} [GeV];'),
    'NJetsCentral' :('NJets_central',linspace(0, 20, 21).tolist(),';central jet multiplicity;'),
    'NJetsForward' :('NJets_forward',linspace(0, 20, 21).tolist(),';forward jet multiplicity;'),
    'NBJets':('NJets_DeepFlavM',linspace(0, 10, 11).tolist(),';b tag multiplicity;'),
    'MET'   :('MET_pt',linspace(0, 1500, 51).tolist(),';#slash{E}_{T} [GeV];'),
    'HT':('Jet_HT',linspace(0, 5000, 51).tolist(),';H_{T} (GeV);'),
    'ST':('Jet_ST',linspace(0, 5000, 51).tolist(),';S_{T} (GeV);'),

    'FatJetEta':('gcFatJet_eta',linspace(-4, 4, 41).tolist(),';AK8 jet #eta;'),
    'FatJetPt' :('gcFatJet_pt',linspace(0, 1500, 51).tolist(),';AK8 jet p_{T} [GeV];'),
    'Tau21'  :('tau21',linspace(0, 1, 51).tolist(),';AK8 Jet #tau_{2}/#tau_{1};'),
    'SoftDrop' :('FatJet_sdMass',linspace(0, 500, 51).tolist(),';AK8 soft drop mass [GeV];'),
    'probj':('pNet_J',linspace(0,1,51).tolist(),';particleNet J score'),
    'probt':('pNet_T',linspace(0,1,51).tolist(),';particleNet t score'),
    'probw':('pNet_W',linspace(0,1,51).tolist(),';particleNet W score'),
    'deeptag':('pNet_tag',linspace(0,10,11).tolist(),';particleNet tag (0 = J, 1 = t, 2 = W)'),
    'NFatJets':('NFatJets',linspace(0, 10, 11).tolist(),';AK8 jet multiplicity;'),
    'nT':('nT_pNet',linspace(0,5,6).tolist(),';N particleNet t-tagged jets'),
    'nW':('nW_pNet',linspace(0,5,6).tolist(),';N particleNet W-tagged jets'),
    'minDR_twoAK8s':('minDR_leadAK8otherAK8',linspace(0,5,51).tolist(),';min #Delta R(leading AK8 jet, other AK8 jet) [GeV]'),

    'tmass':('t_mass',linspace(0,500,51).tolist(),';M(t) [GeV]'),
    'tpt':('t_pt',linspace(0,1000,51).tolist(),';p_{T}(t) [GeV]'),
    'Wdrlep':('DR_W_lep',linspace(0,5,51).tolist(),';leptonic W, #DeltaR(W,lepton)'),
    'tdrWb':('DR_W_b',linspace(0,6.3,51).tolist(),';leptonic t, #DeltaR(W,b)'),
    'isLepW':('leptonicParticle',linspace(0,2,3).tolist(),';lepton from W'),
    'minMlj':('minM_lep_Jet',linspace(0,1000,51).tolist(),';min[M(l,jet)] [GeV], 0 b tags'),
    'PtRel':('ptRel_lep_Jet',linspace(0,500,51).tolist(),';p_{T,rel}(l, closest jet) [GeV]'),
    'PtRelAK8':('ptRel_lep_FatJet',linspace(0,500,51).tolist(),';p_{T,rel}(l, closest AK8 jet) [GeV]'),
    'minDR':('minDR_lep_Jet',linspace(0,5,51).tolist(),';#Delta R(l, closest jet) [GeV]'),
    'minDRAK8':('minDR_lep_FatJet',linspace(0,5,51).tolist(),';#Delta R(l, closest AK8 jet) [GeV]'),

    'BpMass':('Bprime_mass',linspace(0,4000,51).tolist(),';M(B) [GeV]'),
    #'BpMass':('Bprime_mass',51,0,4000,';M(B) [GeV]'),
    'BpPt':('Bprime_pt',linspace(0,3000,51).tolist(),';B quark p_{T} [GeV]'),
    'BpEta':('Bprime_eta',linspace(-5,5,51).tolist(),';B quark #eta'),
    'BpPhi':('Bprime_phi',linspace(-3.14,3.14,51).tolist(),';B quark #phi'),
    'BpDeltaR':('Bprime_DR',linspace(0,5,51).tolist(),';#DeltaR(B quark product jets)'),
    'BpPtBal':('Bprime_ptbal',linspace(0,3,51).tolist(),';B quark t/W p_{T} ratio'),
    'BpChi2':('Bprime_chi2',linspace(0,1000,51).tolist(),';B quark reconstruction #chi^{2}'), # CHECK ME, what range?

    #'DnnTprime':('dnnAll_Tprime',linspace(0,1,51).tolist(),';DNN T score'), # later, these should exist
    #'DnnTTbar':('dnnAll_ttbar',linspace(0,1,51).tolist(),';DNN-T t#bar{t} score'),
    #'DnnWJets':('dnnAll_WJets',linspace(0,1,51).tolist(),';DNN-T W+jets score'),
}

histos = {}
start_time = time.time()
print 'Building histograms...'
for sample in samples.keys():
    print '\t Sample =',sample
    df = RDataFrame('Events',"root://cmseos.fnal.gov/"+indir+"/"+samples[sample])
    sigregion = df.Filter("NJets_forward > 0 && Bprime_mass > 0")\
                  .Define("weight","Generator_weight*{}*{}/({}*abs(Generator_weight))".format(lumi,xsec[sample],nRun[sample]))
    
    for tag in tags.keys():
        print '\t\t tag =',tag
        df_tag = sigregion.Filter(tags[tag])

        for iPlot in plotList.keys():
            print '\t\t\t iPlot =',iPlot

            histo = df_tag.Histo1D((iPlot+'_'+sample+'_'+tag, ';'+plotList[iPlot][-1]+';Events / bin',len(plotList[iPlot][1])-1,array('d',plotList[iPlot][1])),plotList[iPlot][0],"weight")
            #histo.Write()

            histos[iPlot+'_'+sample+'_'+tag] = histo

#outputFile.Close();
#print histos

hists = {}
if doPlotting:

    print 'Getting histograms from pointers...'

    ptr_time = time.time()
    for iPtr in histos.keys():
        hists[iPtr] = histos[iPtr].GetValue()

    #print hists

    plot_time = time.time()
    for tag in tags.keys():
        print '\t tag =',tag
        for iPlot in plotList.keys():

            histPrefix = iPlot+'_'+tag
        
            # add the backgrounds together
            bkgHT = hists[iPlot+'_ttbar_'+tag].Clone()
            for proc in bkgProcList:
                if proc == 'ttbar': continue
                bkgHT.Add(hists[iPlot+'_'+proc+'_'+tag])

            if bkgHT.Integral() < 1: 
                print 'SOMETHING WRONG HERE, no background!',tag,iPlot
                continue
        
            # set a larger uncertainty
            for ibin in range(1,bkgHT.GetNbinsX()+1):
                bkgHT.SetBinError(ibin,bkgHT.GetBinContent(ibin)*lumiSys)

        
            hsig800 = hists[iPlot+'_Bp800_'+tag]
            #hsig1400 = hists['Bp1400_'+tag].Clone()
            hsig2000 = hists[iPlot+'_Bp2000_'+tag]
            if hsig2000.Integral()/bkgHT.Integral() > 0.1: sigScaleFact = 1
            hsig800.Scale(sigScaleFact)
            #hsig1400.Scale(sigScaleFact)
            hsig2000.Scale(sigScaleFact)
        
            ############################################################
            ############## Making Plots of e+jets, mu+jets and e/mu+jets 
            ############################################################
		
            stackbkgHT = THStack("stackbkgHT","")
            # set each background a different color
            for proc in bkgProcList:
                stackbkgHT.Add(hists[iPlot+'_'+proc+'_'+tag])
                hists[iPlot+'_'+proc+'_'+tag].SetLineColor(bkgHistColors[proc])
                hists[iPlot+'_'+proc+'_'+tag].SetFillColor(bkgHistColors[proc])
                hists[iPlot+'_'+proc+'_'+tag].SetLineWidth(2)

            sig1Color= kBlack
            sig2Color= kBlack
            sig3Color= kBlack
			
            hsig800.SetLineColor(sig1Color)
            hsig800.SetFillStyle(0)
            hsig800.SetLineWidth(3)
            hsig2000.SetLineColor(sig2Color)
            hsig2000.SetLineStyle(7)#5)
            hsig2000.SetFillStyle(0)
            hsig2000.SetLineWidth(3)
    
            bkgHT.SetFillStyle(3004)
            bkgHT.SetFillColor(kBlack)
    
            gStyle.SetOptStat(0)
            c1 = TCanvas("c1","c1",1200,1000)
    
            stackbkgHT.Draw("HIST")
            hsig800.Draw("SAME HIST")
            hsig2000.Draw("SAME HIST")
            #hsig3.Draw("SAME HIST")
            bkgHT.Draw("SAME E2")
    
            leg = TLegend(0.45,0.64,0.95,0.89)
            leg.SetShadowColor(0)
            leg.SetFillColor(0)
            leg.SetFillStyle(0)
            leg.SetLineColor(0)
            leg.SetLineStyle(0)
            leg.SetBorderSize(0) 
            leg.SetNColumns(2)
            leg.SetTextFont(62)#42)
            scaleFact1Str = ' x'+str(sigScaleFact)
            scaleFact2Str = ' x'+str(sigScaleFact)
    
            leg.AddEntry(hsig800,sig1leg+scaleFact1Str,"l") #left
            leg.AddEntry(hists[iPlot+'_ttbar_'+tag],"t#bar{t}","f") #right
            leg.AddEntry(hsig2000,sig2leg+scaleFact2Str,"l") #left
            leg.AddEntry(hists[iPlot+'_wjets_'+tag],"W+jets","f") #right
            leg.AddEntry(hists[iPlot+'_singleT_'+tag],"single t","f") #right
            #leg.AddEntry(hsig3,sig3leg+scaleFact3Str,"l") #left
            leg.AddEntry(bkgHT,"Bkg. uncert.","f") #right
            
            leg.Draw("same")
            
            prelimTex2=TLatex()
            prelimTex2.SetNDC()
            prelimTex2.SetTextFont(61)
            prelimTex2.SetLineWidth(2)
            prelimTex2.SetTextSize(0.08)
            if blind: prelimTex2.SetTextSize(0.08)
            prelimTex2.DrawLatex(0.12,0.93,"CMS")
            
            prelimTex3=TLatex()
            prelimTex3.SetNDC()
            prelimTex3.SetTextAlign(12)
            prelimTex3.SetTextFont(52)
            prelimTex3.SetTextSize(0.055)
            prelimTex3.SetLineWidth(2)
            if blind: prelimTex3.DrawLatex(0.26,0.945,"Simulation work in progress") #"Preliminary")
    
            # making up names of the image files, making the folder exists, and doing SaveAs
            savePrefix = outdir+histPrefix
            if yLog: savePrefix+='_logy'
            #c1.SaveAs(savePrefix+".pdf")
            c1.SaveAs(savePrefix+".png")

print("--- %s s for booking, %s min to access pointers, %s s for plots ---" % (round(ptr_time - start_time, 2), round(plot_time - ptr_time,2)/60, round(time.time() - plot_time, 2)))
    



