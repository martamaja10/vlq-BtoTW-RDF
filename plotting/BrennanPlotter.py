#!/usr/bin/python

import os,sys,time,math,pickle,itertools
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
from ROOT import *
#from weights import *
#from modSyst import *
#from utils import *

gROOT.SetBatch(1)
start_time = time.time()

lumi=137 #for plots #56.1 #
lumiInTemplates= '137'#str(targetlumi/1000).replace('.','p') # 1/fb

saveKey = '' # tag for plot names

# labels for the legend
sig1leg='B (0.8 TeV)'
sig2leg='B (1.4 TeV)'
sig3leg='B (2.0 TeV)'

# zoom in/out on signal
sigScaleFact = 400
print 'Scale factor = ',sigScaleFact

bkgProcList = ['ewk','ttbar','st']
bkgHistColors = {'ttbar':kAzure+8,'ewk':kMagenta-2,'st':kGreen-3} #TT

yLog  = True

lumiSys = 0.20 # 20% uncertainty on background

tagList = ['T','W','Wb']
for tag in tagList: # could think about looping over your T, W, Wb

        histPrefix='BprimeMass_'
        catStr='is'+tag
        histPrefix+=catStr
	
        # add the backgrounds together
        bkgHT = histoTTbar_certaincat.Clone()
        bkgHT.Add(histoWjets_certaincat)

        # set a larger uncertainty
        for ibin in range(1,bkgHT.GetNbinsX()+1):
                bkgHT.SetBinError(ibin,bkgHT.GetBinContent(ibin)*lumiSys)

        hsig800.Scale(sigScaleFact)
	hsig1400.Scale(sigScaleFact)
	hsig2000.Scale(sigScaleFact)

        ############################################################
        ############## Making Plots of e+jets, mu+jets and e/mu+jets 
        ############################################################
		

        stackbkgHT = THStack("stackbkgHT","")
        stackbkgHT.Add(ttbarHisto)
        stackbkgHT.Add(wjetsHisto)
        # set each background a different color
        for proc in bkgProcList:
                bkghists[proc+catStr].SetLineColor(bkgHistColors[proc])
                bkghists[proc+catStr].SetFillColor(bkgHistColors[proc])
                bkghists[proc+catStr].SetLineWidth(2)


        sig1Color= kBlack
        sig2Color= kBlack
        sig3Color= kBlack
			
	hsig1.SetLineColor(sig1Color)
	hsig1.SetFillStyle(0)
	hsig1.SetLineWidth(3)
	hsig2.SetLineColor(sig2Color)
	hsig2.SetLineStyle(7)#5)
	hsig2.SetFillStyle(0)
	hsig2.SetLineWidth(3)
	
        bkgHT.SetFillStyle(3004)
        bkgHT.SetFillColor(kBlack)

        gStyle.SetOptStat(0)
        c1 = TCanvas("c1","c1",1200,1000)

        stackbkgHT.Draw("HIST")
        hsig1.Draw("SAME HIST")
        hsig2.Draw("SAME HIST")
        hsig3.Draw("SAME HIST")
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
        scaleFact1Str = ' x'+str(scaleFact1)
        scaleFact2Str = ' x'+str(scaleFact2)

        leg.AddEntry(hsig1,sig1leg+scaleFact1Str,"l") #left
        leg.AddEntry(bkghists['ttbar'+catStr],"t#bar{t}","f") #right
        leg.AddEntry(hsig2,sig2leg+scaleFact2Str,"l") #left
        leg.AddEntry(bkghists['ewk'+catStr],"EW","f") #right
        leg.AddEntry(hsig3,sig3leg+scaleFact3Str,"l") #left
        leg.AddEntry(bkgHTgerr,"Bkg. uncert.","f") #right

        
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
        savePrefix = templateDir+templateDir.split('/')[-2]+'plots/'
        if not os.path.exists(savePrefix): os.system('mkdir '+savePrefix)
        savePrefix+=histPrefix+isRebinned.replace('_rebinned_stat1p1','')+saveKey
        if yLog: savePrefix+='_logy'
        c1.SaveAs(savePrefix+".pdf")
        c1.SaveAs(savePrefix+".png")

print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))


