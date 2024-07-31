import os
from ROOT import *

gROOT.SetBatch(1)

  #   .Histo2D<float,float>({"TrigEff_Dptbins_jet185_El",";electron p_T (GeV);electron #eta",5,ptbins,10,etabinsEl},"probe_pt","probe_eta");
  # auto e04 = JetVars.Filter("isTagMuProbeEl == 1 && passElTrig == 1 && leadjetpt > 185")
  #   .Histo2D<float,float>({"TrigEff_Nptbins_jet185_El",";electron p_T (GeV);electron #eta",5,ptbins,10,etabinsEl},"probe_pt","probe_eta");
  # auto e05 = JetVars.Filter("isTagMuProbeEl == 1 && goodjet_HT > 250")
  #   .Histo2D<float,float>({"TrigEff_Dptbins_ht250_El",";electron p_T (GeV);electron #eta",5,ptbins,10,etabinsEl},"probe_pt","probe_eta");
  # auto e06 = JetVars.Filter("isTagMuProbeEl == 1 && passElTrig == 1 && goodjet_HT > 250")
  #   .Histo2D<float,float>({"TrigEff_Nptbins_ht250_El",";electron p_T (GeV);electron #eta",5,ptbins,10,etabinsEl},"probe_pt","probe_eta");
  # auto e07 = JetVars.Filter("isTagMuProbeEl == 1 && goodjet_HT > 250 && NJets_DeepFlavM > 1")
  #   .Histo2D<float,float>({"TrigEff_Dptbins_ht250b2_El",";electron p_T (GeV);electron #eta",5,ptbins,10,etabinsEl},"probe_pt","probe_eta");
  # auto e08 = JetVars.Filter("isTagMuProbeEl == 1 && passElTrig == 1 && goodjet_HT > 250 && NJets_DeepFlavM > 1")
  #   .Histo2D<float,float>({"TrigEff_Nptbins_ht250b2_El",";electron p_T (GeV);electron #eta",5,ptbins,10,etabinsEl},"probe_pt","probe_eta");
  # auto e09 = JetVars.Filter("isTagMuProbeEl == 1 && leadjetpt > 185 && goodjet_HT > 250")
  #   .Histo2D<float,float>({"TrigEff_Dptbins_jet185ht250_El",";electron p_T (GeV);electron #eta",5,ptbins,10,etabinsEl},"probe_pt","probe_eta");
  # auto e10 = JetVars.Filter("isTagMuProbeEl == 1 && passElTrig == 1 && leadjetpt > 185 && goodjet_HT > 250")
  #   .Histo2D<float,float>({"TrigEff_Nptbins_jet185ht250_El",";electron p_T (GeV);electron #eta",5,ptbins,10,etabinsEl},"probe_pt","probe_eta");
  # auto e11 = JetVars.Filter("isTagMuProbeEl == 1 && leadjetpt > 185 && goodjet_HT > 250 && NJets_DeepFlavM > 1")
  #   .Histo2D<float,float>({"TrigEff_Dptbins_jet185ht250b2_El",";electron p_T (GeV);electron #eta",5,ptbins,10,etabinsEl},"probe_pt","probe_eta");
  # auto e12 = JetVars.Filter("isTagMuProbeEl == 1 && passElTrig == 1 && leadjetpt > 185 && goodjet_HT > 250 && NJets_DeepFlavM > 1")
  #   .Histo2D<float,float>({"TrigEff_Nptbins_jet185ht250b2_El",";electron p_T (GeV);electron #eta",5,ptbins,10,etabinsEl},"probe_pt","probe_eta");

haddsamples = False
printSamples = True
plotSamples = False
labelsEl = ['jet185ht250_']#,'jet185_','ht250_']#,'jet185ht250b2_','ht250b2_']
labelsMu = ['ht250_','ht250b2_','ht250_','ht250b2_','', ]
#plotYears = False

output = 'Jun2024Effs_elIDsNew'
if not os.path.exists(output): os.system('mkdir '+output)
samples = ['TTTo2L2Nu','SingleMuon']#,'SingleElectron']#
years = ['2016APV','2016','2017','2018']
#probes = ['Mu']
probes = ['El']#,'El2','El3','El4']
targetlumi = {'2016APV':19.5, '2016':16.8, '2017':41.5, '2018':59.8, 'all':138}

if haddsamples:
    for year in years: 
        for sample in samples:
            grepsample = sample
            if year == '2018' and sample == 'SingleElectron':
                grepsample = 'EGamma'
            haddcommand = 'hadd root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Jun2024_trigEffs_elIDsNew_hadds/RDF_'+sample+'_'+year+'_TrigEff.root '
            haddcommand += '`xrdfs root://cmseos.fnal.gov/ ls -u /store/user/jmanagan/BtoTW_Jun2024_trigEffs_elIDs/ | grep "'+grepsample+'" | grep "_'+year+'_"` '
            os.system(haddcommand)

#for year in years:
#    for flav in probes:
#        print('std::vector<std::vector<float>> hltsf_'+flav+'_'+year+';')

for ilabel in range(len(labelsEl)):
    labelEl = labelsEl[ilabel]
    labelMu = labelsMu[ilabel]
    for year in years:
        print('Working on '+year+'...')
        trighists = {}

        for probe in probes:
            for sample in samples:
                ifile = TFile.Open("root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Jun2024_trigEffs_elIDsNew_hadds/RDF_"+sample+"_"+year+"_TrigEff.root");

                print('\t\t Working on '+sample+'_'+probe+'...')

                if 'El' in probe:
                    trighists["den_"+sample+"_"+probe] = ifile.Get("TrigEff_Dptbins_"+labelEl+probe)
                    trighists["num_"+sample+"_"+probe] = ifile.Get("TrigEff_Nptbins_"+labelEl+probe)
                else:
                    trighists["den_"+sample+"_"+probe] = ifile.Get("TrigEff_Dptbins_"+labelMu+probe)
                    trighists["num_"+sample+"_"+probe] = ifile.Get("TrigEff_Nptbins_"+labelMu+probe)
                trighists["den_"+sample+"_"+probe].SetDirectory(0)
                trighists["num_"+sample+"_"+probe].SetDirectory(0)
                trighists["num_"+sample+"_"+probe].Divide(trighists["num_"+sample+"_"+probe],trighists["den_"+sample+"_"+probe],1.0,1.0,"B")
                ifile.Close()                

            dataprobe = 'SingleMuon'
            if probe == 'Mu': dataprobe = 'SingleElectron'
            trighists["ratio_"+probe] = trighists["num_TTTo2L2Nu_"+probe].Clone("ratio_"+probe)
            trighists["ratio_"+probe].Divide(trighists["num_"+dataprobe+"_"+probe],trighists["num_TTTo2L2Nu_"+probe],1.0,1.0,"B")

            if(printSamples):
                print('hltsf_'+probe+'_'+year+' = {')
                for jbin in range(1,trighists["ratio_"+probe].GetNbinsY()+1): # rows are eta
                    toprint = '\t{'
                    for ibin in range(1,trighists["ratio_"+probe].GetNbinsX()+1): # columns are pt
                        if ibin < trighists["ratio_"+probe].GetNbinsX(): 
                            toprint += str(round(trighists["ratio_"+probe].GetBinContent(ibin,jbin),5))+', '
                        else: 
                            toprint += str(round(trighists["ratio_"+probe].GetBinContent(ibin,jbin),5))
                    print(toprint+'},')
                print('};')
                print('hltsfUnc_'+probe+'_'+year+' = {')
                for jbin in range(1,trighists["ratio_"+probe].GetNbinsY()+1):
                    toprint = '\t{'
                    for ibin in range(1,trighists["ratio_"+probe].GetNbinsX()+1):
                        if ibin < trighists["ratio_"+probe].GetNbinsX(): 
                            toprint += str(round(trighists["ratio_"+probe].GetBinError(ibin,jbin),5))+', '
                        else: 
                            toprint += str(round(trighists["ratio_"+probe].GetBinError(ibin,jbin),5))
                    print(toprint+'},')
                print('};')

        if not plotSamples: continue
            
        c1 = TCanvas("c1","c1",1600,600)
        gStyle.SetOptStat(0)
        gStyle.SetPaintTextFormat("1.2f")
        c1.SetLogx(True)

        for probe in probes:
            lepstr = "electron"
            dataprobe = 'SingleMuon'
            tagstr = labelEl
            if probe == 'Mu':
                dataprobe = 'SingleElectron'
                lepstr = "muon"
                tagstr = labelMu

            trighists["ratio_"+probe].GetXaxis().SetTitle(lepstr+" p_{T} [GeV]")
            trighists["ratio_"+probe].GetYaxis().SetTitle(lepstr+" #eta")
            trighists["ratio_"+probe].GetZaxis().SetRangeUser(0.0,1.30)
            trighists["ratio_"+probe].Draw("colz texte")
            
            prelimTex=TLatex()
            prelimTex.SetNDC()
            prelimTex.SetTextAlign(31) # align right
            prelimTex.SetTextFont(42)
            prelimTex.SetTextSize(0.05)
            prelimTex.SetLineWidth(2)
            prelimTex.DrawLatex(0.95,0.94,str(targetlumi[year])+" fb^{-1} (13 TeV)")
            
            prelimTex3=TLatex()
            prelimTex3.SetNDC()
            prelimTex3.SetTextAlign(12)
            prelimTex3.SetTextSize(0.035)
            prelimTex3.SetLineWidth(2)
            prelimTex3.DrawLatex(0.15,0.945,"Private work (CMS Data & Simulation)") #"Preliminary")
            
            c1.SaveAs(output+"/TrigEff_SF_"+year+"_"+probe+tagstr+".png");
            c1.SaveAs(output+"/TrigEff_SF_"+year+"_"+probe+tagstr+".pdf");
            c1.SaveAs(output+"/TrigEff_SF_"+year+"_"+probe+tagstr+".root");
            
            trighists["num_"+dataprobe+"_"+probe].GetXaxis().SetTitle(lepstr+" p_{T} [GeV]")
            trighists["num_"+dataprobe+"_"+probe].GetYaxis().SetTitle(lepstr+" #eta");
            trighists["num_"+dataprobe+"_"+probe].GetZaxis().SetRangeUser(0.0,1.30)
            trighists["num_"+dataprobe+"_"+probe].Draw("colz texte")
            
            prelimTex=TLatex()
            prelimTex.SetNDC()
            prelimTex.SetTextAlign(31) # align right
            prelimTex.SetTextFont(42)
            prelimTex.SetTextSize(0.05)
            prelimTex.SetLineWidth(2)
            prelimTex.DrawLatex(0.95,0.94,str(targetlumi[year])+" fb^{-1} (13 TeV)")
            
            prelimTex3=TLatex()
            prelimTex3.SetNDC()
            prelimTex3.SetTextAlign(12)
            prelimTex3.SetTextSize(0.035)
            prelimTex3.SetLineWidth(2)
            prelimTex3.DrawLatex(0.15,0.945,"Private work (CMS Data)") #"Preliminary")
            
            c1.SaveAs(output+"/TrigEff_DataEff_"+year+"_"+probe+tagstr+".png");
            c1.SaveAs(output+"/TrigEff_DataEff_"+year+"_"+probe+tagstr+".pdf");
            c1.SaveAs(output+"/TrigEff_DatEff_"+year+"_"+probe+tagstr+".root");
            
            trighists["num_TTTo2L2Nu_"+probe].GetXaxis().SetTitle(lepstr+" p_{T} [GeV]")
            trighists["num_TTTo2L2Nu_"+probe].GetYaxis().SetTitle(lepstr+" #eta");
            trighists["num_TTTo2L2Nu_"+probe].GetZaxis().SetRangeUser(0.0,1.30)
            trighists["num_TTTo2L2Nu_"+probe].Draw("colz texte")
            
            prelimTex=TLatex()
            prelimTex.SetNDC()
            prelimTex.SetTextAlign(31) # align right
            prelimTex.SetTextFont(42)
            prelimTex.SetTextSize(0.05)
            prelimTex.SetLineWidth(2)
            prelimTex.DrawLatex(0.95,0.94,str(targetlumi[year])+" fb^{-1} (13 TeV)")
            
            prelimTex3=TLatex()
            prelimTex3.SetNDC()
            prelimTex3.SetTextAlign(12)
            prelimTex3.SetTextSize(0.035)
            prelimTex3.SetLineWidth(2)
            prelimTex3.DrawLatex(0.15,0.945,"Private work (CMS Data & Simulation)") #"Preliminary")
            
            c1.SaveAs(output+"/TrigEff_MCEff_"+year+"_"+probe+tagstr+".png");
            c1.SaveAs(output+"/TrigEff_MCEff_"+year+"_"+probe+tagstr+".pdf");
            c1.SaveAs(output+"/TrigEff_MCEff_"+year+"_"+probe+tagstr+".root");

