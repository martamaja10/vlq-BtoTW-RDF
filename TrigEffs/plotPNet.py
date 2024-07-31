import os
from ROOT import *

gROOT.SetBatch(1)

# // Preference here would be t efficiency with 3 masses -- one plot per year? Yes, could implement that way  
# // t-as-t, t-as-W, t-as-J, W-as-t, W-as-W, W-as-J, J-as-t, J-as-W, J-as-J

# // auto h01 = Efficiency.Histo1D({"PNEff_Dptbins_t",";Jet p_T (GeV);",9,ptbins},"t_match","weight");
# // auto h02 = Efficiency.Histo1D({"PNEff_Dptbins_W",";Jet p_T (GeV);",9,ptbins},"W_match","weight");
# // auto h03 = Efficiency.Histo1D({"PNEff_Dptbins_J",";Jet p_T (GeV);",9,ptbins},"J_match","weight");
# // auto h04 = Efficiency.Histo1D({"PNEff_Nptbins_t_t",";Jet p_T (GeV);t efficiency",9,ptbins},"t_match_t_tag","weight");
# // auto h05 = Efficiency.Histo1D({"PNEff_Nptbins_t_W",";Jet p_T (GeV);t-as-W rate",9,ptbins},"t_match_W_tag","weight");
# // auto h06 = Efficiency.Histo1D({"PNEff_Nptbins_t_J",";Jet p_T (GeV);t-as-QCD rate",9,ptbins},"t_match_J_tag","weight");
# // auto h07 = Efficiency.Histo1D({"PNEff_Nptbins_W_t",";Jet p_T (GeV);W-as-t rate",9,ptbins},"W_match_t_tag","weight");
# // auto h08 = Efficiency.Histo1D({"PNEff_Nptbins_W_W",";Jet p_T (GeV);W efficiency",9,ptbins},"W_match_W_tag","weight");
# // auto h09 = Efficiency.Histo1D({"PNEff_Nptbins_W_J",";Jet p_T (GeV);W-as-QCD rate",9,ptbins},"W_match_J_tag","weight");
# // auto h10 = Efficiency.Histo1D({"PNEff_Nptbins_J_t",";Jet p_T (GeV);QCD-as-t rate",9,ptbins},"J_match_t_tag","weight");
# // auto h11 = Efficiency.Histo1D({"PNEff_Nptbins_J_W",";Jet p_T (GeV);QCD-as-W rate",9,ptbins},"J_match_W_tag","weight");
# // auto h12 = Efficiency.Histo1D({"PNEff_Nptbins_J_J",";Jet p_T (GeV);QCD efficiency",9,ptbins},"J_match_J_tag","weight");

haddbkgs = False
plotSamples = False
plotYears = False

output = 'Apr2024Effs_tWseparate'
if not os.path.exists(output): os.system('mkdir '+output)
masses = ['800','1000','1200','1300','1400','1500','1600','1700','1800','2000','2200']
bkgs = ['ttbar','single-t','ttbarVH','Wjets','Zjets','VV','qcd']
bkgstrs = {'ttbar':['TTTo','TT_Mtt'],'single-t':['ST'],'ttbarVH':['TTW','TTZ','ttH'],'Wjets':['WJets'],'Zjets':['DY'],'VV':['WW','WZ','ZZ'],'qcd':['QCD']}
years = ['2016APV','2016','2017','2018','all']
matches = ['t','W','J']
tags = ['t','W']
targetlumi = {'2016APV':19.5, '2016':16.8, '2017':41.5, '2018':59.8, 'all':138}

if haddbkgs:
    for year in years: 
        if year == 'all': continue
        for bkg in bkgs:
            haddcommand = 'hadd root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Apr2024_pNetEffs_tWseparate/RDF_'+bkg+'_'+year+'_PNetEff.root '
            for bkgstr in bkgstrs[bkg]:
                haddcommand += '`xrdfs root://cmseos.fnal.gov/ ls -u /store/user/jmanagan/BtoTW_Apr2024_pNetEffs_tWseparate/ | grep "RDF_'+bkgstr+'" | grep "'+year+'_"` '
            os.system(haddcommand)
    for bkg in bkgs:
        haddcommand = 'hadd root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Apr2024_pNetEffs_tWseparate/RDF_'+bkg+'_all_PNetEff.root '
        for bkgstr in bkgstrs[bkg]:
            haddcommand += '`xrdfs root://cmseos.fnal.gov/ ls -u /store/user/jmanagan/BtoTW_Apr2024_pNetEffs_tWseparate/ | grep "RDF_'+bkgstr+'"` '
        os.system(haddcommand)
    for mass in masses:
        haddcommand = 'hadd root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Apr2024_pNetEffs_tWseparate/RDF_BprimeBtoTW_M-'+mass+'_NWALO_TuneCP5_13TeV-madgraph-pythia8_all_PNetEff.root '
        for bkgstr in bkgstrs[bkg]:
            haddcommand += '`xrdfs root://cmseos.fnal.gov/ ls -u /store/user/jmanagan/BtoTW_Apr2024_pNetEffs_tWseparate/ | grep "RDF_BprimeBtoTW_M-'+mass+'"` '
        os.system(haddcommand)


for match in matches:
    for tag in tags:
        print('std::vector<std::vector<float>> pnet_'+match+'_'+tag+';')

for year in years:
    #    print('Working on '+year+'...')
    sighists = {}
    bkghists = {}

    for bkg in bkgs:
        #print('\t Working on '+bkg+'...')
        ifile = TFile.Open("root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Apr2024_pNetEffs_tWseparate/RDF_"+bkg+"_"+year+"_PNetEff.root");
        for match in matches:
            #print('\t\t Working on '+match+'...')
            bkghists["den_"+bkg+"_"+match] = ifile.Get("PNEff_Dptbins_"+match)
            bkghists["den_"+bkg+"_"+match].SetDirectory(0)
            for tag in tags:
                #print('\t\t\t Working on '+tag+'...')
                bkghists["num_"+bkg+"_"+match+"_"+tag] = ifile.Get("PNEff_Nptbins_"+match+"_"+tag)
                bkghists["num_"+bkg+"_"+match+"_"+tag].SetDirectory(0)
                bkghists["num_"+bkg+"_"+match+"_"+tag].Divide(bkghists["num_"+bkg+"_"+match+"_"+tag],bkghists["den_"+bkg+"_"+match],1.0,1.0,"B")
        ifile.Close()                

    for mass in masses:
        #print('\t Working on '+mass+'...')
        ifile = TFile.Open("root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Apr2024_pNetEffs_tWseparate/RDF_BprimeBtoTW_M-"+mass+"_NWALO_TuneCP5_13TeV-madgraph-pythia8_"+year+"_PNetEff.root");
        for match in matches:
            #print('\t\t Working on '+match+'...')
            sighists["den_"+mass+"_"+match] = ifile.Get("PNEff_Dptbins_"+match)
            sighists["den_"+mass+"_"+match].SetDirectory(0)
            for tag in tags:
                #print('\t\t\t Working on '+tag+'...')
                sighists["num_"+mass+"_"+match+"_"+tag] = ifile.Get("PNEff_Nptbins_"+match+"_"+tag)
                sighists["num_"+mass+"_"+match+"_"+tag].SetDirectory(0)
                sighists["num_"+mass+"_"+match+"_"+tag].Divide(sighists["num_"+mass+"_"+match+"_"+tag],sighists["den_"+mass+"_"+match],1.0,1.0,"B")
        ifile.Close()

    if year == 'all':
        for match in matches:
            for tag in tags:
                print('pnet_'+match+'_'+tag+' = {')
                for bkg in bkgs:
                    toprint = '\t{'
                    for ibin in range(bkghists["num_"+bkg+"_"+match+"_"+tag].GetNbinsX()+1):
                        if ibin < bkghists["num_"+bkg+"_"+match+"_"+tag].GetNbinsX(): 
                            toprint += str(round(bkghists["num_"+bkg+"_"+match+"_"+tag].GetBinContent(ibin),5))+', '
                        else: 
                            toprint += str(round(bkghists["num_"+bkg+"_"+match+"_"+tag].GetBinContent(ibin),5))
                    toprint += '}, // '+bkg
                    print(toprint)
                for mass in masses:
                    toprint = '\t{'
                    for ibin in range(sighists["num_"+mass+"_"+match+"_"+tag].GetNbinsX()+1):
                        if ibin < sighists["num_"+mass+"_"+match+"_"+tag].GetNbinsX(): 
                            toprint += str(round(sighists["num_"+mass+"_"+match+"_"+tag].GetBinContent(ibin),5))+', '
                        else: 
                            toprint += str(round(sighists["num_"+mass+"_"+match+"_"+tag].GetBinContent(ibin),5))
                    toprint += '}, // BprimeM'+mass
                    print(toprint)
                print('};')




    if not plotSamples: continue
    c1 = TCanvas("c1","c1",800,600)
    gStyle.SetOptStat(0)

    for match in matches:
        for tag in tags:
            icolor = 1
            for mass in masses:
                if mass == '800':
                    sighists["num_"+mass+"_"+match+"_"+tag].GetYaxis().SetRangeUser(0.0,1.0)
                    sighists["num_"+mass+"_"+match+"_"+tag].GetYaxis().SetTitle(year+" PNet: "+match+"-match with "+tag+"-tag rate")
                    sighists["num_"+mass+"_"+match+"_"+tag].GetXaxis().SetTitle("AK8 jet p_{T} [GeV]");
                sighists["num_"+mass+"_"+match+"_"+tag].SetLineWidth(2)
                sighists["num_"+mass+"_"+match+"_"+tag].SetMarkerStyle(20)
                sighists["num_"+mass+"_"+match+"_"+tag].SetMarkerColor(icolor)
                sighists["num_"+mass+"_"+match+"_"+tag].SetLineColor(icolor)
                icolor += 1

                if mass == '800': 
                    sighists["num_"+mass+"_"+match+"_"+tag].Draw()
                else:
                    sighists["num_"+mass+"_"+match+"_"+tag].Draw("same")

            leg = TLegend(0.8,0.3,0.98,0.7);
            leg.SetHeader(year.replace('all','full Run 2'))
            for mass in masses:
                leg.AddEntry(sighists["num_"+mass+"_"+match+"_"+tag],"B "+mass+" GeV","pl");
            leg.Draw();

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
            prelimTex3.DrawLatex(0.15,0.945,"Private work (CMS Simulation)") #"Preliminary")

            c1.SaveAs(output+"/PNetEff_Signal_"+year+"_"+match+"_"+tag+".png");
            c1.SaveAs(output+"/PNetEff_Signal_"+year+"_"+match+"_"+tag+".pdf");
            c1.SaveAs(output+"/PNetEff_Signal_"+year+"_"+match+"_"+tag+".root");

    c2 = TCanvas("c2","c2",800,600)
    gStyle.SetOptStat(0)

    for match in matches:
        for tag in tags:
            icolor = 1
            for bkg in bkgs:
                if bkg == 'ttbar':
                    bkghists["num_"+bkg+"_"+match+"_"+tag].GetYaxis().SetRangeUser(0.0,1.0)
                    bkghists["num_"+bkg+"_"+match+"_"+tag].GetYaxis().SetTitle(year+" PNet: "+match+"-match with "+tag+"-tag rate")
                    bkghists["num_"+bkg+"_"+match+"_"+tag].GetXaxis().SetTitle("AK8 jet p_{T} [GeV]");
                bkghists["num_"+bkg+"_"+match+"_"+tag].SetLineWidth(2)
                bkghists["num_"+bkg+"_"+match+"_"+tag].SetMarkerStyle(20)
                bkghists["num_"+bkg+"_"+match+"_"+tag].SetMarkerColor(icolor)
                bkghists["num_"+bkg+"_"+match+"_"+tag].SetLineColor(icolor)
                icolor += 1

                if bkg == 'ttbar': 
                    bkghists["num_"+bkg+"_"+match+"_"+tag].Draw()
                else:
                    bkghists["num_"+bkg+"_"+match+"_"+tag].Draw("same")

            leg = TLegend(0.8,0.3,0.98,0.7);
            leg.SetHeader(year.replace('all','full Run 2'))
            for bkg in bkgs:
                leg.AddEntry(bkghists["num_"+bkg+"_"+match+"_"+tag],bkg.replace('tbar','#bar{t}').replace('single-t','single t').replace('VH','+V/H').replace('jets','+jets').replace('qcd','QCD'),"pl");
            leg.Draw();

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
            prelimTex3.DrawLatex(0.15,0.945,"Private work (CMS Simulation)") #"Preliminary")

            c2.SaveAs(output+"/PNetEff_Bkgs_"+year+"_"+match+"_"+tag+".png");
            c2.SaveAs(output+"/PNetEff_Bkgs_"+year+"_"+match+"_"+tag+".pdf");
            c2.SaveAs(output+"/PNetEff_Bkgs_"+year+"_"+match+"_"+tag+".root");


            #    print('Working on '+year+'...')

if not plotYears: exit(0)

for bkg in bkgs:
    bkghists = {}

    for year in years:
        print('\t Working on '+bkg+'...')
        ifile = TFile.Open("root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Apr2024_pNetEffs_tWseparate/RDF_"+bkg+"_"+year+"_PNetEff.root");
        for match in matches:
            #print('\t\t Working on '+match+'...')
            bkghists["den_"+year+"_"+match] = ifile.Get("PNEff_Dptbins_"+match)
            bkghists["den_"+year+"_"+match].SetDirectory(0)
            for tag in tags:
                #print('\t\t\t Working on '+tag+'...')
                bkghists["num_"+year+"_"+match+"_"+tag] = ifile.Get("PNEff_Nptbins_"+match+"_"+tag)
                bkghists["num_"+year+"_"+match+"_"+tag].SetDirectory(0)
                bkghists["num_"+year+"_"+match+"_"+tag].Divide(bkghists["num_"+year+"_"+match+"_"+tag],bkghists["den_"+year+"_"+match],1.0,1.0,"B")
                # can add a print here if desired
        ifile.Close()

    c3 = TCanvas("c3","c3",800,600)
    gStyle.SetOptStat(0)

    for match in matches:
        for tag in tags:
            icolor = 1
            for year in years:
                if year == '2016APV':
                    bkghists["num_"+year+"_"+match+"_"+tag].GetYaxis().SetRangeUser(0.0,1.0)
                    bkghists["num_"+year+"_"+match+"_"+tag].GetYaxis().SetTitle(year+" PNet: "+match+"-match with "+tag+"-tag rate")
                    bkghists["num_"+year+"_"+match+"_"+tag].GetXaxis().SetTitle("AK8 jet p_{T} [GeV]");
                bkghists["num_"+year+"_"+match+"_"+tag].SetLineWidth(2)
                bkghists["num_"+year+"_"+match+"_"+tag].SetMarkerStyle(20)
                bkghists["num_"+year+"_"+match+"_"+tag].SetMarkerColor(icolor)
                bkghists["num_"+year+"_"+match+"_"+tag].SetLineColor(icolor)
                icolor += 1

                if year == '2016APV': 
                    bkghists["num_"+year+"_"+match+"_"+tag].Draw()
                else:
                    bkghists["num_"+year+"_"+match+"_"+tag].Draw("same")

            leg = TLegend(0.8,0.3,0.98,0.7);
            leg.SetHeader(bkg.replace('tbar','#bar{t}').replace('single-t','single t').replace('VH','+V/H').replace('jets','+jets').replace('qcd','QCD'))
            for year in years:
                leg.AddEntry(bkghists["num_"+year+"_"+match+"_"+tag],year.replace('all','full Run 2'),"pl");
            leg.Draw();

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
            prelimTex3.DrawLatex(0.15,0.945,"Private work (CMS Simulation)") #"Preliminary")

            c3.SaveAs(output+"/PNetEff_Years_"+bkg+"_"+match+"_"+tag+".png");
            c3.SaveAs(output+"/PNetEff_Years_"+bkg+"_"+match+"_"+tag+".pdf");
            c3.SaveAs(output+"/PNetEff_Years_"+bkg+"_"+match+"_"+tag+".root");


        
for mass in masses:
    sighists = {}

    for year in years:
        #print('\t Working on '+mass+'...')
        ifile = TFile.Open("root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Apr2024_pNetEffs_tWseparate/RDF_BprimeBtoTW_M-"+mass+"_NWALO_TuneCP5_13TeV-madgraph-pythia8_"+year+"_PNetEff.root");
        for match in matches:
            #print('\t\t Working on '+match+'...')
            sighists["den_"+year+"_"+match] = ifile.Get("PNEff_Dptbins_"+match)
            sighists["den_"+year+"_"+match].SetDirectory(0)
            for tag in tags:
                #print('\t\t\t Working on '+tag+'...')
                sighists["num_"+year+"_"+match+"_"+tag] = ifile.Get("PNEff_Nptbins_"+match+"_"+tag)
                sighists["num_"+year+"_"+match+"_"+tag].SetDirectory(0)
                sighists["num_"+year+"_"+match+"_"+tag].Divide(sighists["num_"+year+"_"+match+"_"+tag],sighists["den_"+year+"_"+match],1.0,1.0,"B")
                # can add a print here if desired
        ifile.Close()

    c4 = TCanvas("c4","c4",800,600)
    gStyle.SetOptStat(0)

    for match in matches:
        for tag in tags:
            icolor = 1
            for year in years:
                if year == '2016APV':
                    sighists["num_"+year+"_"+match+"_"+tag].GetYaxis().SetRangeUser(0.0,1.0)
                    sighists["num_"+year+"_"+match+"_"+tag].GetYaxis().SetTitle(year+" PNet: "+match+"-match with "+tag+"-tag rate")
                    sighists["num_"+year+"_"+match+"_"+tag].GetXaxis().SetTitle("AK8 jet p_{T} [GeV]");
                sighists["num_"+year+"_"+match+"_"+tag].SetLineWidth(2)
                sighists["num_"+year+"_"+match+"_"+tag].SetMarkerStyle(20)
                sighists["num_"+year+"_"+match+"_"+tag].SetMarkerColor(icolor)
                sighists["num_"+year+"_"+match+"_"+tag].SetLineColor(icolor)
                icolor += 1

                if year == '2016APV': 
                    sighists["num_"+year+"_"+match+"_"+tag].Draw()
                else:
                    sighists["num_"+year+"_"+match+"_"+tag].Draw("same")

            leg = TLegend(0.8,0.3,0.98,0.7);
            leg.SetHeader("B "+mass+" GeV");
            for year in years:
                leg.AddEntry(sighists["num_"+year+"_"+match+"_"+tag],year.replace('all','full Run 2'),"pl");
            leg.Draw();

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
            prelimTex3.DrawLatex(0.15,0.945,"Private work (CMS Simulation)") #"Preliminary")

            c4.SaveAs(output+"/PNetEff_Years_"+mass+"_"+match+"_"+tag+".png");
            c4.SaveAs(output+"/PNetEff_Years_"+mass+"_"+match+"_"+tag+".pdf");
            c4.SaveAs(output+"/PNetEff_Years_"+mass+"_"+match+"_"+tag+".root");

