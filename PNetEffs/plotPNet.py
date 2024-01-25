from ROOT import *

gROOT.SetBatch(1)

# // Preference here would be t efficiency with 3 masses -- one plot per year? Yes, could implement that way  
# // t-as-t, t-as-W, t-as-J, W-as-t, W-as-W, W-as-J, J-as-t, J-as-J

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
# // auto h11 = Efficiency.Histo1D({"PNEff_Nptbins_J_J",";Jet p_T (GeV);QCD efficiency",9,ptbins},"J_match_J_tag","weight");

masses = ['800','1000','1200','1400','1600','1800','2000','2200']
bkgs = ['TTToSemiLeptonic'] # add this later when they are done running
years = ['2016APV','2016','2017','2018']
matches = ['t','W','J']
tags = ['t','W','J']

sighists = {}
bkghists = {}
for year in years:
    for mass in masses:
        ifile = TFile.Open("root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Jan2024_pNetEffs/RDF_BprimeBtoTW_M-"+mass+"_NWALO_TuneCP5_13TeV-madgraph-pythia8_"+year+"_PNetEff.root");
        for match in matches:
            sighists["den_"+mass+"_"+match] = ifile.Get("PNEff_Dptbins_"+match)
            sighists["den_"+mass+"_"+match].SetDirectory(0)
            for tag in tags:
                if tag == 'W' and match == 'J': continue
                sighists["num_"+mass+"_"+match+"_"+tag] = ifile.Get("PNEff_Nptbins_"+match+"_"+tag)
                sighists["num_"+mass+"_"+match+"_"+tag].SetDirectory(0)
                sighists["num_"+mass+"_"+match+"_"+tag].Divide(sighists["num_"+mass+"_"+match+"_"+tag],sighists["den_"+mass+"_"+match],1.0,1.0,"B")
                # can add a print here if desired
        ifile.Close()

    c1 = TCanvas("c1","c1",800,600)
    gStyle.SetOptStat(0)

    for match in matches:
        for tag in tags:
            if tag == 'W' and match == 'J': continue

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
            for mass in masses:
                leg.AddEntry(sighists["num_"+mass+"_"+match+"_"+tag],"B "+mass+" TeV","pl");
            leg.Draw();
            c1.SaveAs("PNetEff_Signal_"+year+"_"+match+"_"+tag+".png");
            c1.SaveAs("PNetEff_Signal_"+year+"_"+match+"_"+tag+".pdf");
            c1.SaveAs("PNetEff_Signal_"+year+"_"+match+"_"+tag+".root");

