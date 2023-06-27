import ROOT
import os 

class sample:
    def __init__(self, prefix, textlist, samplename): #self, color, style, fill, leglabel, label, name=""):
        self.prefix = prefix
        self.textlist = textlist
        self.samplename = samplename
        #self.color = color # I don't know what the attributes I should put in here are
        #self.style = style
        #self.fill = fill
        #self.leglabel = leglabel
        #self.label = label
        #if name == "":
        #    self.name = label
        #else:
        #    self.name = name
            
singleTb = sample("singleTb", "singleTbNanoList.txt", 
                  "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8") 
singleT = sample("singleT", "singleTNanoList.txt", 
                  "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8") 
ttbar = sample("ttbar", "ttbarNanoList.txt", 
                  "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8") 
ttjetsTb = sample("ttjetsTb", "ttjetsTbNanoList.txt", 
                  "TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8") 
ttjetsT = sample("ttjetsT", "ttjetsTNanoList.txt", 
                  "TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets200 = sample("WJets200", "WJets200NanoList.txt", 
                  "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets400 = sample("WJets400", "WJets400NanoList.txt", 
                  "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets600 = sample("WJets600", "WJets600NanoList.txt", 
                  "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets800 = sample("WJets800", "WJets800NanoList.txt", 
                  "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets1200 = sample("WJets1200", "WJets1200NanoList.txt", 
                  "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8") 
WJets2500 = sample("WJets2500", "WJets2500NanoList.txt", 
                  "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8") 
QCD200 = sample("QCD200", "QCD200NanoList.txt", 
                "QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD300 = sample("QCD300", "QCD300NanoList.txt", 
                "QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD500 = sample("QCD500", "QCD500NanoList.txt", 
                "QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD700 = sample("QCD700", "QCD700NanoList.txt", 
                "QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD1000 = sample("QCD1000", "QCD1000NanoList.txt", 
                "QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD1500 = sample("QCD1500", "QCD1500NanoList.txt", 
                "QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8")
QCD2000 = sample("QCD2000", "QCD2000NanoList.txt", 
                "QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8")

samples={
    "singleTb":singleTb,
    "singleT":singleT,
    "ttbar":ttbar,
    "ttjetsTb":ttjetsTb,
    "ttjetsT":ttjetsT,
    "WJets200":WJets200,
    "WJets400":WJets400,
    "WJets600":WJets600,
    "WJets800":WJets800,
    "WJets1200":WJets1200,
    "WJets2500":WJets2500,
    "QCD200":QCD200,
    "QCD300":QCD300,
    "QCD500":QCD500,
    "QCD700":QCD700,
    "QCD1000":QCD1000,
    "QCD1500":QCD1500,
    "QCD2000":QCD2000,

}