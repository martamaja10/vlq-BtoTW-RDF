import ROOT as gbl
import correctionlib
correctionlib.register_pyroot_binding()

f = gbl.TFile.Open("root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_RDF/BpM800_hadd.root")

df = gbl.ROOT.RDataFrame(f.Get("Events"))

gbl.gInterpreter.Declare('auto csetEl = correction::CorrectionSet::from_file("electron.json");')
gbl.gInterpreter.Declare('auto csetEl_2016preID = csetEl->at("UL-Electron-ID-SF");')

#sf
df = df.Filter("isEl >0").Define(
        "leadElSF",
        ('csetEl_2016preID->evaluate({"2018", "sf", "wp90noiso", '
            'std::abs(lepton_eta), lepton_phi})'))

#sfdown
df = df.Filter("isEl >0").Define(
        "leadElSFDown",
        ('csetEl_2016preID->evaluate({"2018", "sfdown", "wp90noiso", '
            'std::abs(lepton_eta), lepton_phi})'))

#sfup
df = df.Filter("isEl >0").Define(
        "leadElSFUp",
        ('csetEl_2016preID->evaluate({"2018", "sfup", "wp90noiso", '
            'std::abs(lepton_eta), lepton_phi})'))

h = df.Histo1D("leadElSF")
h.Draw()

#example of how to save code as "snapshot" from Hogan
#but it's cpp?
std::cout << "-------------------------------------------------" << std::endl << ">>> Saving " << sample << " Snapshot..." << std::endl;
  TString finalFile = "RDF_"+sample+"_finalsel_"+testNum+".root";
  const char* stdfinalFile = finalFile;
  postPresel.Snapshot("Events", stdfinalFile);
