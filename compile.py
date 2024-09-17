import ROOT


ROOT.gSystem.CompileMacro("functions.cpp", "kO")


ROOT.gInterpreter.Load("functions_cpp.so")
