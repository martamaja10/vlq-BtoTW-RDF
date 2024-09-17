#include "ROOT/RVec.hxx"
#include "TF1.h"
//#include "analyzer_RDF.h"

using namespace ROOT::VecOps;

// Pass the TF1 objects as function arguments
/*auto wjetHTpoly = [](TF1* poly2, TF1* poly2U, TF1* poly2D, float &LHE_HT) {
    RVec<double> sf = {poly2->Eval(LHE_HT), poly2U->Eval(LHE_HT), poly2D->Eval(LHE_HT)};
    return sf;
};

auto topHTpoly = [](TF1* polyHT, TF1* polyHTU, TF1* polyHTD, float &AK4HT) {
    RVec<double> sf = {polyHT->Eval(AK4HT), polyHTU->Eval(AK4HT), polyHTD->Eval(AK4HT)};
    return sf;
};*/
//polyHT = ROOT.TF1("polyHT", "min(1.0,max([0] + [1]*x,[2]))", 700, 5000)
//polyHTU = ROOT.TF1("polyHTU", "min(1.0,max([0] + [1]*x + sqrt([3] + 2*x*[4] + x*x*[5]),[2] + [6]))", 700, 5000)
//polyHTD = ROOT.TF1("polyHTD", "min(1.0,max([0] + [1]*x - sqrt([3] + 2*x*[4] + x*x*[5]),[2] - [6]))", 700, 5000)

RVec<double> topHTpoly(TF1* polyHT, TF1* polyHTU, TF1* polyHTD, float &AK4HT) {
    RVec<double> sf = {polyHT->Eval(AK4HT), polyHTU->Eval(AK4HT), polyHTD->Eval(AK4HT)};
    return sf;
}
