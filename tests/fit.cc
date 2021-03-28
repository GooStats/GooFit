#include "fit.h"
#include "goofit/Variable.h"
#include "goofit/BinnedDataSet.h"
#include "TRandom.h"
#include "goofit/PDFs/ExpPdf.hh"
#include "goofit/PDFs/SumPdf.h"
#include "goofit/FitManager.h"
namespace GooFit {
  void fit(double &a,double &b) {
    Variable *energy = new Variable("energy",1,100); // nbins = 100
    BinnedDataSet *data = new BinnedDataSet(energy);
    gRandom->SetSeed(1);
    for(int i = 0;i<100;++i) {
      data->setBinContent(i,gRandom->Poisson(1000*exp(-(i+0.5)/370)));
    }
    Variable *lambda = new Variable("lambda",-1./400,-1./40,-1./4000);
    GooPdf *pdf = new ExpPdf("exp",energy,lambda);
    std::vector<Variable*> Ns;
    Ns.push_back(new Variable("N",1000,1,0,10000));
    std::vector<PdfBase*> pdfs;
    pdfs.push_back(pdf);
    SumPdf *likelihood = new SumPdf("likelihood",1,Ns,pdfs,energy);
    likelihood->setData(data);

    FitManager *fit = new FitManager(likelihood);
    fit->fit();
    fit->getMinuitValues();

    a = lambda->value;
    b = Ns.front()->value;
  }
}
