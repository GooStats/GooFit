#include "goofit/FitManagerMinuit1.h"
#include "TMinuit.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/PdfBase.h"
#include "goofit/Variable.h"

PdfBase *pdfPointer;
int numPars = 0;
std::vector<Variable *> FitManager::vars;
bool FitManager::minim_conv = false;
bool FitManager::hesse_conv = false;

FitManager::FitManager(PdfBase *dat) : minuit(0), overrideCallLimit(-1) { pdfPointer = dat; }

FitManager::~FitManager() {
  if (minuit)
    delete minuit;
}

void FitManager::setupMinuit() {
  vars.clear();
  pdfPointer->getParameters(vars);

  numPars = vars.size();
  if (minuit)
    delete minuit;
  minuit = new TMinuit(numPars);
  int maxIndex = 0;
  int counter = 0;
  for (std::vector<Variable *>::iterator i = vars.begin(); i != vars.end(); ++i) {
    minuit->DefineParameter(counter, (*i)->name.c_str(), (*i)->value, (*i)->error, (*i)->lowerlimit, (*i)->upperlimit);
    if ((*i)->fixed)
      minuit->FixParameter(counter);
    counter++;
    if (maxIndex < (*i)->getIndex())
      maxIndex = (*i)->getIndex();
  }

  numPars = maxIndex + 1;
  pdfPointer->copyParams();
  minuit->SetFCN(FitFun);
  if (!static_cast<GooPdf *>(pdfPointer)->IsChisquareFit())
    minuit->SetErrorDef(0.5);
}

void FitManager::fit() {
  setupMinuit();
  runMigrad();
}

void FitManager::runMigrad() {
  assert(minuit);
  //  int ierflg_simplex;
  int ierflg_minimize;
  //  int ierflg_migrad;
  int ierflg_hesse;
  //  minuit->mnexcm("SIMPLEX", 0,0,ierflg_simplex);
  if (0 < overrideCallLimit) {
    //std::cout << "Calling MIGRAD with call limit " << overrideCallLimit << std::endl;
    //    std::cout << "Calling MINIMIZE with call limit " << overrideCallLimit << std::endl;
    double plist[1];
    plist[0] = overrideCallLimit;
    minuit->mnexcm("MINIMIZE", plist, 1, ierflg_minimize);
    minim_conv = ierflg_minimize == 0;
  } else
    minuit->Migrad();
  //  minuit->mnexcm("SEEK", 0,0,ierflg_hesse);
  //  minuit->mnexcm("SIMPLEX", 0,0,ierflg_simplex);
  minuit->mnexcm("HESSE", 0, 0, ierflg_hesse);
  hesse_conv = ierflg_hesse == 0;
  //  int ierflg;
  //  minuit->mnexcm("MINOS", 0,0,ierflg);
  if (minim_conv)
    ;  //std::cout<<"MINIMIZE exit successfully"<<std::endl;
  else
    std::cout << "MINIMIZE did not exit normally" << std::endl;
  if (hesse_conv)
    ;  //std::cout<<"HESSE exit successfully"<<std::endl;
  else
    std::cout << "HESSE did not exit normally" << std::endl;
}

void FitManager::getMinuitValues() const {
  int counter = 0;
  for (std::vector<Variable *>::iterator i = vars.begin(); i != vars.end(); ++i) {
    minuit->GetParameter(counter++, (*i)->value, (*i)->error);
  }
}

void FitFun(int &__attribute__((__unused__)) npar,
            double *__attribute__((__unused__)) gin,
            double &fun,
            double *fp,
            int __attribute__((__unused__)) iflag) {
  std::vector<double> pars;
  // Notice that npar is number of variable parameters, not total.
  pars.resize(numPars);
  int counter = 0;
  for (auto i : FitManager::vars) {
    if (std::isnan(fp[counter]))
      std::cout << "Variable " << i->name << " " << i->index << " is NaN\n";
    pars[i->getIndex()] = fp[counter++] + i->blind;
  }

  pdfPointer->copyParams(pars);
  fun = pdfPointer->calculateNLL();
}
