#include "goofit/PdfBase.h"
#include "goofit/GlobalCudaDefines.h"
#include "goofit/Variable.h"
#include <algorithm>
#include <map>
#include <set>
#include <utility>

fptype *dev_event_array[maxDataSet];
fptype host_normalisation[maxIndicies];
fptype host_params[maxParams];
unsigned int host_indices[maxIndicies];
unsigned int totalParams = 0;
unsigned int totalConstants = 1;  // First constant is reserved for number of events.
std::map<Variable *, std::set<PdfBase *>> variableRegistry;
std::map<PdfBase *, int> pdfIdMap;
std::string pdfName[maxIndicies];

PdfBase::PdfBase(Variable *x, std::string n)
    : name(std::move(n)),
      numEvents(0),
      numEntries(0),
      normRanges(nullptr),
      fitControl(nullptr),
      integrationBins(-1),
      specialMask(0),
      cachedParams(nullptr),
      properlyInitialised(true)  // Special-case PDFs should set to false.
      ,
      pdfId(-1) {
  if (x != nullptr) {
    registerObservable(x);
  }
}

__host__ void PdfBase::checkInitStatus(std::vector<std::string> &unInited) const {
  if (!properlyInitialised) {
    unInited.push_back(getName());
  }
  for (auto *component : components) {
    component->checkInitStatus(unInited);
  }
}

__host__ void PdfBase::recursiveSetNormalisation(fptype norm) const {
  host_normalisation[parameters] = norm;
  for (auto *component : components) {
    component->recursiveSetNormalisation(norm);
  }
}

int pdfId_global = -1;
__host__ int PdfBase::registerPdf() {
  if (pdfIdMap.find(this) == pdfIdMap.end()) {
    pdfIdMap.insert(std::make_pair(this, ++pdfId_global));
    pdfId = pdfId_global;
  }
  return pdfIdMap.at(this);
}
__host__ unsigned int PdfBase::registerParameter(Variable *var) {
  if (var == nullptr) {
    std::cout << "Error: Attempt to register null Variable with " << getName() << ", aborting.\n";
    assert(var);
    exit(1);
  }

  if (std::find(parameterList.begin(), parameterList.end(), var) != parameterList.end()) {
    return static_cast<unsigned int>(var->getIndex());
  }

  parameterList.push_back(var);
  variableRegistry[var].insert(this);
  if (0 > var->getIndex()) {
    unsigned int unusedIndex = 0;
    while (true) {
      bool canUse = true;
      for (auto &p : variableRegistry) {
        if (static_cast<int>(unusedIndex) != p.first->index) {
          continue;
        }
        canUse = false;
        break;
      }
      if (canUse) {
        break;
      }
      unusedIndex++;
    }

    var->index = unusedIndex;
  }
  return static_cast<unsigned int>(var->getIndex());
}

__host__ void PdfBase::unregisterParameter(Variable *var) {
  if (var == nullptr) {
    return;
  }
  auto pos = std::find(parameterList.begin(), parameterList.end(), var);
  if (pos != parameterList.end()) {
    parameterList.erase(pos);
  }
  variableRegistry[var].erase(this);
  if (variableRegistry[var].empty()) {
    var->index = -1;
  }
  for (auto &component : components) {
    component->unregisterParameter(var);
  }
}

__host__ void PdfBase::getParameters(parCont &ret) const {
  for (auto *p : parameterList) {
    if (std::find(ret.begin(), ret.end(), p) != ret.end()) {
      continue;
    }
    ret.push_back(p);
  }
  for (auto *component : components) {
    component->getParameters(ret);
  }
}

__host__ Variable *PdfBase::getParameterByName(std::string n) const {
  for (auto *p : parameterList) {
    if (p->name == n) {
      return p;
    }
  }
  for (auto *component : components) {
    Variable *cand = component->getParameterByName(n);
    if (cand != nullptr) {
      return cand;
    }
  }
  return nullptr;
}

__host__ void PdfBase::getObservables(std::vector<Variable *> &ret) const {
  for (auto p = obsCBegin(); p != obsCEnd(); ++p) {
    if (std::find(ret.begin(), ret.end(), *p) != ret.end()) {
      continue;
    }
    ret.push_back(*p);
  }
  for (auto *component : components) {
    component->getObservables(ret);
  }
}

__host__ unsigned int PdfBase::registerConstants(unsigned int amount) {
  assert(totalConstants + amount < maxConsts);
  cIndex = totalConstants;
  totalConstants += amount;
  return cIndex;
}

void PdfBase::registerObservable(Variable *obs) {
  if (obs == nullptr) {
    return;
  }
  if (find(observables.begin(), observables.end(), obs) != observables.end()) {
    return;
  }
  //#ifdef DEBUG_GPU
  //  cout<<getName()<<" register observable: "<<obs->name<<" add "<<obs<<" index "<<obs->index<<" val "<<obs->value<<endl;
  //#endif
  observables.push_back(obs);
}

__host__ void PdfBase::setIntegrationFineness(int i) {
  integrationBins = i;
  generateNormRange();
}

__host__ bool PdfBase::parametersChanged() const {
  if (cachedParams == nullptr) {
    return true;
  }

  parCont params;
  getParameters(params);
  int counter = 0;
  for (auto &param : params) {
    if (cachedParams[counter++] != host_params[param->index]) {
      return true;
    }
  }

  return false;
}

__host__ void PdfBase::storeParameters() const {
  parCont params;
  getParameters(params);
  if (cachedParams == nullptr) {
    cachedParams = new fptype[params.size()];
  }

  int counter = 0;
  for (auto &param : params) {
    cachedParams[counter++] = host_params[param->index];
  }
}

void dummySynch() {}
