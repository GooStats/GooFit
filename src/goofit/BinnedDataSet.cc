#include "goofit/BinnedDataSet.h"

#include "goofit/Variable.h"

// Special constructor for one variable
BinnedDataSet::BinnedDataSet(Variable *var, std::string n) : DataSet(var, n) {
  cacheNumBins();
  binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(std::vector<Variable *> &vars, std::string n) : DataSet(vars, n) {
  cacheNumBins();
  binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(std::set<Variable *> &vars, std::string n) : DataSet(vars, n) {
  cacheNumBins();
  binvalues.resize(getNumBins());
}

BinnedDataSet::~BinnedDataSet() = default;

void BinnedDataSet::addEventVector(std::vector<fptype> &vals, fptype weight) {
  numEventsAdded++;
  std::vector<unsigned int> localBins = convertValuesToBins(vals);
  unsigned int bin = localToGlobal(localBins);
  if (bin >= binvalues.size()) {
    std::cout << "Bad bin number " << bin << " / " << binvalues.size() << " from input values ";
    for (double val : vals) {
      std::cout << val << " ";
    }
    std::cout << "\n";
    assert(bin < binvalues.size());
  }
  binvalues[bin] += weight;
}

void BinnedDataSet::cacheNumBins() {
  for (auto v = varsBegin(); v != varsEnd(); ++v) {
    cachedNumBins[*v] = (*v)->numbins;
  }
}

unsigned int BinnedDataSet::getBinNumber() const {
  std::vector<fptype> vals = getCurrentValues();
  std::vector<unsigned int> locals = convertValuesToBins(vals);
  return localToGlobal(locals);
}

unsigned int BinnedDataSet::localToGlobal(std::vector<unsigned int> &locals) const {
  unsigned int priorMatrixSize = 1;
  unsigned int ret = 0;
  for (auto v = varsBegin(); v != varsEnd(); ++v) {
    unsigned int localBin = locals[indexOfVariable(*v)];
    ret += localBin * priorMatrixSize;
    priorMatrixSize *= cachedNumBins.at(*v);  // Use 'at' to preserve const-ness.
  }
  return ret;
}

void BinnedDataSet::globalToLocal(std::vector<unsigned int> &locals, unsigned int global) const {
  locals.clear();

  // To convert global bin number to (x,y,z...) coordinates: For each dimension, take the mod
  // with the number of bins in that dimension. Then divide by the number of bins, in effect
  // collapsing so the grid has one fewer dimension. Rinse and repeat.
  for (auto v = varsBegin(); v != varsEnd(); ++v) {
    int localBin = global % cachedNumBins.at(*v);
    locals.push_back(localBin);
    global /= cachedNumBins.at(*v);
  }
}

fptype BinnedDataSet::getBinCenter(Variable *var, unsigned int bin) const {
  std::vector<unsigned int> locals;
  globalToLocal(locals, bin);
  unsigned int varIndex = indexOfVariable(var);
  unsigned int localBin = locals[varIndex];
  fptype ret = var->upperlimit;
  ret -= var->lowerlimit;
  ret /= cachedNumBins.at(var);
  ret *= (localBin + 0.5);
  ret += var->lowerlimit;
  return ret;
}

fptype BinnedDataSet::getBinVolume(unsigned int __attribute__((__unused__)) bin) const {
  fptype ret = 1;
  for (auto v = varsBegin(); v != varsEnd(); ++v) {
    fptype step = (*v)->upperlimit;
    step -= (*v)->lowerlimit;
    step /= cachedNumBins.at(*v);
    ret *= step;
  }
  return ret;
}

fptype BinnedDataSet::getBinError(unsigned int bin) const {
  if (binerrors.empty()) {
    return SQRT(binvalues[bin]);
  }
  assert(bin < binerrors.size());
  return binerrors[bin];
}

void BinnedDataSet::setBinError(unsigned int bin, fptype error) {
  if (binerrors.empty()) {
    binerrors.resize(binvalues.size());
  }
  assert(bin < binerrors.size());
  binerrors[bin] = error;
}

unsigned int BinnedDataSet::getNumBins() const {
  unsigned int ret = 1;
  for (auto v = varsBegin(); v != varsEnd(); ++v) {
    ret *= cachedNumBins.at(*v);
  }
  return ret;
}

fptype BinnedDataSet::getNumEvents() const {
  fptype ret = 0;
  for (double binvalue : binvalues) {
    ret += binvalue;
  }
  return ret;
}

std::vector<unsigned int> BinnedDataSet::convertValuesToBins(const std::vector<fptype> &vals) const {
  std::vector<unsigned int> localBins;
  auto currVar = varsBegin();
  for (double currVal : vals) {
    assert(currVar != varsEnd());
    if (currVal < (*currVar)->lowerlimit) {
      std::cout << "Warning: Value " << currVal << " less than minimum " << (*currVar)->lowerlimit << " for "
                << (*currVar)->name << "; clamping to minimum.\n";
      currVal = (*currVar)->lowerlimit;
    }
    if (currVal > (*currVar)->upperlimit) {
      std::cout << "Warning: Value " << currVal << " more than maximum " << (*currVar)->upperlimit << " for "
                << (*currVar)->name << "; clamping to maximum.\n";
      currVal = (*currVar)->upperlimit;
    }

    fptype step = (*currVar)->upperlimit;
    step -= (*currVar)->lowerlimit;
    step /= cachedNumBins.at(*currVar);

    fptype curr = currVal;
    curr -= (*currVar)->lowerlimit;
    curr /= step;
    localBins.push_back(static_cast<unsigned int>(floor(curr)));

    ++currVar;
  }
  return localBins;
}
