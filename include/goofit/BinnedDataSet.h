#ifndef BINNED_DATASET_HH
#define BINNED_DATASET_HH

#include "goofit/DataSet.h"
#include <map>

class BinnedDataSet : public DataSet {
  // Class for rectangularly binned datasets - every bin the same size.

public:
  BinnedDataSet(Variable *var, std::string n = "");
  BinnedDataSet(std::vector<Variable *> &vars, std::string n = "");
  BinnedDataSet(std::set<Variable *> &vars, std::string n = "");
  ~BinnedDataSet();

  virtual void addEventVector(std::vector<fptype> &vals, fptype weight);

  fptype getBinContent(unsigned int bin) const { return binvalues[bin]; }
  fptype getBinCenter(Variable *var, unsigned int bin) const;
  unsigned int getBinNumber() const;
  fptype getBinVolume(unsigned int bin) const;
  fptype getBinError(unsigned int bin) const;
  unsigned int getNumBins() const;
  fptype getNumEvents() const;

  void setBinContent(unsigned int bin, fptype value) { binvalues[bin] = value; }
  void setBinError(unsigned int bin, fptype error);

private:
  void cacheNumBins();
  std::vector<unsigned int> convertValuesToBins(const std::vector<fptype> &vals) const;
  unsigned int localToGlobal(std::vector<unsigned int> &locals) const;
  void globalToLocal(std::vector<unsigned int> &locals, unsigned int global) const;

  std::vector<fptype> binvalues;
  std::vector<fptype> binerrors;
  std::map<Variable *, int>
      cachedNumBins;  // Store these numbers in case they change on the user end - vast confusion possible.
};

#endif
