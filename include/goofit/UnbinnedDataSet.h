#ifndef UNBINNED_DATASET_HH
#define UNBINNED_DATASET_HH

#include "goofit/DataSet.h"

class UnbinnedDataSet : public DataSet {
  // Class for unbinned datasets. 

public: 
  UnbinnedDataSet (Variable* var, std::string n = "");
  UnbinnedDataSet (std::vector<Variable*>& vars, std::string n = "");
  UnbinnedDataSet (std::set<Variable*>& vars, std::string n = "");
  ~UnbinnedDataSet (); 

  virtual void addEventVector (std::vector<fptype>& vals, fptype weight = 1); 
  int getNumEvents () const {return data.size();} 
  fptype getValue (Variable* var, int idx) const; 
  void loadEvent (int idx); 
  void setValueForAllEvents (Variable* var); 

private:
  std::vector<std::map<Variable*, fptype> > data; 

};

#endif
