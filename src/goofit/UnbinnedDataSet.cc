#include "goofit/UnbinnedDataSet.h"

#include "goofit/Variable.h"

// Special constructor for one variable
UnbinnedDataSet::UnbinnedDataSet(Variable *var, std::string n) : DataSet(var, n) {}

UnbinnedDataSet::UnbinnedDataSet(std::vector<Variable *> &vars, std::string n) : DataSet(vars, n) {}

UnbinnedDataSet::UnbinnedDataSet(std::set<Variable *> &vars, std::string n) : DataSet(vars, n) {}

UnbinnedDataSet::~UnbinnedDataSet() = default;

fptype UnbinnedDataSet::getValue(Variable *var, int idx) const {
  if (idx >= getNumEvents()) {
    std::cout << "UnbinnedDataSet error: Attempt to find " << var->name << " in event " << idx << " when only "
              << getNumEvents() << " events exist.\n";
    return -1;
  }
  if (0 > idx) {
    std::cout << "UnbinnedDataSet error: Cannot return value of " << var->name << " for index " << idx << std::endl;
    return -1;
  }

  auto event = data[idx].find(var);
  if (event == data[idx].end()) {
    static std::map<std::string, bool> printed;
    if (!printed[var->name]) {
      std::cout << "UnbinnedDataSet error: Could not find variable " << var->name << " in event " << idx
                << " of data set " << getName() << ". Returning -1. This message will only be printed once.\n";
      printed[var->name] = true;
    }
    return -1;
  }

  return (*event).second;
}

void UnbinnedDataSet::addEventVector(std::vector<fptype> &vals, fptype /*weight*/) {
  // NB, unbinned data set ignores weights.
  numEventsAdded++;

  auto currVar = varsBegin();
  std::map<Variable *, fptype> currEvent;
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

    currEvent[*currVar] = currVal;
    ++currVar;
  }
  assert(currVar == varsEnd());
  data.push_back(currEvent);
}

void UnbinnedDataSet::loadEvent(int idx) {
  if (idx >= getNumEvents()) {
    std::cout << "UnbinnedDataSet error: Attempt to load event " << idx << " when only " << getNumEvents()
              << " events exist.\n";
    return;
  }
  if (0 > idx) {
    std::cout << "UnbinnedDataSet error: Cannot load event " << idx << std::endl;
    return;
  }

  for (auto &v : data[idx]) {
    v.first->value = v.second;
  }
}

void UnbinnedDataSet::setValueForAllEvents(Variable *var) {
  for (auto &i : data) {
    i[var] = var->value;
  }
}
