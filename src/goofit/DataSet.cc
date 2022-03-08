#include "goofit/DataSet.h"

#include "goofit/Variable.h"

DataSet::DataSet(Variable *var, std::string n) : numEventsAdded(0), name(n) {
  variables.push_back(var);
  if (n.empty()) {
    generateName();
  }
}

DataSet::DataSet(std::vector<Variable *> &vars, std::string n) : numEventsAdded(0), name(n) {
  for (auto &var : vars) {
    variables.push_back(var);
  }
  if (n.empty()) {
    generateName();
  }
}

DataSet::DataSet(std::set<Variable *> &vars, std::string n) : numEventsAdded(0), name(n) {
  variables.resize(vars.size());

  for (auto *var : vars) {
    variables[var->index] = var;
  }
  if (n.empty()) {
    generateName();
  }
}

DataSet::~DataSet() { variables.clear(); }

void DataSet::addEvent() {
  std::vector<fptype> vals = getCurrentValues();
  addEventVector(vals);
}

void DataSet::addWeightedEvent(fptype weight) {
  std::vector<fptype> vals = getCurrentValues();
  addEventVector(vals, weight);
}

void DataSet::addEvent(fptype val) {
  // Helper method to avoid the user having to wrap
  // every one-dimensional event in a std::vector.
  assert(1 == variables.size());

  std::vector<fptype> helper;
  helper.push_back(val);
  addEventVector(helper);
}

std::vector<fptype> DataSet::getCurrentValues() const {
  std::vector<fptype> values;
  for (auto v = varsBegin(); v != varsEnd(); ++v) {
    values.push_back((*v)->value);
  }
  return values;
}

void DataSet::getVariables(std::vector<Variable *> &vars) {
  for (auto &variable : variables) {
    vars.push_back(variable);
  }
}

unsigned int DataSet::indexOfVariable(Variable *var) const {
  for (unsigned int i = 0; i < variables.size(); ++i) {
    if (var != variables[i]) {
      continue;
    }
    return i;
  }

  std::cout << "Error: Attempted to get index of variable " << var->name << " in DataSet of ";
  for (auto v = varsBegin(); v != varsEnd(); ++v) {
    std::cout << "\n  " << (*v)->name << std::endl;
  }
  std::cout << "\nAborting." << std::endl;
  abort();
}

void DataSet::generateName() {
  // Create default name as list of variables.
  if (!name.empty()) {
    return;
  }
  for (auto v = varsBegin(); v != varsEnd(); ++v) {
    name += (*v)->name;
    auto next = v;
    next++;
    if (next == varsEnd()) {
      continue;
    }
    name += ", ";
  }
}
