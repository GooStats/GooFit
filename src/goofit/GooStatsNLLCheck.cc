#include "GooStatsNLLCheck.h"
#include "TFile.h"
ClassImp(GooStatsNLLCheck)  // NOLINT
    GooStatsNLLCheck *GooStatsNLLCheck::get() {
  static GooStatsNLLCheck *me = nullptr;
  if (!me)
    me = new GooStatsNLLCheck();
  return me;
}
void GooStatsNLLCheck::init(const std::string &fname, const std::string &name) {
  file = TFile::Open(fname.c_str(), "RECREATE");
  this->SetName(name.c_str());
  this->SetTitle(name.c_str());
}
void GooStatsNLLCheck::new_LL(double _totLL) {
  totLL.push_back(_totLL - _s_totLL);
  _s_totLL = _totLL;
  results.push_back(result);
  result.clear();
}
void GooStatsNLLCheck::new_LL_single(double _singleLL) {
  totLL.push_back(_singleLL);
  _s_totLL += _singleLL;
  results.push_back(result);
  result.clear();
}
void GooStatsNLLCheck::record_LL(int bin, double E, double M, double T, double LL) {
  result[bin].E = E;
  result[bin].M = M;
  result[bin].T = T;
  result[bin].LL = LL;
}
void GooStatsNLLCheck::record_species(int bin, const std::string &name, double T) {
  result[bin].compositions[name] = T;
}
void GooStatsNLLCheck::save() {
  file->cd();
  this->Write();
  file->Close();
  delete file;
  file = nullptr;
}
#include <iostream>
template <typename... Args>
std::string string_format(const std::string &format, Args... args) {
  int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) + 1;  // Extra space for '\0'
  if (size_s <= 0) {
    throw std::runtime_error("Error during formatting.");
  }
  auto size = static_cast<size_t>(size_s);
  auto buf = std::unique_ptr<char[]>(new char[size]);
  std::snprintf(buf.get(), size, format.c_str(), args...);
  return {buf.get(), buf.get() + size - 1};  // We don't want the '\0' inside
}
void GooStatsNLLCheck::print() const {
  for (auto ele : get_results()) {
    for (auto bin : ele) {
      std::cout << string_format(
                       "log(L) %.12le b %lf M %lf tot %.12le", bin.second.LL, bin.second.E, bin.second.M, bin.second.T)
                << std::endl;
      for (auto spc : bin.second.compositions) {
        std::cout << string_format(" %s %.12le", spc.first.c_str(), spc.second);
      }
      std::cout << std::endl;
    }
  }
  for (auto LL : get_totLL()) {
    std::cout << string_format("log(L) %.12le", LL) << std::endl;
  }
  std::cout << string_format("final log(L) %.12le", finalLL) << std::endl;
}
