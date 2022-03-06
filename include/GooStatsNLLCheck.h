/*****************************************************************************/
// Author: Xuefeng Ding <xuefeng.ding.physics@gmail.com>
// Institute: Gran Sasso Science Institute, L'Aquila, 67100, Italy
// Date: 2018 April 7th
// Version: v1.0
// Description: GooStats, a statistical analysis toolkit that runs on GPU.
//
// All rights reserved. 2018 copyrighted.
/*****************************************************************************/
#ifndef GooStatsNLLCheck_H
#define GooStatsNLLCheck_H
#include "TNamed.h"
#include <map>
#include <string>
class TFile;
class GooStatsNLLCheck : public TNamed {
public:
  GooStatsNLLCheck() : _s_totLL(0) {}
  virtual ~GooStatsNLLCheck() {}
  static GooStatsNLLCheck *get();
  void init(const std::string &fname, const std::string &name);
  void new_LL(double _totLL);
  void new_LL_single(double _singleLL);
  void record_LL(int bin, double E, double M, double T, double LL);
  void record_species(int bin, const std::string &name, double T);
  void record_finalLL(double LL) { finalLL = LL; }
  void save();
  void print() const;
  struct Info {
    double E = 0;
    double M = 0;
    double T = 0;
    double LL = 0;
    std::map<std::string, double> compositions;
  };
  const std::vector<std::map<int, GooStatsNLLCheck::Info>> &get_results() const { return results; }
  const std::vector<double> &get_totLL() const { return totLL; }
  double get_finalLL() const { return finalLL; }

private:
  std::vector<double> totLL;
  std::vector<std::map<int, GooStatsNLLCheck::Info>> results;
  double finalLL;

  TFile *file;                                   //! do not save it
  std::map<int, GooStatsNLLCheck::Info> result;  //! temporary
  double _s_totLL;                               //! temporal
  ClassDef(GooStatsNLLCheck, 1)
};
#endif
