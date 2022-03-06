#include <cmath>

#include "GooStatsNLLCheck.h"
#include "TFile.h"
#include "fit.h"
#include "gtest/gtest.h"

TEST (GooFit, NLLTest) {
  double a = NAN;
  double b = NAN;
  GooStatsNLLCheck::get()->init("NLL_CHECK_test.root","test");
  GooFit::fit(a,b);
  GooStatsNLLCheck::get()->new_LL_single(0);
  GooStatsNLLCheck::get()->save();
  auto results = GooStatsNLLCheck::get()->get_results();
  TFile *reference_f = TFile::Open("NLL_CHECK_reference.root");
  auto *ref_obj = dynamic_cast<GooStatsNLLCheck*>(reference_f->Get("test"));
  auto reference = (ref_obj->get_results())[0];
  for(auto x: results[0]) {
    EXPECT_NEAR(x.second.LL, reference[x.first].LL,x.second.LL*5e-11);
  }
}
