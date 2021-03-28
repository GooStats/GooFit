#include "gtest/gtest.h"
#include "fit.h"

TEST (GooFit, testFitResult) {
  double a,b;
  GooFit::fit(a,b);
  EXPECT_NEAR(-1./370, a,1./370*2e-2);
  EXPECT_NEAR(1000, b,1000*2e-2);
}
