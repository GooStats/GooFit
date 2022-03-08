#include "fit.h"
#include "gtest/gtest.h"

TEST(GooFit, testFitResult) {
  double a;
  double b;
  GooFit::fit(a, b);
  EXPECT_NEAR(-1. / 370, a, 1. / 370 * 2e-2);
  EXPECT_NEAR(1000, b, 1000 * 2e-2);
  int c[3];
  std::cout << c[3] << std::endl;
}
