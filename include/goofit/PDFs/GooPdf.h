#ifndef THRUST_PDF_FUNCTOR_HH
#define THRUST_PDF_FUNCTOR_HH
#include "goofit/PdfBase.h"
#include <cmath>

#define CALLS_TO_PRINT 10

EXEC_TARGET int dev_powi(int base, int exp);  // Implemented in SmoothHistogramPdf.
void *getMetricPointer(std::string name);

typedef fptype (*device_function_ptr)(fptype *,
                                      fptype *,
                                      unsigned int *);  // Pass event, parameters, index into parameters.
typedef fptype (*device_metric_ptr)(fptype, fptype *, unsigned int);

extern void *host_fcn_ptr;

EXEC_TARGET fptype callFunction(fptype *eventAddress, unsigned int functionIdx, unsigned int paramIdx);

class MetricTaker;
class BinnedMetricTaker;

class GooPdf : public PdfBase {
public:
  GooPdf(Variable *x, std::string n);
  __host__ virtual double calculateNLL() const;
  __host__ double evaluateAtPoints(const std::vector<double> &v);
  __host__ void evaluateAtPoints(Variable *var, std::vector<fptype> &res);

  __host__ virtual fptype normalise() const;
  __host__ virtual fptype integrate(__attribute__((unused)) fptype lo, __attribute__((unused)) fptype hi) const {
    return 0;
  }

  __host__ virtual bool hasAnalyticIntegral() const { return false; }

  __host__ fptype getValue();
  __host__ void getCompProbsAtDataPoints(std::vector<std::vector<fptype>> &values);
  __host__ void initialise(std::vector<unsigned int> pindices, void *dev_functionPtr = host_fcn_ptr);
  __host__ void scan(Variable *var, std::vector<fptype> &values);
  __host__ virtual void setFitControl(FitControl *const fc, bool takeOwnerShip);
  __host__ bool IsChisquareFit() const;
  __host__ virtual void setMetrics();
  __host__ void setParameterConstantness(bool constant = true);

  __host__ virtual void transformGrid(fptype *host_output);
  static __host__ int findFunctionIdx(void *dev_functionPtr);

  __host__ void debug() const;

protected:
  __host__ virtual double sumOfNll(int numVars) const;
  MetricTaker *logger;              // for calculating Nll
  BinnedMetricTaker *binnedlogger;  // for calculating Nll
};

class MetricTaker : public thrust::unary_function<thrust::tuple<int, fptype *, int>, fptype> {
public:
  MetricTaker(PdfBase *dat, void *dev_functionPtr);
  EXEC_TARGET fptype operator()(thrust::tuple<int, fptype *, int> t)
      const;  // Event number, dev_event_array(pass this way for nvcc reasons), event size

private:
  unsigned int metricIndex;  // Function-pointer index of processing function, eg logarithm, chi-square, other metric.
  unsigned int functionIdx;  // Function-pointer index of actual PDF
  unsigned int parameters;
};
class BinnedMetricTaker : public thrust::unary_function<thrust::tuple<int, int, fptype *>, fptype> {
public:
  BinnedMetricTaker(PdfBase *dat, void *dev_functionPtr);
  EXEC_TARGET fptype operator()(thrust::tuple<int, int, fptype *> t) const;  // Event number, event size, low/up/n

private:
  unsigned int metricIndex;  // Function-pointer index of processing function, eg logarithm, chi-square, other metric.
  unsigned int functionIdx;  // Function-pointer index of actual PDF
  unsigned int parameters;
};

#endif
