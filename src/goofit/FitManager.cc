#include "goofit/FitManager.h"

#if MINUIT_VERSION == 1
#include "FitManagerMinuit1.cc"
#elif MINUIT_VERSION == 2
#include "FitManagerMinuit2.cc"
#else 
#include "FitManagerMinuit3.cc"
#endif 
