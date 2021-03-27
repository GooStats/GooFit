#ifndef EXP_PDF_HH
#define EXP_PDF_HH

#include "goofit/PDFs/GooPdf.h" 

class ExpPdf : public GooPdf {
public:
  ExpPdf (std::string n, Variable* _x, Variable* alpha, Variable* offset = 0); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return (1 == host_indices[parameters]);} 



private:

};

#endif
