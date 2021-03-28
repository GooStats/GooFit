#ifndef EXP_PDF_HH
#define EXP_PDF_HH

#include "goofit/PDFs/GooPdf.h" 

class ExpPdf : public GooPdf {
public:
  ExpPdf (std::string n, Variable* _x, Variable* alpha, Variable* offset = 0); 


private:

};

#endif
