#include "goofit/PDFs/ExpPdf.hh"

EXEC_TARGET fptype device_Exp (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x = evt[indices[2 + indices[0]]]; 
  fptype alpha = p[indices[1]];

  fptype ret = EXP(alpha*x); 
  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_Exp = device_Exp; 

__host__ ExpPdf::ExpPdf (std::string n, Variable* _x, Variable* alpha, Variable* )
  : GooPdf(_x, n) 
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(alpha));
  GET_FUNCTION_ADDR(ptr_to_Exp);
  initialise(pindices); 
}
